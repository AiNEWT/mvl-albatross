# HGLM Run ####

observeEvent(input$re_hglm_g_run, {
    g_resetresult <<- TRUE
})

re_hglm_results <- eventReactive(input$re_hglm_g_run, {
    # ¦§ Warning & Notify ####
    # if user doesn't input some essential part, we wil return null value and show notification
    if (input$re_hglm_m_resp == "")  {
        showNotification("Please choose response", type="warning")
        return()
    }
    
    if(input$re_hglm_m_check_offset && input$re_hglm_m_offset == "") {
        showNotification("Please choose offset variable", type="warning")
        return()
    }
    
    if(!is.null(input$re_hglm_m_check_binomd) && input$re_hglm_m_check_binomd && input$re_hglm_m_binomd == "") {
        showNotification("Please choose binomial denominator", type="warning")
        return()
    }
    
    # Start Progressbar
    withProgress(message = 'HGLM', style = "notification", value = 0.1, {
    removeTab(inputId = "re_hglm_resulttabset", target = "Comparison Model")
    
    # ¦§ Variable Declaration 1 ####
    # variable declaration for which is used widely in re_hglm_results
    dataLength = nrow(data2)
    
    # variable declaration for which is used in making MM(MeanModel), DM(DispersionModel)
    meanLink = input$re_hglm_m_link
    phiLink = "log"
    lambdaLink = NULL
    # In the past, they don't use input for phiLink and lambdaLink. why? 
    if(!is.null(input$re_hglm_p_check_phi) && input$re_hglm_p_check_phi)
        phiLink = input$re_hglm_p_link
    if(!is.null(input$re_hglm_l_check_lambda) && input$re_hglm_l_check_lambda)
        lambdaLink = input$re_hglm_l_link
    
    meanFormula = input$re_hglm_m_model
    phiFormula = "constant"               
    lambdaFormula = "lambda ~ 1"  
    
    if(!is.null(input$re_hglm_p_check_phi) && input$re_hglm_p_check_phi)
        phiFormula = input$re_hglm_p_model
    if(!is.null(input$re_hglm_l_check_lambda) && input$re_hglm_l_check_lambda)
        lambdaFormula = input$re_hglm_l_model
    
    meanNumberRandom = length(strsplit(meanFormula, "|", fixed = TRUE)[[1]]) - 1
    phiNumberRandom = 0
    lambdaNumberRandom = 0
    
    meanRandFamily = NULL
    phiRandFamily = NULL
    lambdaRandFamily = "gaussian"
    
    if(meanNumberRandom > 0)
        for (i in 1:meanNumberRandom) 
            meanRandFamily = c(meanRandFamily, input[[paste0("re_hglm_m_randfamily_", i)]])

    meanOffsetVariable = NULL
    phiOffsetVariable = NULL
    # dhglmfit doesn't provide offset variable for lambda. 
    
    if(input$re_hglm_m_check_offset)
        meanOffsetVariable = data2[[input$re_hglm_m_offset]]
    if(!is.null(input$re_hglm_p_check_offset) && input$re_hglm_p_check_offset)
        phiOffsetVariable = data2[[input$re_hglm_p_offset]]
    
    Lmatrix = NULL
    spatial = NULL
    longitude = NULL
    latitude = NULL
    spline = NULL
    adjacent = NULL # variable for neighborhood
    # LinPredRandVariance = NULL      # same as lambdaFormula
    # RandDistRandVariance = "gaussian" # same as lambdaRandFamily
    
    phiCorr = NULL # for GARCH and AR(1) model
    
    # if IAR and MRF, we set spatial,adjacent and revise meanRandFamily
    if (meanNumberRandom > 0) {
        if (meanRandFamily[1] == "IAR" || meanRandFamily[1] == "MRF") {
            spatial = meanRandFamily[1]
            meanRandFamily = rep("gaussian", meanNumberRandom)
            
            neighborfile = re_hglm_m_neighborfile()
            adjacent <- matrix(0, ncol = dataLength, nrow = dataLength)
            for (i in 1:dataLength) {
                neigh = as.numeric(strsplit(as.character(neighborfile[i, 2]), ",", fixed = TRUE)[[1]])
                neighLength = length(neigh)
                for (j in 1:neighLength) {
                    location <- neigh[j]
                    adjacent[i, location] <- 1
                }
            }
        }
    } 
    
    #if Matern, we set spatial, longitude, latitude and revise meanRandFamily.
    if (meanNumberRandom > 0) {
        if (meanRandFamily[1] == "Matern") {
            spatial = meanRandFamily[1]
            meanRandFamily <- rep("gaussian", meanNumberRandom)
            longitude = (data2 %>% dplyr::select("longitude"))[[1]]
            latitude = (data2 %>% dplyr::select("latitude"))[[1]]
        }
    }
    
    #if random work or local linear trend of quarterly seasonal, we set Lmatrix and revise meanRandFamily.
    if (meanNumberRandom > 0) {
        if (meanRandFamily[1] == "random work" ||
            meanRandFamily[1] == "local linear trend" ||
            meanRandFamily[1] == "quarterly seasonal") {
            
            Lmatrix_RandomWork <- matrix(0, ncol = dataLength, nrow = dataLength)
            Lmatrix_LocalLinearTrend <- matrix(0, ncol = dataLength, nrow = dataLength)
            wt <- matrix(0, ncol = dataLength, nrow = dataLength)
            Lmatrix_QuarterlySeasonal <- matrix(0, ncol = dataLength, nrow = dataLength)
            
            for (i in 1:dataLength) 
                for (j in 1:i) 
                    Lmatrix_RandomWork[i, j] <- 1
            for (i in 1:dataLength) 
                for (j in 1:i) 
                    Lmatrix_LocalLinearTrend[i, j] <- (i - j + 1)
            for (i in 1:dataLength) { #wt : weighted matrix for quarterly seasonal
                if (i == 1)
                    wt[i, 1] <- 1
                if (i == 2) {
                    wt[i, 1] <- 1
                    wt[i, 2] <- 1
                }
                if (i == 3) {
                    wt[i, 1] <- 1
                    wt[i, 2] <- 1
                    wt[i, 3] <- 1
                }
                if (i > 3) {
                    for (j in 1:4) {
                        ii <- i - j + 1
                        wt[i, ii] <- 1
                    }
                }
            }
            Lmatrix_QuarterlySeasonal = solve(wt)
            
            if (meanRandFamily[1] == "random work") Lmatrix = list(Lmatrix_RandomWork)
            if (meanRandFamily[1] == "local linear trend") Lmatrix = list(Lmatrix_LocalLinearTrend)
            if (meanRandFamily[1] == "quarterly seasonal") Lmatrix = list(Lmatrix_QuarterlySeasonal)
            
            
            #bad algorithm
            #algorithm is not good, but, I can't revise because I don't know dhglmfit_Lmatrix specifically
            if (meanNumberRandom >= 2) {
                for (i in 2:meanNumberRandom) {
                    if (meanRandFamily[i] == "random work") Lmatrix1 <- Lmatrix_RandomWork
                    if (meanRandFamily[i] == "local linear trend") Lmatrix1 <- Lmatrix_LocalLinearTrend
                    if (meanRandFamily[i] == "quarterly seasonal") Lmatrix1 <- Lmatrix_QuarterlySeasonal
                    Lmatrix = list(Lmatrix[[1]], Lmatrix1)
                }
            }
            meanRandFamily = rep("gaussian", meanNumberRandom)
        }
    }
    
    if (phiNumberRandom > 0) {
        if (phiRandFamily[1] == "GARCH") {
            phiCorr = phiRandFamily[1]
            phiRandFamily[1] = "gaussian"
        }
        if (phiRandFamily[1] == "AR(1)") {
            phiCorr = "AR"
            phiRandFamily[1] = "gaussian"
        }
    }
    
    # ¦§ MeanModel (MM) ####
    if (meanNumberRandom == 0 || lambdaFormula == "lambda ~ 1") 
        MM = DHGLMMODELING(
            Model = "mean",
            Link = meanLink,
            LinPred = as.formula(meanFormula),
            RandDist = meanRandFamily,
            Offset = meanOffsetVariable,
            LMatrix = Lmatrix,
            spatial = spatial,
            longitude = longitude,
            latitude = latitude,
            spline = spline,
            Neighbor = adjacent
        )
    else         
        MM = DHGLMMODELING(
            Model = "mean",
            Link = meanLink,
            LinPred = as.formula(meanFormula),
            RandDist = meanRandFamily,
            Offset = meanOffsetVariable,
            LMatrix = Lmatrix,
            LinkRandVariance = lambdaLink, #in the past : "log"
            LinPredRandVariance = as.formula(lambdaFormula), 
            RandDistRandVariance = lambdaRandFamily,    #in the past : RandDistRandVariance
            spatial = spatial,
            longitude = longitude,
            latitude = latitude,
            spline = spline,
            Neighbor = adjacent
        )
    
    # ¦§ DispersionModel (DM) ####
    
    if(phiFormula == "constant") 
        DM = DHGLMMODELING(
            Model = "dispersion",
            Link = "log"
        )
    else
        DM = DHGLMMODELING(
            Model = "dispersion",
            Link = phiLink,
            LinPred = as.formula(phiFormula),
            RandDist = phiRandFamily,
            Offset = phiOffsetVariable,
            corr = phiCorr
        )
    
    incProgress(0.4, detail = paste("Loading..."))
    
    # ¦§ Variable Declaration 2 ####
    # Variable Declaration for dhglmfit
    if (is.null(input$re_hglm_a_mord))
        mord = 1
    else
        mord = as.numeric(input$re_hglm_a_mord)
    
    if (is.null(input$re_hglm_a_dord))
        dord = 1
    else
        dord = as.numeric(input$re_hglm_a_dord)
    
    if (is.null(input$re_hglm_a_dispersion))
        Dmethod = "deviance"
    else
        Dmethod = input$re_hglm_a_dispersion
    
    if (is.null(input$re_hglm_a_reml))
        reml = TRUE
    else 
        reml = FALSE
    
    BinomialDen = NULL
    if (!is.null(input$re_hglm_m_check_binomd) && input$re_hglm_m_check_binomd) 
        BinomialDen = data2[[input$re_hglm_m_binomd]]
    
    # ¦§ Fitted Model ####
    
    fittedModel <<-
        dhglmfit(
            RespDist = isolate(input$re_hglm_m_dist),
            BinomialDen = BinomialDen,
            DataMain = data2,
            MeanModel = MM,
            DispersionModel = DM,
            Dmethod = Dmethod,
            mord = mord,
            dord = dord,
            REML = reml
        )
    
    incProgress(0.3, detail = paste("Loading..."))

    
    fittedModel$option = list(
        RespDist = isolate(input$re_hglm_m_dist),
        BinomialDen = BinomialDen,
        DataMain = data2,
        MeanModel = MM,
        DispersionModel = DM,
        Dmethod = Dmethod,
        mord = mord,
        dord = dord,
        REML = reml
    ) 
    
    # ¦§ Description 1 ####
    
    nMrand <- length(strsplit(input$re_hglm_m_model, "|", fixed = TRUE)[[1]]) - 1
    mRandstring = NULL
    pRandstring = NA
    lRandstring = NA
    if (nMrand >= 1)
        for (i in 1:nMrand)
            mRandstring = paste0(c(mRandstring, input[[paste0("re_hglm_m_randfamily_", i)]]), collapse = ", ")
    else
        mRandstring = NA
    
    modelDesc = matrix(c(
        input$re_hglm_m_model,
        "constant",
        "lambda ~ 1",
        input$re_hglm_m_link,
        phiLink,
        "log",
        input$re_hglm_m_dist,
        "gaussian",
        "gaussian",
        mRandstring,
        NA,
        NA
    ),nrow=3, ncol=4)
    if (!is.null(input$re_hglm_p_check_phi) && input$re_hglm_p_check_phi) {
        modelDesc[2, 1] <- input$re_hglm_p_model
        modelDesc[2, 2] <- input$re_hglm_p_link
    }
    if (!is.null(input$re_hglm_l_check_lambda) && input$re_hglm_l_check_lambda) {
        modelDesc[3, 1] <- input$re_hglm_l_model
        modelDesc[3, 2] <- input$re_hglm_l_link
    }
    colnames(modelDesc) <- c("Model", "Link", "Dist", "Rand")
    rownames(modelDesc) <- c("Mean", "Phi", "Lambda")
    fittedModel$modelDesc <- modelDesc
    
    # ¦§ Likelihood ####
    
    if(!is.null(fittedModel$likeli_coeff) && nrow(fittedModel$likeli_coeff) ==5 ) {
        likeli_coeff <- matrix(fittedModel$likeli_coeff, ncol = 5)
        colnames(likeli_coeff) <- c("-2ML", "-2RL", "cAIC", "Scaled Deviance", "df")
        fittedModel$likeli_coeff = likeli_coeff
    }
    
    if(!is.null(spatial) && (spatial == "IAR" || spatial == "MRF")) {
        beta_h = fittedModel$beta_h
        se_beta = fittedModel$se_beta
        z_beta<-beta_h/se_beta
        ##        pval <- 2 * pnorm(abs(z_beta), lower.tail = FALSE)
        beta_coeff<-cbind(matrix(beta_h),matrix(se_beta),matrix(z_beta))
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- NULL
        fittedModel$beta_coeff = beta_coeff
        
        likeli_coeff <- matrix(
            c(
                fittedModel$like_value[-1]
            ),
            ncol = 4
        )
        colnames(likeli_coeff)<-c("-2ML", "-2RL", "cAIC", "df")
        fittedModel$likeli_coeff = likeli_coeff
    }
    
    if(!is.null(spatial) && spatial == "Matern") {
        likeli_coeff <- matrix(fittedModel$likeli_coeff, ncol = 3)
        colnames(likeli_coeff)<-c("-2ML", "-2RL", "cAIC")
        fittedModel$likeli_coeff = likeli_coeff
    }
    
    if(!is.null(phiCorr) && phiCorr=="AR") {
        likeli_coeff <- matrix(
            c(
                fittedModel$ml, 
                fittedModel$rl, 
                fittedModel$caic,
                fittedModel$scaled_dv,
                fittedModel$df
            ),
            ncol = 5
        )
        colnames(likeli_coeff) <- c("-2ML", "-2RL", "cAIC", "Scaled Deviance", "df")
        fittedModel$likeli_coeff = likeli_coeff
    }
    
    if(!is.null(phiCorr) && phiCorr=="GARCH") {
        likeli_coeff <- matrix(
            c(
                fittedModel$ml, 
                fittedModel$rl, 
                fittedModel$caic,
                fittedModel$df
            ),
            ncol = 4
        )
        colnames(likeli_coeff) <- c("-2ML", "-2RL", "cAIC", "df")
        fittedModel$likeli_coeff = likeli_coeff
    }
    
    # ¦§ Description 2 (Exp) ####
    
    if (length(colnames(fittedModel$beta_coeff)) == 3) {
        pValue <- 2 * pnorm(abs(fittedModel$beta_coeff[, 3]), lower.tail = FALSE)
        llValue <- fittedModel$beta_coeff[, 1] - 1.96 * fittedModel$beta_coeff[, 2]
        ulValue <- fittedModel$beta_coeff[, 1] + 1.96 * fittedModel$beta_coeff[, 2]
        
        fittedModel$beta_coeff <- cbind(fittedModel$beta_coeff, pValue, llValue, ulValue)
        colnames(fittedModel$beta_coeff) <- c("Estimate", "Std. Error", "t-value", "p_val", "LL", "UL")
    }
    
    if (!is.null(fittedModel$beta_coeff) && !is.null(input$re_hglm_m_check_exp) && input$re_hglm_m_check_exp) {
        fittedModel$beta_coeff[, 1] <- exp(fittedModel$beta_coeff[, 1])
        colnames(fittedModel$beta_coeff)[1] <- "exp(Estimate)"
        
        checkColnames <- colnames(fittedModel$beta_coeff)
        if (checkColnames[5] == "LL" && checkColnames[6] == "UL") {
            fittedModel$beta_coeff[, 5] <- exp(fittedModel$beta_coeff[, 5])
            colnames(fittedModel$beta_coeff)[5] <- "exp(LL)"
            fittedModel$beta_coeff[, 6] <- exp(fittedModel$beta_coeff[, 6])
            colnames(fittedModel$beta_coeff)[6] <- "exp(UL)"
        }

        fittedModel$beta_coeff <- subset(fittedModel$beta_coeff, select=-2)
    }
    
    if (!is.null(fittedModel$phi_coeff)) {
        fittedModel$phi_coeff <- cbind(fittedModel$phi_coeff, "exp(Estimate)" = exp(fittedModel$phi_coeff[, 1]))
        fittedModel$phi_coeff <- subset(fittedModel$phi_coeff, select=c(1,4,2,3))
    }

    if (!is.null(fittedModel$lambda_coeff)) {
        fittedModel$lambda_coeff <- cbind(fittedModel$lambda_coeff, "exp(Estimate)" = exp(fittedModel$lambda_coeff[, 1]))
        fittedModel$lambda_coeff <- subset(fittedModel$lambda_coeff, select=c(1,4,2,3))
    }
    
    # ¦§ Comparison Model ####
    
    hglmComparisonmodel1 = NULL
    hglmComparisonmodel2 = NULL
    hglmComparisonmodel3 = NULL
    
    if (!is.null(input$re_hglm_m_check_comparison) && input$re_hglm_m_check_comparison) {
        spatialVariable = NA
        if (!is.null(spatial))
            spatialVariable = spatial
        
        hglmComparisonmodel1 <- matrix(
            c(
                input$re_hglm_m_model,
                input$re_hglm_m_link,
                input$re_hglm_m_dist,
                mRandstring,
                spatialVariable
            ),
            ncol = 5
        )
        
        hglmComparisonmodel2 <- matrix(
            c(
                phiFormula,
                lambdaFormula
            ),
            ncol = 2
        )
        
        hglmComparisonmodel3 <- matrix(
            c(
                ifelse(is.null(fittedModel$likeli_coeff), fittedModel$ml, fittedModel$likeli_coeff[1,"-2ML"]),
                ifelse(is.null(fittedModel$likeli_coeff), fittedModel$rl, fittedModel$likeli_coeff[1,"-2RL"]),
                ifelse(is.null(fittedModel$likeli_coeff), fittedModel$caic, fittedModel$likeli_coeff[1,"cAIC"]),
                ifelse(is.null(fittedModel$likeli_coeff) || spatial == "Matern", fittedModel$df, fittedModel$likeli_coeff[1,"df"]),
                ifelse(is.null(fittedModel$n.mean), NA, round(fittedModel$n.mean, digits=0)),
                ifelse(is.null(fittedModel$n.disp), NA, round(fittedModel$n.disp, digits=0))
            ),
            ncol = 6
        )
        hglmComparisonmodel1 <- rbind(g_re_hglm_r_comparisonmodel_1, hglmComparisonmodel1)
        hglmComparisonmodel2 <- rbind(g_re_hglm_r_comparisonmodel_2, hglmComparisonmodel2)
        hglmComparisonmodel3 <- rbind(g_re_hglm_r_comparisonmodel_3, hglmComparisonmodel3)
        
        colnames(hglmComparisonmodel1) <- c("Model", "Link", "Dist", "Rand", "Spatial")
        colnames(hglmComparisonmodel2) <- c("Phi Model", "Lambda Model")
        colnames(hglmComparisonmodel3) <- c("-2ML", "-2RL", "cAIC", "df", "n.mean", "n.disp")
        
        g_re_hglm_r_comparisonmodel_1 <<- hglmComparisonmodel1
        g_re_hglm_r_comparisonmodel_2 <<- hglmComparisonmodel2
        g_re_hglm_r_comparisonmodel_3 <<- hglmComparisonmodel3
        
        fittedModel$option$Comparisonmodel1 = hglmComparisonmodel1
        fittedModel$option$Comparisonmodel2 = hglmComparisonmodel2
        fittedModel$option$Comparisonmodel3 = hglmComparisonmodel3
        
        appendTab(inputId = "re_hglm_resulttabset",
            tabPanel(
                "Comparison Model",
                br(),
                tableOutput("re_hglm_r_comparisonmodel1"),
                tableOutput("re_hglm_r_comparisonmodel2"),
                tableOutput("re_hglm_r_comparisonmodel3")
            )
        )
    }
    
    # ¦¦ R Codes ####
    
    if (!is.null(input$re_hglm_m_check_rcodes) && input$re_hglm_m_check_rcodes) {
        Rcode1 <- NULL
        Rcode2 <- NULL
        Rcode3 <- NULL
        offsetRcode = NULL
        if (!is.null(meanOffsetVariable))
            offsetRcode <- paste0(", Offset = data[[\"", input$re_hglm_m_offset, "\"]]")
        mRandRcode = NULL
        if (nMrand >= 1)
            for (i in 1:nMrand)
                mRandRcode = paste0(c(mRandRcode, input[[paste0("re_hglm_m_randfamily_", i)]]), collapse = "\", \"")
        
        if (lambdaFormula == "lambda ~ 1") {
            Rcode1 <- paste0(
                "MM <- DHGLMMODELING(Model = \"mean\", Link = \"",
                input$re_hglm_m_link,
                "\", LinPred = ",
                meanFormula,
                ", RandDist = c(\"",
                mRandRcode,
                "\")",
                offsetRcode,
                ")"
            )
        } else {
            Rcode1 <- paste0(
                "MM <- DHGLMMODELING(Model = \"mean\", Link = \"",
                input$re_hglm_m_link,
                "\", LinPred = ",
                meanFormula,
                ", RandDist = c(\"",
                mRandRcode,
                "\"), LinPredRandVariance = ",
                lambdaFormula,
                ", LinkRandVariance = \"log\"",
                offsetRcode,
                ")"
            )
        }
        if(phiFormula == "constant") 
            Rcode2 <- "DM <- DHGLMMODELING(Model = \"dispersion\", Link = \"log\")"
        else
            Rcode2 <- paste0(
                "DM <- DHGLMMODELING(Model = \"dispersion\", Link = \"", phiLink, "\", LinPred = ",
                phiFormula,
                ")"
            )
        
        Rcode3 <- paste0(
            "dhglmfit(RespDist = \"", 
            input$re_hglm_m_dist,
            "\", DataMain = data, MeanModel = MM, DispersionModel = DM)"
        )
        
        Rcodes <- matrix(c(Rcode1, Rcode2, Rcode3), ncol = 1)
        colnames(Rcodes) <- "Call"
        rownames(Rcodes) <- c("MM", "DM", "dhglmfit")
        fittedModel$option$Rcodes = Rcodes
    }
            
    setProgress(1, detail = "Finish")
    return(fittedModel)
    })
})

# HGLM UI (Accordion) ####

output$re_hglm_r_accordion <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    hglmAccordion <- bs_accordion_sidebar(
        id = "re_hglm_r_accordion",
        spec_side = c(width = 3, offset = 0),
        spec_main = c(width = 9, offset = 0)    
    )
    
    # ¦§ Mean ####
    
    hglmAccordion <- hglmAccordion %>%
        bs_append(
            title_side = "Mean",
            content_side = NULL,
            content_main = div(
                fluidRow(
                    column(
                        9, 
                        h3(strong("Model for Mean"))
                    ),
                    column(
                        3,
                        style = "text-align:right; padding:15px"
                    )
                ),
                uiOutput("re_hglm_m_resp"), 
                uiOutput("re_hglm_m_variable"), 
                fluidRow(
                    column(
                        8,
                        uiOutput("re_hglm_m_interaction") 
                    ),
                    column(
                        4,
                        uiOutput("re_hglm_m_interactionappend"),
                        style = "text-align:right; padding:15px"
                    )
                ),
                uiOutput("re_hglm_m_rand"),
                uiOutput("re_hglm_m_check_slope"),
                uiOutput("re_hglm_m_check_withoutslope"),
                uiOutput("re_hglm_m_slope"),
                uiOutput("re_hglm_m_check_randinteraction"),
                fluidRow(
                    column(
                        8,
                        uiOutput("re_hglm_m_randinteraction")
                    ),
                    column(
                        4,
                        uiOutput("re_hglm_m_randinteractionappend"),
                        style = "text-align:right; padding:15px"
                    )
                ), #new components
                
                uiOutput("re_hglm_m_dist"),
                uiOutput("re_hglm_m_check_binomd"),
                uiOutput("re_hglm_m_binomd"),
                uiOutput("re_hglm_m_link"),
                
                hr(style = "border-color: #2C3E50;"),
                uiOutput("re_hglm_m_randfamily"),
                uiOutput("re_hglm_m_upload_neighbor"),
                
                uiOutput("re_hglm_m_check_nointercept"),
                uiOutput("re_hglm_m_check_offset"),
                uiOutput("re_hglm_m_offset"),
                uiOutput("re_hglm_m_check_exp"),

                hr(style = "border-color: #2C3E50;"),
                uiOutput("re_hglm_m_check_comparison"),
                uiOutput("re_hglm_m_check_rcodes")
            )
        )
    
    # ¦§ Phi ####
    
    hglmAccordion <- hglmAccordion %>%
    bs_append(
        title_side = "Phi",
        content_side = uiOutput("re_hglm_p_check_phi"),
        content_main = div(
            h3(strong("Model for Phi")),
            uiOutput("re_hglm_p_resp"),
            uiOutput("re_hglm_p_variable"),
            fluidRow(
                column(
                    8,
                    uiOutput("re_hglm_p_interaction")
                ),
                column(
                    4,
                    uiOutput("re_hglm_p_interactionappend"),
                    style = "text-align:right; padding:15px"
                )
            ),

            uiOutput("re_hglm_p_link"),
            uiOutput("re_hglm_p_randfamily"),
            
            uiOutput("re_hglm_p_check_nointercept"),
            uiOutput("re_hglm_p_check_offset"),
            uiOutput("re_hglm_p_offset")
        )
    )
    
    # ¦§ Lambda ####
    
    hglmAccordion <- hglmAccordion %>%
    bs_append(
        title_side = "Lambda",
        content_side = uiOutput("re_hglm_l_check_lambda"),
        content_main = div(
            h3(strong("Model for Lambda")),
            uiOutput("re_hglm_l_resp"),
            uiOutput("re_hglm_l_variable"),
            fluidRow(
                column(
                    8,
                    uiOutput("re_hglm_l_interaction")
                ),
                column(
                    4,
                    uiOutput("re_hglm_l_interactionappend"),
                    style = "text-align:right; padding:15px"
                )
            ),

            uiOutput("re_hglm_l_link"),
            uiOutput("re_hglm_l_randfamily"),
            
            uiOutput("re_hglm_l_check_nointercept")
            # uiOutput("re_hglm_l_check_offset"),
            # uiOutput("re_hglm_l_offset")
            # DHGLM doesn't support offset variable in lambda
        )
    )
    
    # ¦¦ Additional Settings ####
    
    hglmAccordion <- hglmAccordion %>%
    bs_append(
        title_side = "Setting",
        content_side = NULL,
        content_main = div(
            h3(strong("Additional Settings")),
            uiOutput("re_hglm_a_dispersion"),
            uiOutput("re_hglm_a_reml"),
            # uiOutput("re_hglm_a_check_orthogonal"),
            uiOutput("re_hglm_a_mord"),
            uiOutput("re_hglm_a_dord"),
            h5(helpText("Order of Laplace Approximation for Likelihood(mean) and Restricted Likelihood(Dispersion)"))
        )
    )
    
    div(
        hglmAccordion,
        use_bs_tooltip(),
        use_bs_accordion_sidebar() # needs to be at end, for some reason
    )
})

# HGLM Components####
# ¦§ Mean ####

output$re_hglm_m_model<-renderUI({ #algorithm not good
    input$re_hglm_m_resp
    input$re_hglm_m_variable$right
    input$re_hglm_m_rand
    input$re_hglm_m_check_nointercept
    input$re_hglm_m_binomd
    input$re_hglm_m_slope
    
    if(any(is.null(input$re_hglm_m_resp), is.null(input$re_hglm_m_variable)))
        return()
    
    responsePart = input$re_hglm_m_resp
    variablePart = paste0(input$re_hglm_m_variable$right, collapse=" + ")
    randomPart = NULL
    meanModel = NULL
    
    # in DHGLM, even if we use binomial denominator, model equation doesn't change. 
    # if (all(!is.null(input$re_hglm_m_check_binomd), input$re_hglm_m_check_binomd, 
    #         !is.null(input$re_hglm_m_binomd), input$re_hglm_m_binomd != ""))
    #     responsePart = paste0("cbind(", responsePart, ",", input$re_hglm_m_binomd,"-", responsePart, ")")
    
    if(!is.null(input$re_hglm_m_rand) && length(input$re_hglm_m_rand) > 0) {
        if (is.null(input$re_hglm_m_check_withoutslope) || !input$re_hglm_m_check_withoutslope) {
            randomPart = paste0(input$re_hglm_m_rand, collapse = ") + (1|")
            randomPart = paste0(" + (1|", randomPart, ")")
        }
        
        #random slope (a|b) form
        if(all(!is.null(input$re_hglm_m_check_slope), input$re_hglm_m_check_slope, length(input$re_hglm_m_slope)>0 )) 
            for(i in 1:length(input$re_hglm_m_rand)) 
                for(j in 1:length(input$re_hglm_m_slope)) 
                    randomPart = paste0(randomPart," + (", input$re_hglm_m_slope[j], "|", input$re_hglm_m_rand[i], ")")
    }
    
    if(!is.null(input$re_hglm_m_check_nointercept) && !input$re_hglm_m_check_nointercept) {
        if(variablePart == "")
            meanModel = paste(responsePart, "~", "1", randomPart)
        else
            meanModel = paste(responsePart, "~", variablePart, randomPart)
    }
    
    if(!is.null(input$re_hglm_m_check_nointercept) && input$re_hglm_m_check_nointercept) {
        if(variablePart == "")
            meanModel = paste(responsePart, "~", "-1", randomPart)
        else
            meanModel = paste(responsePart, "~", "-1 +", variablePart, randomPart)
    }
    
    textAreaInput("re_hglm_m_model", "Model for Mean", value = meanModel, height = "60px")
})

output$re_hglm_m_resp <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue = names(select_if(data2, is.numeric))
    
    div(
        selectInput(
            "re_hglm_m_resp",
            "Response Variable",
            choices = as.list(c("", nameValue)), 
            selected = NULL, 
            multiple = FALSE
        ),
        bsTooltip(
            "re_hglm_m_resp", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$re_hglm_m_variable <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$re_hglm_m_interactionappend
    
    nameValue <- c(names(data2), g_re_hglm_m_interaction)
    
    chooserInput("re_hglm_m_variable", "Variable", "Selected", nameValue, c(), size = 15, multiple = TRUE, width = 125)
})

output$re_hglm_m_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue <- names(data2)
    
    selectInput(
        "re_hglm_m_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$re_hglm_m_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    actionButton("re_hglm_m_interactionappend", "Append")
})

observeEvent(input$re_hglm_m_interactionappend, {
    if(length(input$re_hglm_m_interaction) > 1)
        g_re_hglm_m_interaction <<- c(g_re_hglm_m_interaction, paste(input$re_hglm_m_interaction, collapse=":"))
})


output$re_hglm_m_rand <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$re_hglm_m_randinteractionappend
    
    nameValue <- c(names(data2), g_re_hglm_m_randinteraction)
    
    selectInput(
        "re_hglm_m_rand",
        "Random Effects",
        choices = as.list(nameValue), 
        multiple = TRUE
    )
})

output$re_hglm_m_check_slope <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_hglm_m_check_slope", "Random Slope Model", value = FALSE)
})

output$re_hglm_m_check_withoutslope <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if (is.null(input$re_hglm_m_check_slope) || !input$re_hglm_m_check_slope)
        return()
    
    checkboxInput("re_hglm_m_check_withoutslope", "Without Random Intercept", value = FALSE)
})

output$re_hglm_m_slope <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if (is.null(input$re_hglm_m_check_slope) || !input$re_hglm_m_check_slope)
        return()
    
    nameValue = names(data2)
    
    selectInput(
        "re_hglm_m_slope",
        "Random slope",
        choices = as.list(nameValue), 
        multiple = TRUE
    )
})

output$re_hglm_m_check_randinteraction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_hglm_m_check_randinteraction", "Interaction in the Random Effect", value = FALSE)
})


output$re_hglm_m_randinteraction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$re_hglm_m_check_randinteraction) || !input$re_hglm_m_check_randinteraction)
        return()
    
    nameValue <- names(data2)
    
    selectInput(
        "re_hglm_m_randinteraction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$re_hglm_m_randinteractionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_hglm_m_check_randinteraction) || !input$re_hglm_m_check_randinteraction)
        return()
    
    actionButton("re_hglm_m_randinteractionappend", "Append")
})

observeEvent(input$re_hglm_m_randinteractionappend, {
    
    if(length(input$re_hglm_m_randinteraction) > 1)
        g_re_hglm_m_randinteraction <<- c(g_re_hglm_m_randinteraction, paste(input$re_hglm_m_randinteraction, collapse=":"))
})

output$re_hglm_m_dist<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    dist = c(
        "gaussian" = "gaussian",
        "binomial" = "binomial",
        "poisson" = "poisson",
        "gamma" = "gamma"
    )
    
    selectInput(
        "re_hglm_m_dist",
        "Distribution", 
        choices = dist, 
        multiple = FALSE
    )
})

output$re_hglm_m_check_binomd<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(any(is.null(input$re_hglm_m_dist), input$re_hglm_m_dist != "binomial"))
        return()
    
    checkboxInput("re_hglm_m_check_binomd", "Binomial Denominator", value = FALSE)
})

output$re_hglm_m_binomd<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$re_hglm_m_check_binomd) || !(input$re_hglm_m_check_binomd))
        return()
    
    nameValue = names(select_if(data2, is.numeric))
    
    selectInput(
        "re_hglm_m_binomd",
        "Binomial Denominator",
        choices = as.list(c("", nameValue)),
        multiple = FALSE
    )
})

output$re_hglm_m_link <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_hglm_m_dist))
        return()
    
    selection = "identity"
    
    if (input$re_hglm_m_dist == "binomial")
        selection = "logit"
    else if (input$re_hglm_m_dist == "poisson")
        selection = "log"
    else if (input$re_hglm_m_dist == "gamma")
        selection = "log"
    
    link = c(
        "identity" = "identity",
        "log" = "log",
        "logit" = "logit",
        "probit" = "probit",
        "cloglog" = "cloglog",
        "inverse" = "inverse"
    )
    
    selectInput(
        "re_hglm_m_link",
        "Link Function",
        choices = link,
        selected = selection,
        multiple = FALSE
    )
})

observe({ #algorithm is not good
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_hglm_m_model) || is.null(input$re_hglm_m_rand) || length(input$re_hglm_m_rand) == 0) {
        output$re_hglm_m_randfamily <- renderUI({
            return()
        })
    }
    
    if(length(input$re_hglm_m_rand) > 0) {
        numberRandom <- length(strsplit(input$re_hglm_m_model, "|", fixed = TRUE)[[1]]) - 1
        for (i in 1:numberRandom) {
            output$re_hglm_m_randfamily <- renderUI({
                re_hglm_m_randfamily_list <- lapply(1:numberRandom, function(i) {
                    selectInput(
                        paste0("re_hglm_m_randfamily_", i),
                        paste("Distribution for Random effects", i),
                        choices = c(
                            "gaussian",
                            "gamma",
                            "beta",
                            "inverse-gamma",
                            "random work",
                            "local linear trend",
                            "quarterly seasonal",
                            # "cubic spline", # it is dummy
                            "IAR",
                            "MRF",
                            "Matern",
                            "AR(1)" 
                        )
                    )
                })
                do.call(tagList, re_hglm_m_randfamily_list)
            })
        }
    }
    else {
        output$re_hglm_m_randfamily <- renderUI({
            return()
        })
    }
})

# Neighbor for IAR, MRF 
# It only reacts with randfamily 1. 
output$re_hglm_m_upload_neighbor <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_hglm_m_randfamily_1) || (input$re_hglm_m_randfamily_1 != "IAR" && input$re_hglm_m_randfamily_1 != "MRF"))
        return()
    
    fileInput("re_hglm_m_neighborfile", "Upload Neighborhood File", multiple = TRUE)
})

re_hglm_m_neighborfile <- reactive({
    re_hglm_m_neighborfile <- input$re_hglm_m_neighborfile
    if (is.null(re_hglm_m_neighborfile)) {
        return()
    }
    read.table(
        file = re_hglm_m_neighborfile$datapath,
        sep = input$sep,
        header = input$header,
        stringsAsFactors = TRUE
    )
})

output$re_hglm_m_check_nointercept<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_hglm_m_check_nointercept", "No Intercept Model", value = FALSE)
})

output$re_hglm_m_check_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_hglm_m_check_offset", "Offset Variable", value = FALSE)
})

output$re_hglm_m_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue = names(select_if(data2, is.numeric))
    
    if(is.null(input$re_hglm_m_check_offset) || !input$re_hglm_m_check_offset)
        return()
    
    div(
        selectInput(
            "re_hglm_m_offset",
            "Offset Variable",
            choices = as.list(c("", nameValue)),
            selected = NULL,
            multiple = FALSE
        ),
        bsTooltip(
            "re_hglm_m_offset", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$re_hglm_m_check_exp<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_hglm_m_check_exp", "Exponential Scale", value = FALSE)
})

output$re_hglm_m_check_vif<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_hglm_m_check_vif", "VIF", value = FALSE)
})

output$re_hglm_m_check_robustse<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_hglm_m_check_robustse", "Robust Standard Errors", value = FALSE)
})

output$re_hglm_m_check_confint<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_hglm_m_check_confint", "Confidence Intervals for Coefficients", value = FALSE)
})

output$re_hglm_m_check_exp <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_m_check_confint
    
    checkboxInput("re_hglm_m_check_exp", "Exponential scale", value = FALSE)
})

output$re_hglm_m_check_comparison<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_hglm_m_check_comparison", "model comparison", value = FALSE)
})

output$re_hglm_m_check_rcodes<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_hglm_m_check_rcodes", "R Codes", value = FALSE)
})

# ¦§ Phi ####

output$re_hglm_p_check_phi<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_hglm_p_check_phi", strong("Use"), value = FALSE)
})

output$re_hglm_p_model<-renderUI({ #algorithm not good
    input$re_hglm_p_resp
    input$re_hglm_p_variable$right
    input$re_hglm_p_check_nointercept
    
    if(is.null(input$re_hglm_p_check_phi) || !input$re_hglm_p_check_phi)
        return()
    
    if(any(is.null(input$re_hglm_p_resp), is.null(input$re_hglm_p_variable)))
        return()
    
    responsePart = input$re_hglm_p_resp
    variablePart = paste0(input$re_hglm_p_variable$right, collapse=" + ")
    phiModel = NULL
    
    if(!is.null(input$re_hglm_p_check_nointercept) && !input$re_hglm_p_check_nointercept) {
        if(variablePart == "")
            phiModel = paste(responsePart, "~", "1")
        else
            phiModel = paste(responsePart, "~", variablePart)
    }
    
    if(!is.null(input$re_hglm_p_check_nointercept) && input$re_hglm_p_check_nointercept) {
        if(variablePart == "")
            phiModel = paste(responsePart, "~", "-1")
        else
            phiModel = paste(responsePart, "~", "-1 +", variablePart)
    }
    
    textAreaInput("re_hglm_p_model", "Model for Phi", value = phiModel, height = "60px")
})

output$re_hglm_p_resp <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$re_hglm_p_check_phi) || !input$re_hglm_p_check_phi)
        return()
    
    nameValue = c("phi" = "phi")
    
    selectInput(
        "re_hglm_p_resp",
        "Residual Variance",
        choices = nameValue,
        selected = NULL, 
        multiple = FALSE
    )
})

output$re_hglm_p_variable <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$re_hglm_p_interactionappend
    
    if(is.null(input$re_hglm_p_check_phi) || !input$re_hglm_p_check_phi)
        return()
    
    nameValue <- c(names(data2), g_re_hglm_p_interaction)
    
    chooserInput("re_hglm_p_variable", "Variable", "Selected", nameValue, c(), size = 15, multiple = TRUE, width = 125)
})

output$re_hglm_p_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$re_hglm_p_check_phi) || !input$re_hglm_p_check_phi)
        return()
    
    nameValue <- names(data2)
    
    selectInput(
        "re_hglm_p_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$re_hglm_p_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_hglm_p_check_phi) || !input$re_hglm_p_check_phi)
        return()
    
    actionButton("re_hglm_p_interactionappend", "Append")
})

observeEvent(input$re_hglm_p_interactionappend, {
    if(length(input$re_hglm_p_interaction) > 1)
        g_re_hglm_p_interaction <<- c(g_re_hglm_p_interaction, paste(input$re_hglm_p_interaction, collapse=":"))
})

output$re_hglm_p_link <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_hglm_p_check_phi) || !input$re_hglm_p_check_phi)
        return()
    
    link = c(
        "log" = "log",
        "inverse" = "inverse",
        "identity" = "identity")
    
    selectInput(
        "re_hglm_p_link",
        "Link Function",
        choices = link,
        multiple = FALSE
    )
})

output$re_hglm_p_check_nointercept<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_p_check_phi
    
    if(is.null(input$re_hglm_p_check_phi) || !input$re_hglm_p_check_phi)
        return()
    
    checkboxInput("re_hglm_p_check_nointercept", "No Intercept Model", value = FALSE)
})

output$re_hglm_p_check_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_p_check_phi
    
    if(is.null(input$re_hglm_p_check_phi) || !input$re_hglm_p_check_phi)
        return()
    
    checkboxInput("re_hglm_p_check_offset", "Offset Variable", value = FALSE)
})

output$re_hglm_p_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$re_hglm_p_check_phi
    
    if(is.null(input$re_hglm_p_check_phi) || !input$re_hglm_p_check_phi)
        return()
    
    nameValue = names(select_if(data2, is.numeric))
    
    if(all(!is.null(input$re_hglm_p_check_offset), input$re_hglm_p_check_offset)) {
        selectInput(
            "re_hglm_p_offset",
            "Offset Variable (Numeric Only)",
            choices = as.list(c("", nameValue)),
            selected = NULL,
            multiple = FALSE
        )
    }
})

# ¦§ Lambda ####

output$re_hglm_l_check_lambda<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_hglm_l_check_lambda", strong("Use"), value = FALSE)
})

output$re_hglm_l_model<-renderUI({ #algorithm not good
    input$re_hglm_l_resp
    input$re_hglm_l_variable$right
    input$re_hglm_l_check_nointercept
    input$re_hglm_l_check_lambda
    
    if(is.null(input$re_hglm_l_check_lambda) || !input$re_hglm_l_check_lambda)
        return()
    
    if(any(is.null(input$re_hglm_l_resp), is.null(input$re_hglm_l_variable)))
        return()
    
    responsePart = input$re_hglm_l_resp
    variablePart = paste0(input$re_hglm_l_variable$right, collapse=" + ")
    lambdaModel = NULL
    
    
    if(!is.null(input$re_hglm_l_check_nointercept) && !input$re_hglm_l_check_nointercept) {
        if(variablePart == "")
            lambdaModel = paste(responsePart, "~", "1")
        else
            lambdaModel = paste(responsePart, "~", variablePart)
    }
    
    if(!is.null(input$re_hglm_l_check_nointercept) && input$re_hglm_l_check_nointercept) {
        if(variablePart == "")
            lambdaModel = paste(responsePart, "~", "-1")
        else
            lambdaModel = paste(responsePart, "~", "-1 +", variablePart)
    }
    
    textAreaInput("re_hglm_l_model", "Model for Lambda", value = lambdaModel, height = "60px")
})

output$re_hglm_l_resp <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$re_hglm_l_check_lambda
    
    if(is.null(input$re_hglm_l_check_lambda) || !input$re_hglm_l_check_lambda)
        return()
    
    nameValue = c("lambda" = "lambda")
    
    selectInput(
        "re_hglm_l_resp",
        "Variance of Random Effect",
        choices = nameValue,
        selected = NULL, 
        multiple = FALSE
    )
})

output$re_hglm_l_variable <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$re_hglm_l_interactionappend
    input$re_hglm_l_check_lambda
    
    if(is.null(input$re_hglm_l_check_lambda) || !input$re_hglm_l_check_lambda)
        return()
    
    nameValue <- c(names(data2), g_re_hglm_l_interaction)
    
    chooserInput("re_hglm_l_variable", "Variable", "Selected", nameValue, c(), size = 15, multiple = TRUE, width = 125)
})

output$re_hglm_l_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$re_hglm_l_check_lambda
    
    if(is.null(input$re_hglm_l_check_lambda) || !input$re_hglm_l_check_lambda)
        return()
    
    nameValue <- names(data2)
    
    selectInput(
        "re_hglm_l_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$re_hglm_l_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_l_check_lambda
    
    if(is.null(input$re_hglm_l_check_lambda) || !input$re_hglm_l_check_lambda)
        return()
    
    actionButton("re_hglm_l_interactionappend", "Append")
})

observeEvent(input$re_hglm_l_interactionappend, {
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(length(input$re_hglm_l_interaction) > 1)
        g_re_hglm_l_interaction <<- c(g_re_hglm_l_interaction, paste(input$re_hglm_l_interaction, collapse=":"))
})

output$re_hglm_l_link <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_l_check_lambda
    
    if(is.null(input$re_hglm_l_check_lambda) || !input$re_hglm_l_check_lambda)
        return()
    
    link = c("log" = "log")
    
    selectInput(
        "re_hglm_l_link",
        "Link Function",
        choices = link,
        multiple = FALSE
    )
})

output$re_hglm_l_check_nointercept<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_l_check_lambda
    
    if(is.null(input$re_hglm_l_check_lambda) || !input$re_hglm_l_check_lambda)
        return()
    
    checkboxInput("re_hglm_l_check_nointercept", "No Intercept Model", value = FALSE)
})

# ¦¦ Additional Settings ####

output$re_hglm_a_dispersion<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    dispersionmethods = c(
        "deviance" = "deviance", 
        "Pearson" = "Pearson"
    )
    
    selectInput(
        "re_hglm_a_dispersion",
        "Method of Fitting Dispersion Model", 
        choices = dispersionmethods, 
        multiple = FALSE
    )
})

output$re_hglm_a_reml <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    selectInput(
        "re_hglm_a_reml",
        "REML or ML",
        choice = c("REML" = TRUE, "ML" = FALSE),
        multiple = FALSE
    )
})

# output$re_hglm_a_check_orthogonal <- renderUI({
#     input$file
#     input$resetData
#     
#     selectInput(
#         "re_hglm_m_check_orthogonal",
#         "SE with orthogonality or not",
#         choice = c(TRUE, FALSE), 
#         selected = FALSE,
#         multiple = FALSE
#     )
# })

output$re_hglm_a_mord <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    selectInput(
        "re_hglm_a_mord",
        "Order for Mean Model",
        choice = c(1, 0),
        multiple = FALSE
    )
})

output$re_hglm_a_dord <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    selectInput(
        "re_hglm_a_dord",
        "Order for Dispersion Model",
        choice = c(1, 2),
        multiple = FALSE
    )
})

# HGLM Results ####

output$re_hglm_r_mainsummary <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_hglm_results()$modelDesc))
        return()
    
    re_hglm_results()$modelDesc
}, rownames = TRUE, bordered = TRUE, caption = "Model Description", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_hglm_m_coeff <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_hglm_results()$beta_coeff))
        return()
    
    re_hglm_results()$beta_coeff
}, rownames = TRUE, bordered = TRUE, caption = "Estimate from Mean Model", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_hglm_l_coeff <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_hglm_results()$lambda_coeff))
        return()
    
    re_hglm_results()$lambda_coeff
}, rownames = TRUE, bordered = TRUE, caption = "Estimate for log(Lambda)", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_hglm_p_coeff <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_hglm_results()$phi_coeff))
        return()
    
    re_hglm_results()$phi_coeff
}, rownames = TRUE, bordered = TRUE, caption = "Estimate from Dispersion Model", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_hglm_r_likelihood <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_hglm_results()$likeli_coeff))
        return()
    
    re_hglm_results()$likeli_coeff
}, rownames = FALSE, bordered = TRUE, caption = "Likelihood", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_hglm_r_comparisonmodel1<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_hglm_results()$option$Comparisonmodel1))
        return()
    
    re_hglm_results()$option$Comparisonmodel1
}, rownames = TRUE, bordered = TRUE, caption = "Mean", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_hglm_r_comparisonmodel2<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_hglm_results()$option$Comparisonmodel2))
        return()
    
    re_hglm_results()$option$Comparisonmodel2
}, rownames = TRUE, bordered = TRUE, caption = "Phi & Lambda", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_hglm_r_comparisonmodel3<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_hglm_results()$option$Comparisonmodel3))
        return()
    
    re_hglm_results()$option$Comparisonmodel3
}, rownames = TRUE, bordered = TRUE, caption = "Likelihood", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_hglm_r_rcodes <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_hglm_results()$option$Rcodes))
        return()
  
    re_hglm_results()$option$Rcodes
}, rownames = TRUE, bordered = TRUE, caption = "R Codes", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

# HGLM Plot ####
# ¦¦ Mean #####

output$re_hglm_m_selectplot1<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_hglm_results()))
        return()
    
    meanNumberRand = length(strsplit(re_hglm_results()$modelDesc[1,1], "|", fixed = TRUE)[[1]]) - 1
    
    choiceName = c("Residuals")
    choiceValue = c("mean")
    if (meanNumberRand == 1) {
        choiceName = c(choiceName, "Random Effect", "95% CI for Random Effect")
        choiceValue = c(choiceValue, "v", "ci")
    }
    else if(meanNumberRand >1) {
        for (i in 1:meanNumberRand) {
            choiceName = c(choiceName, paste0("Random Effect", i))
            choiceValue = c(choiceValue, paste0("v", i))
        }
        for (i in 1:meanNumberRand) {
            choiceName = c(choiceName, paste0("95% CI for Random Effect", i))
            choiceValue = c(choiceValue, paste0("ci", i))
        }
    }
    radioGroupButtons(
        inputId = "re_hglm_m_selectplot1",
        label = "Select Plot", 
        choiceNames = choiceName,
        choiceValues = choiceValue, 
        selected = "mean",
        direction = "vertical"
    )
})

output$re_hglm_m_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if(g_resetresult == FALSE || is.null(re_hglm_results()))
        return()
    
    local({
        output$re_hglm_m_downloadplot1 <- downloadPlot({
            if(is.null(input$re_hglm_m_selectplot1))
                plotType = "mean"
            else
                plotType = input$re_hglm_m_selectplot1
            res<-re_hglm_results()
            spatial = ""
            
            if(!is.null(res$model_mu[[13]]) && res$model_mu[[13]]=="Matern")
                spatial = res$model_mu[[13]]
            else if(!is.null(res$option$MeanModel[[13]]) && (res$option$MeanModel[[13]] == "IAR" || res$option$MeanModel[[13]] == "MRF"))
                spatial = res$option$MeanModel[[13]]
            
            ggplotdhglm(res, type = plotType, spatial)
        })
    })
    
    if(is.null(input$re_hglm_m_selectplot1))
        plotType = "mean"
    else
        plotType = input$re_hglm_m_selectplot1
    
    res<-re_hglm_results()
    spatial = ""
    
    if(!is.null(res$model_mu[[13]]) && res$model_mu[[13]]=="Matern")
        spatial = res$model_mu[[13]]
    else if(!is.null(res$option$MeanModel[[13]]) && (res$option$MeanModel[[13]] == "IAR" || res$option$MeanModel[[13]] == "MRF"))
        spatial = res$option$MeanModel[[13]]
    
    ggplotdhglm(res, type = plotType, spatial)
})

output$re_hglm_m_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run

    if(g_resetresult == FALSE || is.null(re_hglm_results()))
        return()
    
    div(
        h4("Model Checking Plots for Mean"),
        div(
            style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
            withSpinner(
                plotOutput("re_hglm_m_plot1", height="600px"),
                type = 1,
                color = "#2c3e50",
                size = 1.2
            ),
            div(
                style = "position: absolute; left: 0.5em; bottom: 0.5em;",
                dropdown(
                    uiOutput("re_hglm_m_selectplot1"),
                    size = "xs",
                    icon = icon("chart-bar", class = "opt"), 
                    up = TRUE
                )
            ),
            div(
                style = "position: absolute; left:4em; bottom: 0.5em;",
                dropdown(
                    downloadButton(outputId = "re_hglm_m_downloadplot1", label = "Download Plot"),
                    size = "xs",
                    icon = icon("download", class = "opt"),
                    up = TRUE
                )
            )
        ),
        br()
    )
})

# ¦§ Phi #####

output$re_hglm_p_selectplot1<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if(g_resetresult == FALSE || is.null(re_hglm_results()$resid_phi))
        return()

    choiceName = c("Residuals")
    choiceValue = c("phi")
    
    radioGroupButtons(
        inputId = "re_hglm_p_selectplot1",
        label = "Select Plot", 
        choiceNames = choiceName,
        choiceValues = choiceValue, 
        selected = "phi",
        direction = "vertical"
    )
})

output$re_hglm_p_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if(g_resetresult == FALSE || is.null(re_hglm_results()$resid_phi))
        return()
    
    local({
        output$re_hglm_p_downloadplot1 <- downloadPlot({
            if(is.null(input$re_hglm_p_selectplot1))
                plotType = "phi"
            else
                plotType = input$re_hglm_p_selectplot1
            
            res<-re_hglm_results()
            
            ggplotdhglm(res, type = plotType)
        })
    })
    
    if(is.null(input$re_hglm_p_selectplot1))
        plotType = "phi"
    else
        plotType = input$re_hglm_p_selectplot1
    
    res<-re_hglm_results()
    
    ggplotdhglm(res, type = plotType)
})

output$re_hglm_p_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run

    if(g_resetresult == FALSE || is.null(re_hglm_results()$resid_phi))
        return()
    
    div(
        h4("Model Checking Plots for Phi"),
        div(
            style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
            
            withSpinner(
                plotOutput("re_hglm_p_plot1", height="600px"),
                type = 1,
                color = "#2c3e50",
                size = 1.2
            ),
            div(
                style = "position: absolute; left: 0.5em; bottom: 0.5em;",
                dropdown(
                    uiOutput("re_hglm_p_selectplot1"),
                    size = "xs",
                    icon = icon("chart-bar", class = "opt"), 
                    up = TRUE
                )
            ),
            div(
                style = "position: absolute; left:4em; bottom: 0.5em;",
                dropdown(
                    downloadButton(outputId = "re_hglm_p_downloadplot1", label = "Download Plot"),
                    size = "xs",
                    icon = icon("download", class = "opt"),
                    up = TRUE
                )
            )
        ),
        br()
    )
})



# ¦¦ Lambda ####

output$re_hglm_l_selectplot1<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if(g_resetresult == FALSE || is.null(re_hglm_results()$resid_lambda))
        return()

    choiceName = c("Residuals")
    choiceValue = c("lambda")
    
    radioGroupButtons(
        inputId = "re_hglm_l_selectplot1",
        label = "Select Plot", 
        choiceNames = choiceName,
        choiceValues = choiceValue, 
        selected = "lambda",
        direction = "vertical"
    )
})

output$re_hglm_l_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if(g_resetresult == FALSE || is.null(re_hglm_results()$resid_lambda))
        return()
    
    local({
        output$re_hglm_l_downloadplot1 <- downloadPlot({
            if(is.null(input$re_hglm_l_selectplot1))
                plotType = "lambda"
            else
                plotType = input$re_hglm_l_selectplot1
            
            res<-re_hglm_results()
            
            ggplotdhglm(res, type = plotType)
        })
    })
    
    if(is.null(input$re_hglm_l_selectplot1))
        plotType = "lambda"
    else
        plotType = input$re_hglm_l_selectplot1
    
    res<-re_hglm_results()
    
    ggplotdhglm(res, type = plotType)
})

output$re_hglm_l_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run

    if(g_resetresult == FALSE || is.null(re_hglm_results()$resid_lambda))
        return()
    
    div(
        h4("Model Checking Plots for Lambda"),
        div(
            style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
            
            withSpinner(
                plotOutput("re_hglm_l_plot1", height="600px"),
                type = 1,
                color = "#2c3e50",
                size = 1.2
            ),
            div(
                style = "position: absolute; left: 0.5em; bottom: 0.5em;",
                dropdown(
                    uiOutput("re_hglm_l_selectplot1"),
                    size = "xs",
                    icon = icon("chart-bar", class = "opt"), 
                    up = TRUE
                )
            ),
            div(
                style = "position: absolute; left:4em; bottom: 0.5em;",
                dropdown(
                    downloadButton(outputId = "re_hglm_l_downloadplot1", label = "Download Plot"),
                    size = "xs",
                    icon = icon("download", class = "opt"),
                    up = TRUE
                )
            )
        ),
        br()
    )
})

# HGLM Prediction ####

output$re_hglm_r_check_95mu <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    checkboxInput("re_hglm_r_check_95mu", "95% Confidence Interval for mu", value = FALSE)
})

output$re_hglm_r_box_95mu <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if(is.null(input$re_hglm_r_check_95mu) || !input$re_hglm_r_check_95mu)
        return()
    
    muLength = nrow(re_hglm_results()$beta_coeff)
    if(rownames(re_hglm_results()$beta_coeff)[1] == "(Intercept)")
        muLength = muLength - 1 
    
    textAreaInput("re_hglm_r_box_95mu", "Type mean Covariates", value = paste(rep(1, muLength), collapse = ","))
})

output$re_hglm_r_check_95phi <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if(is.null(re_hglm_results()$modelDesc) || re_hglm_results()$modelDesc[2,1] == "phi ~ 1"
       || re_hglm_results()$modelDesc[2,1] == "phi~1")
        return()
    
    checkboxInput("re_hglm_r_check_95phi", "95% Confidence Interval for phi", value = FALSE)
})

output$re_hglm_r_box_95phi <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if(is.null(input$re_hglm_r_check_95phi) || !input$re_hglm_r_check_95phi)
        return()
    
    phiLength = nrow(re_hglm_results()$phi_coeff)
    if(rownames(re_hglm_results()$phi_coeff)[1] == "(Intercept)")
        phiLength = phiLength - 1 
    
    textAreaInput("re_hglm_r_box_95phi", "Type phi Covariates", value = paste(rep(1, phiLength), collapse = ","))
})

output$re_hglm_r_prediction <- renderPrint({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if(!is.null(input$re_hglm_r_check_95mu) && input$re_hglm_r_check_95mu && !is.null(re_hglm_results()$beta_coeff)) {
        
        muLength = nrow(re_hglm_results()$beta_coeff)
        if(rownames(re_hglm_results()$beta_coeff)[1] == "(Intercept)")
            muLength = muLength - 1 
        
        muCovariates <- as.numeric(strsplit(input$re_hglm_r_box_95mu, ",", fixed = TRUE)[[1]])
        if(muLength >= length(muCovariates))
            muCovariates <- c(muCovariates, rep(0,muLength - length(muCovariates)))
        else muCovariates <- muCovariates[1:muLength]
        
        if(rownames(re_hglm_results()$beta_coeff)[1] == "(Intercept)")
            muCovariates<-c(1,muCovariates)
        
        muCovMatrix <- matrix(data = muCovariates, ncol=1)
        rownames(muCovMatrix) <- rownames(re_hglm_results()$beta_coeff)
        print("Covariates for mean Predictor")
        print(t(muCovMatrix))
        muPrediction <- sum(as.vector(muCovMatrix) * re_hglm_results()$beta_coeff[,1])
        muSE <- sqrt(t(muCovMatrix) %*% re_hglm_results()$cov_mu %*% muCovMatrix)
        LL <- muPrediction - 1.96 * muSE
        UL <- muPrediction + 1.96 * muSE
        print(paste0("Lower Limit for Confidence Intervel of mean predictor : ", LL))
        print(paste0("Upper Limit for Confidence Intervel of mean predictor : ", UL))
    }
    
    if(!is.null(input$re_hglm_r_check_95phi) && input$re_hglm_r_check_95phi 
       && !is.null(re_hglm_results()$phi_coeff)) {
        phiLength = nrow(re_hglm_results()$phi_coeff)
        if(rownames(re_hglm_results()$phi_coeff)[1] == "(Intercept)")
            phiLength = phiLength - 1
        
        phiCovariates <- as.numeric(strsplit(input$re_hglm_r_box_95phi, ",", fixed = TRUE)[[1]])
        if(phiLength >= length(phiCovariates))
            phiCovariates <- c(phiCovariates, rep(0,phiLength - length(phiCovariates)))
        else phiCovariates <- phiCovariates[1:phiLength]
        
        if(rownames(re_hglm_results()$phi_coeff)[1] == "(Intercept)")
            phiCovariates<-c(1,phiCovariates)
        
        phiCovMatrix <- matrix(data = phiCovariates, ncol=1)
        rownames(phiCovMatrix) <- rownames(re_hglm_results()$phi_coeff)
        print("Covariates for phi Predictor")
        print(t(phiCovMatrix))
        phiPrediction <- sum(as.vector(phiCovMatrix) * re_hglm_results()$phi_coeff[,1])
        phiSE <- sqrt(t(phiCovMatrix) %*% re_hglm_results()$cov_phi %*% phiCovMatrix)
        LL <- phiPrediction - 1.96 * phiSE
        UL <- phiPrediction + 1.96 * phiSE
        print(paste0("Lower Limit for Confidence Intervel of phi predictor : ", LL))
        print(paste0("Upper Limit for Confidence Intervel of phi predictor : ", UL))
    }
})


output$re_hglm_r_check_calculator <-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    checkboxInput("re_hglm_r_check_calculator", "Calculator for data", value = FALSE)
})

output$re_hglm_r_box_calculator <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if(is.null(input$re_hglm_r_check_calculator) || !input$re_hglm_r_check_calculator)
        return()
    
    textAreaInput("re_hglm_r_box_calculator", "Type Arithmetic Form (e.g. sum(log(y))", value = paste0("sum(", input$re_hglm_m_resp, ")"))
})

output$re_hglm_r_calculator <- renderPrint({
    input$file
    input$resetData
    input$di_option_run
    input$re_hglm_g_run
    
    if(is.null(input$re_hglm_r_check_calculator) || !input$re_hglm_r_check_calculator 
       || is.null(input$re_hglm_r_box_calculator))
        return()
    
    calcEquation <- input$re_hglm_r_box_calculator
    
    quo_var <- rlang::parse_expr(quo_name(enquo(calcEquation)))
    calcResults <- dplyr::mutate(data2, new1 = !!quo_var)
    
    print("Calculator Results")
    print(calcResults[1,ncol(calcResults)])
})
