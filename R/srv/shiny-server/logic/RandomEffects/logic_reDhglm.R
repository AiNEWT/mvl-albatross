# DHGLM Run ####

observeEvent(input$re_dhglm_g_run, {
    g_resetresult <<- TRUE
})

re_dhglm_results <- eventReactive(input$re_dhglm_g_run, {
    # ¦§ Warning & Notify ####
    # if user doesn't input some essential part, we wil return null value and show notification
    if (input$re_dhglm_m_resp == "")  {
        showNotification("Please choose response", type="warning")
        return()
    }
    
    if(input$re_dhglm_m_check_offset && input$re_dhglm_m_offset == "") {
        showNotification("Please choose offset variable", type="warning")
        return()
    }
    
    if(!is.null(input$re_dhglm_m_check_binomd) && input$re_dhglm_m_check_binomd && input$re_dhglm_m_binomd == "") {
        showNotification("Please choose binomial denominator", type="warning")
        return()
    }

    # Start Progressbar
    withProgress(message = 'DHGLM', style = "notification", value = 0.1, {    
    removeTab(inputId = "re_dhglm_resulttabset", target = "Model Comparison")
    
    # ¦§ Variable Declaration 1 ####
    # variable declaration for which is used widely in re_dhglm_results
    dataLength = nrow(data2)
    
    # variable declaration for which is used in making MM(MeanModel), DM(DispersionModel)
    meanLink = input$re_dhglm_m_link
    phiLink = "log"
    lambdaLink = NULL
    # In the past, they don't use input for phiLink and lambdaLink. why? 
    if(!is.null(input$re_dhglm_p_check_phi) && input$re_dhglm_p_check_phi)
        phiLink = input$re_dhglm_p_link
    if(!is.null(input$re_dhglm_l_check_lambda) && input$re_dhglm_l_check_lambda)
        lambdaLink = input$re_dhglm_l_link
    
    meanFormula = input$re_dhglm_m_model
    phiFormula = "constant"               
    lambdaFormula = "lambda ~ 1"  
    
    if(!is.null(input$re_dhglm_p_check_phi) && input$re_dhglm_p_check_phi)
        phiFormula = input$re_dhglm_p_model
    if(!is.null(input$re_dhglm_l_check_lambda) && input$re_dhglm_l_check_lambda)
        lambdaFormula = input$re_dhglm_l_model
    
    meanNumberRandom = length(strsplit(meanFormula, "|", fixed = TRUE)[[1]]) - 1
    phiNumberRandom = length(strsplit(phiFormula, "|", fixed = TRUE)[[1]]) - 1
    lambdaNumberRandom = length(strsplit(lambdaFormula, "|", fixed = TRUE)[[1]]) - 1
    
    meanRandFamily = NULL
    phiRandFamily = NULL
    lambdaRandFamily = NULL #in the past : NULL
    
    if(meanNumberRandom > 0)
        for (i in 1:meanNumberRandom) 
            meanRandFamily = c(meanRandFamily, input[[paste0("re_dhglm_m_randfamily_", i)]])
    if(phiNumberRandom > 0)
        for (i in 1:phiNumberRandom) 
            phiRandFamily = c(phiRandFamily, input[[paste0("re_dhglm_p_randfamily_", i)]])
    if(lambdaNumberRandom > 0)
        for (i in 1:lambdaNumberRandom) 
            lambdaRandFamily = c(lambdaRandFamily, input[[paste0("re_dhglm_l_randfamily_", i)]])
    else lambdaRandFamily = "gaussian"
    # In fact, dhglmfit doesn't support multi random effects in phi or lambda except "salamander" data. 
    # So this code should be revised in the future. 
    
    meanOffsetVariable = NULL
    phiOffsetVariable = NULL
    # dhglmfit doesn't provide offset variable for lambda. 
    
    if(input$re_dhglm_m_check_offset)
        meanOffsetVariable = data2[[input$re_dhglm_m_offset]]
    if(!is.null(input$re_dhglm_p_check_offset) && input$re_dhglm_p_check_offset)
        phiOffsetVariable = data2[[input$re_dhglm_p_offset]]
    
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

            neighborfile = re_dhglm_m_neighborfile()
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
            
            responsePart = input$re_dhglm_p_resp
            variablePart = paste0(input$re_dhglm_p_variable$right, collapse=" + ")

            
            
            if(input$re_dhglm_p_check_nointercept) {
                if(variablePart == "")
                    phiFormula = paste(responsePart, "~", "-1")
                else
                    phiFormula = paste(responsePart, "~", "-1 +", variablePart)
            }
            else {
                if(variablePart == "")
                    phiFormula = paste(responsePart, "~", "1")
                else
                    phiFormula = paste(responsePart, "~", variablePart)
            }
            
            
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
    if (is.null(input$re_dhglm_a_mord))
        mord = 1
    else
        mord = as.numeric(input$re_dhglm_a_mord)
    
    if (is.null(input$re_dhglm_a_dord))
        dord = 1
    else
        dord = as.numeric(input$re_dhglm_a_dord)
    
    if (is.null(input$re_dhglm_a_dispersion))
        Dmethod = "deviance"
    else
        Dmethod = input$re_dhglm_a_dispersion
    
    if (is.null(input$re_dhglm_a_reml))
        reml = TRUE
    else 
        reml = FALSE
    
    BinomialDen = NULL
    if (!is.null(input$re_dhglm_m_check_binomd) && input$re_dhglm_m_check_binomd) 
        BinomialDen = data2[[input$re_dhglm_m_binomd]]
    
    betafix = NULL
    if (!is.null(input$re_dhglm_a_check_betafix) && input$re_dhglm_a_check_betafix)
        betafix = as.numeric(input$re_dhglm_a_betafix)
    
    # ¦§ Fitted Model ####
    
    fittedModel <<-
        dhglmfit(
            RespDist = isolate(input$re_dhglm_m_dist),
            BinomialDen = BinomialDen,
            DataMain = data2,
            MeanModel = MM,
            DispersionModel = DM,
            Dmethod = Dmethod,
            mord = mord,
            dord = dord,
            REML = reml,
            BetaFix = betafix
        )

    incProgress(0.3, detail = paste("Loading..."))
    
    # VIF, RobuseSE, ConfInt, Exp
    # This code is for VIF, Robuse SE, Conf Int, Exp. 
    # However, it doesn't work and it isn'd correlated with DHGLM, but glm. 
    # How to revise and how to return it?
    # if (input$re_dhglm_m_check_vif) {
    #     responsePart = input$re_dhglm_m_resp
    #     variablePart = paste0(input$re_dhglm_m_variable$right, collapse = "+")
    #     meanModel = paste0(responsePart, " ~ ", variablePart)
    #     reg <- glm(as.formula(meanModel), family = input$re_dhglm_m_dist, data = data2)
    # }
    # if (input$re_dhglm_m_check_robustse) {
    #     responsePart = input$re_dhglm_m_resp
    #     variablePart = paste0(input$re_dhglm_m_variable$right, collapse = "+")
    #     meanModel = paste0(responsePart, " ~ ", variablePart)
    #     reg <- glm(as.formula(meanModel), family = input$re_dhglm_m_dist, data = data2)
    #     print(coeftest(reg, vcov = vcovHC(reg, "HC1")))
    #     print("print(coeftest(reg, vcov =)")
    # }
    # 
    # if (input$re_dhglm_m_check_confint) {
    #     responsePart = input$re_dhglm_m_resp
    #     variablePart = paste0(input$re_dhglm_m_variable$right, collapse = "+")
    #     meanModel = paste0(responsePart, " ~ ", variablePart)
    #     reg <- glm(as.formula(meanModel), family = input$re_dhglm_m_dist, data = data2)
    #     print(confint(reg))
    #     print("print(confint(reg))")
    # }
    # if (input$re_dhglm_m_check_exp) {
    #     responsePart = input$re_dhglm_m_resp
    #     variablePart = paste0(input$re_dhglm_m_variable$right, collapse = "+")
    #     meanModel = paste0(responsePart, " ~ ", variablePart)
    #     reg <- glm(as.formula(meanModel), family = input$re_dhglm_m_dist, data = data2)
    #     print(exp(confint(reg)))
    #     print("print(exp(confint(reg)))")
    # }

    # if(input$dhglm_m_check_cs==TRUE) MM[[16]]="cubic"
    # if(input$dhglm_m_check_gp==TRUE) MM[[16]]="cov_kernel"
    # if(input$dhglm_m_check_svm==TRUE) MM[[16]]="precision_kernel"
    # if(input$dhglm_p_check_cs==TRUE) DM[[16]]="cubic"
    # if(input$dhglm_p_check_gp==TRUE) DM[[16]]="cov_kernel"
    # if(input$dhglm_p_check_svm==TRUE) DM[[16]]="precision_kernel"
    
    fittedModel$option = list(
        RespDist = isolate(input$re_dhglm_m_dist),
        BinomialDen = BinomialDen,
        DataMain = data2,
        MeanModel = MM,
        DispersionModel = DM,
        Dmethod = Dmethod,
        mord = mord,
        dord = dord,
        REML = reml,
        Betafix = betafix
    ) 
    
    # ¦§ Description 1 ####
    
    nMrand <- length(strsplit(input$re_dhglm_m_model, "|", fixed = TRUE)[[1]]) - 1
    mRandstring = NULL
    pRandstring = "NULL"
    lRandstring = "NULL"
    if (nMrand >= 1)
        for (i in 1:nMrand)
            mRandstring = paste0(c(mRandstring, input[[paste0("re_dhglm_m_randfamily_", i)]]), collapse = ", ")
    else
        mRandstring = "NULL"
        
    modelDesc = matrix(c(
        input$re_dhglm_m_model,
        "constant",
        "lambda ~ 1",
        input$re_dhglm_m_link,
        phiLink,
        "log",
        input$re_dhglm_m_dist,
        "gaussian",
        "gaussian",
        mRandstring,
        "NULL",
        "NULL"
    ),nrow=3, ncol=4)
    if(!is.null(input$re_dhglm_p_check_phi) && input$re_dhglm_p_check_phi) {
        nPrand <- length(strsplit(input$re_dhglm_p_model, "|", fixed = TRUE)[[1]]) - 1
        pRandstring = NULL
        if (nPrand >= 1)
            for (i in 1:nPrand)
                pRandstring = paste0(c(pRandstring, input[[paste0("re_dhglm_p_randfamily_", i)]]), collapse = ", ")
        else
            pRandstring = "NULL"
        
        modelDesc[2, 1] <- phiFormula
        modelDesc[2, 2] <- input$re_dhglm_p_link
        modelDesc[2, 4] <- pRandstring
    }
    if(!is.null(input$re_dhglm_l_check_lambda) && input$re_dhglm_l_check_lambda) {
        nLrand <- length(strsplit(input$re_dhglm_l_model, "|", fixed = TRUE)[[1]]) - 1
        lRandstring = NULL
        if (nLrand >= 1)
            for (i in 1:nLrand)
                lRandstring = paste0(c(lRandstring, input[[paste0("re_dhglm_l_randfamily_", i)]]), collapse = ", ")
        else
            lRandstring = "NULL"
        
        modelDesc[3, 1] <- input$re_dhglm_l_model
        modelDesc[3, 2] <- input$re_dhglm_l_link
        modelDesc[3, 4] <- lRandstring
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
        likeli_coeff <- t(fittedModel$likeli_coeff)
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
    
    if (!is.null(fittedModel$beta_coeff) && !is.null(input$re_dhglm_m_check_exp) && input$re_dhglm_m_check_exp) {
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

    if (!is.null(fittedModel$alpha_coeff)) {
        fittedModel$alpha_coeff <- cbind(fittedModel$alpha_coeff, "exp(Estimate)" = exp(fittedModel$alpha_coeff[, 1]))
        fittedModel$alpha_coeff <- subset(fittedModel$alpha_coeff, select=c(1,4,2,3))
    }
    
    if (!is.null(fittedModel$lambda_coeff)) {
        fittedModel$lambda_coeff <- cbind(fittedModel$lambda_coeff, "exp(Estimate)" = exp(fittedModel$lambda_coeff[, 1]))
        fittedModel$lambda_coeff <- subset(fittedModel$lambda_coeff, select=c(1,4,2,3))
    }
    
    if (!is.null(fittedModel$tau_coeff)) {
        fittedModel$tau_coeff <- cbind(fittedModel$tau_coeff, "exp(Estimate)" = exp(fittedModel$tau_coeff[, 1]))
        fittedModel$tau_coeff <- subset(fittedModel$tau_coeff, select=c(1,4,2,3))
    }
    
    # ¦§ Comparison Model ####
    
    dhglmComparisonmodel1 = NULL
    dhglmComparisonmodel2 = NULL
    dhglmComparisonmodel3 = NULL
    
    if (!is.null(input$re_dhglm_m_check_comparison) && input$re_dhglm_m_check_comparison) {
        spatialVariable = "NULL"
        if (!is.null(spatial))
            spatialVariable = spatial
            
        dhglmComparisonmodel1 <- matrix(
            c(
                input$re_dhglm_m_model,
                input$re_dhglm_m_link,
                input$re_dhglm_m_dist,
                mRandstring,
                spatialVariable
            ),
            ncol = 5
        )
        
        dhglmComparisonmodel2 <- matrix(
            c(
                phiFormula,
                pRandstring,
                lambdaFormula,
                lRandstring
            ),
            ncol = 4
        )
        n.params=fittedModel$n.mean+fittedModel$n.disp
        n.disp=fittedModel$n.disp
        dhglmComparisonmodel3 <- matrix(
            c(
                ifelse(is.null(fittedModel$likeli_coeff), fittedModel$ml, fittedModel$likeli_coeff[1]),
                ifelse(is.null(fittedModel$likeli_coeff), fittedModel$ml+2*n.params, fittedModel$likeli_coeff[1]+2*n.params),
                ifelse(is.null(fittedModel$likeli_coeff), fittedModel$rl, fittedModel$likeli_coeff[2]),
                ifelse(is.null(fittedModel$likeli_coeff), fittedModel$rl+2*n.disp, fittedModel$likeli_coeff[2]+2*n.disp),
                ifelse(is.null(fittedModel$likeli_coeff), fittedModel$sd, fittedModel$likeli_coeff[4]),
                ifelse(is.null(fittedModel$likeli_coeff), fittedModel$df, fittedModel$likeli_coeff[5]),
                ifelse(is.null(fittedModel$likeli_coeff), fittedModel$caic, fittedModel$likeli_coeff[3])
            ),
            ncol = 7
        )
        dhglmComparisonmodel1 <- rbind(g_re_dhglm_r_comparisonmodel_1, dhglmComparisonmodel1)
        dhglmComparisonmodel2 <- rbind(g_re_dhglm_r_comparisonmodel_2, dhglmComparisonmodel2)
        dhglmComparisonmodel3 <- rbind(g_re_dhglm_r_comparisonmodel_3, dhglmComparisonmodel3)
        
        colnames(dhglmComparisonmodel1) <- c("Model", "Link", "Dist", "Rand", "Spatial")
        colnames(dhglmComparisonmodel2) <- c("Phi", "Phi Rand", "Lambda", "Lambda Rand")
        colnames(dhglmComparisonmodel3) <- c("-2ML", "AIC.ML", "-2RL", "AIC.RL", "Scaled Deviance","df", "cAIC")
        
        g_re_dhglm_r_comparisonmodel_1 <<- dhglmComparisonmodel1
        g_re_dhglm_r_comparisonmodel_2 <<- dhglmComparisonmodel2
        g_re_dhglm_r_comparisonmodel_3 <<- dhglmComparisonmodel3

        fittedModel$option$Comparisonmodel1 = dhglmComparisonmodel1
        fittedModel$option$Comparisonmodel2 = dhglmComparisonmodel2
        fittedModel$option$Comparisonmodel3 = dhglmComparisonmodel3
        
        appendTab(inputId = "re_dhglm_resulttabset",
            tabPanel(
                "Model Comparison",
                br(),
                h4("History of Models"),
                tableOutput("re_dhglm_r_comparisonmodel1"),
                tableOutput("re_dhglm_r_comparisonmodel2"),
                tableOutput("re_dhglm_r_comparisonmodel3")
            )
        )
    }
    
    # ¦¦ R Codes ####

    if (!is.null(input$re_dhglm_m_check_rcodes) && input$re_dhglm_m_check_rcodes) {
        Rcode1 <- NULL
        Rcode2 <- NULL
        Rcode3 <- NULL
        offsetRcode = NULL
        if (!is.null(meanOffsetVariable))
            offsetRcode <- paste0(", Offset = data[[\"", input$re_dhglm_m_offset, "\"]]")
        mRandRcode = NULL
        if (nMrand >= 1)
            for (i in 1:nMrand)
                mRandRcode = paste0(c(mRandRcode, input[[paste0("re_dhglm_m_randfamily_", i)]]), collapse = "\", \"")
        
        if (lambdaFormula == "lambda ~ 1") {
            Rcode1 <- paste0(
                "MM <- DHGLMMODELING(Model = \"mean\", Link = \"",
                input$re_dhglm_m_link,
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
                input$re_dhglm_m_link,
                "\", LinPred = ",
                meanFormula,
                ", RandDist = c(\"",
                mRandRcode,
                "\"), LinPredRandVariance = ",
                lambdaFormula,
                ", RandDistRandVariance = ",
                lambdaRandFamily,
                ", LinkRandVariance = \"log\"",
                offsetRcode,
                ")"
            )
        }
        
        pRandRcode = NULL
        if (!is.null(phiRandFamily) && phiRandFamily != "")
            pRandRcode <- paste0(", RandDist = \"", phiRandFamily, "\"")
        
        if(phiFormula == "constant") 
            Rcode2 <- "DM <- DHGLMMODELING(Model = \"dispersion\", Link = \"log\")"
        else
            Rcode2 <- paste0(
                "DM <- DHGLMMODELING(Model = \"dispersion\", Link = \"", phiLink, "\", LinPred = ",
                phiFormula,
                pRandRcode,
                ")"
            )
        
        Rcode3 <- paste0(
            "dhglmfit(RespDist = \"", 
            input$re_dhglm_m_dist,
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

# DHGLM UI (Accordion) ####

output$re_dhglm_r_accordion<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    dhglmAccordion <- bs_accordion_sidebar(
        id = "re_dhglm_r_accordion",
        spec_side = c(width = 3, offset = 0),
        spec_main = c(width = 9, offset = 0)    
    )
    
    # ¦§ Mean ####
    
    dhglmAccordion <- dhglmAccordion %>%
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
                uiOutput("re_dhglm_m_resp"), 
                uiOutput("re_dhglm_m_variable"), 
                fluidRow(
                    column(
                        8,
                        uiOutput("re_dhglm_m_interaction") 
                    ),
                    column(
                        4,
                        uiOutput("re_dhglm_m_interactionappend"),
                        style = "text-align:right; padding:15px"
                    )
                ),
                uiOutput("re_dhglm_m_rand"),
                uiOutput("re_dhglm_m_check_slope"),
                uiOutput("re_dhglm_m_check_withoutslope"),
                uiOutput("re_dhglm_m_slope"),
                uiOutput("re_dhglm_m_check_randinteraction"),
                fluidRow(
                    column(
                        8,
                        uiOutput("re_dhglm_m_randinteraction")
                    ),
                    column(
                        4,
                        uiOutput("re_dhglm_m_randinteractionappend"),
                        style = "text-align:right; padding:15px"
                    )
                ), #new components
                
                uiOutput("re_dhglm_m_dist"),
                uiOutput("re_dhglm_m_check_binomd"),
                uiOutput("re_dhglm_m_binomd"),
                uiOutput("re_dhglm_m_link"),
                
                hr(style = "border-color: #2C3E50;"),
                uiOutput("re_dhglm_m_randfamily"),
                uiOutput("re_dhglm_m_upload_neighbor"),

                uiOutput("re_dhglm_m_check_nointercept"),
                uiOutput("re_dhglm_m_check_offset"),
                uiOutput("re_dhglm_m_offset"),
                
                # They don't work
                # uiOutput("re_dhglm_m_check_vif"),
                # uiOutput("re_dhglm_m_check_robustse"),
                # uiOutput("re_dhglm_m_check_confint"),
                uiOutput("re_dhglm_m_check_exp"),
                
                hr(style = "border-color: #2C3E50;"),
                uiOutput("re_dhglm_m_check_comparison"),
                uiOutput("re_dhglm_m_check_rcodes")
            )
        )
    
    # ¦§ Phi ####
    
    dhglmAccordion <- dhglmAccordion %>%
    bs_append(
        title_side = "Phi",
        content_side = uiOutput("re_dhglm_p_check_phi"),
        content_main = div(
            h3(strong("Model for Phi")),
            uiOutput("re_dhglm_p_resp"),
            uiOutput("re_dhglm_p_variable"),
            fluidRow(
                column(
                    8,
                    uiOutput("re_dhglm_p_interaction")
                ),
                column(
                    4,
                    uiOutput("re_dhglm_p_interactionappend"),
                    style = "text-align:right; padding:15px"
                )
            ),
            uiOutput("re_dhglm_p_rand"),
            uiOutput("re_dhglm_p_check_slope"),
            uiOutput("re_dhglm_p_check_withoutslope"),
            uiOutput("re_dhglm_p_slope"),
            uiOutput("re_dhglm_p_check_randinteraction"),
            fluidRow(
                column(
                    8,
                    uiOutput("re_dhglm_p_randinteraction")
                ),
                column(
                    4,
                    uiOutput("re_dhglm_p_randinteractionappend"),
                    style = "text-align:right; padding:15px"
                )
            ), #new components

            uiOutput("re_dhglm_p_link"),
            uiOutput("re_dhglm_p_randfamily"),

            uiOutput("re_dhglm_p_check_nointercept"),
            uiOutput("re_dhglm_p_check_offset"),
            uiOutput("re_dhglm_p_offset")
        )
    )
    
    # ¦§ Lambda ####
    
    dhglmAccordion <- dhglmAccordion %>%
    bs_append(
        title_side = "Lambda",
        content_side = uiOutput("re_dhglm_l_check_lambda"),
        content_main = div(
            h3(strong("Model for Lambda")),
            uiOutput("re_dhglm_l_resp"),
            uiOutput("re_dhglm_l_variable"),
            fluidRow(
                column(
                    8,
                    uiOutput("re_dhglm_l_interaction")
                ),
                column(
                    4,
                    uiOutput("re_dhglm_l_interactionappend"),
                    style = "text-align:right; padding:15px"
                )
            ),
            uiOutput("re_dhglm_l_rand"),
            uiOutput("re_dhglm_l_check_randinteraction"),
            fluidRow(
                column(
                    8,
                    uiOutput("re_dhglm_l_randinteraction")
                ),
                column(
                    4,
                    uiOutput("re_dhglm_l_randinteractionappend"),
                    style = "text-align:right; padding:15px"
                )
            ), 
    
            uiOutput("re_dhglm_l_link"),
            uiOutput("re_dhglm_l_randfamily"),
    
            uiOutput("re_dhglm_l_check_nointercept")
            # uiOutput("re_dhglm_l_check_offset"),
            # uiOutput("re_dhglm_l_offset")
            # DHGLM doesn't support offset variable in lambda
        )
    )
    
    # ¦¦ Additional Settings ####
    
    dhglmAccordion <- dhglmAccordion %>%
    bs_append(
        title_side = "Setting",
        content_side = NULL,
        content_main = div(
            h3(strong("Additional Settings")),
            uiOutput("re_dhglm_a_dispersion"),
            uiOutput("re_dhglm_a_reml"),
            # uiOutput("re_dhglm_a_check_orthogonal"),
#            uiOutput("re_dhglm_a_check_betafix"),
            uiOutput("re_dhglm_a_betafix"),
            uiOutput("re_dhglm_a_mord"),
            uiOutput("re_dhglm_a_dord"),
            h5(helpText("Order of Laplace Approximation for Likelihood(mean) and Restricted Likelihood(Dispersion)")),
#            uiOutput("re_dhglm_a_check_correlation"),
#           uiOutput("re_dhglm_accordion_c")
        )
    )
    


    # ¦¦ Correlation Structure ####

  dhglmAccordion <- dhglmAccordion %>%
    bs_append(
     title_side = "Correlation Structure",
      content_side = NULL,
      content_main = div(
        uiOutput("re_dhglm_corr_variable"),
        fluidRow(
          column(
            6,
            uiOutput("re_dhglm_corr_structure_1")
          ),
          column(
            6,
            uiOutput("re_dhglm_corr_structure_2")
          )
        ),
        fluidRow(
          column(
            6,
            splitLayout(style = "padding: 6px;"),
            actionButton("re_dhglm_corr_structure_1_run", "Add", icon = icon("play")),
            style = "text-align:center;"
          ) ,
          column(
            6,
            splitLayout(style = "padding: 6px;"),
            actionButton("re_dhglm_corr_structure_2_run", "Add", icon = icon("play")),
            style = "text-align:center;"
          ) 
        ),
        fluidRow(
          column(
            6,
            uiOutput("re_dhglm_corr_structure_11")
          ),
          column(
            6,
            uiOutput("re_dhglm_corr_structure_21")
          )
        )
        
      )
    )
    div(
      dhglmAccordion,
      use_bs_tooltip(),
      use_bs_accordion_sidebar() # needs to be at end, for some reason
    )
})


# DHGLM Components ####
# ¦§ Mean ####

output$re_dhglm_m_model<-renderUI({ #algorithm not good
    input$re_dhglm_m_resp
    input$re_dhglm_m_variable$right
    input$re_dhglm_m_rand
    input$re_dhglm_m_check_nointercept
    input$re_dhglm_m_binomd
    input$re_dhglm_m_slope
    
    if(any(is.null(input$re_dhglm_m_resp), is.null(input$re_dhglm_m_variable)))
        return()
    
    responsePart = input$re_dhglm_m_resp
    variablePart = paste0(input$re_dhglm_m_variable$right, collapse=" + ")
    randomPart = NULL
    meanModel = NULL
    
    # in DHGLM, even if we use binomial denominator, model equation doesn't change. 
    # if (all(!is.null(input$re_dhglm_m_check_binomd), input$re_dhglm_m_check_binomd, 
    #         !is.null(input$re_dhglm_m_binomd), input$re_dhglm_m_binomd != ""))
    #     responsePart = paste0("cbind(", responsePart, ",", input$re_dhglm_m_binomd,"-", responsePart, ")")
    
    if(!is.null(input$re_dhglm_m_rand) && length(input$re_dhglm_m_rand) > 0) {
        if (is.null(input$re_dhglm_m_check_withoutslope) || !input$re_dhglm_m_check_withoutslope) {
            randomPart = paste0(input$re_dhglm_m_rand, collapse = ") + (1|")
            randomPart = paste0("+ (1|", randomPart, ")")
        }
        
        #random slope (a|b) form
        if(all(!is.null(input$re_dhglm_m_check_slope), input$re_dhglm_m_check_slope, length(input$re_dhglm_m_slope) > 0)) 
            for(i in 1:length(input$re_dhglm_m_rand)) 
                for(j in 1:length(input$re_dhglm_m_slope)) 
                    randomPart = paste0(randomPart," + (", input$re_dhglm_m_slope[j], "|", input$re_dhglm_m_rand[i], ")")
    }
    
    if(!is.null(input$re_dhglm_m_check_nointercept) && !input$re_dhglm_m_check_nointercept) {
        if(variablePart == "")
            meanModel = paste(responsePart, "~", "1", randomPart)
        else
            meanModel = paste(responsePart, "~", variablePart, randomPart)
    }
    
    if(!is.null(input$re_dhglm_m_check_nointercept) && input$re_dhglm_m_check_nointercept) {
        if(variablePart == "")
            meanModel = paste(responsePart, "~", "-1", randomPart)
        else
            meanModel = paste(responsePart, "~", "-1 +", variablePart, randomPart)
    }
   
    textAreaInput("re_dhglm_m_model", "Model for Mean", value = meanModel, height = "60px")
})

output$re_dhglm_m_resp <- renderUI({
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
            "re_dhglm_m_resp",
            "Response Variable",
            choices = as.list(c("", nameValue)), 
            selected = NULL, 
            multiple = FALSE
        ),
        bsTooltip(
            "re_dhglm_m_resp", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$re_dhglm_m_variable <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$re_dhglm_m_interactionappend
    
    nameValue <- c(names(data2), g_re_dhglm_m_interaction)
    
    chooserInput("re_dhglm_m_variable", "Variable", "Selected", nameValue, c(), size = 15, multiple = TRUE, width = 125)
})

output$re_dhglm_m_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue <- names(data2)
    
    selectInput(
        "re_dhglm_m_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$re_dhglm_m_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    actionButton("re_dhglm_m_interactionappend", "Append")
})

observeEvent(input$re_dhglm_m_interactionappend, {
    if(length(input$re_dhglm_m_interaction) > 1)
        g_re_dhglm_m_interaction <<- c(g_re_dhglm_m_interaction, paste(input$re_dhglm_m_interaction, collapse=":"))
})


output$re_dhglm_m_rand <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$re_dhglm_m_randinteractionappend
    
    nameValue <- c(names(data2), g_re_dhglm_m_randinteraction)
    
    selectInput(
        "re_dhglm_m_rand",
        "Random Effects",
        choices = as.list(nameValue), 
        multiple = TRUE
    )
})

output$re_dhglm_m_check_slope <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_dhglm_m_check_slope", "Random Slope Model", value = FALSE)
})

output$re_dhglm_m_check_withoutslope <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if (is.null(input$re_dhglm_m_check_slope) || !input$re_dhglm_m_check_slope)
        return()
    
    checkboxInput("re_dhglm_m_check_withoutslope", "Without Random Intercept", value = FALSE)
})

output$re_dhglm_m_slope <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if (is.null(input$re_dhglm_m_check_slope) || !input$re_dhglm_m_check_slope)
        return()
    
    nameValue = names(data2)
    
    selectInput(
        "re_dhglm_m_slope",
        "Random slope",
        choices = as.list(nameValue), 
        multiple = TRUE
    )
})

output$re_dhglm_m_check_randinteraction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_dhglm_m_check_randinteraction", "Interaction in the Random Effect", value = FALSE)
})

output$re_dhglm_m_randinteraction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$re_dhglm_m_check_randinteraction) || !input$re_dhglm_m_check_randinteraction)
        return()
    
    nameValue <- names(data2)
    
    selectInput(
        "re_dhglm_m_randinteraction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$re_dhglm_m_randinteractionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_dhglm_m_check_randinteraction) || !input$re_dhglm_m_check_randinteraction)
        return()
    
    actionButton("re_dhglm_m_randinteractionappend", "Append")
})

observeEvent(input$re_dhglm_m_randinteractionappend, {
    
    if(length(input$re_dhglm_m_randinteraction) > 1)
        g_re_dhglm_m_randinteraction <<- c(g_re_dhglm_m_randinteraction, paste(input$re_dhglm_m_randinteraction, collapse=":"))
})

output$re_dhglm_m_dist<-renderUI({
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
        "re_dhglm_m_dist",
        "Distribution", 
        choices = dist, 
        multiple = FALSE
    )
})

output$re_dhglm_m_check_binomd<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(any(is.null(input$re_dhglm_m_dist), input$re_dhglm_m_dist != "binomial"))
        return()
        
    checkboxInput("re_dhglm_m_check_binomd", "Binomial Denominator", value = FALSE)
})

output$re_dhglm_m_binomd<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$re_dhglm_m_check_binomd) || !(input$re_dhglm_m_check_binomd))
        return()
        
    nameValue = names(select_if(data2, is.numeric))
    
    selectInput(
        "re_dhglm_m_binomd",
        "Binomial Denominator",
        choices = as.list(c("", nameValue)),
        multiple = FALSE
    )
})

output$re_dhglm_m_link <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_dhglm_m_dist))
        return()
    
    selection = "identity"
    
    if (input$re_dhglm_m_dist == "binomial")
        selection = "logit"
    else if (input$re_dhglm_m_dist == "poisson")
        selection = "log"
    else if (input$re_dhglm_m_dist == "gamma")
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
        "re_dhglm_m_link",
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
    
    if(is.null(input$re_dhglm_m_model) || is.null(input$re_dhglm_m_rand) || length(input$re_dhglm_m_rand) == 0) {
        output$re_dhglm_m_randfamily <- renderUI({
            return()
        })
    }
    
    if(length(input$re_dhglm_m_rand) > 0) {
        numberRandom <- length(strsplit(input$re_dhglm_m_model, "|", fixed = TRUE)[[1]]) - 1
        for (i in 1:numberRandom) {
            output$re_dhglm_m_randfamily <- renderUI({
                re_dhglm_m_randfamily_list <- lapply(1:numberRandom, function(i) {
                    selectInput(
                        paste0("re_dhglm_m_randfamily_", i),
                        paste("Distribution for Random effects", i),
                        choices = c(
                            "gaussian",
                            "gamma",
                            "beta",
                            "inverse-gamma",
                            "saturated",
                            "shared",
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
                do.call(tagList, re_dhglm_m_randfamily_list)
            })
        }
    }
    else {
        output$re_dhglm_m_randfamily <- renderUI({
            return()
        })
    }
})

# Neighbor for IAR, MRF
# It only reacts with randfamily 1. 
output$re_dhglm_m_upload_neighbor <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_dhglm_m_randfamily_1) || (input$re_dhglm_m_randfamily_1 != "IAR" && input$re_dhglm_m_randfamily_1 != "MRF"))
        return()
    
    fileInput("re_dhglm_m_neighborfile", "Upload Neighborhood File", multiple = TRUE)
})

re_dhglm_m_neighborfile <- reactive({
    re_dhglm_m_neighborfile <- input$re_dhglm_m_neighborfile
    if (is.null(re_dhglm_m_neighborfile)) {
        return()
    }
    read.table(
        file = re_dhglm_m_neighborfile$datapath,
        sep = input$sep,
        header = input$header,
        stringsAsFactors = TRUE
    )
})

output$re_dhglm_m_check_nointercept<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_dhglm_m_check_nointercept", "No Intercept Model", value = FALSE)
})

output$re_dhglm_m_check_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_dhglm_m_check_offset", "Offset Variable", value = FALSE)
})

output$re_dhglm_m_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue = names(select_if(data2, is.numeric))
    
    if(is.null(input$re_dhglm_m_check_offset) || !input$re_dhglm_m_check_offset)
        return()
       
    div(
        selectInput(
            "re_dhglm_m_offset",
            "Offset Variable",
            choices = as.list(c("", nameValue)),
            selected = NULL,
            multiple = FALSE
        ),
        bsTooltip(
            "re_dhglm_m_offset", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$re_dhglm_m_check_vif<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_dhglm_m_check_vif", "VIF", value = FALSE)
})

output$re_dhglm_m_check_robustse<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_dhglm_m_check_robustse", "Robust Standard Errors", value = FALSE)
})

output$re_dhglm_m_check_confint<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_dhglm_m_check_confint", "Confidence Intervals for Coefficients", value = FALSE)
})

output$re_dhglm_m_check_exp <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_dhglm_m_check_exp", "Exponential scale", value = FALSE)
})

output$re_dhglm_m_check_comparison<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_dhglm_m_check_comparison", "Model Comparison", value = FALSE)
})

output$re_dhglm_m_check_rcodes<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_dhglm_m_check_rcodes", "R Codes", value = FALSE)
})

# ¦§ Phi ####

output$re_dhglm_p_check_phi<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_dhglm_p_check_phi", strong("Use"), value = FALSE)
})

output$re_dhglm_p_model<-renderUI({ #algorithm not good
    input$re_dhglm_p_resp
    input$re_dhglm_p_variable$right
    input$re_dhglm_p_rand
    input$re_dhglm_p_check_nointercept
    input$re_dhglm_p_slope
    
    if(is.null(input$re_dhglm_p_check_phi) || !input$re_dhglm_p_check_phi)
        return()
    
    if(any(is.null(input$re_dhglm_p_resp), is.null(input$re_dhglm_p_variable)))
        return()
    
    responsePart = input$re_dhglm_p_resp
    variablePart = paste0(input$re_dhglm_p_variable$right, collapse=" + ")
    randomPart = NULL
    phiModel = NULL
    
    if(!is.null(input$re_dhglm_p_rand) && length(input$re_dhglm_p_rand) > 0) {
        if (is.null(input$re_dhglm_p_check_withoutslope) || !input$re_dhglm_p_check_withoutslope) {
            randomPart = paste0(input$re_dhglm_p_rand, collapse = ") + (1|")
            randomPart = paste0("+ (1|", randomPart, ")")
        }
        
        #random slope (a|b) form
        if(all(!is.null(input$re_dhglm_p_check_slope), input$re_dhglm_p_check_slope, length(input$re_dhglm_p_slope) > 0)) 
            for(i in 1:length(input$re_dhglm_p_rand)) 
                for(j in 1:length(input$re_dhglm_p_slope)) 
                    randomPart = paste0(randomPart," + (", input$re_dhglm_p_slope[j], "|", input$re_dhglm_p_rand[i], ")")
    }
    
    if(!is.null(input$re_dhglm_p_check_nointercept) && !input$re_dhglm_p_check_nointercept) {
        if(variablePart == "")
            phiModel = paste(responsePart, "~", "1", randomPart)
        else
            phiModel = paste(responsePart, "~", variablePart, randomPart)
    }
    
    if(!is.null(input$re_dhglm_p_check_nointercept) && input$re_dhglm_p_check_nointercept) {
        if(variablePart == "")
            phiModel = paste(responsePart, "~", "-1", randomPart)
        else
            phiModel = paste(responsePart, "~", "-1 +", variablePart, randomPart)
    }
    
    textAreaInput("re_dhglm_p_model", "Model for Phi", value = phiModel, height = "60px")
})

output$re_dhglm_p_resp <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$re_dhglm_p_check_phi) || !input$re_dhglm_p_check_phi)
        return()
    
    nameValue = c("phi" = "phi")
    
    selectInput(
        "re_dhglm_p_resp",
        "Residual Variance",
        choices = nameValue,
        selected = NULL, 
        multiple = FALSE
    )
})

output$re_dhglm_p_variable <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$re_dhglm_p_interactionappend
    
    if(is.null(input$re_dhglm_p_check_phi) || !input$re_dhglm_p_check_phi)
        return()
    
    nameValue <- c(names(data2), g_re_dhglm_p_interaction)
    
    chooserInput("re_dhglm_p_variable", "Variable", "Selected", nameValue, c(), size = 15, multiple = TRUE, width = 150)
})

output$re_dhglm_p_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$re_dhglm_p_check_phi) || !input$re_dhglm_p_check_phi)
        return()
    
    nameValue <- names(data2)
    
    selectInput(
        "re_dhglm_p_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$re_dhglm_p_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_dhglm_p_check_phi) || !input$re_dhglm_p_check_phi)
        return()
    
    actionButton("re_dhglm_p_interactionappend", "Append")
})

observeEvent(input$re_dhglm_p_interactionappend, {
    if(length(input$re_dhglm_p_interaction) > 1)
        g_re_dhglm_p_interaction <<- c(g_re_dhglm_p_interaction, paste(input$re_dhglm_p_interaction, collapse=":"))
})


output$re_dhglm_p_rand <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$re_dhglm_p_randinteractionappend
    
    if(is.null(input$re_dhglm_p_check_phi) || !input$re_dhglm_p_check_phi)
        return()
    
    nameValue <- c(names(data2), g_re_dhglm_p_randinteraction)
    
    selectInput(
        "re_dhglm_p_rand",
        "Random Effects",
        choices = as.list(nameValue), 
        multiple = TRUE
    )
})

output$re_dhglm_p_check_slope <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_dhglm_p_check_phi) || !input$re_dhglm_p_check_phi)
        return()
    
    checkboxInput("re_dhglm_p_check_slope", "Random Slope Model", value = FALSE)
})

output$re_dhglm_p_check_withoutslope <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_dhglm_p_check_phi) || !input$re_dhglm_p_check_phi)
        return()
    
    if (is.null(input$re_dhglm_p_check_slope) || !input$re_dhglm_p_check_slope)
        return()
    
    checkboxInput("re_dhglm_p_check_withoutslope", "Without Random Intercept", value = FALSE)
})

output$re_dhglm_p_slope <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$re_dhglm_p_check_phi) || !input$re_dhglm_p_check_phi)
        return()
    
    nameValue = names(data2)
    
    if (!is.null(input$re_dhglm_p_check_slope) && input$re_dhglm_p_check_slope) {
        selectInput(
            "re_dhglm_p_slope",
            "Random Slope",
            choices = as.list(nameValue), #maybe have to revise?
            multiple = TRUE
        )
    }
})

output$re_dhglm_p_check_randinteraction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_dhglm_p_check_phi) || !input$re_dhglm_p_check_phi)
        return()
    
    checkboxInput("re_dhglm_p_check_randinteraction", "Interaction in the Random Effect", value = FALSE)
})

output$re_dhglm_p_randinteraction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$re_dhglm_p_check_phi) || !input$re_dhglm_p_check_phi)
        return()
    
    nameValue <- names(data2)
    
    if(!is.null(input$re_dhglm_p_check_randinteraction) && input$re_dhglm_p_check_randinteraction)
        selectInput(
            "re_dhglm_p_randinteraction",
            "Make Interaction Variable",
            choices = as.list(c("", nameValue)), 
            selected = NULL, 
            multiple = TRUE
        )
})

output$re_dhglm_p_randinteractionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_dhglm_p_check_phi) || !input$re_dhglm_p_check_phi)
        return()
    
    if(!is.null(input$re_dhglm_p_check_randinteraction) && input$re_dhglm_p_check_randinteraction)
        actionButton("re_dhglm_p_randinteractionappend", "Append")
})

observeEvent(input$re_dhglm_p_randinteractionappend, {
    if(length(input$re_dhglm_p_randinteraction) > 1)
        g_re_dhglm_p_randinteraction <<- c(g_re_dhglm_p_randinteraction, paste(input$re_dhglm_p_randinteraction, collapse=":"))
})

output$re_dhglm_p_link <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_dhglm_p_check_phi) || !input$re_dhglm_p_check_phi)
        return()

    link = c(
        "log" = "log",
        "inverse" = "inverse",
        "identity" = "identity")
    
    selectInput(
        "re_dhglm_p_link",
        "Link Function",
        choices = link,
        multiple = FALSE
    )
})

observe({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_p_model
    input$re_dhglm_p_rand
    input$re_dhglm_p_check_phi
    
    if(is.null(input$re_dhglm_p_model) || is.null(input$re_dhglm_p_rand) || length(input$re_dhglm_p_rand) == 0
       || is.null(input$re_dhglm_p_check_phi) || !input$re_dhglm_p_check_phi ) {
        output$re_dhglm_p_randfamily <- renderUI({
            return()
        })
    }
    
    if(!is.null(input$re_dhglm_p_check_phi) && input$re_dhglm_p_check_phi && length(input$re_dhglm_p_rand) > 0) {
        numberRandom <- length(strsplit(input$re_dhglm_p_model, "|", fixed = TRUE)[[1]]) - 1
        for (i in 1:numberRandom) {
            output$re_dhglm_p_randfamily <- renderUI({
                re_dhglm_p_randfamily_list <- lapply(1:numberRandom, function(i) {
                    selectInput(
                        paste0("re_dhglm_p_randfamily_", i),
                        paste("Distribution for Random effects", i),
                        choices = c(
                            "gaussian",
                            "gamma",
                            "inverse-gamma",
                            "AR(1)",
                            "GARCH",
                            "cubic spline"
                        )
                    )
                })
                do.call(tagList, re_dhglm_p_randfamily_list)
            })
        }
    }
    else {
        output$re_dhglm_p_randfamily <- renderUI({
            return()
        })
    }
})

output$re_dhglm_p_check_nointercept<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_p_check_phi
    
    if(is.null(input$re_dhglm_p_check_phi) || !input$re_dhglm_p_check_phi)
        return()
    
    checkboxInput("re_dhglm_p_check_nointercept", "No Intercept Model", value = FALSE)
})

output$re_dhglm_p_check_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_p_check_phi
    
    if(is.null(input$re_dhglm_p_check_phi) || !input$re_dhglm_p_check_phi)
        return()
    
    checkboxInput("re_dhglm_p_check_offset", "Offset Variable", value = FALSE)
})

output$re_dhglm_p_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$re_dhglm_p_check_phi
    
    if(is.null(input$re_dhglm_p_check_phi) || !input$re_dhglm_p_check_phi)
        return()
    
    if(is.null(input$re_dhglm_p_check_offset) || !input$re_dhglm_p_check_offset)
        return()
    
    nameValue = names(select_if(data2, is.numeric))
    
    div(
        selectInput(
            "re_dhglm_p_offset",
            "Offset Variable",
            choices = as.list(c("", nameValue)),
            selected = NULL,
            multiple = FALSE
        ),
        bsTooltip(
            "re_dhglm_p_offset", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

# ¦§ Lambda ####

output$re_dhglm_l_check_lambda<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_dhglm_l_check_lambda", strong("Use"), value = FALSE)
})

output$re_dhglm_l_model<-renderUI({ #algorithm not good
    input$re_dhglm_l_resp
    input$re_dhglm_l_variable$right
    input$re_dhglm_l_rand
    input$re_dhglm_l_check_nointercept
    input$re_dhglm_l_slope
    input$re_dhglm_l_check_lambda
    
    if(is.null(input$re_dhglm_l_check_lambda) || !input$re_dhglm_l_check_lambda)
        return()
    
    if(any(is.null(input$re_dhglm_l_resp), is.null(input$re_dhglm_l_variable)))
        return()
    
    responsePart = input$re_dhglm_l_resp
    variablePart = paste0(input$re_dhglm_l_variable$right, collapse=" + ")
    randomPart = NULL
    lambdaModel = NULL
    
    if(!is.null(input$re_dhglm_l_rand) && length(input$re_dhglm_l_rand) > 0) {
        randomPart = paste0(input$re_dhglm_l_rand, collapse = ") + (1|")
        randomPart = paste0("+ (1|", randomPart, ")")
    }
    
    if(!is.null(input$re_dhglm_l_check_nointercept) && !input$re_dhglm_l_check_nointercept) {
        if(variablePart == "")
            lambdaModel = paste(responsePart, "~", "1", randomPart)
        else
            lambdaModel = paste(responsePart, "~", variablePart, randomPart)
    }
    
    if(!is.null(input$re_dhglm_l_check_nointercept) && input$re_dhglm_l_check_nointercept) {
        if(variablePart == "")
            lambdaModel = paste(responsePart, "~", "-1", randomPart)
        else
            lambdaModel = paste(responsePart, "~", "-1 +", variablePart, randomPart)
    }
    
    textAreaInput("re_dhglm_l_model", "Model for Lambda", value = lambdaModel, height = "60px")
})

output$re_dhglm_l_resp <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$re_dhglm_l_check_lambda
    
    if(is.null(input$re_dhglm_l_check_lambda) || !input$re_dhglm_l_check_lambda)
        return()
    
    nameValue = c("lambda" = "lambda")
    
    selectInput(
        "re_dhglm_l_resp",
        "Variance of Random Effect",
        choices = nameValue,
        selected = NULL, 
        multiple = FALSE
    )
})

output$re_dhglm_l_variable <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$re_dhglm_l_interactionappend
    input$re_dhglm_l_check_lambda
    
    if(is.null(input$re_dhglm_l_check_lambda) || !input$re_dhglm_l_check_lambda)
        return()
    
    nameValue <- c(names(data2), g_re_dhglm_l_interaction)
    
    chooserInput("re_dhglm_l_variable", "Variable", "Selected", nameValue, c(), size = 15, multiple = TRUE, width = 150)
})

output$re_dhglm_l_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$re_dhglm_l_check_lambda
    
    if(is.null(input$re_dhglm_l_check_lambda) || !input$re_dhglm_l_check_lambda)
        return()
    
    nameValue <- names(data2)
    
    selectInput(
        "re_dhglm_l_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$re_dhglm_l_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_l_check_lambda
    
    if(is.null(input$re_dhglm_l_check_lambda) || !input$re_dhglm_l_check_lambda)
        return()
    
    actionButton("re_dhglm_l_interactionappend", "Append")
})

observeEvent(input$re_dhglm_l_interactionappend, {
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(length(input$re_dhglm_l_interaction) > 1)
        g_re_dhglm_l_interaction <<- c(g_re_dhglm_l_interaction, paste(input$re_dhglm_l_interaction, collapse=":"))
})


output$re_dhglm_l_rand <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$re_dhglm_l_randinteractionappend
    input$re_dhglm_l_check_lambda
    
    if(is.null(input$re_dhglm_l_check_lambda) || !input$re_dhglm_l_check_lambda)
        return()
    
    nameValue <- c(names(data2), g_re_dhglm_l_randinteraction)
    
    selectInput(
        "re_dhglm_l_rand",
        "Random Effects",
        choices = as.list(nameValue), 
        multiple = TRUE
    )
})

output$re_dhglm_l_check_randinteraction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_l_check_lambda
    
    if(is.null(input$re_dhglm_l_check_lambda) || !input$re_dhglm_l_check_lambda)
        return()
    
    checkboxInput("re_dhglm_l_check_randinteraction", "Interaction in the Random Effect", value = FALSE)
})


output$re_dhglm_l_randinteraction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$re_dhglm_l_check_randinteraction
    input$re_dhglm_l_check_lambda
    
    if(is.null(input$re_dhglm_l_check_lambda) || !input$re_dhglm_l_check_lambda)
        return()
    
    nameValue <- names(data2)
    
    if(!is.null(input$re_dhglm_l_check_randinteraction) && input$re_dhglm_l_check_randinteraction)
        selectInput(
            "re_dhglm_l_randinteraction",
            "Make Interaction Variable",
            choices = as.list(c("", nameValue)), 
            selected = NULL, 
            multiple = TRUE
        )
})

output$re_dhglm_l_randinteractionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_l_check_randinteraction
    input$re_dhglm_l_check_lambda
    
    if(is.null(input$re_dhglm_l_check_lambda) || !input$re_dhglm_l_check_lambda)
        return()
    
    if(!is.null(input$re_dhglm_l_check_randinteraction) && input$re_dhglm_l_check_randinteraction)
        actionButton("re_dhglm_l_randinteractionappend", "Append")
})

observeEvent(input$re_dhglm_l_randinteractionappend, {
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(length(input$re_dhglm_l_randinteraction) > 1)
        g_re_dhglm_l_randinteraction <<- c(g_re_dhglm_l_randinteraction, paste(input$re_dhglm_l_randinteraction, collapse=":"))
})

output$re_dhglm_l_link <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_l_check_lambda
    
    if(is.null(input$re_dhglm_l_check_lambda) || !input$re_dhglm_l_check_lambda)
        return()
    
    link = c("log" = "log")
    
    selectInput(
        "re_dhglm_l_link",
        "Link Function",
        choices = link,
        multiple = FALSE
    )
})

observe({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_l_model
    input$re_dhglm_l_rand
    input$re_dhglm_l_check_lambda
    
    if(is.null(input$re_dhglm_l_model) || is.null(input$re_dhglm_l_rand) || length(input$re_dhglm_l_rand) == 0
       || is.null(input$re_dhglm_l_check_lambda) || !input$re_dhglm_l_check_lambda ) {
        output$re_dhglm_l_randfamily <- renderUI({
            return()
        })
    }
    
    
    if(!is.null(input$re_dhglm_l_check_lambda) && input$re_dhglm_l_check_lambda && length(input$re_dhglm_l_rand) > 0) {
        numberRandom <- length(strsplit(input$re_dhglm_l_model, "|", fixed = TRUE)[[1]]) - 1
        for (i in 1:numberRandom) {
            output$re_dhglm_l_randfamily <- renderUI({
                re_dhglm_l_randfamily_list <- lapply(1:numberRandom, function(i) {
                    # radioButtons(
                    #     inputId = paste0("mVariable", i),
                    #     label = paste0("mVariable", i),
                    #     choices = c("A", "B", "C")
                    # )
                    selectInput(
                        paste0("re_dhglm_l_randfamily_", i),
                        paste("Distribution for Random effects", i),
                        choices = c(
                            "gaussian",
                            "gamma",
                            "inverse-gamma"
                        )
                    )
                })
                do.call(tagList, re_dhglm_l_randfamily_list)
            })
        }
    }
    else {
        output$re_dhglm_l_randfamily <- renderUI({
            return()
        })
    }
})

output$re_dhglm_l_check_nointercept<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_l_check_lambda
    
    if(is.null(input$re_dhglm_l_check_lambda) || !input$re_dhglm_l_check_lambda)
        return()
    
    checkboxInput("re_dhglm_l_check_nointercept", "No Intercept Model", value = FALSE)
})

# DHGLMfit doesn't support offset variable in lambda
# output$re_dhglm_l_check_offset<-renderUI({
#     input$file
#     input$resetData
#     input$re_dhglm_l_check_lambda
#     
#     if(is.null(input$re_dhglm_l_check_lambda) || !input$re_dhglm_l_check_lambda)
#         return()
#     
#     checkboxInput("re_dhglm_l_check_offset", "Offset Variable", value = FALSE)
# })
# 
# output$re_dhglm_l_offset<-renderUI({
#     input$file
#     input$resetData
#     input$dm_mv_run
#     input$dm_ct_run
#     input$dm_sr_run
#     input$dm_md_run
#     input$re_dhglm_l_check_lambda
#     
#     if(is.null(input$re_dhglm_l_check_lambda) || !input$re_dhglm_l_check_lambda)
#         return()
#     
#     nameValue = names(select_if(data2, is.numeric))
#     
#     if(all(!is.null(input$re_dhglm_l_check_offset), input$re_dhglm_l_check_offset)) {
#         selectInput(
#             "re_dhglm_l_offset",
#             "Offset Variable (Numeric Only)",
#             choices = as.list(c("", nameValue)),
#             selected = NULL,
#             multiple = FALSE
#         )
#     }
# })

# ¦¦ Additional Settings ####

output$re_dhglm_a_dispersion<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    dispersionmethods = c(
        "deviance" = "deviance", 
        "Pearson" = "Pearson"
    )
    
    selectInput(
        "re_dhglm_a_dispersion",
        "Method of Fitting Dispersion Model", 
        choices = dispersionmethods, 
        multiple = FALSE
    )
})

output$re_dhglm_a_reml <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    selectInput(
        "re_dhglm_a_reml",
        "REML or ML",
        choice = c("REML" = TRUE,"ML" = FALSE),
        multiple = FALSE
    )
})

output$re_dhglm_a_check_betafix <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput(
        "re_dhglm_a_check_betafix",
        "Fix Beta Values",
        value = FALSE
    )
})

output$re_dhglm_a_betafix <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_a_check_betafix
    
    if(is.null(input$re_dhglm_a_check_betafix) || !input$re_dhglm_a_check_betafix)
        return()
    
    textAreaInput("re_dhglm_a_betafix", label = NULL, value = "0")
})

output$re_dhglm_a_check_correlation <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput(
        "re_dhglm_a_check_correlation",
        "Correlation Structure",
        value = FALSE
    )
})

output$re_dhglm_corr_structure_1<-renderUI({
  names=c("1:Female","2:Female","3:Female","1:Male","2:Male","3:Male")
  selectInput("re_dhglm_corr_structure_1","Shared ",choices =as.list(names), multiple = TRUE)
})

output$re_dhglm_corr_variable<-renderUI({
  names=names(data2)
  selectInput("re_dhglm_corr_variable","Variable for Correlation Structure",choices =as.list(names), multiple = FALSE)
})

output$re_dhglm_corr_structure_2<-renderUI({
  result=re_dhglm_output_mean_1()
  names=c("1:Female","2:Female","3:Female","1:Male","2:Male","3:Male")
  names1=paste0(names,"")
  names2=names
  if (!is.null(result)==TRUE) {
    txt1=strsplit(re_dhglm_Model_Sel_1,"\n")[[1]]
    length1=length(txt1)
    if(length1>=1) {
      for (i in 1:length1) {
        txt2=strsplit(txt1[[i]],",")[[1]]
        length2=length(txt2)
        for (j in 2:length2) {
          names2<-gsub(txt2[j],"",names1)
          names1<-paste0(names2,"")
        }
      }
    }
  } else {
    names2=names
  }
  selectInput("re_dhglm_corr_structure_2","Saturated",choices =as.list(names2), multiple = TRUE)
})

re_dhglm_output_mean_1 <-eventReactive(input$re_dhglm_corr_structure_1_run, {
  select1=paste(input$re_dhglm_corr_structure_1, collapse=",")
  if (!is.null(select1)) {
    if (is.null(re_dhglm_Model_Sel_1)) re_dhglm_Model_Sel_1 <<- select1
    else re_dhglm_Model_Sel_1 <<- paste0(re_dhglm_Model_Sel_1,"\n",select1)
  } else {
    re_dhglm_Model_Sel_1 <<- re_dhglm_Model_Sel_1
  }
  return(re_dhglm_Model_Sel_1)
})

re_dhglm_output_mean_2 <-eventReactive(input$re_dhglm_corr_structure_2_run, {
  select2=paste(input$re_dhglm_corr_structure_2, collapse=",")
  if (is.null(re_dhglm_Model_Sel_2)) re_dhglm_Model_Sel_2 <<- select2
  else re_dhglm_Model_Sel_2 <<- paste0(re_dhglm_Model_Sel_2,"\n",select2)
  return(re_dhglm_Model_Sel_2)
})

output$re_dhglm_corr_structure_11 <- renderUI({
  Height = paste0(80, "px")
  value=re_dhglm_output_mean_1()
  textAreaInput("re_dhglm_corr_structure_11", "", value = value, height = Height)
})

output$re_dhglm_corr_structure_21 <- renderUI({
  Height = paste0(80, "px")
  value=re_dhglm_output_mean_2()
  textAreaInput("re_dhglm_corr_structure_21", "", value = value, height = Height)
})


# output$re_dhglm_a_check_orthogonal <- renderUI({
#     input$file
#     input$resetData
#     
#     selectInput(
#         "re_dhglm_m_check_orthogonal",
#         "SE with orthogonality or not",
#         choice = c(TRUE, FALSE), 
#         selected = FALSE,
#         multiple = FALSE
#     )
# })

output$re_dhglm_a_mord <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    selectInput(
        "re_dhglm_a_mord",
        "Order for Mean Model",
        choice = c(1, 0),
        multiple = FALSE
    )
})

output$re_dhglm_a_dord <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    selectInput(
        "re_dhglm_a_dord",
        "Order for Dispersion Model",
        choice = c(1, 2),
        multiple = FALSE
    )
})

# DHGLM Results ####

output$re_dhglm_r_mainsummary <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_dhglm_results()$modelDesc))
        return()
  
    re_dhglm_results()$modelDesc
}, rownames = TRUE, bordered = TRUE, caption = "Model Description", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_dhglm_m_coeff <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_dhglm_results()$beta_coeff))
        return()
  
    re_dhglm_results()$beta_coeff
}, rownames = TRUE, bordered = TRUE, caption = "Estimate from Mean Model", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_dhglm_l_coeff <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_dhglm_results()$lambda_coeff))
        return()
  
    if (sum(abs(re_dhglm_results()$special))<0.00001) re_dhglm_results()$lambda_coeff
    else re_dhglm_results()$special
}, rownames = TRUE, bordered = TRUE, caption = "Estimate for log(Lambda)", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_dhglm_l_taucoeff <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_dhglm_results()$alpha_coeff))
        return()
  
    re_dhglm_results()$alpha_coeff
}, rownames = TRUE, bordered = TRUE, caption = "Estimate for log(Tau)", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_dhglm_p_coeff <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_dhglm_results()$phi_coeff))
        return()
  
    re_dhglm_results()$phi_coeff
}, rownames = TRUE, bordered = TRUE, caption = "Estimate from Dispersion Model", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_dhglm_p_alphacoeff <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_dhglm_results()$tau_coeff))
        return()
  
    re_dhglm_results()$tau_coeff
}, rownames = TRUE, bordered = TRUE, caption = "Estimate for log(Alpha)", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_dhglm_r_likelihood <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_dhglm_results()$likeli_coeff))
        return()
  
    re_dhglm_results()$likeli_coeff
}, rownames = FALSE, bordered = TRUE, caption = "Likelihood", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_dhglm_r_comparisonmodel1<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_dhglm_results()$option$Comparisonmodel1))
        return()
  
    re_dhglm_results()$option$Comparisonmodel1
}, rownames = TRUE, bordered = TRUE, caption = "Mean", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_dhglm_r_comparisonmodel2<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_dhglm_results()$option$Comparisonmodel2))
        return()
  
    re_dhglm_results()$option$Comparisonmodel2
}, rownames = TRUE, bordered = TRUE, caption = "Phi & Lambda", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_dhglm_r_comparisonmodel3<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_dhglm_results()$option$Comparisonmodel3))
        return()
  
    re_dhglm_results()$option$Comparisonmodel3
}, rownames = TRUE, bordered = TRUE, caption = "Likelihood", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_dhglm_r_rcodes <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_dhglm_results()$option$Rcodes))
        return()
  
    re_dhglm_results()$option$Rcodes
}, rownames = TRUE, bordered = TRUE, caption = "R Codes", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
# DHGLM Plots ####
# ¦§ Mean #####

output$re_dhglm_m_selectplot1<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if (g_resetresult == FALSE || is.null(re_dhglm_results()))
        return()
    
    meanNumberRand = length(strsplit(re_dhglm_results()$modelDesc[1,1], "|", fixed = TRUE)[[1]]) - 1
    
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
        inputId = "re_dhglm_m_selectplot1",
        label = "Select Plot", 
        choiceNames = choiceName,
        choiceValues = choiceValue, 
        selected = "mean",
        direction = "vertical"
    )
})

output$re_dhglm_m_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    local({
        output$re_dhglm_m_downloadplot1 <- downloadPlot({
            if(is.null(input$re_dhglm_m_selectplot1))
                plotType = "mean"
            else
                plotType = input$re_dhglm_m_selectplot1
            res<-re_dhglm_results()
            spatial = ""
            
            if(!is.null(res$model_mu[[13]]) && res$model_mu[[13]]=="Matern")
                spatial = res$model_mu[[13]]
            else if(!is.null(res$option$MeanModel[[13]]) && (res$option$MeanModel[[13]] == "IAR" || res$option$MeanModel[[13]] == "MRF"))
                spatial = res$option$MeanModel[[13]]
            
            ggplotdhglm(res, type = plotType, spatial)
        })
    })
        
    if(g_resetresult == FALSE || is.null(re_dhglm_results()))
        return()
    
    if(is.null(input$re_dhglm_m_selectplot1))
        plotType = "mean"
    else
        plotType = input$re_dhglm_m_selectplot1
    
    res<-re_dhglm_results()
    spatial = ""
    
    if(!is.null(res$model_mu[[13]]) && res$model_mu[[13]]=="Matern")
        spatial = res$model_mu[[13]]
    else if(!is.null(res$option$MeanModel[[13]]) && (res$option$MeanModel[[13]] == "IAR" || res$option$MeanModel[[13]] == "MRF"))
        spatial = res$option$MeanModel[[13]]
    
    ggplotdhglm(res, type = plotType, spatial)
})

output$re_dhglm_m_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run

    if(g_resetresult == FALSE || is.null(re_dhglm_results()))
        return()
    
    div(
        h4("Model Checking Plots for Mean"),
        div(
            style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
            withSpinner(
                plotOutput("re_dhglm_m_plot1", height="600px"),
                type = 1,
                color = "#2c3e50",
                size = 1.2
            ),
            div(
                style = "position: absolute; left: 0.5em; bottom: 0.5em;",
                dropdown(
                    uiOutput("re_dhglm_m_selectplot1"),
                    size = "xs",
                    icon = icon("chart-bar", class = "opt"), 
                    up = TRUE
                )
            ),
            div(
                style = "position: absolute; left:4em; bottom: 0.5em;",
                dropdown(
                    downloadButton(outputId = "re_dhglm_m_downloadplot1", label = "Download Plot"),
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

output$re_dhglm_p_selectplot1<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if(g_resetresult == FALSE || is.null(re_dhglm_results()$resid_phi))
        return()

    choiceName = c("Residuals")
    choiceValue = c("phi")
    
    if(!is.null(re_dhglm_results()$phi_sv_h)) {
        choiceName = c(choiceName, "Random Effect")
        choiceValue = c(choiceValue, "phiv")
    }
    if(!is.null(re_dhglm_results()$phi_v_h)){
        choiceName = c(choiceName, "95% CI for Random Effect")
        choiceValue = c(choiceValue, "phici")
    }
    
    radioGroupButtons(
        inputId = "re_dhglm_p_selectplot1",
        label = "Select Plot", 
        choiceNames = choiceName,
        choiceValues = choiceValue, 
        selected = "phi",
        direction = "vertical"
    )
})

output$re_dhglm_p_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    local({
        output$re_dhglm_p_downloadplot1 <- downloadPlot({
            if(is.null(input$re_dhglm_p_selectplot1))
                plotType = "phi"
            else
                plotType = input$re_dhglm_p_selectplot1
            
            res<-re_dhglm_results()
            
            ggplotdhglm(res, type = plotType)
        })  
    })
        
    if(g_resetresult == FALSE || is.null(re_dhglm_results()$resid_phi))
        return()
    
    if(is.null(input$re_dhglm_p_selectplot1))
        plotType = "phi"
    else
        plotType = input$re_dhglm_p_selectplot1
    
    res<-re_dhglm_results()
    
    ggplotdhglm(res, type = plotType)
})

output$re_dhglm_p_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run

    if(g_resetresult == FALSE || is.null(re_dhglm_results()$resid_phi))
        return()
    
    div(
        h4("Model Checking Plots for Phi"),
        div(
            style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
            
            withSpinner(
                plotOutput("re_dhglm_p_plot1", height="600px"),
                type = 1,
                color = "#2c3e50",
                size = 1.2
            ),
            div(
                style = "position: absolute; left: 0.5em; bottom: 0.5em;",
                dropdown(
                    uiOutput("re_dhglm_p_selectplot1"),
                    size = "xs",
                    icon = icon("chart-bar", class = "opt"), 
                    up = TRUE
                )
            ),
            div(
                style = "position: absolute; left:4em; bottom: 0.5em;",
                dropdown(
                    downloadButton(outputId = "re_dhglm_p_downloadplot1", label = "Download Plot"),
                    size = "xs",
                    icon = icon("download", class = "opt"),
                    up = TRUE
                )
            )
        ),
        br()
    )
})


# ¦¦ Lambda #####

output$re_dhglm_l_selectplot1<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if(g_resetresult == FALSE || is.null(re_dhglm_results()$resid_lambda))
        return()
    
    choiceName = c("Residuals")
    choiceValue = c("lambda")
    
    if (!is.null(re_dhglm_results()$lambda_sv_h)) {
        choiceName = c(choiceName, "Random Effect")
        choiceValue = c(choiceValue, "lambdav")
    }
    if (!is.null(re_dhglm_results()$lambda_v_h)) {
        choiceName = c(choiceName, "95% CI for Random Effect")
        choiceValue = c(choiceValue, "lambdaci")
    }
    
    radioGroupButtons(
        inputId = "re_dhglm_l_selectplot1",
        label = "Select Plot", 
        choiceNames = choiceName,
        choiceValues = choiceValue, 
        selected = "lambda",
        direction = "vertical"
    )
})

output$re_dhglm_l_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    local({
        output$re_dhglm_l_downloadplot1 <- downloadPlot({
            if(is.null(input$re_dhglm_l_selectplot1))
                plotType = "lambda"
            else
                plotType = input$re_dhglm_l_selectplot1
            
            res<-re_dhglm_results()
            
            ggplotdhglm(res, type = plotType)
        })
    })
    
    if(g_resetresult == FALSE || is.null(re_dhglm_results()$resid_lambda))
        return()
    
    if(is.null(input$re_dhglm_l_selectplot1))
        plotType = "lambda"
    else
        plotType = input$re_dhglm_l_selectplot1
    
    res<-re_dhglm_results()
    
    ggplotdhglm(res, type = plotType)
})

output$re_dhglm_l_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run

    if(g_resetresult == FALSE || is.null(re_dhglm_results()$resid_lambda))
        return()
    
    div(
        h4("Model Checking Plots for Lambda"),
        div(
            style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
            
            withSpinner(
                plotOutput("re_dhglm_l_plot1", height="600px"),
                type = 1,
                color = "#2c3e50",
                size = 1.2
            ),
            div(
                style = "position: absolute; left: 0.5em; bottom: 0.5em;",
                dropdown(
                    uiOutput("re_dhglm_l_selectplot1"),
                    size = "xs",
                    icon = icon("chart-bar", class = "opt"), 
                    up = TRUE
                )
            ),
            div(
                style = "position: absolute; left:4em; bottom: 0.5em;",
                dropdown(
                    downloadButton(outputId = "re_dhglm_l_downloadplot1", label = "Download Plot"),
                    size = "xs",
                    icon = icon("download", class = "opt"),
                    up = TRUE
                )
            )
        ),
        br()
    )
})

# DHGLM Prediction ####

output$re_dhglm_r_check_95mu <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    checkboxInput("re_dhglm_r_check_95mu", "95% Confidence Interval for mu", value = FALSE)
})

output$re_dhglm_r_box_95mu <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if(is.null(input$re_dhglm_r_check_95mu) || !input$re_dhglm_r_check_95mu)
        return()
    
    muLength = nrow(re_dhglm_results()$beta_coeff)
    if(rownames(re_dhglm_results()$beta_coeff)[1] == "(Intercept)")
        muLength = muLength - 1 
    
    textAreaInput("re_dhglm_r_box_95mu", "Type mean Covariates", value = paste(rep(1, muLength), collapse = ","))
})

output$re_dhglm_r_check_95phi <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if(is.null(re_dhglm_results()$modelDesc) || re_dhglm_results()$modelDesc[2,1] == "phi ~ 1"
       || re_dhglm_results()$modelDesc[2,1] == "phi~1")
        return()
    
    checkboxInput("re_dhglm_r_check_95phi", "95% Confidence Interval for phi", value = FALSE)
})

output$re_dhglm_r_box_95phi <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if(is.null(input$re_dhglm_r_check_95phi) || !input$re_dhglm_r_check_95phi)
        return()
    
    phiLength = nrow(re_dhglm_results()$phi_coeff)
    if(rownames(re_dhglm_results()$phi_coeff)[1] == "(Intercept)")
        phiLength = phiLength - 1 
    
    textAreaInput("re_dhglm_r_box_95phi", "Type phi Covariates", value = paste(rep(1, phiLength), collapse = ","))
})

output$re_dhglm_r_prediction <- renderPrint({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if(!is.null(input$re_dhglm_r_check_95mu) && input$re_dhglm_r_check_95mu && !is.null(re_dhglm_results()$beta_coeff)) {
        
        muLength = nrow(re_dhglm_results()$beta_coeff)
        if(rownames(re_dhglm_results()$beta_coeff)[1] == "(Intercept)")
            muLength = muLength - 1 
        
        muCovariates <- as.numeric(strsplit(input$re_dhglm_r_box_95mu, ",", fixed = TRUE)[[1]])
        if(muLength >= length(muCovariates))
            muCovariates <- c(muCovariates, rep(0,muLength - length(muCovariates)))
        else muCovariates <- muCovariates[1:muLength]
        
        if(rownames(re_dhglm_results()$beta_coeff)[1] == "(Intercept)")
            muCovariates<-c(1,muCovariates)
            
        muCovMatrix <- matrix(data = muCovariates, ncol=1)
        rownames(muCovMatrix) <- rownames(re_dhglm_results()$beta_coeff)
        print("Covariates for mean Predictor")
        print(t(muCovMatrix))
        muPrediction <- sum(as.vector(muCovMatrix) * re_dhglm_results()$beta_coeff[,1])
        muSE <- sqrt(t(muCovMatrix) %*% re_dhglm_results()$cov_mu %*% muCovMatrix)
        LL <- muPrediction - 1.96 * muSE
        UL <- muPrediction + 1.96 * muSE
        print(paste0("Lower Limit for Confidence Intervel of mean predictor : ", LL))
        print(paste0("Upper Limit for Confidence Intervel of mean predictor : ", UL))
    }
    
    if(!is.null(input$re_dhglm_r_check_95phi) && input$re_dhglm_r_check_95phi 
       && !is.null(re_dhglm_results()$phi_coeff)) {
        phiLength = nrow(re_dhglm_results()$phi_coeff)
        if(rownames(re_dhglm_results()$phi_coeff)[1] == "(Intercept)")
            phiLength = phiLength - 1
        
        phiCovariates <- as.numeric(strsplit(input$re_dhglm_r_box_95phi, ",", fixed = TRUE)[[1]])
        if(phiLength >= length(phiCovariates))
            phiCovariates <- c(phiCovariates, rep(0,phiLength - length(phiCovariates)))
        else phiCovariates <- phiCovariates[1:phiLength]
        
        if(rownames(re_dhglm_results()$phi_coeff)[1] == "(Intercept)")
            phiCovariates<-c(1,phiCovariates)
        
        phiCovMatrix <- matrix(data = phiCovariates, ncol=1)
        rownames(phiCovMatrix) <- rownames(re_dhglm_results()$phi_coeff)
        print("Covariates for phi Predictor")
        print(t(phiCovMatrix))
        phiPrediction <- sum(as.vector(phiCovMatrix) * re_dhglm_results()$phi_coeff[,1])
        phiSE <- sqrt(t(phiCovMatrix) %*% re_dhglm_results()$cov_phi %*% phiCovMatrix)
        LL <- phiPrediction - 1.96 * phiSE
        UL <- phiPrediction + 1.96 * phiSE
        print(paste0("Lower Limit for Confidence Intervel of phi predictor : ", LL))
        print(paste0("Upper Limit for Confidence Intervel of phi predictor : ", UL))
    }
})


output$re_dhglm_r_check_calculator <-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    checkboxInput("re_dhglm_r_check_calculator", "Calculator for data", value = FALSE)
})

output$re_dhglm_r_box_calculator <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if(is.null(input$re_dhglm_r_check_calculator) || !input$re_dhglm_r_check_calculator)
        return()
    
    textAreaInput("re_dhglm_r_box_calculator", "Type Arithmetic Form (e.g. sum(log(y))", value = paste0("sum(", input$re_dhglm_m_resp, ")"))
})

output$re_dhglm_r_calculator <- renderPrint({
    input$file
    input$resetData
    input$di_option_run
    input$re_dhglm_g_run
    
    if(is.null(input$re_dhglm_r_check_calculator) || !input$re_dhglm_r_check_calculator 
       || is.null(input$re_dhglm_r_box_calculator))
        return()
    
    calcEquation <- input$re_dhglm_r_box_calculator
    
    quo_var <- rlang::parse_expr(quo_name(enquo(calcEquation)))
    calcResults <- dplyr::mutate(data2, new1 = !!quo_var)
    
    print("Calculator Results")
    print(calcResults[1,ncol(calcResults)])
})
