# GLMM Run ####

observeEvent(input$re_glmm_g_run, {
    g_resetresult <<- TRUE
})

re_glmm_results <- eventReactive(input$re_glmm_g_run, {
    # ¦§ Warning & Notify ####
    # if user doesn't input some essential part, we wil return null value and show notification
    if (input$re_glmm_m_resp == "")  {
        showNotification("Please choose response", type="warning")
        return()
    }
    
    if(input$re_glmm_m_check_offset && input$re_glmm_m_offset == "") {
        showNotification("Please choose offset variable", type="warning")
        return()
    }
    
    if(!is.null(input$re_glmm_m_check_binomd) && input$re_glmm_m_check_binomd && input$re_glmm_m_binomd == "") {
        showNotification("Please choose binomial denominator", type="warning")
        return()
    }
    
    # Start Progressbar
    withProgress(message = 'GLMM', style = "notification", value = 0.1, {
        
    # ¦§ Variable Declaration 1 ####
    # variable declaration for which is used widely in re_glmm_results
    dataLength = nrow(data2)
    
    # variable declaration for which is used in making MM(MeanModel), DM(DispersionModel)
    meanLink = input$re_glmm_m_link
    
    meanFormula = input$re_glmm_m_model
    
    meanNumberRandom = length(strsplit(meanFormula, "|", fixed = TRUE)[[1]]) - 1
    
    meanRandFamily = NULL
    
    if(meanNumberRandom > 0)
        for (i in 1:meanNumberRandom) 
            meanRandFamily = c(meanRandFamily, input[[paste0("re_glmm_m_randfamily_", i)]])
    
    meanOffsetVariable = NULL
    
    if(input$re_glmm_m_check_offset)
        meanOffsetVariable = data2[[input$re_glmm_m_offset]]
    
    # ¦§ MeanModel (MM) ####
    
    MM = DHGLMMODELING(
        Model = "mean",
        Link = meanLink,
        LinPred = as.formula(meanFormula),
        RandDist = meanRandFamily,
        Offset = meanOffsetVariable
    )
    
    # ¦§ DispersionModel (DM) ####
    
    DM = DHGLMMODELING(
        Model = "dispersion",
        Link = "log"
    )
    
    incProgress(0.4, detail = paste("Loading..."))
    
    # ¦§ Variable Declaration 2 ####
    # Variable Declaration for dhglmfit
    if (is.null(input$re_glmm_a_mord))
        mord = 1
    else
        mord = as.numeric(input$re_glmm_a_mord)
    
    if (is.null(input$re_glmm_a_dord))
        dord = 1
    else
        dord = as.numeric(input$re_glmm_a_dord)
    
    if (is.null(input$re_glmm_a_dispersion))
        Dmethod = "deviance"
    else
        Dmethod = input$re_glmm_a_dispersion
    
    if (is.null(input$re_glmm_a_reml))
        reml = TRUE
    else 
        reml = FALSE
    
    se_orhogonal_dhglm  = input$re_glmm_a_check_orthgonal #wrong word, but...
    
    BinomialDen = NULL
    if (!is.null(input$re_glmm_m_check_binomd) && input$re_glmm_m_check_binomd) 
        BinomialDen = data2[[input$re_glmm_m_binomd]]
    
    # ¦§ Fitted Model ####
    
    fittedModel <<-
        dhglmfit(
            RespDist = isolate(input$re_glmm_m_dist),
            BinomialDen = BinomialDen,
            DataMain = data2,
            MeanModel = MM,
            DispersionModel = DM,
            Dmethod = Dmethod,
            mord = mord,
            dord = dord,
            REML = reml,
            se_orhogonal_dhglm = se_orhogonal_dhglm
        )
    
    incProgress(0.3, detail = paste("Loading..."))
    
    fittedModel$option = list(
        RespDist = isolate(input$re_glmm_m_dist),
        BinomialDen = BinomialDen,
        DataMain = data2,
        MeanModel = MM,
        DispersionModel = DM,
        Dmethod = Dmethod,
        mord = mord,
        dord = dord,
        REML = reml,
        se_orhogonal_dhglm = se_orhogonal_dhglm
    ) 
        
    # ¦§ Description 1 ####
    
    nMrand <- length(strsplit(input$re_glmm_m_model, "|", fixed = TRUE)[[1]]) - 1
    mRandstring = NULL
    if (nMrand >= 1)
        for (i in 1:nMrand)
            mRandstring = paste0(c(mRandstring, input[[paste0("re_glmm_m_randfamily_", i)]]), collapse = ", ")
    else
        mRandstring = NA
    
    modelDesc = matrix(c(
        input$re_glmm_m_model,
        input$re_glmm_m_link,
        input$re_glmm_m_dist,
        mRandstring
    ),nrow=1, ncol=4)
    
    colnames(modelDesc) <- c("Model", "Link", "Dist", "Rand")
    rownames(modelDesc) <- c("Mean")
    fittedModel$modelDesc <- modelDesc
    
    # ¦§ Likelihood ####
    
    likeli_coeff <- matrix(fittedModel$likeli_coeff, ncol = 5)
    colnames(likeli_coeff) <- c("-2ML", "-2RL", "cAIC", "Scaled Deviance", "df")
    fittedModel$likeli_coeff = likeli_coeff
    
    # ¦§ Description 2 (Exp) ####
    
    if (length(colnames(fittedModel$beta_coeff)) == 3) {
        pValue <- 2 * pnorm(abs(fittedModel$beta_coeff[, 3]), lower.tail = FALSE)
        llValue <- fittedModel$beta_coeff[, 1] - 1.96 * fittedModel$beta_coeff[, 2]
        ulValue <- fittedModel$beta_coeff[, 1] + 1.96 * fittedModel$beta_coeff[, 2]
        
        fittedModel$beta_coeff <- cbind(fittedModel$beta_coeff, pValue, llValue, ulValue)
        colnames(fittedModel$beta_coeff) <- c("Estimate", "Std. Error", "t-value", "p_val", "LL", "UL")
    }
    
    if (!is.null(fittedModel$beta_coeff) && !is.null(input$re_glmm_m_check_exp) && input$re_glmm_m_check_exp) {
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
    
    glmmComparisonmodel1 = NULL
    glmmComparisonmodel2 = NULL
    
    if (!is.null(input$re_glmm_m_check_comparison) && input$re_glmm_m_check_comparison) {
        
        glmmComparisonmodel1 <- matrix(
            c(
                input$re_glmm_m_model,
                input$re_glmm_m_link,
                input$re_glmm_m_dist,
                mRandstring
            ),
            ncol = 4
        )
        
        glmmComparisonmodel2 <- matrix(
            c(
                fittedModel$likeli_coeff[1,1],
                fittedModel$likeli_coeff[1,2],
                fittedModel$likeli_coeff[1,3],
                fittedModel$likeli_coeff[1,5],
                round(fittedModel$n.mean, digits=0),
                round(fittedModel$n.disp, digits=0)
            ),
            ncol = 6
        )
        glmmComparisonmodel1 <- rbind(g_re_glmm_r_comparisonmodel_1, glmmComparisonmodel1)
        glmmComparisonmodel2 <- rbind(g_re_glmm_r_comparisonmodel_2, glmmComparisonmodel2)
        
        colnames(glmmComparisonmodel1) <- c("Model", "Link", "Dist", "Rand")
        colnames(glmmComparisonmodel2) <- c("-2ML", "-2RL", "cAIC", "df", "n.mean", "n.disp")
        
        g_re_glmm_r_comparisonmodel_1 <<- glmmComparisonmodel1
        g_re_glmm_r_comparisonmodel_2 <<- glmmComparisonmodel2
        
        fittedModel$option$Comparisonmodel1 = glmmComparisonmodel1
        fittedModel$option$Comparisonmodel2 = glmmComparisonmodel2
    }
    
    # ¦¦ R Codes ####
    
    if (!is.null(input$re_glmm_m_check_rcodes) && input$re_glmm_m_check_rcodes) {
        Rcode1 <- NULL
        Rcode2 <- NULL
        Rcode3 <- NULL
        offsetRcode = NULL
        if (!is.null(meanOffsetVariable))
            offsetRcode <- paste0(", Offset = data[[\"", input$re_glmm_m_offset, "\"]]")
        mRandRcode = NULL
        if (nMrand >= 1) {
            for (i in 1:nMrand)
                mRandRcode = paste0(c(mRandRcode, input[[paste0("re_glmm_m_randfamily_", i)]]), collapse = "\", \"")
            
            mRandRcode = paste0(", RandDist = c(\"", mRandRcode, "\")")
        }
        
        Rcode1 <- paste0(
            "MM <- DHGLMMODELING(Model = \"mean\", Link = \"",
            input$re_glmm_m_link,
            "\", LinPred = ",
            meanFormula,
            mRandRcode,
            offsetRcode,
            ")"
        )
        
        Rcode2 <- paste0(
            "DM <- DHGLMMODELING(Model = \"dispersion\", Link = \"log\")"
        )
        
        Rcode3 <- paste0(
            "dhglmfit(RespDist = \"", 
            input$re_glmm_m_dist,
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

# GLMM UI (Accordion) ####

output$re_glmm_r_accordion<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    glmmAccordion <- bs_accordion_sidebar(
        id = "re_glmm_r_accordion",
        spec_side = c(width = 3, offset = 0),
        spec_main = c(width = 9, offset = 0)    
    )
    
    # ¦§ Mean ####
    
    glmmAccordion <- glmmAccordion %>%
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
            uiOutput("re_glmm_m_resp"), 
            uiOutput("re_glmm_m_variable"), 
            fluidRow(
                column(
                    8,
                    uiOutput("re_glmm_m_interaction") 
                ),
                column(
                    4,
                    uiOutput("re_glmm_m_interactionappend"),
                    style = "text-align:right; padding:15px"
                )
            ),
            uiOutput("re_glmm_m_rand"),
            uiOutput("re_glmm_m_check_slope"),
            uiOutput("re_glmm_m_check_withoutslope"),
            uiOutput("re_glmm_m_slope"),
            uiOutput("re_glmm_m_check_randinteraction"),
            fluidRow(
                column(
                    8,
                    uiOutput("re_glmm_m_randinteraction")
                ),
                column(
                    4,
                    uiOutput("re_glmm_m_randinteractionappend"),
                    style = "text-align:right; padding:15px"
                )
            ), #new components
            
            uiOutput("re_glmm_m_dist"),
            uiOutput("re_glmm_m_check_binomd"),
            uiOutput("re_glmm_m_binomd"),
            uiOutput("re_glmm_m_link"),
            
            hr(style = "border-color: #2C3E50;"),
            uiOutput("re_glmm_m_randfamily"),
            uiOutput("re_glmm_m_check_nointercept"),
            uiOutput("re_glmm_m_check_offset"),
            uiOutput("re_glmm_m_offset"),
            uiOutput("re_glmm_m_check_exp"),
            
            hr(style = "border-color: #2C3E50;"),
            uiOutput("re_glmm_m_check_comparison"),
            uiOutput("re_glmm_m_check_rcodes")
        )
    )
    
    # ¦¦ Additional Settings ####
    
    glmmAccordion <- glmmAccordion %>%
    bs_append(
        title_side = "Setting",
        content_side = NULL,
        content_main = div(
            h3(strong("Additional Settings")),
            uiOutput("re_glmm_a_dispersion"),
            uiOutput("re_glmm_a_reml"),
            uiOutput("re_glmm_a_mord"),
            uiOutput("re_glmm_a_dord"),
            h5(helpText("Order of Laplace Approximation for Likelihood(mean) and Restricted Likelihood(Dispersion)"))
        )
    )
    
    div(
        glmmAccordion,
        use_bs_tooltip(),
        use_bs_accordion_sidebar() # needs to be at end, for some reason
    )
})

# GLMM Components####
# ¦§ Mean ####

output$re_glmm_m_model<-renderUI({ #algorithm not good
    input$re_glmm_m_resp
    input$re_glmm_m_variable$right
    input$re_glmm_m_rand
    input$re_glmm_m_check_nointercept
    input$re_glmm_m_binomd
    input$re_glmm_m_slope
    
    if(any(is.null(input$re_glmm_m_resp), is.null(input$re_glmm_m_variable)))
        return()
    
    responsePart = input$re_glmm_m_resp
    variablePart = paste0(input$re_glmm_m_variable$right, collapse=" + ")
    randomPart = NULL
    meanModel = NULL
    
    # in DHGLM, even if we use binomial denominator, model equation doesn't change. 
    # if (all(!is.null(input$re_glmm_m_check_binomd), input$re_glmm_m_check_binomd, 
    #         !is.null(input$re_glmm_m_binomd), input$re_glmm_m_binomd != ""))
    #     responsePart = paste0("cbind(", responsePart, ",", input$re_glmm_m_binomd,"-", responsePart, ")")
    
    if(!is.null(input$re_glmm_m_rand) && length(input$re_glmm_m_rand) > 0) {
        if (is.null(input$re_glmm_m_check_withoutslope) || !input$re_glmm_m_check_withoutslope) {
            randomPart = paste0(input$re_glmm_m_rand, collapse = ") + (1|")
            randomPart = paste0(" + (1|", randomPart, ")")
        }
        
        #random slope (a|b) form
        if(all(!is.null(input$re_glmm_m_check_slope), input$re_glmm_m_check_slope, length(input$re_glmm_m_slope)>0 )) 
            for(i in 1:length(input$re_glmm_m_rand)) 
                for(j in 1:length(input$re_glmm_m_slope)) 
                    randomPart = paste0(randomPart," + (", input$re_glmm_m_slope[j], "|", input$re_glmm_m_rand[i], ")")
    }
    
    if(!is.null(input$re_glmm_m_check_nointercept) && !input$re_glmm_m_check_nointercept) {
        if(variablePart == "")
            meanModel = paste(responsePart, "~", "1", randomPart)
        else
            meanModel = paste(responsePart, "~", variablePart, randomPart)
    }
    
    if(!is.null(input$re_glmm_m_check_nointercept) && input$re_glmm_m_check_nointercept) {
        if(variablePart == "")
            meanModel = paste(responsePart, "~", "-1", randomPart)
        else
            meanModel = paste(responsePart, "~", "-1 +", variablePart, randomPart)
    }
    
    textAreaInput("re_glmm_m_model", "Model for Mean", value = meanModel, height = "60px")
})

output$re_glmm_m_resp <- renderUI({
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
            "re_glmm_m_resp",
            "Response Variable",
            choices = as.list(c("", nameValue)), 
            selected = NULL, 
            multiple = FALSE
        ),
        bsTooltip(
            "re_glmm_m_resp", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$re_glmm_m_variable <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$re_glmm_m_interactionappend
    
    nameValue <- c(names(data2), g_re_glmm_m_interaction)
    
    chooserInput("re_glmm_m_variable", "Variable", "Selected", nameValue, c(), size = 15, multiple = TRUE, width = 150)
})

output$re_glmm_m_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue <- names(data2)
    
    selectInput(
        "re_glmm_m_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$re_glmm_m_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    actionButton("re_glmm_m_interactionappend", "Append")
})

observeEvent(input$re_glmm_m_interactionappend, {
    if(length(input$re_glmm_m_interaction) > 1)
        g_re_glmm_m_interaction <<- c(g_re_glmm_m_interaction, paste(input$re_glmm_m_interaction, collapse=":"))
})


output$re_glmm_m_rand <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$re_glmm_m_randinteractionappend
    
    nameValue <- c(names(data2), g_re_glmm_m_randinteraction)
    
    selectInput(
        "re_glmm_m_rand",
        "Random Effects",
        choices = as.list(nameValue), 
        multiple = TRUE
    )
})

output$re_glmm_m_check_slope <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_glmm_m_check_slope", "Random Slope Model", value = FALSE)
})

output$re_glmm_m_check_withoutslope <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if (is.null(input$re_glmm_m_check_slope) || !input$re_glmm_m_check_slope)
        return()
    
    checkboxInput("re_glmm_m_check_withoutslope", "Without Random Intercept", value = FALSE)
})

output$re_glmm_m_slope <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if (is.null(input$re_glmm_m_check_slope) || !input$re_glmm_m_check_slope)
        return()
    
    nameValue = names(data2)
    
    selectInput(
        "re_glmm_m_slope",
        "Random slope",
        choices = as.list(nameValue), 
        multiple = TRUE
    )
})

output$re_glmm_m_check_randinteraction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_glmm_m_check_randinteraction", "Interaction in the Random Effect", value = FALSE)
})


output$re_glmm_m_randinteraction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$re_glmm_m_check_randinteraction) || !input$re_glmm_m_check_randinteraction)
        return()
    
    nameValue <- names(data2)
    
    selectInput(
        "re_glmm_m_randinteraction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$re_glmm_m_randinteractionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_glmm_m_check_randinteraction) || !input$re_glmm_m_check_randinteraction)
        return()
    
    actionButton("re_glmm_m_randinteractionappend", "Append")
})

observeEvent(input$re_glmm_m_randinteractionappend, {
    
    if(length(input$re_glmm_m_randinteraction) > 1)
        g_re_glmm_m_randinteraction <<- c(g_re_glmm_m_randinteraction, paste(input$re_glmm_m_randinteraction, collapse=":"))
})

output$re_glmm_m_dist<-renderUI({
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
        "re_glmm_m_dist",
        "Distribution", 
        choices = dist, 
        multiple = FALSE
    )
})

output$re_glmm_m_check_binomd<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(any(is.null(input$re_glmm_m_dist), input$re_glmm_m_dist != "binomial"))
        return()
    
    checkboxInput("re_glmm_m_check_binomd", "Binomial Denominator", value = FALSE)
})

output$re_glmm_m_binomd<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$re_glmm_m_check_binomd) || !(input$re_glmm_m_check_binomd))
        return()
    
    nameValue = names(select_if(data2, is.numeric))
    
    selectInput(
        "re_glmm_m_binomd",
        "Binomial Denominator",
        choices = as.list(c("", nameValue)),
        multiple = FALSE
    )
})

output$re_glmm_m_link <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_glmm_m_dist))
        return()
    
    selection = "identity"
    
    if (input$re_glmm_m_dist == "binomial")
        selection = "logit"
    else if (input$re_glmm_m_dist == "poisson")
        selection = "log"
    else if (input$re_glmm_m_dist == "gamma")
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
        "re_glmm_m_link",
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
    
    if(is.null(input$re_glmm_m_model) || is.null(input$re_glmm_m_rand) || length(input$re_glmm_m_rand) == 0) {
        output$re_glmm_m_randfamily <- renderUI({
            return()
        })
    }
    
    if(length(input$re_glmm_m_rand) > 0) {
        numberRandom <- length(strsplit(input$re_glmm_m_model, "|", fixed = TRUE)[[1]]) - 1
        for (i in 1:numberRandom) {
            output$re_glmm_m_randfamily <- renderUI({
                re_glmm_m_randfamily_list <- lapply(1:numberRandom, function(i) {
                    selectInput(
                        paste0("re_glmm_m_randfamily_", i),
                        paste("Distribution for Random effects", i),
                        choices = c("gaussian")
                    )
                })
                do.call(tagList, re_glmm_m_randfamily_list)
            })
        }
    }
    else {
        output$re_glmm_m_randfamily <- renderUI({
            return()
        })
    }
})

output$re_glmm_m_check_nointercept<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_glmm_m_check_nointercept", "No Intercept Model", value = FALSE)
})

output$re_glmm_m_check_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_glmm_m_check_offset", "Offset Variable", value = FALSE)
})

output$re_glmm_m_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue = names(select_if(data2, is.numeric))
    
    if(is.null(input$re_glmm_m_check_offset) || !input$re_glmm_m_check_offset)
        return()
    
    div(
        selectInput(
            "re_glmm_m_offset",
            "Offset Variable",
            choices = as.list(c("", nameValue)),
            selected = NULL,
            multiple = FALSE
        ),
        bsTooltip(
            "re_glmm_m_offset", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$re_glmm_m_check_exp <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_glmm_m_check_exp", "Exponential scale", value = FALSE)
})

output$re_glmm_m_check_comparison<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_glmm_m_check_comparison", "model comparison", value = FALSE)
})

output$re_glmm_m_check_rcodes<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_glmm_m_check_rcodes", "R Codes", value = FALSE)
})

# ¦¦ Additional Settings ####

output$re_glmm_a_dispersion<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    dispersionmethods = c(
        "deviance" = "deviance", 
        "Pearson" = "Pearson"
    )
    
    selectInput(
        "re_glmm_a_dispersion",
        "Method of Fitting Dispersion Model", 
        choices = dispersionmethods, 
        multiple = FALSE
    )
})

output$re_glmm_a_reml <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    selectInput(
        "re_glmm_a_reml",
        "REML or ML",
        choice = c("REML" = TRUE, "ML" = FALSE),
        multiple = FALSE
    )
})

output$re_glmm_m_check_orthogonal <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    selectInput(
        "re_glmm_m_check_orthogonal",
        "SE with orthogonality or not",
        choice = c(TRUE, FALSE), 
        selected = FALSE,
        multiple = FALSE
    )
})

output$re_glmm_a_mord <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    selectInput(
        "re_glmm_a_mord",
        "Order for Mean Model",
        choice = c(1, 0),
        multiple = FALSE
    )
})

output$re_glmm_a_dord <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    selectInput(
        "re_glmm_a_dord",
        "Order for Dispersion Model",
        choice = c(1, 2),
        multiple = FALSE
    )
})

# GLMM Results ####

output$re_glmm_r_mainsummary <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_glmm_g_run
    
    if (g_resetresult == FALSE || is.null(re_glmm_results()$modelDesc))
        return()
    
    re_glmm_results()$modelDesc
}, rownames = TRUE, bordered = TRUE, caption = "Model Description", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_glmm_m_coeff <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_glmm_g_run
    
    if (g_resetresult == FALSE || is.null(re_glmm_results()$beta_coeff))
        return()
    
    re_glmm_results()$beta_coeff
}, rownames = TRUE, bordered = TRUE, caption = "Estimate from Mean Model", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_glmm_l_coeff <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_glmm_g_run
    
    if (g_resetresult == FALSE || is.null(re_glmm_results()$lambda_coeff))
        return()
    
    re_glmm_results()$lambda_coeff
}, rownames = TRUE, bordered = TRUE, caption = "Estimate for log(Lambda)", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_glmm_p_coeff <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_glmm_g_run
    
    if (g_resetresult == FALSE || is.null(re_glmm_results()$phi_coeff))
        return()
    
    re_glmm_results()$phi_coeff
}, rownames = TRUE, bordered = TRUE, caption = "Estimate from Dispersion Model", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_glmm_r_likelihood <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_glmm_g_run
    
    if (g_resetresult == FALSE || is.null(re_glmm_results()$likeli_coeff))
        return()
    
    re_glmm_results()$likeli_coeff
}, rownames = FALSE, bordered = TRUE, caption = "Likelihood", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_glmm_r_comparisonmodel1<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_glmm_g_run
    
    if (g_resetresult == FALSE || is.null(re_glmm_results()$option$Comparisonmodel1))
        return()
    
    re_glmm_results()$option$Comparisonmodel1
}, rownames = TRUE, bordered = TRUE, caption = "Comparison of Mean Model", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_glmm_r_comparisonmodel2<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_glmm_g_run
    
    if (g_resetresult == FALSE || is.null(re_glmm_results()$option$Comparisonmodel2))
        return()
    
    re_glmm_results()$option$Comparisonmodel2
}, rownames = TRUE, bordered = TRUE, caption = "Comparison of Likelihood", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$re_glmm_r_rcodes <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$re_glmm_g_run
    
    if (g_resetresult == FALSE || is.null(re_glmm_results()$option$Rcodes))
        return()
  
    re_glmm_results()$option$Rcodes
}, rownames = TRUE, bordered = TRUE, caption = "R Codes", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

# GLMM Plot ####
# ¦¦ Mean #####

output$re_glmm_m_selectplot1<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_glmm_g_run
    
    if (g_resetresult == FALSE || is.null(re_glmm_results()))
        return()
    
    meanNumberRand = length(strsplit(re_glmm_results()$modelDesc[1,1], "|", fixed = TRUE)[[1]]) - 1
    
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
        inputId = "re_glmm_m_selectplot1",
        label = "Select Plot", 
        choiceNames = choiceName,
        choiceValues = choiceValue, 
        selected = "mean",
        direction = "vertical"
    )
})

output$re_glmm_m_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$re_glmm_g_run
    
    if(g_resetresult == FALSE || is.null(re_glmm_results()))
        return()
    
    local({
        output$re_glmm_m_downloadplot1 <- downloadPlot({
            if(is.null(input$re_glmm_m_selectplot1))
                plotType = "mean"
            else
                plotType = input$re_glmm_m_selectplot1
            res<-re_glmm_results()
            spatial = ""
            
            ggplotdhglm(res, type = plotType, spatial)
        })
    })
    
    if(is.null(input$re_glmm_m_selectplot1))
        plotType = "mean"
    else
        plotType = input$re_glmm_m_selectplot1
    
    res<-re_glmm_results()
    spatial = ""
    
    ggplotdhglm(res, type = plotType, spatial)
})

output$re_glmm_m_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_glmm_g_run

    if(g_resetresult == FALSE || is.null(re_glmm_results()))
        return()
    
    div(
        h4("Model Checking Plots for Mean"),
        div(
            style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
            withSpinner(
                plotOutput("re_glmm_m_plot1", height="600px"),
                type = 1,
                color = "#2c3e50",
                size = 1.2
            ),
            div(
                style = "position: absolute; left: 0.5em; bottom: 0.5em;",
                dropdown(
                    uiOutput("re_glmm_m_selectplot1"),
                    size = "xs",
                    icon = icon("chart-bar", class = "opt"), 
                    up = TRUE
                )
            ),
            div(
                style = "position: absolute; left:4em; bottom: 0.5em;",
                dropdown(
                    downloadButton(outputId = "re_glmm_m_downloadplot1", label = "Download Plot"),
                    size = "xs",
                    icon = icon("download", class = "opt"),
                    up = TRUE
                )
            )
        ),
        br()
    )
})

# GLMM Prediction ####

output$re_glmm_r_check_95mu <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_glmm_g_run
    
    checkboxInput("re_glmm_r_check_95mu", "95% Confidence Interval for mu", value = FALSE)
})

output$re_glmm_r_box_95mu <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_glmm_g_run
    
    if(is.null(input$re_glmm_r_check_95mu) || !input$re_glmm_r_check_95mu)
        return()
    
    muLength = nrow(re_glmm_results()$beta_coeff)
    if(rownames(re_glmm_results()$beta_coeff)[1] == "(Intercept)")
        muLength = muLength - 1 
    
    textAreaInput("re_glmm_r_box_95mu", "Type mean Covariates", value = paste(rep(1, muLength), collapse = ","))
})

output$re_glmm_r_prediction <- renderPrint({
    input$file
    input$resetData
    input$di_option_run
    input$re_glmm_g_run
    
    if(!is.null(input$re_glmm_r_check_95mu) && input$re_glmm_r_check_95mu && !is.null(re_glmm_results()$beta_coeff)) {
        
        muLength = nrow(re_glmm_results()$beta_coeff)
        if(rownames(re_glmm_results()$beta_coeff)[1] == "(Intercept)")
            muLength = muLength - 1 
        
        muCovariates <- as.numeric(strsplit(input$re_glmm_r_box_95mu, ",", fixed = TRUE)[[1]])
        if(muLength >= length(muCovariates))
            muCovariates <- c(muCovariates, rep(0,muLength - length(muCovariates)))
        else muCovariates <- muCovariates[1:muLength]
        
        if(rownames(re_glmm_results()$beta_coeff)[1] == "(Intercept)")
            muCovariates<-c(1,muCovariates)
        
        muCovMatrix <- matrix(data = muCovariates, ncol=1)
        rownames(muCovMatrix) <- rownames(re_glmm_results()$beta_coeff)
        print("Covariates for mean Predictor")
        print(t(muCovMatrix))
        muPrediction <- sum(as.vector(muCovMatrix) * re_glmm_results()$beta_coeff[,1])
        muSE <- sqrt(t(muCovMatrix) %*% re_glmm_results()$cov_mu %*% muCovMatrix)
        LL <- muPrediction - 1.96 * muSE
        UL <- muPrediction + 1.96 * muSE
        print(paste0("Lower Limit for Confidence Intervel of mean predictor : ", LL))
        print(paste0("Upper Limit for Confidence Intervel of mean predictor : ", UL))
    }
})


output$re_glmm_r_check_calculator <-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_glmm_g_run
    
    checkboxInput("re_glmm_r_check_calculator", "Calculator for data", value = FALSE)
})

output$re_glmm_r_box_calculator <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_glmm_g_run
    
    if(is.null(input$re_glmm_r_check_calculator) || !input$re_glmm_r_check_calculator)
        return()
    
    textAreaInput("re_glmm_r_box_calculator", "Type Arithmetic Form (e.g. sum(log(y))", value = paste0("sum(", input$re_glmm_m_resp, ")"))
})

output$re_glmm_r_calculator <- renderPrint({
    input$file
    input$resetData
    input$di_option_run
    input$re_glmm_g_run
    
    if(is.null(input$re_glmm_r_check_calculator) || !input$re_glmm_r_check_calculator 
       || is.null(input$re_glmm_r_box_calculator))
        return()
    
    calcEquation <- input$re_glmm_r_box_calculator
    
    quo_var <- rlang::parse_expr(quo_name(enquo(calcEquation)))
    calcResults <- dplyr::mutate(data2, new1 = !!quo_var)
    
    print("Calculator Results")
    print(calcResults[1,ncol(calcResults)])
})
