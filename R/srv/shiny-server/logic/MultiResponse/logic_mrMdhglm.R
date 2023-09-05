# MDHGLM Components ####

output$mr_mdhglm_g_respslider <- renderUI({
    sliderInput(
        inputId = sprintf("mr_mdhglm_g_resp_count"),
        label = "Number of Response Variables",
        value = 2,
        min = 2,
        max = 6
    )
})

output$mr_mdhglm_g_corr_structure<-renderUI({
    selectInput(
        "mr_mdhglm_g_corr_structure",
        "Correlation Structure",
        choices = c("correlated", "independent", "selection", "shared"),
        selected = "correlated",
        multiple = FALSE
    )
})

mr_mdhglm_g_resp_value <- reactive({
    if (is.null(input$mr_mdhglm_g_resp_count))
        return(2)
    else
        return(input$mr_mdhglm_g_resp_count)
})


# ¦§ Mean, Phi, Lambda ####
source('logic/MultiResponse/MDHGLM/logic_mrMdhglm_m.R', local = T)
source('logic/MultiResponse/MDHGLM/logic_mrMdhglm_p.R', local = T)
source('logic/MultiResponse/MDHGLM/logic_mrMdhglm_l.R', local = T)

# ¦¦ Additional Settings ####
output$mr_mdhglm_a_mord<-renderUI({
    selectInput(
        "mr_mdhglm_a_mord",
        "Order for Mean Model",
        choices=c("0"="0","1"="1"),
        multiple=FALSE
    )
})

output$mr_mdhglm_a_dord<-renderUI({
    selectInput(
        "mr_mdhglm_a_dord",
        "Order for Dispersion Model",
        choices=c("1"="1","2"="2"),
        multiple=FALSE
    )
})

output$mr_mdhglm_a_order<-renderUI({
    selectInput(
        "mr_mdhglm_a_order",
        "Order for Laplace's Approximation",
        choices=c("1"="1","2"="2"),
        multiple=FALSE
    )
})

output$mr_mdhglm_a_reml<-renderUI({
    selectInput(
        "mr_mdhglm_a_reml",
        "REML or ML",
        choice = c("REML" = TRUE,"ML" = FALSE),
        multiple=FALSE
    )
})

# MDHGLM Run ####

observeEvent(input$mr_mdhglm_g_run, {
    g_resetresult <<- TRUE
})

mr_mdhglm_results <- eventReactive(input$mr_mdhglm_g_run, {
    # Start Progressbar  
    withProgress(message = 'MDHGLM', style = "notification", value = 0.1, {
    
    removeTab(inputId = "mr_mdhglm_resulttabset", target = "Model Checking plots")
    removeTab(inputId = "mr_mdhglm_resulttabset", target = "Comparison Model")
    
    # ¦§ Variable Declaration 1 ####
    MM_list <- vector(mode = "list", length = input$mr_mdhglm_g_resp_count)
    DM_list <- vector(mode = "list", length = input$mr_mdhglm_g_resp_count)
    Rcodes <- vector(mode = "list", length = input$mr_mdhglm_g_resp_count)
    f <- data2
    
    for(k in 1:input$mr_mdhglm_g_resp_count) {
        m_form <- input[[paste0("mr_mdhglm_m_model_", k)]]
        p_form <- NULL
        l_form <- NULL

        if (!is.null(input[[paste0("mr_mdhglm_p_check_phi_", k)]]) && input[[paste0("mr_mdhglm_p_check_phi_", k)]])
            p_form <- input[[paste0("mr_mdhglm_p_model_", k)]]
        else
            p_form <- phi~1

        if (!is.null(input[[paste0("mr_mdhglm_l_check_lambda_", k)]]) && input[[paste0("mr_mdhglm_l_check_lambda_", k)]])
            l_form <- input[[paste0("mr_mdhglm_l_model_", k)]]
        else
            l_form <- lambda~1

        res <- m_form
        m_RandDistM <- NULL
        p_RandDistM <- NULL
        l_RandDistM <- NULL
        m_nRand <- length(input[[paste0("mr_mdhglm_m_rand_", k)]])
        p_nRand <- length(input[[paste0("mr_mdhglm_p_rand_", k)]])
        l_nRand <- length(input[[paste0("mr_mdhglm_l_rand_", k)]])

        if (m_nRand > 0) {
            for (i in 1:m_nRand) {
                m_RandDistM <- c(m_RandDistM,input[[paste0("mr_mdhglm_m_rf_",i,"_", k)]])
            }
        }
        if (!is.null(input[[paste0("mr_mdhglm_p_check_phi_", k)]]) && input[[paste0("mr_mdhglm_p_check_phi_", k)]]) {
            if (p_nRand > 0) {
                for (i in 1:p_nRand) {
                    p_RandDistM <- c(p_RandDistM,input[[paste0("mr_mdhglm_p_rf_",i,"_", k)]])
                }
            }
        }
        LinPredRandVariance=NULL
        if (!is.null(input[[paste0("mr_mdhglm_l_check_lambda_", k)]]) && input[[paste0("mr_mdhglm_l_check_lambda_", k)]]) {
            if (l_nRand > 0) {
                for (i in 1:l_nRand) {
                    l_RandDistM <- c(l_RandDistM,input[[paste0("mr_mdhglm_l_rf_",i,"_", k)]])
                }
            }
            LinPredRandVariance=l_form
        }

        MM = NULL
        DM = NULL
        
        # ¦§ MeanModel (MM) ####
        
        MM <-
            DHGLMMODELING(
                Model = "mean",
                Link = isolate(input[[paste0("mr_mdhglm_m_link_", k)]]),
                LinPred = as.formula(m_form),
                RandDist = m_RandDistM,
                LinPredRandVariance = as.formula(LinPredRandVariance),
                RandDistRandVariance = l_RandDistM,
                LinkRandVariance = "log"
            )

        # ¦§ DispersionModel (DM) ####
        
        DM <-
            DHGLMMODELING(
                Model = "dispersion",
                Link = "log",
                LinPred = as.formula(p_form),
                RandDist = p_RandDistM
            )
        if (p_form=="phi ~ 1") {
            DM <- DHGLMMODELING(Model = "dispersion", Link = "log")
        }
        if (l_form=="lambda ~ 1") {
            MM <-
                DHGLMMODELING(
                  Model = "mean",
                  Link = isolate(input[[paste0("mr_mdhglm_m_link_", k)]]),
                  LinPred = as.formula(m_form),
                  RandDist = m_RandDistM
                )
        }

        if (is.null(DM[[4]])) DM[[4]]<-"gaussian"
        
        # ¦§ R Codes ####
        
        if (input[[paste0('mr_mdhglm_m_check_rcodes_', k)]]) {
            Rcode1 <- NULL
            Rcode2 <- NULL
            Rcode3 <- NULL
            if (l_form == "lambda ~ 1") {
                Rcode1 <- paste0(
                    "MM <- DHGLMMODELING(Model = \"mean\", Link = \"",
                    input[[paste0("mr_mdhglm_m_link_", k)]],
                    "\", LinPred = ",
                    m_form,
                    ", RandDist = \"",
                    m_RandDistM,
                    ")"
                )
            } else {
                Rcode1 <- paste0(
                    "MM <- DHGLMMODELING(Model = \"mean\", Link = \"",
                    input[[paste0("mr_mdhglm_m_link_", k)]],
                    "\", LinPred = ",
                    m_form,
                    ", RandDist = \"",
                    m_RandDistM,
                    ", LinPredRandVariance = ",
                    LinPredRandVariance,
                    ", RandDistRandVariance = ",
                    l_RandDistM,
                    ", LinkRandVariance = \"log\")"
                )
            }
            
            if (p_form == "phi ~ 1") {
                Rcode2 <- "DM <- DHGLMMODELING(Model = \"dispersion\", Link = \"log\")"
            } else {
                Rcode2 <- paste0(
                    "DM <- DHGLMMODELING(Model = \"dispersion\", Link = \"log\", LinPred = ",
                    p_form,
                    " RandDist = ",
                    p_RandDistM,
                    "\")"
                )
            }
            
            Rcode3 <- paste0(
                "jointfit(RespDist = \"", 
                input[[paste0("mr_mdhglm_m_dist_", k)]],
                "\", DataMain = data, structure = \"",
                input$mr_mdhglm_g_corr_structure,
                "\", MeanModel = MM_list, DispersionModel = DM_list)"
            )
            
            Rcodes[[k]] <- matrix(c(Rcode1, Rcode2, Rcode3), ncol = 1)
            colnames(Rcodes[[k]]) <- "Call"
            rownames(Rcodes[[k]]) <- c("MM", "DM", "dhglmfit")
        }
        
        MM_list[[k]] <- MM
        DM_list[[k]] <- DM
    }
    
    # Increase Progressbar
    incProgress(0.4, detail = paste("Loading..."))
    input_datamain <- vector(mode = "list", length = input$mr_mdhglm_g_resp_count)
    input_respdist = NULL
    input_factor = NULL

    for(k in 1:input$mr_mdhglm_g_resp_count) {
        if (input[['mr_mdhglm_g_corr_structure']] == "factor") {
            input_factor <- as.numeric(c(input_factor, input[[paste0("mr_mdhglm_m_factor_", k)]]))
        }
        input_respdist <- c(input_respdist, input[[paste0("mr_mdhglm_m_dist_", k)]])
        input_datamain[[k]] <- f
    }
    
    # ¦§ Variable Declaration 2 ####
    
    mord = NULL
    dord = NULL
    orderLaplace = NULL
    reml = NULL
    if (is.null(input$mr_mdhglm_a_mord))
        mord = 0
    else
        mord = as.numeric(input$re_dhglm_a_mord)
    
    if (is.null(input$mr_mdhglm_a_dord))
        dord = 1
    else
        dord = as.numeric(input$mr_mdhglm_a_dord)
    
    if (is.null(input$mr_mdhglm_a_order))
        orderLaplace = 1
    else
        orderLaplace = as.numeric(input$mr_mdhglm_a_order)
    
    if (is.null(input$mr_mdhglm_a_reml))
        reml = TRUE
    else 
        reml = FALSE
    
    corrStructure = input$mr_mdhglm_g_corr_structure
    
    # ¦§ Fitted Model ####
    
    fittedModel <- jointfit(
        RespDist = input_respdist,
        DataMain = input_datamain,
        structure = corrStructure,
        MeanModel = MM_list,
        DispersionModel = DM_list,
        mord = mord,
        dord = dord,
        order = orderLaplace,
        REML = reml,
        factor = input_factor
    )
    
    # Increase Progressbar
    incProgress(0.3, detail = paste("Loading..."))
    
    # ¦§ Description 1 ####
    
    modelDesclist <- vector(mode = "list", length = 6)
    for(k in 1:input$mr_mdhglm_g_resp_count) {
        nMrand <- length(strsplit(input[[paste0('mr_mdhglm_m_model_',k)]], "|", fixed = TRUE)[[1]]) - 1
        mRandstring = NULL
        if (nMrand >= 1)
            for (i in 1:nMrand)
                mRandstring = paste0(c(mRandstring, input[[paste0("mr_mdhglm_m_rf_", i, "_", k)]]), collapse = ", ")
        else
            mRandstring = "NULL"
        
        modelDesc = matrix(c(
            input[[paste0('mr_mdhglm_m_model_',k)]],
            "phi ~ 1",
            "lambda ~ 1",
            input[[paste0('mr_mdhglm_m_link_',k)]],
            "log",
            "log",
            input[[paste0('mr_mdhglm_m_dist_',k)]],
            "gaussian",
            "gaussian",
            mRandstring,
            "NULL",
            "NULL"
        ),nrow=3, ncol=4)
        if (!is.null(input[[paste0("mr_mdhglm_p_check_phi_", k)]]) && input[[paste0("mr_mdhglm_p_check_phi_", k)]]) {
            nPrand <- length(strsplit(input[[paste0('mr_mdhglm_p_model_',k)]], "|", fixed = TRUE)[[1]]) - 1
            pRandstring = NULL
            if (nPrand >= 1)
                for (i in 1:nPrand)
                    pRandstring = paste0(c(pRandstring, input[[paste0("mr_mdhglm_p_rf_", i, "_", k)]]), collapse = ", ")
            else
                pRandstring = "NULL"
        
            modelDesc[2, 1] <- input[[paste0('mr_mdhglm_p_model_',k)]]
            modelDesc[2, 2] <- input[[paste0('mr_mdhglm_p_link_',k)]]
            modelDesc[2, 4] <- pRandstring
        }
        if (!is.null(input[[paste0("mr_mdhglm_l_check_lambda_", k)]]) && input[[paste0("mr_mdhglm_l_check_lambda_", k)]]) {
            nLrand <- length(strsplit(input[[paste0('mr_mdhglm_l_model_',k)]], "|", fixed = TRUE)[[1]]) - 1
            lRandstring = NULL
            if (nLrand >= 1)
                for (i in 1:nLrand)
                    lRandstring = paste0(c(lRandstring, input[[paste0("mr_mdhglm_l_rf_", i, "_", k)]]), collapse = ", ")
            else
                pRandstring = "NULL"
            
            modelDesc[3, 1] <- input[[paste0('mr_mdhglm_l_model_',k)]]
            modelDesc[3, 2] <- input[[paste0('mr_mdhglm_l_link_',k)]]
            modelDesc[3, 4] <- lRandstring
        }
        colnames(modelDesc) <- c("Model", "Link", "Dist", "Rand")
        rownames(modelDesc) <- c("Mean", "Phi", "Lambda")
        modelDesclist[[k]] <- modelDesc
    }
    
    # * except Shared ####
    if (input[['mr_mdhglm_g_corr_structure']] != "shared") {
        # ¦§ Likelihood ####
        
        Likelihoodlist <- vector(mode = "list", length = input$mr_mdhglm_g_resp_count)
        for(k in 1:input$mr_mdhglm_g_resp_count) {
            likelihood = matrix(
                c(
                    fittedModel[[k]]$ml, 
                    fittedModel[[k]]$rl, 
                    fittedModel[[k]]$caic,
                    fittedModel[[k]]$scaled_dv,
                    fittedModel[[k]]$df
                ), 
                ncol = 5
            )
            colnames(likelihood) <- c("-2ML", "-2RL", "cAIC", "Scaled Deviance", "df")
            Likelihoodlist[[k]] <- likelihood
        }
        
        # ¦§ Description 2 (Exp) ####
        # betaCoefflist <- vector(mode = "list", length = input$mr_mdhglm_g_resp_count)
        # phiCoefflist <- vector(mode = "list", length = input$mr_mdhglm_g_resp_count)
        # lambdaCoefflist <- vector(mode = "list", length = input$mr_mdhglm_g_resp_count)
        # 
        for(k in 1:input$mr_mdhglm_g_resp_count) {
            if (length(colnames(fittedModel[[k]]$beta_coeff)) == 3) {
                pValue <- 2 * pnorm(abs(fittedModel[[k]]$beta_coeff[, 3]), lower.tail = FALSE)
                llValue <- fittedModel[[k]]$beta_coeff[, 1] - 1.96 * fittedModel[[k]]$beta_coeff[, 2]
                ulValue <- fittedModel[[k]]$beta_coeff[, 1] + 1.96 * fittedModel[[k]]$beta_coeff[, 2]
    
                fittedModel[[k]]$beta_coeff <- cbind(fittedModel[[k]]$beta_coeff, pValue, llValue, ulValue)
                colnames(fittedModel[[k]]$beta_coeff) <- c("Estimate", "Std. Error", "t-value", "p_val", "LL", "UL")
            }
    
            if (length(colnames(fittedModel[[k]]$beta_coeff)) == 4) {
                llValue <- fittedModel[[k]]$beta_coeff[, 1] - 1.96 * fittedModel[[k]]$beta_coeff[, 2]
                ulValue <- fittedModel[[k]]$beta_coeff[, 1] + 1.96 * fittedModel[[k]]$beta_coeff[, 2]
    
                fittedModel[[k]]$beta_coeff <- cbind(fittedModel[[k]]$beta_coeff, llValue, ulValue)
                colnames(fittedModel[[k]]$beta_coeff) <- c("Estimate", "Std. Error", "t-value", "p_val", "LL", "UL")
            }
    
            if (!is.null(fittedModel[[k]]$beta_coeff) && !is.null(input[[paste0('mr_mdhglm_m_check_exp_', k)]]) && input[[paste0('mr_mdhglm_m_check_exp_', k)]]) {
                fittedModel[[k]]$beta_coeff[, 1] <- exp(fittedModel[[k]]$beta_coeff[, 1])
                colnames(fittedModel[[k]]$beta_coeff)[1] <- "Estimate(Exp)"
                fittedModel[[k]]$beta_coeff[, 5] <- exp(fittedModel[[k]]$beta_coeff[, 5])
                colnames(fittedModel[[k]]$beta_coeff)[5] <- "LL(Exp)"
                fittedModel[[k]]$beta_coeff[, 6] <- exp(fittedModel[[k]]$beta_coeff[, 6])
                colnames(fittedModel[[k]]$beta_coeff)[6] <- "UL(Exp)"
                fittedModel[[k]]$beta_coeff <- subset(fittedModel[[k]]$beta_coeff, select=-2)
            }
    
            if (!is.null(fittedModel[[k]]$phi_coeff)) {
                fittedModel[[k]]$phi_coeff <- cbind(fittedModel[[k]]$phi_coeff, "Estimate(Exp)" = exp(fittedModel[[k]]$phi_coeff[, 1]))
                fittedModel[[k]]$phi_coeff <- subset(fittedModel[[k]]$phi_coeff, select=c(1,4,2,3))
            }
          
            if (!is.null(fittedModel[[k]]$alpha_coeff)) {
                fittedModel[[k]]$alpha_coeff <- cbind(fittedModel[[k]]$alpha_coeff, "Estimate(Exp)" = exp(fittedModel[[k]]$alpha_coeff[, 1]))
                fittedModel[[k]]$alpha_coeff <- subset(fittedModel[[k]]$alpha_coeff, select=c(1,4,2,3))
            }
          
            if (!is.null(fittedModel[[k]]$lambda_coeff)) {
                fittedModel[[k]]$lambda_coeff <- cbind(fittedModel[[k]]$lambda_coeff, "Estimate(Exp)" = exp(fittedModel[[k]]$lambda_coeff[, 1]))
                fittedModel[[k]]$lambda_coeff <- subset(fittedModel[[k]]$lambda_coeff, select=c(1,4,2,3))
            }
          
            if (!is.null(fittedModel[[k]]$tau_coeff)) {
                fittedModel[[k]]$tau_coeff <- cbind(fittedModel[[k]]$tau_coeff, "Estimate(Exp)" = exp(fittedModel[[k]]$tau_coeff[, 1]))
                fittedModel[[k]]$tau_coeff <- subset(fittedModel[[k]]$tau_coeff, select=c(1,4,2,3))
            }
        }
    
        # ¦§ Output Likelihood ####
        outputLikelihooddf = 0
        for(k in 1:input$mr_mdhglm_g_resp_count) {
            outputLikelihooddf <- outputLikelihooddf + fittedModel[[k]]$df
        }
        
        outputLikelihoodlist <- matrix(c(0, 0, 0, outputLikelihooddf), nrow = 1)
        colnames(outputLikelihoodlist) <- c("-2ML", "-2RL", "cAIC", "df")
        for(k in 1:input$mr_mdhglm_g_resp_count) {
            outputLikelihoodlist[1] <- outputLikelihoodlist[1] + fittedModel[[k]]$ml
            outputLikelihoodlist[2] <- outputLikelihoodlist[2] + fittedModel[[k]]$rl
            outputLikelihoodlist[3] <- outputLikelihoodlist[3] + fittedModel[[k]]$caic
        }
        
        # ¦§ Option ####
        fittedModel$option <-list(
            RespDist = input_respdist,
            DataMain = input_datamain,
            structure = input$mr_mdhglm_g_corr_structure,
            MeanModel = MM_list,
            DispersionModel = DM_list,
            mord = mord,
            dord = dord,
            order = orderLaplace,
            REML = reml,
            factor = input_factor,
            ModelDesc = modelDesclist,
            outputModelDesc = modelDesclist,
            Likelihood = Likelihoodlist,
            outputLikelihood = outputLikelihoodlist,
            RespCount = input$mr_mdhglm_g_resp_count,
            Rcodes = Rcodes
        )
        
        appendTab(inputId = "mr_mdhglm_resulttabset",
            tabPanel(
                "Model Checking plots",
                br(),
                mainPanel(
                    uiOutput("mr_mdhglm_g_mcplots"),
                    width = 12
                )
            )
        )
        
        # ¦§ Comparison Model ####
        mdhglmComparisonmodel = NULL
        
        if (!is.null(input$mr_mdhglm_m_check_comparison_1) && input$mr_mdhglm_m_check_comparison_1) {
            g_mr_mdhglm_r_comparisonmodelcount <<- g_mr_mdhglm_r_comparisonmodelcount + 1
            mdhglmComparisonmodel <- matrix(
                c(
                    as.character(round(outputLikelihoodlist[1], digits=5)),
                    as.character(round(outputLikelihoodlist[2], digits=5)),
                    as.character(round(outputLikelihoodlist[3], digits=5)),
                    as.character(round(outputLikelihoodlist[4], digits=5))
                ),
                ncol = 4
            )
            
            mdhglmComparisonmodel <- rbind(g_mr_mdhglm_r_comparisonmodel_1, mdhglmComparisonmodel)
            colnames(mdhglmComparisonmodel) <- c("-2ML", "-2RL", "cAIC", "df")
            g_mr_mdhglm_r_comparisonmodel_1 <<- mdhglmComparisonmodel
            fittedModel$option$Comparisonmodel = mdhglmComparisonmodel
            
            g_mr_mdhglm_r_comparisonmodel_2[[g_mr_mdhglm_r_comparisonmodelcount]] <<- modelDesclist
            fittedModel$option$ComparisonmodelDesc = g_mr_mdhglm_r_comparisonmodel_2
            fittedModel$option$Comparisonmodelcount = g_mr_mdhglm_r_comparisonmodelcount
            
            g_mr_mdhglm_r_comparisonmodel_3[[g_mr_mdhglm_r_comparisonmodelcount]] <<- matrix(input$mr_mdhglm_g_corr_structure)
            colnames(g_mr_mdhglm_r_comparisonmodel_3[[g_mr_mdhglm_r_comparisonmodelcount]]) <<- "Correlation Structure"
            fittedModel$option$ComparisonmodelCorr = g_mr_mdhglm_r_comparisonmodel_3
            
            appendTab(inputId = "mr_mdhglm_resulttabset",
                tabPanel(
                    "Comparison Model",
                    br(),
                    uiOutput("mr_mdhglm_r_comparisonmodeldesc"),
                    tableOutput("mr_mdhglm_r_comparisonmodel")
                )
            )
        }
        
        # REMOVE ####
        toto <<- fittedModel
        
        # Set Progressbar
        setProgress(1, detail = "Finish")
        return(fittedModel)
    
    } else {
        # * Shared ####
        
        sfittedModel <- list()
        # ¦§ Beta_Coefficient ####
        betaCoeff <- cbind(fittedModel$Beta_h, fittedModel$SE_Beta, fittedModel$t_Beta, fittedModel$p_Beta)
        colnames(betaCoeff) <- c("Estimate", "Std. Error", "t_value", "p.value")
        count = 0
        sbetaCoeff <- vector(mode = "list", length = input$mr_mdhglm_g_resp_count)
        for(k in 1:input$mr_mdhglm_g_resp_count) {
            lengthCount = length(input[[paste0('mr_mdhglm_m_variable_', k)]]$right) + 1
            sbetaCoeff[[k]] <- betaCoeff[(1+count):(count+lengthCount),]
            rownames(sbetaCoeff[[k]]) <- c("(Intercept)", input[[paste0('mr_mdhglm_m_variable_', k)]]$right)
            count <- count + lengthCount
        }
        sfittedModel$sbeta <- sbetaCoeff
        
        # ¦§ Phi_Coefficient ####
        phiCoeff <- cbind(fittedModel$Log_Phi, fittedModel$SE_Log_Phi)
        colnames(phiCoeff) <- c("Estimate", "Std. Error")
        sphiCoeff <- vector(mode = "list", length = input$mr_mdhglm_g_resp_count)
        for(k in 1:input$mr_mdhglm_g_resp_count) {
            sphiCoeff[[k]] <- phiCoeff[k,]  
            sphiCoeff[[k]] <- matrix(sphiCoeff[[k]], nrow = 1)
            colnames(sphiCoeff[[k]]) <- c("Estimate", "Std. Error")
        }
        sfittedModel$sphi <- sphiCoeff
        
        # ¦§ Lambda_Coefficient ####
        lambdaCoeff <- cbind(fittedModel$Log_Lambda, fittedModel$SE_Log_Lambda)
        colnames(lambdaCoeff) <- c("Estimate", "Std. Error")
        slambdaCoeff <- vector(mode = "list", length = input$mr_mdhglm_g_resp_count)
        slambdaCoeff[[1]] <- lambdaCoeff
        
        sfittedModel$slambda <- slambdaCoeff
        
        # ¦§ Shared Parameter ####
        sharedParameter <- cbind(fittedModel$Shared, fittedModel$SE_Shared)
        colnames(sharedParameter) <- c("Estimate", "Std. Error")
        
        # ¦§ Likelihood ####
        likelihood = matrix(fittedModel$CAIC)
        colnames(likelihood) <- c("cAIC")
        
        # ¦§ Option ####
        sfittedModel$option<-list(
            ModelDesc = modelDesclist,
            structure = input$mr_mdhglm_g_corr_structure,
            RespCount = input$mr_mdhglm_g_resp_count,
            Rcodes = Rcodes,
            sharedParameter = sharedParameter,
            sLikelihood = likelihood
        )
        # Set Progressbar
        
        # ¦¦ Comparison Model ####
        mdhglmComparisonmodel = NULL
        
        if (!is.null(input$mr_mdhglm_m_check_comparison_1) && input$mr_mdhglm_m_check_comparison_1) {
            g_mr_mdhglm_r_comparisonmodelcount <<- g_mr_mdhglm_r_comparisonmodelcount + 1
            mdhglmComparisonmodel <- matrix(
                c(
                    "NULL",
                    "NULL",
                    as.character(round(fittedModel$CAIC, digits=5)),
                    "NULL"
                ),
                ncol = 4
            )
            
            mdhglmComparisonmodel <- rbind(g_mr_mdhglm_r_comparisonmodel_1, mdhglmComparisonmodel)
            
            colnames(mdhglmComparisonmodel) <- c("-2ML", "-2RL", "cAIC", "df")
            
            g_mr_mdhglm_r_comparisonmodel_1 <<- mdhglmComparisonmodel
            
            sfittedModel$option$Comparisonmodel = mdhglmComparisonmodel
            
            g_mr_mdhglm_r_comparisonmodel_2[[g_mr_mdhglm_r_comparisonmodelcount]] <<- modelDesclist
            sfittedModel$option$ComparisonmodelDesc = g_mr_mdhglm_r_comparisonmodel_2
            sfittedModel$option$Comparisonmodelcount = g_mr_mdhglm_r_comparisonmodelcount
            
            g_mr_mdhglm_r_comparisonmodel_3[[g_mr_mdhglm_r_comparisonmodelcount]] <<- matrix(input$mr_mdhglm_g_corr_structure)
            colnames(g_mr_mdhglm_r_comparisonmodel_3[[g_mr_mdhglm_r_comparisonmodelcount]]) <<- "Correlation Structure"
            sfittedModel$option$ComparisonmodelCorr = g_mr_mdhglm_r_comparisonmodel_3
            
            appendTab(inputId = "mr_mdhglm_resulttabset",
                tabPanel(
                    "Comparison Model",
                    br(),
                    uiOutput("mr_mdhglm_r_comparisonmodeldesc"),
                    tableOutput("mr_mdhglm_r_comparisonmodel")
                )
            )
        }
        
        setProgress(1, detail = "Finish")
        return(sfittedModel)
    }
    
    }) # End Progressbar
})

# MDHGLM tabpanel ####
output$mr_mdhglm_g_tabpanel <- renderUI({
    do.call(tabsetPanel, c(id = 'mdhglm_tabp', lapply(1:(mr_mdhglm_g_resp_value() + 1), function(i) {
        if (i <= mr_mdhglm_g_resp_value()) {
            tabPanel(
                title = paste0('Response', i),
                br(),
                uiOutput(paste0("mr_mdhglm_m_model_", i)),
                uiOutput(paste0("mr_mdhglm_p_model_", i)),
                uiOutput(paste0("mr_mdhglm_l_model_", i)),
                uiOutput(paste0("mr_mdhglm_r_accordion_", i))
            )
        } else {
            tabPanel(
                title = paste0('Settings'),
                h3(strong("Additional Settings")),
                br(),
                uiOutput("mr_mdhglm_a_mord"),
                uiOutput("mr_mdhglm_a_dord"),
                uiOutput("mr_mdhglm_a_order"),
                uiOutput("mr_mdhglm_a_reml")
            )
        }

    })))
})

# MDHGLM Accordion ####

lapply(1:6, function(k) {  
    output[[paste0('mr_mdhglm_r_accordion_',k)]]<-renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        mdhglmAccordion <- bs_accordion_sidebar(
            id = paste0("mr_mdhglm_r_accordion_", k),
            spec_side = c(width = 3, offset = 0),
            spec_main = c(width = 9, offset = 0)    
        )
        
        mdhglmAccordion <- mdhglmAccordion %>%
            bs_append(
                title_side = "Mean",
                content_side = NULL,
                content_main = div(
                    h3(strong("Model for Mean")),
                    uiOutput(paste0("mr_mdhglm_m_resp_", k)),
                    uiOutput(paste0("mr_mdhglm_m_variable_", k)),
                    fluidRow(
                        column(
                            8,
                            uiOutput(paste0("mr_mdhglm_m_interaction_", k))
                        ),
                        column(
                            4,
                            uiOutput(paste0("mr_mdhglm_m_interactionappend_", k)),
                            style = "text-align:right; padding:15px"
                        )
                    ),
                    uiOutput(paste0("mr_mdhglm_m_rand_", k)),
                    uiOutput(paste0("mr_mdhglm_m_check_slope_", k)),
                    uiOutput(paste0("mr_mdhglm_m_check_withoutslope_", k)),
                    uiOutput(paste0("mr_mdhglm_m_slope_", k)),
                    uiOutput(paste0("mr_mdhglm_m_check_randinteraction_", k)),
                    fluidRow(
                        column(
                            8,
                            uiOutput(paste0("mr_mdhglm_m_randinteraction_", k))
                        ),
                        column(
                            4,
                            uiOutput(paste0("mr_mdhglm_m_randinteractionappend_", k)),
                            style = "text-align:right; padding:15px"
                        )
                    ), #new components
                    uiOutput(paste0("mr_mdhglm_m_randfamily_", k)),
                    hr(style = "border-color: #2C3E50;"),
                    uiOutput(paste0("mr_mdhglm_m_dist_", k)),
                    # uiOutput(paste0("mr_mdhglm_m_check_binomd_", k)),
                    # uiOutput(paste0("mr_mdhglm_m_binomd_", k)),
                    uiOutput(paste0("mr_mdhglm_m_link_", k)),
                    uiOutput(paste0("mr_mdhglm_m_check_nointercept_", k)),
                    uiOutput(paste0("mr_mdhglm_m_check_offset_", k)),
                    uiOutput(paste0("mr_mdhglm_m_offset_", k)),
                    uiOutput(paste0("mr_mdhglm_m_factor_", k)),
                    hr(style = "border-color: #2C3E50;"),
                    uiOutput(paste0("mr_mdhglm_m_check_rcodes_", k)),
                    uiOutput(paste0("mr_mdhglm_m_check_exp_", k))
                )
            )
        
        mdhglmAccordion <- mdhglmAccordion %>%
        bs_append(
            title_side = "Phi",
            content_side = uiOutput(paste0("mr_mdhglm_p_check_phi_", k)),
            content_main = div(
                h3(strong("Model for phi")),
                
                uiOutput(paste0("mr_mdhglm_p_resp_", k)),
                uiOutput(paste0("mr_mdhglm_p_variable_", k)),
                fluidRow(
                    column(
                        8,
                        uiOutput(paste0("mr_mdhglm_p_interaction_", k))
                    ),
                    column(
                        4,
                        uiOutput(paste0("mr_mdhglm_p_interactionappend_", k)),
                        style = "text-align:right; padding:15px"
                    )
                ),
                uiOutput(paste0("mr_mdhglm_p_rand_", k)),
                uiOutput(paste0("mr_mdhglm_p_check_slope_", k)),
                uiOutput(paste0("mr_mdhglm_p_check_withoutslope_", k)),
                uiOutput(paste0("mr_mdhglm_p_slope_", k)),
                uiOutput(paste0("mr_mdhglm_p_check_randinteraction_", k)),
                fluidRow(
                    column(
                        8,
                        uiOutput(paste0("mr_mdhglm_p_randinteraction_", k))
                    ),
                    column(
                        4,
                        uiOutput(paste0("mr_mdhglm_p_randinteractionappend_", k)),
                        style = "text-align:right; padding:15px"
                    )
                ), #new components
                uiOutput(paste0("mr_mdhglm_p_randfamily_", k)),
                uiOutput(paste0("mr_mdhglm_p_link_", k)),
                uiOutput(paste0("mr_mdhglm_p_check_offset_", k)),
                uiOutput(paste0("mr_mdhglm_p_offset_", k))
            )
        )
        
        mdhglmAccordion <- mdhglmAccordion %>%
        bs_append(
            title_side = "Lambda",
            content_side = uiOutput(paste0("mr_mdhglm_l_check_lambda_", k)),
            content_main = div(
                h3(strong("Model for Lambda")),
                
                uiOutput(paste0("mr_mdhglm_l_resp_", k)),
                uiOutput(paste0("mr_mdhglm_l_variable_", k)),
                fluidRow(
                    column(
                        8,
                        uiOutput(paste0("mr_mdhglm_l_interaction_", k))
                    ),
                    column(
                        4,
                        uiOutput(paste0("mr_mdhglm_l_interactionappend_", k)),
                        style = "text-align:right; padding:15px"
                    )
                ),
                uiOutput(paste0("mr_mdhglm_l_rand_", k)),
                uiOutput(paste0("mr_mdhglm_l_check_randinteraction_", k)),
                fluidRow(
                    column(
                        8,
                        uiOutput(paste0("mr_mdhglm_l_randinteraction_", k))
                    ),
                    column(
                        4,
                        uiOutput(paste0("mr_mdhglm_l_randinteractionappend_", k)),
                        style = "text-align:right; padding:15px"
                    )
                ), #new components
                uiOutput(paste0("mr_mdhglm_l_randfamily_", k)),
                uiOutput(paste0("mr_mdhglm_l_link_", k))
            )
        )
        
        div(
            mdhglmAccordion,
            use_bs_tooltip(),
            use_bs_accordion_sidebar() # needs to be at end, for some reason
        )
    })
})

# MDHGLM output ####
lapply(1:6, function(k) {
    output[[paste0('mr_mdhglm_r_outputmainsummary_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()$option$outputModelDesc[[k]]))
            return()
      
        mr_mdhglm_results()$option$outputModelDesc[[k]]
    }, rownames = TRUE, bordered = TRUE, caption = paste0("Model Description ", k), spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    output[[paste0('mr_mdhglm_r_mainsummary_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()$option$ModelDesc[[k]]))
            return()
      
        mr_mdhglm_results()$option$ModelDesc[[k]]
    }, rownames = TRUE, bordered = TRUE, caption = "Model Description", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    output[[paste0('mr_mdhglm_m_coeff_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()[[k]]$beta_coeff))
            return()
      
        mr_mdhglm_results()[[k]]$beta_coeff
    }, rownames = TRUE, bordered = TRUE, caption = "Estimate from Mean Model", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    output[[paste0('mr_mdhglm_l_coeff_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()[[k]]$lambda_coeff))
            return()
      
        mr_mdhglm_results()[[k]]$lambda_coeff
    }, rownames = TRUE, bordered = TRUE, caption = "Estimate for log(Lambda)", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    output[[paste0('mr_mdhglm_l_taucoeff_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()[[k]]$alpha_coeff))
            return()
      
        mr_mdhglm_results()[[k]]$alpha_coeff
    }, rownames = TRUE, bordered = TRUE, caption = "Estimate for log(Tau)", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)    

    output[[paste0('mr_mdhglm_p_coeff_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()[[k]]$phi_coeff))
            return()
      
        mr_mdhglm_results()[[k]]$phi_coeff
    }, rownames = TRUE, bordered = TRUE, caption = "Estimate from Dispersion Model", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    output[[paste0('mr_mdhglm_p_alphacoeff_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()[[k]]$tau_coeff))
            return()
      
        mr_mdhglm_results()[[k]]$tau_coeff
    }, rownames = TRUE, bordered = TRUE, caption = "Estimate for log(Alpha)", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    output[[paste0('mr_mdhglm_r_likelihood_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()$option$Likelihood[[k]]))
            return()
      
        mr_mdhglm_results()$option$Likelihood[[k]]
    }, rownames = FALSE, bordered = TRUE, caption = "Likelihood", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    
    output[[paste0('mr_mdhglm_r_sbeta_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()$sbeta[[k]]))
            return()
      
        mr_mdhglm_results()$sbeta[[k]]
    }, rownames = TRUE, bordered = TRUE, caption = "Beta Coefficient", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    output[[paste0('mr_mdhglm_r_sphi_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()$sphi[[k]]))
            return()
      
        mr_mdhglm_results()$sphi[[k]]
    }, rownames = FALSE, bordered = TRUE, caption = "Phi Coefficient", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    output[[paste0('mr_mdhglm_r_slambda_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()$slambda[[k]]))
            return()
      
        mr_mdhglm_results()$slambda[[k]]
    }, rownames = FALSE, bordered = TRUE, caption = "Lambda Coefficient", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    output[[paste0('mr_mdhglm_r_rcodes_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()$option$Rcodes[[k]]))
            return()
      
        mr_mdhglm_results()$option$Rcodes[[k]]
    }, rownames = TRUE, bordered = TRUE, caption = "R Codes", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
})

output$mr_mdhglm_r_cor <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$mr_mdhglm_g_run
    
    if (g_resetresult == FALSE || is.null(mr_mdhglm_results()$cor))
        return()
  
    mr_mdhglm_results()$cor
}, rownames = TRUE, colnames = TRUE, bordered = TRUE, caption = "Estimate for Correlation", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_mdhglm_r_outputlikelihood <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$mr_mdhglm_g_run
    
    if (g_resetresult == FALSE || is.null(mr_mdhglm_results()$option$outputLikelihood))
        return()
  
    mr_mdhglm_results()$option$outputLikelihood
}, rownames = FALSE, bordered = TRUE, caption = "Likelihood", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_mdhglm_r_sharedparameter <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$mr_mdhglm_g_run
    
    if (g_resetresult == FALSE || is.null(mr_mdhglm_results()$option$sharedParameter))
        return()
  
    mr_mdhglm_results()$option$sharedParameter
}, rownames = TRUE, bordered = TRUE, caption = "Shared Parameter", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_mdhglm_r_slikelihood <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$mr_mdhglm_g_run
    
    if (g_resetresult == FALSE || is.null(mr_mdhglm_results()$option$sLikelihood))
        return()
  
    mr_mdhglm_results()$option$sLikelihood
}, rownames = TRUE, bordered = TRUE, caption = "Likelihood", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_mdhglm_r_comparisonmodel <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$mr_mdhglm_g_run
    
    if (g_resetresult == FALSE || is.null(mr_mdhglm_results()$option$Comparisonmodel))
        return()
  
    mr_mdhglm_results()$option$Comparisonmodel
}, rownames = TRUE, bordered = TRUE, caption = "Likelihood", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

# MDHGLM Plots ####
lapply(1:6, function(k) {
    # ¦§ Mean Plot #####
    output[[paste0('mr_mdhglm_m_selectplot1_',k)]]<-renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()[[k]]$nrand))
            return()
        
        meanNumberRand = mr_mdhglm_results()[[k]]$nrand
        
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
            inputId = sprintf("mr_mdhglm_m_selectplot1_%s", k),
            label = "Select Plot", 
            choiceNames = choiceName,
            choiceValues = choiceValue, 
            selected = "mean",
            direction = "vertical"
        )
    })
    
    output[[paste0('mr_mdhglm_m_plot1_',k)]] <- renderPlot({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()))
            return()
        
        local({
            output[[paste0('mr_mdhglm_m_downloadplot1_',k)]] <- downloadPlot({
                if(is.null(input[[paste0("mr_mdhglm_m_selectplot1_", k)]]))
                    plotType = "mean"
                else
                    plotType = input[[paste0("mr_mdhglm_m_selectplot1_", k)]]
                
                res <- mr_mdhglm_results()[[k]]
                ggplotdhglm(res, type = plotType)
            })
        })
        
        if(is.null(input[[paste0("mr_mdhglm_m_selectplot1_", k)]]))
            plotType = "mean"
        else
            plotType = input[[paste0("mr_mdhglm_m_selectplot1_", k)]]
        
        res <- mr_mdhglm_results()[[k]]
        ggplotdhglm(res, type = plotType)
    })
    
    output[[paste0('mr_mdhglm_m_showplot1_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()))
            return()
        
        div(
            h4("Model Checking Plots for Mean"),
            div(
                style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
                withSpinner(
                    plotOutput(paste0('mr_mdhglm_m_plot1_', k), height = "600px"),
                    type = 1,
                    color = "#2c3e50",
                    size = 1.2
                ),
                div(
                    style = "position: absolute; left: 0.5em; bottom: 0.5em;",
                    dropdown(
                        uiOutput(paste0('mr_mdhglm_m_selectplot1_', k)),
                        size = "xs",
                        icon = icon("chart-bar", class = "opt"), 
                        up = TRUE
                    )
                ),
                div(
                    style = "position: absolute; left:4em; bottom: 0.5em;",
                    dropdown(
                        downloadButton(outputId = paste0('mr_mdhglm_m_downloadplot1_', k), label = "Download Plot"),
                        size = "xs",
                        icon = icon("download", class = "opt"),
                        up = TRUE
                    )
                )
            ),
            br()
        )
    })
    
    # ¦§ Phi Plot #####
    output[[paste0('mr_mdhglm_p_selectplot1_',k)]]<-renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()[[k]]$resid_phi))
            return()
        
        choiceName = c("Residuals")
        choiceValue = c("phi")
    
        radioGroupButtons(
            inputId = sprintf("mr_mdhglm_p_selectplot1_%s", k),
            label = "Select Plot", 
            choiceNames = choiceName,
            choiceValues = choiceValue, 
            selected = "phi",
            direction = "vertical"
        )
    })
    
    output[[paste0('mr_mdhglm_p_plot1_',k)]] <- renderPlot({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        local({
            output[[paste0('mr_mdhglm_p_downloadplot1_',k)]] <- downloadPlot({
                if(is.null(input[[paste0("mr_mdhglm_p_selectplot1_", k)]]))
                    plotType = "phi"
                else
                    plotType = input[[paste0("mr_mdhglm_p_selectplot1_", k)]]
                
                res <- mr_mdhglm_results()[[k]]
                ggplotdhglm(res, type = plotType)
            })
        })
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()[[k]]$resid_phi))
            return()
        
        if(is.null(input[[paste0("mr_mdhglm_p_selectplot1_", k)]]))
            plotType = "phi"
        else
            plotType = input[[paste0("mr_mdhglm_p_selectplot1_", k)]]
        
        res <- mr_mdhglm_results()[[k]]
        ggplotdhglm(res, type = plotType)
    })
    
    output[[paste0('mr_mdhglm_p_showplot1_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()[[k]]$resid_phi))
            return()
        
        div(
            h4("Model Checking Plots for Phi"),
            div(
                style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
                withSpinner(
                    plotOutput(paste0('mr_mdhglm_p_plot1_', k), height = "600px"),
                    type = 1,
                    color = "#2c3e50",
                    size = 1.2
                ),
                div(
                    style = "position: absolute; left: 0.5em; bottom: 0.5em;",
                    dropdown(
                        uiOutput(paste0('mr_mdhglm_p_selectplot1_', k)),
                        size = "xs",
                        icon = icon("chart-bar", class = "opt"), 
                        up = TRUE
                    )
                ),
                div(
                    style = "position: absolute; left:4em; bottom: 0.5em;",
                    dropdown(
                        downloadButton(outputId = paste0('mr_mdhglm_p_downloadplot1_', k), label = "Download Plot"),
                        size = "xs",
                        icon = icon("download", class = "opt"),
                        up = TRUE
                    )
                )
            ),
            br()
        )
    })
    
    # ¦¦ Lambda Plot #####
    output[[paste0('mr_mdhglm_l_selectplot1_',k)]]<-renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()[[k]]$resid_lambda))
            return()
        
        choiceName = c("Residuals")
        choiceValue = c("lambda")
    
        radioGroupButtons(
            inputId = sprintf("mr_mdhglm_l_selectplot1_%s", k),
            label = "Select Plot", 
            choiceNames = choiceName,
            choiceValues = choiceValue, 
            selected = "lambda",
            direction = "vertical"
        )
    })
    
    output[[paste0('mr_mdhglm_l_plot1_',k)]] <- renderPlot({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        local({
            output[[paste0('mr_mdhglm_l_downloadplot1_',k)]] <- downloadPlot({
                if(is.null(input[[paste0("mr_mdhglm_l_selectplot1_", k)]]))
                    plotType = "lambda"
                else
                    plotType = input[[paste0("mr_mdhglm_l_selectplot1_", k)]]
                
                res <- mr_mdhglm_results()[[k]]
                ggplotdhglm(res, type = plotType)
            })
        })
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()[[k]]$resid_lambda))
            return()
        
        if(is.null(input[[paste0("mr_mdhglm_l_selectplot1_", k)]]))
            plotType = "lambda"
        else
            plotType = input[[paste0("mr_mdhglm_l_selectplot1_", k)]]
        
        res <- mr_mdhglm_results()[[k]]
        ggplotdhglm(res, type = plotType)
    })
    
    output[[paste0('mr_mdhglm_l_showplot1_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()[[k]]$resid_lambda))
            return()
        
        div(
            h4("Model Checking Plots for Lambda"),
            div(
                style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
                withSpinner(
                    plotOutput(paste0('mr_mdhglm_l_plot1_', k), height = "600px"),
                    type = 1,
                    color = "#2c3e50",
                    size = 1.2
                ),
                div(
                    style = "position: absolute; left: 0.5em; bottom: 0.5em;",
                    dropdown(
                        uiOutput(paste0('mr_mdhglm_l_selectplot1_', k)),
                        size = "xs",
                        icon = icon("chart-bar", class = "opt"), 
                        up = TRUE
                    )
                ),
                div(
                    style = "position: absolute; left:4em; bottom: 0.5em;",
                    dropdown(
                        downloadButton(outputId = paste0('mr_mdhglm_l_downloadplot1_', k), label = "Download Plot"),
                        size = "xs",
                        icon = icon("download", class = "opt"),
                        up = TRUE
                    )
                )
            ),
            br()
        )
    })
})
    
# MDHGLM Comparison Model ####

lapply(1:10, function(k) {
    output[[paste0('mr_mdhglm_r_comparisoncorr_', k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_mdhglm_g_run
        
        if (g_resetresult == FALSE || is.null(mr_mdhglm_results()$option$ComparisonmodelCorr[[k]]))
            return()
      
        mr_mdhglm_results()$option$ComparisonmodelCorr[[k]]
    }, rownames = FALSE, bordered = TRUE, caption = "", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
  
    lapply(1:6, function(j) {
        output[[paste0('mr_mdhglm_r_comparisonsummary_', k, "_", j)]] <- renderTable({
            input$file
            input$resetData
            input$di_option_run
            input$mr_mdhglm_g_run
            
            if (g_resetresult == FALSE || is.null(mr_mdhglm_results()$option$ComparisonmodelDesc[[k]][[j]]))
                return()
          
            mr_mdhglm_results()$option$ComparisonmodelDesc[[k]][[j]]
        }, rownames = TRUE, bordered = TRUE, caption = paste0("Model Description ", j), spacing = "m",
        caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    })
})

output$mr_mdhglm_r_comparisonmodeldesc <- renderUI({
    do.call(tabsetPanel, c(id = 'comparisontabp', lapply(1:mr_mdhglm_results()$option$Comparisonmodelcount, function(i) {
        tabPanel(
            paste0("Result ", i),
            uiOutput(paste0('mr_mdhglm_r_comparisoncorr_', i)),
            uiOutput(paste0('mr_mdhglm_r_comparisonsummary_', i, "_", 1)),
            uiOutput(paste0('mr_mdhglm_r_comparisonsummary_', i, "_", 2)),
            uiOutput(paste0('mr_mdhglm_r_comparisonsummary_', i, "_", 3)),
            uiOutput(paste0('mr_mdhglm_r_comparisonsummary_', i, "_", 4)),
            uiOutput(paste0('mr_mdhglm_r_comparisonsummary_', i, "_", 5)),
            uiOutput(paste0('mr_mdhglm_r_comparisonsummary_', i, "_", 6)),
            hr()
        )
    })))
})

# MDHGLM summary ####
output$mr_mdhglm_g_summary <- renderUI({
    do.call(tabsetPanel, c(id = 'summ', lapply(0:mr_mdhglm_results()$option$RespCount, function(i) {
        if (i == 0) {
            tabPanel(
                paste0("Model Output"),
                
                # except shared
                uiOutput("mr_mdhglm_r_outputmainsummary_1"),
                uiOutput("mr_mdhglm_r_outputmainsummary_2"),
                uiOutput("mr_mdhglm_r_outputmainsummary_3"),
                uiOutput("mr_mdhglm_r_outputmainsummary_4"),
                uiOutput("mr_mdhglm_r_outputmainsummary_5"),
                uiOutput("mr_mdhglm_r_outputmainsummary_6"),
                uiOutput("mr_mdhglm_r_cor"),
                uiOutput("mr_mdhglm_r_outputlikelihood"),
                
                # shared
                uiOutput("mr_mdhglm_r_sharedparameter"),
                uiOutput("mr_mdhglm_r_slikelihood")
            )
        } else {
            tabPanel(
                paste0("Response ", i),
                uiOutput(paste0('mr_mdhglm_r_mainsummary_', i)),
                
                # except Shared
                uiOutput(paste0('mr_mdhglm_m_coeff_', i)),
                uiOutput(paste0('mr_mdhglm_l_coeff_', i)),
                uiOutput(paste0('mr_mdhglm_l_taucoeff_', i)),
                uiOutput(paste0('mr_mdhglm_p_coeff_', i)),
                uiOutput(paste0('mr_mdhglm_p_alphacoeff_', i)),
                uiOutput(paste0('mr_mdhglm_r_likelihood_', i)),
                
                # Shared
                uiOutput(paste0('mr_mdhglm_r_sbeta_', i)),
                uiOutput(paste0('mr_mdhglm_r_sphi_', i)),
                uiOutput(paste0('mr_mdhglm_r_slambda_', i)),
                
                uiOutput(paste0('mr_mdhglm_r_rcodes_', i))
            )
        }

    })))
})

output$mr_mdhglm_g_mcplots <- renderUI({
    do.call(tabsetPanel, c(id = 'plot', lapply(1:mr_mdhglm_results()$option$RespCount, function(i) {
        tabPanel(
            paste0("Response ", i),
            div(style = "padding:10px"),
            uiOutput(paste0('mr_mdhglm_m_showplot1_', i)),
            uiOutput(paste0('mr_mdhglm_p_showplot1_', i)),
            uiOutput(paste0('mr_mdhglm_l_showplot1_', i))
        )
    })))
})