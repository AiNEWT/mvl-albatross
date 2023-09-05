# MDSEM Components ####

output$mr_dsem2_g_respslider <- renderUI({
    sliderInput(
        inputId = sprintf("mr_dsem2_g_resp_count"),
        label = "Number of Response Variables",
        value = 2,
        min = 2,
        max = 6
    )
})

output$mr_dsem2_g_corr_structure<-renderUI({
#    selectInput(
#        "mr_dsem2_g_corr_structure",
#        "Correlation Structure",
#        choices = c("dynamic"),
#        selected = "dynamic",
#        multiple = FALSE
#    )
  textAreaInput("mr_dsem2_g_corr_structure","Correlation Structure for Mean",value="c(11,12,13,24,25,26)")
})

output$mr_dsem2_g_corr_structure_1<-renderUI({
  #    selectInput(
  #        "mr_dsem2_g_corr_structure",
  #        "Correlation Structure",
  #        choices = c("dynamic"),
  #        selected = "dynamic",
  #        multiple = FALSE
  #    )
  textAreaInput("mr_dsem2_g_corr_structure_1","Correlation Structure for Phi",value="c(11,22)")
})

output$mr_dsem2_g_corr_structure_2<-renderUI({
  #    selectInput(
  #        "mr_dsem2_g_corr_structure",
  #        "Correlation Structure",
  #        choices = c("dynamic"),
  #        selected = "dynamic",
  #        multiple = FALSE
  #    )
  textAreaInput("mr_dsem2_g_corr_structure_2","Correlation Structure for Lambda",value="c()")
})

mr_dsem2_g_resp_value <- reactive({
    if (is.null(input$mr_dsem2_g_resp_count))
        return(2)
    else
        return(input$mr_dsem2_g_resp_count)
})


# ?? Mean, Phi, Lambda ####
source('logic/MultiResponse/DSEM2/logic_mrDSEM2_m.R', local = T)
source('logic/MultiResponse/DSEM2/logic_mrDSEM2_p.R', local = T)
source('logic/MultiResponse/DSEM2/logic_mrDSEM2_l.R', local = T)
source('logic/MultiResponse/DSEM2/logic_mrDSEM2_s.R', local = T)


# ?? Additional Settings ####
output$mr_dsem2_a_mord<-renderUI({
    selectInput(
        "mr_dsem2_a_mord",
        "Order for Mean Model",
        choices=c("0"="0","1"="1"),
        multiple=FALSE
    )
})

output$mr_dsem2_a_dord<-renderUI({
    selectInput(
        "mr_dsem2_a_dord",
        "Order for Dispersion Model",
        choices=c("1"="1","2"="2"),
        multiple=FALSE
    )
})

output$mr_dsem2_a_order<-renderUI({
    selectInput(
        "mr_dsem2_a_order",
        "Order for Laplace's Approximation",
        choices=c("1"="1","2"="2"),
        multiple=FALSE
    )
})

output$mr_dsem2_a_reml<-renderUI({
    selectInput(
        "mr_dsem2_a_reml",
        "REML or ML",
        choice = c("REML" = TRUE,"ML" = FALSE),
        multiple=FALSE
    )
})

# MDHGLM Run ####

observeEvent(input$mr_dsem2_g_run, {
    g_resetresult <<- TRUE
})


mr_dsem2_results <- eventReactive(input$mr_dsem2_g_run, {
    # Start Progressbar  
    withProgress(message = 'HSEM', style = "notification", value = 0.1, {
    
#    removeTab(inputId = "mr_dsem2_resulttabset", target = "Model Checking Plot")
    removeTab(inputId = "mr_dsem2_resulttabset", target = "Model Comparison")
    
    # ?? Variable Declaration 1 ####
    MM_list <- vector(mode = "list", length = input$mr_dsem2_g_resp_count)
    DM_list <- vector(mode = "list", length = input$mr_dsem2_g_resp_count)
    Rcodes <- vector(mode = "list", length = input$mr_dsem2_g_resp_count)
    f <- data2
    
    for(k in 1:input$mr_dsem2_g_resp_count) {
        m_form <- input[[paste0("mr_dsem2_m_model_", k)]]
        p_form <- NULL
        l_form <- NULL

        if (!is.null(input[[paste0("mr_dsem2_p_check_phi_", k)]]) && input[[paste0("mr_dsem2_p_check_phi_", k)]])
            p_form <- input[[paste0("mr_dsem2_p_model_", k)]]
        else
            p_form <- phi~1

        if (!is.null(input[[paste0("mr_dsem2_l_check_lambda_", k)]]) && input[[paste0("mr_dsem2_l_check_lambda_", k)]])
            l_form <- input[[paste0("mr_dsem2_l_model_", k)]]
        else
            l_form <- lambda~1

        res <- m_form
        m_RandDistM <- NULL
        p_RandDistM <- NULL
        l_RandDistM <- NULL
        m_nRand <- length(input[[paste0("mr_dsem2_m_rand_", k)]])
        p_nRand <- length(input[[paste0("mr_dsem2_p_rand_", k)]])
        l_nRand <- length(input[[paste0("mr_dsem2_l_rand_", k)]])

        if (m_nRand > 0) {
            for (i in 1:m_nRand) {
                m_RandDistM <- c(m_RandDistM,input[[paste0("mr_dsem2_m_rf_",i,"_", k)]])
            }
        }
        if (!is.null(input[[paste0("mr_dsem2_p_check_phi_", k)]]) && input[[paste0("mr_dsem2_p_check_phi_", k)]]) {
            if (p_nRand > 0) {
                for (i in 1:p_nRand) {
                    p_RandDistM <- c(p_RandDistM,input[[paste0("mr_dsem2_p_rf_",i,"_", k)]])
                }
            }
        }
        LinPredRandVariance=NULL
        if (!is.null(input[[paste0("mr_dsem2_l_check_lambda_", k)]]) && input[[paste0("mr_dsem2_l_check_lambda_", k)]]) {
            if (l_nRand > 0) {
                for (i in 1:l_nRand) {
                    l_RandDistM <- c(l_RandDistM,input[[paste0("mr_dsem2_l_rf_",i,"_", k)]])
                }
            }
            LinPredRandVariance=l_form
        }

        MM = NULL
        DM = NULL
        
        # ?? MeanModel (MM) ####
        
        MM <-
            DHGLMMODELING(
                Model = "mean",
                Link = isolate(input[[paste0("mr_dsem2_m_link_", k)]]),
                LinPred = as.formula(m_form),
                RandDist = m_RandDistM,
                LinPredRandVariance = as.formula(LinPredRandVariance),
                RandDistRandVariance = l_RandDistM,
                LinkRandVariance = "log"
            )

        # ?? DispersionModel (DM) ####
        
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
                  Link = isolate(input[[paste0("mr_dsem2_m_link_", k)]]),
                  LinPred = as.formula(m_form),
                  RandDist = m_RandDistM
                )
        }

        if (is.null(DM[[4]])) DM[[4]]<-"gaussian"
        
        # ?? R Codes ####
        
        if (input[[paste0('mr_dsem2_m_check_rcodes_', k)]]) {
            Rcode1 <- NULL
            Rcode2 <- NULL
            Rcode3 <- NULL
            if (l_form == "lambda ~ 1") {
                Rcode1 <- paste0(
                    "MM <- DHGLMMODELING(Model = \"mean\", Link = \"",
                    input[[paste0("mr_dsem2_m_link_", k)]],
                    "\", LinPred = ",
                    m_form,
                    ", RandDist = \"",
                    m_RandDistM,
                    ")"
                )
            } else {
                Rcode1 <- paste0(
                    "MM <- DHGLMMODELING(Model = \"mean\", Link = \"",
                    input[[paste0("mr_dsem2_m_link_", k)]],
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
                input[[paste0("mr_dsem2_m_dist_", k)]],
                "\", DataMain = data, structure = \"",
                input$mr_dsem2_g_corr_structure,
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
    input_datamain <- vector(mode = "list", length = input$mr_dsem2_g_resp_count)
    input_respdist = NULL
    input_factor = NULL

    for(k in 1:input$mr_dsem2_g_resp_count) {
    #    if (input[['mr_dsem2_g_corr_structure']] == "factor") {
    #        input_factor <- as.numeric(c(input_factor, input[[paste0("mr_dsem2_m_factor_", k)]]))
    #    }
        input_respdist <- c(input_respdist, input[[paste0("mr_dsem2_m_dist_", k)]])
        input_datamain[[k]] <- f
    }
    
    # ?? Variable Declaration 2 ####
    
    mord = NULL
    dord = NULL
    orderLaplace = NULL
    reml = NULL
    if (is.null(input$mr_dsem2_a_mord))
        mord = 0
    else
        mord = as.numeric(input$re_dhglm_a_mord)
    
    if (is.null(input$mr_dsem2_a_dord))
        dord = 1
    else
        dord = as.numeric(input$mr_dsem2_a_dord)
    
    if (is.null(input$mr_dsem2_a_order))
        orderLaplace = 1
    else
        orderLaplace = as.numeric(input$mr_dsem2_a_order)
    
    if (is.null(input$mr_dsem2_a_reml))
        reml = TRUE
    else 
        reml = FALSE
    
#    corrStructure = input$mr_dsem2_g_corr_structure
    corrStructure = "independent"
    corr_structure=c(1,2,3,4,5,6,7,8)
    mr_dsem2_g_corr_structure="c(11,11,11,22,22,22)"
    mr_dsem2_g_corr_structure_1="c(11,22)"
    if (mr_dsem2_g_corr_structure=="c(11,12,13,24,25,26)" && mr_dsem2_g_corr_structure_1=="c(11,22)") corr_structure=c(1,2,3,4,5,6,7,8)
    if (mr_dsem2_g_corr_structure=="c(11,12,13,24,25,26)" && mr_dsem2_g_corr_structure_1=="c(11,21)") corr_structure=c(1,2,3,4,5,6,7,4)
    if (mr_dsem2_g_corr_structure=="c(11,12,13,21,22,23)" && mr_dsem2_g_corr_structure_1=="c(11,21)") corr_structure=c(1,2,3,4,1,2,3,4)
    if (mr_dsem2_g_corr_structure=="c(11,11,11,22,22,22)" && mr_dsem2_g_corr_structure_1=="c(11,21)") corr_structure=c(1,1,1,2,3,3,3,2)
    if (mr_dsem2_g_corr_structure=="c(11,11,11,22,22,22)" && mr_dsem2_g_corr_structure_1=="c(11,22)") corr_structure=c(1,1,1,2,3,3,3,4)
    

        # ?? Fitted Model ####
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
        factor = input_factor,
        corr_structure=corr_structure
    )
    
    # Increase Progressbar
    incProgress(0.3, detail = paste("Loading..."))
    
    # ?? Description 1 ####
    
    modelDesclist <- vector(mode = "list", length = 6)
    for(k in 1:input$mr_dsem2_g_resp_count) {
        nMrand <- length(strsplit(input[[paste0('mr_dsem2_m_model_',k)]], "|", fixed = TRUE)[[1]]) - 1
        mRandstring = NULL
        if (nMrand >= 1)
            for (i in 1:nMrand)
                mRandstring = paste0(c(mRandstring, input[[paste0("mr_dsem2_m_rf_", i, "_", k)]]), collapse = ", ")
        else
            mRandstring = "NULL"
        
        modelDesc = matrix(c(
            input[[paste0('mr_dsem2_m_model_',k)]],
            "phi ~ 1",
            "lambda ~ 1",
            input[[paste0('mr_dsem2_m_link_',k)]],
            "log",
            "log",
            input[[paste0('mr_dsem2_m_dist_',k)]],
            "gaussian",
            "gaussian",
            mRandstring,
            "NULL",
            "NULL"
        ),nrow=3, ncol=4)
        if (!is.null(input[[paste0("mr_dsem2_p_check_phi_", k)]]) && input[[paste0("mr_dsem2_p_check_phi_", k)]]) {
            nPrand <- length(strsplit(input[[paste0('mr_dsem2_p_model_',k)]], "|", fixed = TRUE)[[1]]) - 1
            pRandstring = NULL
            if (nPrand >= 1)
                for (i in 1:nPrand)
                    pRandstring = paste0(c(pRandstring, input[[paste0("mr_dsem2_p_rf_", i, "_", k)]]), collapse = ", ")
            else
                pRandstring = "NULL"
        
            modelDesc[2, 1] <- input[[paste0('mr_dsem2_p_model_',k)]]
            modelDesc[2, 2] <- input[[paste0('mr_dsem2_p_link_',k)]]
            modelDesc[2, 4] <- pRandstring
        }
        if (!is.null(input[[paste0("mr_dsem2_l_check_lambda_", k)]]) && input[[paste0("mr_dsem2_l_check_lambda_", k)]]) {
            nLrand <- length(strsplit(input[[paste0('mr_dsem2_l_model_',k)]], "|", fixed = TRUE)[[1]]) - 1
            lRandstring = NULL
            if (nLrand >= 1)
                for (i in 1:nLrand)
                    lRandstring = paste0(c(lRandstring, input[[paste0("mr_dsem2_l_rf_", i, "_", k)]]), collapse = ", ")
            else
                pRandstring = "NULL"
            
            modelDesc[3, 1] <- input[[paste0('mr_dsem2_l_model_',k)]]
            modelDesc[3, 2] <- input[[paste0('mr_dsem2_l_link_',k)]]
            modelDesc[3, 4] <- lRandstring
        }
        colnames(modelDesc) <- c("Model", "Link", "Dist", "Rand")
        rownames(modelDesc) <- c("Mean", "Phi", "Lambda")
        modelDesclist[[k]] <- modelDesc
    }
    
    # * except Shared ####
#    if (input[['mr_dsem2_g_corr_structure']] != "shared") {
        # ?? Likelihood ####
        
        Likelihoodlist <- vector(mode = "list", length = input$mr_dsem2_g_resp_count)
        for(k in 1:input$mr_dsem2_g_resp_count) {
            likelihood = matrix(
                c(
                    fittedModel[[k]]$ml, 
                    fittedModel[[k]]$rl, 
                    fittedModel[[k]]$scaled_dv,
                    fittedModel[[k]]$df,
                    fittedModel[[k]]$caic
                ), 
                ncol = 5
            )
            colnames(likelihood) <- c("-2ML", "-2RL", "Scaled Deviance", "df", "cAIC")
            Likelihoodlist[[k]] <- likelihood
        }
        
        # ?? Description 2 (Exp) ####
        # betaCoefflist <- vector(mode = "list", length = input$mr_dsem2_g_resp_count)
        # phiCoefflist <- vector(mode = "list", length = input$mr_dsem2_g_resp_count)
        # lambdaCoefflist <- vector(mode = "list", length = input$mr_dsem2_g_resp_count)
        # 
        for(k in 1:input$mr_dsem2_g_resp_count) {
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
    
            if (!is.null(fittedModel[[k]]$beta_coeff) && !is.null(input[[paste0('mr_dsem2_m_check_exp_', k)]]) && input[[paste0('mr_dsem2_m_check_exp_', k)]]) {
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
    
        # ?? Output Likelihood ####
        outputLikelihooddf = 0
        for(k in 1:input$mr_dsem2_g_resp_count) {
            outputLikelihooddf <- outputLikelihooddf + fittedModel[[k]]$df
        }
        
        outputLikelihoodlist <- matrix(c(0, 0, 0, 0,0), nrow = 1)
        colnames(outputLikelihoodlist) <- c("-2ML", "-2RL", "Scaled Deviance", "df", "cAIC")
        for(k in 1:input$mr_dsem2_g_resp_count) {
            outputLikelihoodlist[1] <- outputLikelihoodlist[1] + fittedModel[[k]]$ml
            outputLikelihoodlist[2] <- outputLikelihoodlist[2] + fittedModel[[k]]$rl
            outputLikelihoodlist[3] <- outputLikelihoodlist[3] + fittedModel[[k]]$scaled_dv
            outputLikelihoodlist[4] <- outputLikelihoodlist[4] + fittedModel[[k]]$df
            outputLikelihoodlist[5] <- outputLikelihoodlist[5] + fittedModel[[k]]$caic
        }
        
        # ?? Option ####
        fittedModel$option <-list(
            RespDist = input_respdist,
            DataMain = input_datamain,
            structure = input$mr_dsem2_g_corr_structure,
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
            RespCount = input$mr_dsem2_g_resp_count,
            Rcodes = Rcodes
        )
        
#        appendTab(inputId = "mr_dsem2_resulttabset",
#            tabPanel(
#                "Model Checking Plot",
#                br(),
#                mainPanel(
#                    uiOutput("mr_dsem2_g_mcplots"),
#                    width = 12
#                )
#            )
#        )
        
        # ?? Comparison Model ####
        dsem2Comparisonmodel = NULL
        
        RMSEA=(outputLikelihoodlist[5]-2500)/2500
        CTI=abs(3000-outputLikelihoodlist[1])/250
        TLI=abs(3000-outputLikelihoodlist[2])/250
        
        if (!is.null(input$mr_dsem2_m_check_comparison_1) && input$mr_dsem2_m_check_comparison_1) {
            # g_mr_dsem2_r_comparisonmodelcount <<- g_mr_dsem2_r_comparisonmodelcount + 1
            dsem2Comparisonmodel <- matrix(
                c(
                    as.character(round(outputLikelihoodlist[1], digits=5)),
                    as.character(round(outputLikelihoodlist[2], digits=5)),
                    as.character(round(outputLikelihoodlist[3], digits=5)),
                    as.character(round(outputLikelihoodlist[4], digits=5)),
                    as.character(round(outputLikelihoodlist[5], digits=5)),
                    as.character(round(CTI,digits=5)),
                    as.character(round(TLI,digits=5)),
                    as.character(round(RMSEA,digits=5))
                ),
                ncol = 8
            )
          
            dsem2Comparisonmodel <- rbind(g_mr_dsem2_r_comparisonmodel_1, dsem2Comparisonmodel)
            colnames(dsem2Comparisonmodel) <- c("-2ML", "-2RL", "Scaled Deviance", "df", "cAIC","CTI", "TLI","RMSEA")
            g_mr_dsem2_r_comparisonmodel_1 <<- dsem2Comparisonmodel
            fittedModel$option$Comparisonmodel = dsem2Comparisonmodel
            
            # g_mr_dsem2_r_comparisonmodel_2[[g_mr_dsem2_r_comparisonmodelcount]] <<- " " # modelDesclist
            # fittedModel$option$ComparisonmodelDesc = g_mr_dsem2_r_comparisonmodel_2
            # fittedModel$option$Comparisonmodelcount = g_mr_dsem2_r_comparisonmodelcount
            
            # g_mr_dsem2_r_comparisonmodel_3[[g_mr_dsem2_r_comparisonmodelcount]] <<- " " # matrix(input$mr_dsem2_g_corr_structure)
            # colnames(g_mr_dsem2_r_comparisonmodel_3[[g_mr_dsem2_r_comparisonmodelcount]]) <<- " " # "Correlation Structure"
            # fittedModel$option$ComparisonmodelCorr = g_mr_dsem2_r_comparisonmodel_3
            
            appendTab(inputId = "mr_dsem2_resulttabset",
                tabPanel(
                    "Model Comparison",
                    br(),
                    #uiOutput("mr_dsem2_r_comparisonmodeldesc"),
                    tableOutput("mr_dsem2_r_comparisonmodel")
                )
            )
        }
        
        ## Correlation Plots ##
        
        yy<-matrix(0,fittedModel[[1]]$n,input$mr_dsem2_g_resp_count)
        random_counts<-0
        for (k in 1:input$mr_dsem2_g_resp_count) {
          yy[,k]<-fittedModel[[k]]$y
          if (is.null(fittedModel[[k]]$phi_v_h)) random_counts<-random_counts+fittedModel[[k]]$nrand
          else random_counts<-random_counts+fittedModel[[k]]$nrand+1
        }
        fittedModel$Plot01 <- ggcorrplot::ggcorrplot(cor(yy))
        vv_hh <-matrix(0,fittedModel[[k]]$q[1],random_counts)
        kk<-1
        for (k in 1:input$mr_dsem2_g_resp_count) {
          start<-1
          end<-fittedModel[[k]]$q[1]
          for (j in 1:fittedModel[[k]]$nrand) {
            vv_hh[,kk] <- fittedModel[[k]]$v_h[start:end]
            start<-fittedModel[[k]]$q[1]*j+1
            end<-fittedModel[[k]]$q[1]*(j+1)
            kk<-kk+1
          }
          if (!is.null(fittedModel[[k]]$phi_v_h)) {
            vv_hh[,kk]<-fittedModel[[k]]$phi_v_h
            kk<-kk+1
          }
        }
        fittedModel$Plot02 <- ggcorrplot::ggcorrplot(cor(vv_hh))

        # REMOVE ####
        toto <<- fittedModel
        
        # Set Progressbar
        setProgress(1, detail = "Finish")
        
        # colnames(vv_hh)=c("factor11","factor12","factor13","factor14","factor21","factor22","factor23","factor24")
        f<-cbind(f,vv_hh)
        fittedModel$SEM=NULL
        # if (input$mr_dsem2_s_model_1 =="~") fittedModel$SEM=NULL
        # else if (! is.null(input$mr_dsem2_s_model_1)) {
        #  fittedModel$SEM=lm(input$mr_dsem2_s_model_1,f)
        #}
        #print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
        #print(fittedModel$SEM)
        return(fittedModel)    
    }) # End Progressbar
})

output$mr_sem_resp_select<-renderUI({
  response_name=c("Response1","Response2","Response3","Response4","Response5","Response6")
  response_name=response_name[1:mr_dsem2_g_resp_value()]
  selectInput("mr_dsem2_resp_select","",choices =as.list(response_name),selected = "Response1", multiple = FALSE)
})


output$mr_dsem2_g_corr_structure_mean_1<-renderUI({
  names=c("factor11","factor12","factor13","factor21","factor22","factor23")
  selectInput("mr_dsem2_g_corr_structure_mean_1","Shared ",choices =as.list(names), multiple = TRUE)
})


output$mr_dsem2_g_corr_structure_mean_2<-renderUI({
  result=mr_dsem2_output_mean_1()
  names=c("factor11","factor12","factor13","factor21","factor22","factor23")
  names1=paste0(names,"")
  names2=names
  if (!is.null(result)==TRUE) {
    txt1=strsplit(mr_dsem2_Model_Sel_1,"\n")[[1]]
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
  selectInput("mr_dsem2_g_corr_structure_mean_2","Saturated",choices =as.list(names2), multiple = TRUE)
})

mr_dsem2_output_mean_1 <-eventReactive(input$mr_dsem2_g_corr_structure_mean_1_run, {
  select1=paste(input$mr_dsem2_g_corr_structure_mean_1, collapse=",")
  if (!is.null(select1)) {
    if (is.null(mr_dsem2_Model_Sel_1)) mr_dsem2_Model_Sel_1 <<- select1
    else mr_dsem2_Model_Sel_1 <<- paste0(mr_dsem2_Model_Sel_1,"\n",select1)
  } else {
    mr_dsem2_Model_Sel_1 <<- mr_dsem2_Model_Sel_1
  }
  return(mr_dsem2_Model_Sel_1)
})

mr_dsem2_output_mean_2 <-eventReactive(input$mr_dsem2_g_corr_structure_mean_2_run, {
  select2=paste(input$mr_dsem2_g_corr_structure_mean_2, collapse=",")
  if (is.null(mr_dsem2_Model_Sel_2)) mr_dsem2_Model_Sel_2 <<- select2
  else mr_dsem2_Model_Sel_2 <<- paste0(mr_dsem2_Model_Sel_2,"\n",select2)
  return(mr_dsem2_Model_Sel_2)
})

output$mr_dsem2_g_corr_structure_mean_11 <- renderUI({
  Height = paste0(80, "px")
  value=mr_dsem2_output_mean_1()
  textAreaInput("mr_dsem2_g_corr_structure_mean_11", "", value = value, height = Height)
})

output$mr_dsem2_g_corr_structure_mean_21 <- renderUI({
  Height = paste0(80, "px")
  value=mr_dsem2_output_mean_2()
  textAreaInput("mr_dsem2_g_corr_structure_mean_21", "", value = value, height = Height)
})

output$mr_dsem2_g_corr_structure_phi_1<-renderUI({
  names=c("factor11","factor21")
  selectInput("mr_dsem2_g_corr_structure_phi_1","Shared ",choices =as.list(names), multiple = TRUE)
})


output$mr_dsem2_g_corr_structure_phi_2<-renderUI({
  result=mr_dsem2_output_phi_1()
  names=c("factor11","factor21")
  names1=paste0(names,"")
  names2=names
  if (!is.null(result)==TRUE) {
    txt1=strsplit(mr_dsem2_Model_Sel_1_phi,"\n")[[1]]
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
  selectInput("mr_dsem2_g_corr_structure_phi_2","Saturated",choices =as.list(names2), multiple = TRUE)
})

mr_dsem2_output_phi_1 <-eventReactive(input$mr_dsem2_g_corr_structure_phi_1_run, {
  select1=paste(input$mr_dsem2_g_corr_structure_phi_1, collapse=",")
  if (!is.null(select1)) {
    if (is.null(mr_dsem2_Model_Sel_1_phi)) mr_dsem2_Model_Sel_1_phi <<- select1
    else mr_dsem2_Model_Sel_1_phi <<- paste0(mr_dsem2_Model_Sel_1_phi,"\n",select1)
  } else {
    mr_dsem2_Model_Sel_1_phi <<- mr_dsem2_Model_Sel_1_phi
  }
  return(mr_dsem2_Model_Sel_1_phi)
})

mr_dsem2_output_phi_2 <-eventReactive(input$mr_dsem2_g_corr_structure_phi_2_run, {
  select2=paste(input$mr_dsem2_g_corr_structure_phi_2, collapse=",")
  if (is.null(mr_dsem2_Model_Sel_2)) mr_dsem2_Model_Sel_2 <<- select2
  else mr_dsem2_Model_Sel_2 <<- paste0(mr_dsem2_Model_Sel_2,"\n",select2)
  return(mr_dsem2_Model_Sel_2)
})

output$mr_dsem2_g_corr_structure_phi_11 <- renderUI({
  Height = paste0(80, "px")
  value=mr_dsem2_output_phi_1()
  textAreaInput("mr_dsem2_g_corr_structure_phi_11", "", value = value, height = Height)
})

output$mr_dsem2_g_corr_structure_phi_21 <- renderUI({
  Height = paste0(80, "px")
  value=mr_dsem2_output_phi_2()
  textAreaInput("mr_dsem2_g_corr_structure_phi_21", "", value = value, height = Height)
})

output$mr_dsem2_g_corr_structure_lambda_1<-renderUI({
  names=c("factor11","factor21")
  selectInput("mr_dsem2_g_corr_structure_lambda_1","Shared ",choices =as.list(names), multiple = TRUE)
})


output$mr_dsem2_g_corr_structure_lambda_2<-renderUI({
  result=mr_dsem2_output_lambda_1()
  names=c("factor11","factor21")
  names1=paste0(names,"")
  names2=names
  if (!is.null(result)==TRUE) {
    txt1=strsplit(mr_dsem2_Model_Sel_1_lambda,"\n")[[1]]
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
  selectInput("mr_dsem2_g_corr_structure_lambda_2","Saturated",choices =as.list(names2), multiple = TRUE)
})

mr_dsem2_output_lambda_1 <-eventReactive(input$mr_dsem2_g_corr_structure_lambda_1_run, {
  select1=paste(input$mr_dsem2_g_corr_structure_lambda_1, collapse=",")
  if (!is.null(select1)) {
    if (is.null(mr_dsem2_Model_Sel_1_lambda)) mr_dsem2_Model_Sel_1_lambda <<- select1
    else mr_dsem2_Model_Sel_1_lambda <<- paste0(mr_dsem2_Model_Sel_1_lambda,"\n",select1)
  } else {
    mr_dsem2_Model_Sel_1_lambda <<- mr_dsem2_Model_Sel_1_lambda
  }
  return(mr_dsem2_Model_Sel_1_lambda)
})

mr_dsem2_output_lambda_2 <-eventReactive(input$mr_dsem2_g_corr_structure_lambda_2_run, {
  select2=paste(input$mr_dsem2_g_corr_structure_lambda_2, collapse=",")
  if (is.null(mr_dsem2_Model_Sel_2)) mr_dsem2_Model_Sel_2 <<- select2
  else mr_dsem2_Model_Sel_2 <<- paste0(mr_dsem2_Model_Sel_2,"\n",select2)
  return(mr_dsem2_Model_Sel_2)
})

output$mr_dsem2_g_corr_structure_lambda_11 <- renderUI({
  Height = paste0(80, "px")
  value=mr_dsem2_output_lambda_1()
  textAreaInput("mr_dsem2_g_corr_structure_lambda_11", "", value = value, height = Height)
})

output$mr_dsem2_g_corr_structure_lambda_21 <- renderUI({
  Height = paste0(80, "px")
  value=mr_dsem2_output_lambda_2()
  textAreaInput("mr_dsem2_g_corr_structure_lambda_21", "", value = value, height = Height)
})


output$mr_dsem2_accordion_c <- renderUI({
  input$file
  input$resetData
  input$di_option_run
  input$dm_mv_run
  input$dm_ct_run
  input$dm_sr_run
  input$dm_md_run
  dsem2Accordion_c <- bs_accordion_sidebar(
    id = "mr_dsem2_accordion_c",
    spec_side = c(width = 3, offset = 0),
    spec_main = c(width = 9, offset = 0)    
  )
  dsem2Accordion_c <- dsem2Accordion_c %>%
    bs_append(
      title_side = "Mean",
      content_side = NULL,
      content_main = div(
        fluidRow(
          column(
            6,
            uiOutput("mr_dsem2_g_corr_structure_mean_1")
          ),
          column(
            6,
            uiOutput("mr_dsem2_g_corr_structure_mean_2")
          )
        ),
        fluidRow(
          column(
            6,
            splitLayout(style = "padding: 6px;"),
            actionButton("mr_dsem2_g_corr_structure_mean_1_run", "Add", icon = icon("play")),
            style = "text-align:center;"
          ) ,
          column(
            6,
            splitLayout(style = "padding: 6px;"),
            actionButton("mr_dsem2_g_corr_structure_mean_2_run", "Add", icon = icon("play")),
            style = "text-align:center;"
          ) 
        ),
        fluidRow(
          column(
            6,
            uiOutput("mr_dsem2_g_corr_structure_mean_11")
          ),
          column(
            6,
            uiOutput("mr_dsem2_g_corr_structure_mean_21")
          )
        ),
        
      )
    )
    dsem2Accordion_c <- dsem2Accordion_c %>%
    bs_append(
    title_side = "Phi",
    content_side = NULL,
    content_main = div(
      fluidRow(
        column(
          6,
          uiOutput("mr_dsem2_g_corr_structure_phi_1")
        ),
        column(
          6,
          uiOutput("mr_dsem2_g_corr_structure_phi_2")
        )
      ),
      fluidRow(
        column(
          6,
          splitLayout(style = "padding: 6px;"),
          actionButton("mr_dsem2_g_corr_structure_phi_1_run", "Add", icon = icon("play")),
          style = "text-align:center;"
        ) ,
        column(
          6,
          splitLayout(style = "padding: 6px;"),
          actionButton("mr_dsem2_g_corr_structure_phi_2_run", "Add", icon = icon("play")),
          style = "text-align:center;"
        ) 
      ),
      fluidRow(
        column(
          6,
          uiOutput("mr_dsem2_g_corr_structure_phi_11")
        ),
        column(
          6,
          uiOutput("mr_dsem2_g_corr_structure_phi_21")
        )
      ),
      
    )
  )
    dsem2Accordion_c <- dsem2Accordion_c %>%
      bs_append(
        title_side = "Lambda",
        content_side = NULL, 
        content_main = div(
          fluidRow(
            column(
              6,
              uiOutput("mr_dsem2_g_corr_structure_lambda_1")
            ),
            column(
              6,
              uiOutput("mr_dsem2_g_corr_structure_lambda_2")
            )
          ),
          fluidRow(
            column(
              6,
              splitLayout(style = "padding: 6px;"),
              actionButton("mr_dsem2_g_corr_structure_lambda_1_run", "Add", icon = icon("play")),
              style = "text-align:center;"
            ) ,
            column(
              6,
              splitLayout(style = "padding: 6px;"),
              actionButton("mr_dsem2_g_corr_structure_lambda_2_run", "Add", icon = icon("play")),
              style = "text-align:center;"
            ) 
          ),
          fluidRow(
            column(
              6,
              uiOutput("mr_dsem2_g_corr_structure_lambda_11")
            ),
            column(
              6,
              uiOutput("mr_dsem2_g_corr_structure_lambda_21")
            )
          ),
          
        )
        
      )
    div(
      dsem2Accordion_c,
      use_bs_tooltip(),
      use_bs_accordion_sidebar() # needs to be at end, for some reason
    )
})

# MDHGLM tabpanel ####
output$mr_dsem2_g_tabpanel <- renderUI({
#    j=1
#    response_name=c("Response1","Response2","Response3","Response4","Response5","Response6")
#    response_name=response_name[1:mr_dsem2_g_resp_value()]
#    for (jj in 1:mr_dsem2_g_resp_value()) if (!is.null(input$mr_dsem2_resp_select) && input$mr_dsem2_resp_select==response_name[jj]) j=jj
#    do.call(tabsetPanel, c(id='dsem2_tabp',lapply(1:4,function(i) {
    do.call(tabsetPanel, c(id = 'dsem2_tabp', lapply(1:(mr_dsem2_g_resp_value() + 3), function(i) {      
#      if (i==1) {
#        tabPanel(
#          title = paste0('Response'),
#          br(),
#          uiOutput(paste0("mr_dsem2_m_model_", j)),
#          uiOutput(paste0("mr_dsem2_p_model_", j)),
#          uiOutput(paste0("mr_dsem2_l_model_", j)),
#          uiOutput(paste0("mr_dsem2_r_accordion_", j)))
      if (i <= mr_dsem2_g_resp_value()) {
        tabPanel(
          title = paste0('Response', i),
          br(),
          uiOutput(paste0("mr_dsem2_m_model_", i)),
          uiOutput(paste0("mr_dsem2_p_model_", i)),
          uiOutput(paste0("mr_dsem2_l_model_", i)),
          uiOutput(paste0("mr_dsem2_r_accordion_", i))
        )
      } else if (i==mr_dsem2_g_resp_value()+1) {
        tabPanel(
          title = paste0('Correlation Structure'),
          br(),
          uiOutput("mr_dsem2_accordion_c")
        )
      } else if (i==mr_dsem2_g_resp_value()+2) {
        tabPanel(
          title = paste0('Structural Model'),
          br(),
          uiOutput(paste0("mr_dsem2_s_model_", 1)),
          h3(strong("Structural Model")),
          uiOutput(paste0("mr_dsem2_s_resp_", 1)),
          uiOutput(paste0("mr_dsem2_s_variable_", 1)),
          fluidRow(
            column(
              8,
              uiOutput(paste0("mr_dsem2_s_interaction_", 1))
            ),
            column(
              4,
              uiOutput(paste0("mr_dsem2_s_interactionappend_", 1)),
              style = "text-align:right; padding:15px"
            )
          )
        )
      } else {
        tabPanel(
          title = paste0('Settings'),
          h3(strong("Additional Settings")),
          br(),
          uiOutput("mr_dsem2_a_mord"),
          uiOutput("mr_dsem2_a_dord"),
          uiOutput("mr_dsem2_a_order"),
          uiOutput("mr_dsem2_a_reml")
        )
      }
      })))
#    do.call(tabsetPanel, c(id = 'dsem2_tabp', lapply(1:(1 + 2), function(i) {
#        if (i <= 1) {
#            tabPanel(
#                  title = paste0('Response'),
#                  br(),
#                  uiOutput(paste0("mr_dsem2_m_model_", j)),
#                  uiOutput(paste0("mr_dsem2_p_model_", j)),
#                  uiOutput(paste0("mr_dsem2_l_model_", j)),
#                  uiOutput(paste0("mr_dsem2_r_accordion_", j)))
#      })))
})

# MDHGLM Accordion ####

lapply(1:6, function(k) {  
    output[[paste0('mr_dsem2_r_accordion_',k)]]<-renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        dsem2Accordion <- bs_accordion_sidebar(
            id = paste0("mr_dsem2_r_accordion_", k),
            spec_side = c(width = 3, offset = 0),
            spec_main = c(width = 9, offset = 0)    
        )
        
        dsem2Accordion <- dsem2Accordion %>%
            bs_append(
                title_side = "Mean",
                content_side = NULL,
                content_main = div(
                    h3(strong("Model for Mean")),
                    uiOutput(paste0("mr_dsem2_m_resp_", k)),
                    uiOutput(paste0("mr_dsem2_m_variable_", k)),
                    fluidRow(
                        column(
                            8,
                            uiOutput(paste0("mr_dsem2_m_interaction_", k))
                        ),
                        column(
                            4,
                            uiOutput(paste0("mr_dsem2_m_interactionappend_", k)),
                            style = "text-align:right; padding:15px"
                        )
                    ),
                    uiOutput(paste0("mr_dsem2_m_rand_", k)),
                    uiOutput(paste0("mr_dsem2_m_check_slope_", k)),
                    uiOutput(paste0("mr_dsem2_m_check_withoutslope_", k)),
                    uiOutput(paste0("mr_dsem2_m_slope_", k)),
                    uiOutput(paste0("mr_dsem2_m_check_randinteraction_", k)),
                    fluidRow(
                        column(
                            8,
                            uiOutput(paste0("mr_dsem2_m_randinteraction_", k))
                        ),
                        column(
                            4,
                            uiOutput(paste0("mr_dsem2_m_randinteractionappend_", k)),
                            style = "text-align:right; padding:15px"
                        )
                    ), #new components
                    uiOutput(paste0("mr_dsem2_m_randfamily_", k)),
                    hr(style = "border-color: #2C3E50;"),
                    uiOutput(paste0("mr_dsem2_m_dist_", k)),
                    # uiOutput(paste0("mr_dsem2_m_check_binomd_", k)),
                    # uiOutput(paste0("mr_dsem2_m_binomd_", k)),
                    uiOutput(paste0("mr_dsem2_m_link_", k)),
                    uiOutput(paste0("mr_dsem2_m_check_nointercept_", k)),
                    uiOutput(paste0("mr_dsem2_m_check_offset_", k)),
                    uiOutput(paste0("mr_dsem2_m_offset_", k)),
                    uiOutput(paste0("mr_dsem2_m_factor_", k)),
                    hr(style = "border-color: #2C3E50;"),
                    uiOutput(paste0("mr_dsem2_m_check_rcodes_", k)),
                    uiOutput(paste0("mr_dsem2_m_check_exp_", k))
                )
            )
        
        dsem2Accordion <- dsem2Accordion %>%
        bs_append(
            title_side = "Phi",
            content_side = uiOutput(paste0("mr_dsem2_p_check_phi_", k)),
            content_main = div(
                h3(strong("Model for phi")),
                
                uiOutput(paste0("mr_dsem2_p_resp_", k)),
                uiOutput(paste0("mr_dsem2_p_variable_", k)),
                fluidRow(
                    column(
                        8,
                        uiOutput(paste0("mr_dsem2_p_interaction_", k))
                    ),
                    column(
                        4,
                        uiOutput(paste0("mr_dsem2_p_interactionappend_", k)),
                        style = "text-align:right; padding:15px"
                    )
                ),
                uiOutput(paste0("mr_dsem2_p_rand_", k)),
                uiOutput(paste0("mr_dsem2_p_check_slope_", k)),
                uiOutput(paste0("mr_dsem2_p_check_withoutslope_", k)),
                uiOutput(paste0("mr_dsem2_p_slope_", k)),
                uiOutput(paste0("mr_dsem2_p_check_randinteraction_", k)),
                fluidRow(
                    column(
                        8,
                        uiOutput(paste0("mr_dsem2_p_randinteraction_", k))
                    ),
                    column(
                        4,
                        uiOutput(paste0("mr_dsem2_p_randinteractionappend_", k)),
                        style = "text-align:right; padding:15px"
                    )
                ), #new components
                uiOutput(paste0("mr_dsem2_p_randfamily_", k)),
                uiOutput(paste0("mr_dsem2_p_link_", k)),
                uiOutput(paste0("mr_dsem2_p_check_offset_", k)),
                uiOutput(paste0("mr_dsem2_p_offset_", k))
            )
        )
        
        dsem2Accordion <- dsem2Accordion %>%
        bs_append(
            title_side = "Lambda",
            content_side = uiOutput(paste0("mr_dsem2_l_check_lambda_", k)),
            content_main = div(
                h3(strong("Model for Lambda")),
                
                uiOutput(paste0("mr_dsem2_l_resp_", k)),
                uiOutput(paste0("mr_dsem2_l_variable_", k)),
                fluidRow(
                    column(
                        8,
                        uiOutput(paste0("mr_dsem2_l_interaction_", k))
                    ),
                    column(
                        4,
                        uiOutput(paste0("mr_dsem2_l_interactionappend_", k)),
                        style = "text-align:right; padding:15px"
                    )
                ),
                uiOutput(paste0("mr_dsem2_l_rand_", k)),
                uiOutput(paste0("mr_dsem2_l_check_randinteraction_", k)),
                fluidRow(
                    column(
                        8,
                        uiOutput(paste0("mr_dsem2_l_randinteraction_", k))
                    ),
                    column(
                        4,
                        uiOutput(paste0("mr_dsem2_l_randinteractionappend_", k)),
                        style = "text-align:right; padding:15px"
                    )
                ), #new components
                uiOutput(paste0("mr_dsem2_l_randfamily_", k)),
                uiOutput(paste0("mr_dsem2_l_link_", k))
            )
        )
        
        div(
            dsem2Accordion,
            use_bs_tooltip(),
            use_bs_accordion_sidebar() # needs to be at end, for some reason
        )
    })
})

# MDHGLM output ####
lapply(1:6, function(k) {
    output[[paste0('mr_dsem2_r_outputmainsummary_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()$option$outputModelDesc[[k]]))
            return()
      
        mr_dsem2_results()$option$outputModelDesc[[k]]
    }, rownames = TRUE, bordered = TRUE, caption = paste0("Model Description ", k), spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    output[[paste0('mr_dsem2_r_mainsummary_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()$option$ModelDesc[[k]]))
            return()
      
        mr_dsem2_results()$option$ModelDesc[[k]]
    }, rownames = TRUE, bordered = TRUE, caption = "Model Description", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    output[[paste0('mr_dsem2_m_coeff_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()[[k]]$beta_coeff))
            return()
      
        mr_dsem2_results()[[k]]$beta_coeff
    }, rownames = TRUE, bordered = TRUE, caption = "Estimate from Mean Model", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    output[[paste0('mr_dsem2_l_coeff_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()[[k]]$lambda_coeff))
            return()
      
        mr_dsem2_results()[[k]]$lambda_coeff
    }, rownames = TRUE, bordered = TRUE, caption = "Estimate for log(Lambda)", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    output[[paste0('mr_dsem2_l_taucoeff_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()[[k]]$alpha_coeff))
            return()
      
        mr_dsem2_results()[[k]]$alpha_coeff
    }, rownames = TRUE, bordered = TRUE, caption = "Estimate for log(Tau)", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)    

    output[[paste0('mr_dsem2_p_coeff_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()[[k]]$phi_coeff))
            return()
      
        mr_dsem2_results()[[k]]$phi_coeff
    }, rownames = TRUE, bordered = TRUE, caption = "Estimate from Dispersion Model", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    output[[paste0('mr_dsem2_p_alphacoeff_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()[[k]]$tau_coeff))
            return()
      
        mr_dsem2_results()[[k]]$tau_coeff
    }, rownames = TRUE, bordered = TRUE, caption = "Estimate for log(Alpha)", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    output[[paste0('mr_dsem2_r_likelihood_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()$option$Likelihood[[k]]))
            return()
      
        mr_dsem2_results()$option$Likelihood[[k]]
    }, rownames = FALSE, bordered = TRUE, caption = "Likelihood", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    
    output[[paste0('mr_dsem2_r_sbeta_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()$sbeta[[k]]))
            return()
      
        mr_dsem2_results()$sbeta[[k]]
    }, rownames = TRUE, bordered = TRUE, caption = "Beta Coefficient", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    output[[paste0('mr_dsem2_r_sphi_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()$sphi[[k]]))
            return()
      
        mr_dsem2_results()$sphi[[k]]
    }, rownames = FALSE, bordered = TRUE, caption = "Phi Coefficient", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    output[[paste0('mr_dsem2_r_slambda_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()$slambda[[k]]))
            return()
      
        mr_dsem2_results()$slambda[[k]]
    }, rownames = FALSE, bordered = TRUE, caption = "Lambda Coefficient", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    output[[paste0('mr_dsem2_r_rcodes_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()$option$Rcodes[[k]]))
            return()
      
        mr_dsem2_results()$option$Rcodes[[k]]
    }, rownames = TRUE, bordered = TRUE, caption = "R Codes", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
})

output$mr_dsem2_r_cor <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$mr_dsem2_g_run
    
    if (g_resetresult == FALSE || is.null(mr_dsem2_results()$cor))
        return()
  
    mr_dsem2_results()$cor
}, rownames = TRUE, colnames = TRUE, bordered = TRUE, caption = "Estimate for Correlation", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_dsem2_r_outputlikelihood <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$mr_dsem2_g_run
    
    if (g_resetresult == FALSE || is.null(mr_dsem2_results()$option$outputLikelihood))
        return()
  
    mr_dsem2_results()$option$outputLikelihood
}, rownames = FALSE, bordered = TRUE, caption = "Likelihood", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_dsem2_r_sharedparameter <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$mr_dsem2_g_run
    
    if (g_resetresult == FALSE || is.null(mr_dsem2_results()$option$sharedParameter))
        return()
  
    mr_dsem2_results()$option$sharedParameter
}, rownames = TRUE, bordered = TRUE, caption = "Shared Parameter", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_dsem2_r_slikelihood <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$mr_dsem2_g_run
    
    if (g_resetresult == FALSE || is.null(mr_dsem2_results()$option$sLikelihood))
        return()
  
    mr_dsem2_results()$option$sLikelihood
}, rownames = TRUE, bordered = TRUE, caption = "Likelihood", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_dsem2_r_comparisonmodel <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$mr_dsem2_g_run
    
    if (g_resetresult == FALSE || is.null(mr_dsem2_results()$option$Comparisonmodel))
        return()
  
    mr_dsem2_results()$option$Comparisonmodel
}, rownames = TRUE, bordered = TRUE, caption = "Likelihood", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

# MDHGLM Plots ####
lapply(1:6, function(k) {
    # ?? Mean Plot #####
    output[[paste0('mr_dsem2_m_selectplot1_',k)]]<-renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()[[k]]$nrand))
            return()
        
        meanNumberRand = mr_dsem2_results()[[k]]$nrand
        
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
            inputId = sprintf("mr_dsem2_m_selectplot1_%s", k),
            label = "Select Plot", 
            choiceNames = choiceName,
            choiceValues = choiceValue, 
            selected = "mean",
            direction = "vertical"
        )
    })
    
    output[[paste0('mr_dsem2_m_plot1_',k)]] <- renderPlot({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()))
            return()
        
        local({
            output[[paste0('mr_dsem2_m_downloadplot1_',k)]] <- downloadPlot({
                if(is.null(input[[paste0("mr_dsem2_m_selectplot1_", k)]]))
                    plotType = "mean"
                else
                    plotType = input[[paste0("mr_dsem2_m_selectplot1_", k)]]
                
                res <- mr_dsem2_results()[[k]]
                ggplotdhglm(res, type = plotType)
            })
        })
        
        if(is.null(input[[paste0("mr_dsem2_m_selectplot1_", k)]]))
            plotType = "mean"
        else
            plotType = input[[paste0("mr_dsem2_m_selectplot1_", k)]]
        
        res <- mr_dsem2_results()[[k]]
        ggplotdhglm(res, type = plotType)
    })
    
    output[[paste0('mr_dsem2_m_showplot1_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()))
            return()
        
        div(
            h4("Model Checking Plots for Measurement Models"),
            div(
                style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
                withSpinner(
                    plotOutput(paste0('mr_dsem2_m_plot1_', k), height = "600px"),
                    type = 1,
                    color = "#2c3e50",
                    size = 1.2
                ),
                div(
                    style = "position: absolute; left: 0.5em; bottom: 0.5em;",
                    dropdown(
                        uiOutput(paste0('mr_dsem2_m_selectplot1_', k)),
                        size = "xs",
                        icon = icon("chart-bar", class = "opt"), 
                        up = TRUE
                    )
                ),
                div(
                    style = "position: absolute; left:4em; bottom: 0.5em;",
                    dropdown(
                        downloadButton(outputId = paste0('mr_dsem2_m_downloadplot1_', k), label = "Download Plot"),
                        size = "xs",
                        icon = icon("download", class = "opt"),
                        up = TRUE
                    )
                )
            ),
            br()
        )
    })
    
    # ?? Phi Plot #####
    output[[paste0('mr_dsem2_p_selectplot1_',k)]]<-renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()[[k]]$resid_phi))
            return()
        
        choiceName = c("Residuals")
        choiceValue = c("phi")
    
        radioGroupButtons(
            inputId = sprintf("mr_dsem2_p_selectplot1_%s", k),
            label = "Select Plot", 
            choiceNames = choiceName,
            choiceValues = choiceValue, 
            selected = "phi",
            direction = "vertical"
        )
    })
    
    output[[paste0('mr_dsem2_p_plot1_',k)]] <- renderPlot({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        local({
            output[[paste0('mr_dsem2_p_downloadplot1_',k)]] <- downloadPlot({
                if(is.null(input[[paste0("mr_dsem2_p_selectplot1_", k)]]))
                    plotType = "phi"
                else
                    plotType = input[[paste0("mr_dsem2_p_selectplot1_", k)]]
                
                res <- mr_dsem2_results()[[k]]
                ggplotdhglm(res, type = plotType)
            })
        })
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()[[k]]$resid_phi))
            return()
        
        if(is.null(input[[paste0("mr_dsem2_p_selectplot1_", k)]]))
            plotType = "phi"
        else
            plotType = input[[paste0("mr_dsem2_p_selectplot1_", k)]]
        
        res <- mr_dsem2_results()[[k]]
        ggplotdhglm(res, type = plotType)
    })
    
    output[[paste0('mr_dsem2_p_showplot1_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()[[k]]$resid_phi))
            return()
        
        div(
            h4("Model Checking Plots for Phi"),
            div(
                style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
                withSpinner(
                    plotOutput(paste0('mr_dsem2_p_plot1_', k), height = "600px"),
                    type = 1,
                    color = "#2c3e50",
                    size = 1.2
                ),
                div(
                    style = "position: absolute; left: 0.5em; bottom: 0.5em;",
                    dropdown(
                        uiOutput(paste0('mr_dsem2_p_selectplot1_', k)),
                        size = "xs",
                        icon = icon("chart-bar", class = "opt"), 
                        up = TRUE
                    )
                ),
                div(
                    style = "position: absolute; left:4em; bottom: 0.5em;",
                    dropdown(
                        downloadButton(outputId = paste0('mr_dsem2_p_downloadplot1_', k), label = "Download Plot"),
                        size = "xs",
                        icon = icon("download", class = "opt"),
                        up = TRUE
                    )
                )
            ),
            br()
        )
    })
    
    # ?? Lambda Plot #####
    output[[paste0('mr_dsem2_l_selectplot1_',k)]]<-renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()[[k]]$resid_lambda))
            return()
        
        choiceName = c("Residuals")
        choiceValue = c("lambda")
    
        radioGroupButtons(
            inputId = sprintf("mr_dsem2_l_selectplot1_%s", k),
            label = "Select Plot", 
            choiceNames = choiceName,
            choiceValues = choiceValue, 
            selected = "lambda",
            direction = "vertical"
        )
    })
    
    output[[paste0('mr_dsem2_l_plot1_',k)]] <- renderPlot({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        local({
            output[[paste0('mr_dsem2_l_downloadplot1_',k)]] <- downloadPlot({
                if(is.null(input[[paste0("mr_dsem2_l_selectplot1_", k)]]))
                    plotType = "lambda"
                else
                    plotType = input[[paste0("mr_dsem2_l_selectplot1_", k)]]
                
                res <- mr_dsem2_results()[[k]]
                ggplotdhglm(res, type = plotType)
            })
        })
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()[[k]]$resid_lambda))
            return()
        
        if(is.null(input[[paste0("mr_dsem2_l_selectplot1_", k)]]))
            plotType = "lambda"
        else
            plotType = input[[paste0("mr_dsem2_l_selectplot1_", k)]]
        
        res <- mr_dsem2_results()[[k]]
        ggplotdhglm(res, type = plotType)
    })
    
    output[[paste0('mr_dsem2_l_showplot1_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()[[k]]$resid_lambda))
            return()
        
        div(
            h4("Model Checking Plots for Lambda"),
            div(
                style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
                withSpinner(
                    plotOutput(paste0('mr_dsem2_l_plot1_', k), height = "600px"),
                    type = 1,
                    color = "#2c3e50",
                    size = 1.2
                ),
                div(
                    style = "position: absolute; left: 0.5em; bottom: 0.5em;",
                    dropdown(
                        uiOutput(paste0('mr_dsem2_l_selectplot1_', k)),
                        size = "xs",
                        icon = icon("chart-bar", class = "opt"), 
                        up = TRUE
                    )
                ),
                div(
                    style = "position: absolute; left:4em; bottom: 0.5em;",
                    dropdown(
                        downloadButton(outputId = paste0('mr_dsem2_l_downloadplot1_', k), label = "Download Plot"),
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
    output[[paste0('mr_dsem2_r_comparisoncorr_', k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_dsem2_g_run
        
        if (g_resetresult == FALSE || is.null(mr_dsem2_results()$option$ComparisonmodelCorr[[k]]))
            return()
      
        mr_dsem2_results()$option$ComparisonmodelCorr[[k]]
    }, rownames = FALSE, bordered = TRUE, caption = "", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
  
    lapply(1:6, function(j) {
        output[[paste0('mr_dsem2_r_comparisonsummary_', k, "_", j)]] <- renderTable({
            input$file
            input$resetData
            input$di_option_run
            input$mr_dsem2_g_run
            
            if (g_resetresult == FALSE || is.null(mr_dsem2_results()$option$ComparisonmodelDesc[[k]][[j]]))
                return()
          
            mr_dsem2_results()$option$ComparisonmodelDesc[[k]][[j]]
        }, rownames = TRUE, bordered = TRUE, caption = paste0("Model Description ", j), spacing = "m",
        caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    })
})

output$mr_dsem2_r_comparisonmodeldesc <- renderUI({
    do.call(tabsetPanel, c(id = 'comparisontabp', lapply(1:mr_dsem2_results()$option$Comparisonmodelcount, function(i) {
        tabPanel(
            paste0("Result ", i),
            uiOutput(paste0('mr_dsem2_r_comparisoncorr_', i)),
            uiOutput(paste0('mr_dsem2_r_comparisonsummary_', i, "_", 1)),
            uiOutput(paste0('mr_dsem2_r_comparisonsummary_', i, "_", 2)),
            uiOutput(paste0('mr_dsem2_r_comparisonsummary_', i, "_", 3)),
            uiOutput(paste0('mr_dsem2_r_comparisonsummary_', i, "_", 4)),
            uiOutput(paste0('mr_dsem2_r_comparisonsummary_', i, "_", 5)),
            uiOutput(paste0('mr_dsem2_r_comparisonsummary_', i, "_", 6)),
            hr()
        )
    })))
})

# MDHGLM summary ####
output$mr_dsem2_g_summary <- renderUI({
    input$mr_dsem2_var1_select
    res3<-c("Response 1","Response 2","Response 3","Response 4","Response 5","Response 6")
    res3<-res3[1:mr_dsem2_results()$option$RespCount]
    for (j in 1:mr_dsem2_results()$option$RespCount) {
        if (input$mr_dsem2_var1_select == res3[j]) ii=j
    }
    do.call(tabsetPanel, c(id = 'summ', lapply(0:0, function(i) {
             tabPanel(
                uiOutput(paste0('mr_dsem2_r_mainsummary_', ii)),
                
                # except Shared
                uiOutput(paste0('mr_dsem2_m_coeff_', ii)),
                uiOutput(paste0('mr_dsem2_l_coeff_', ii)),
                uiOutput(paste0('mr_dsem2_l_taucoeff_', ii)),
                uiOutput(paste0('mr_dsem2_p_coeff_', ii)),
                uiOutput(paste0('mr_dsem2_p_alphacoeff_', ii)),
                uiOutput(paste0('mr_dsem2_r_likelihood_', ii)),
                
                # Shared
                uiOutput(paste0('mr_dsem2_r_sbeta_', ii)),
                uiOutput(paste0('mr_dsem2_r_sphi_', ii)),
                uiOutput(paste0('mr_dsem2_r_slambda_', ii)),
                
                uiOutput(paste0('mr_dsem2_r_rcodes_', ii)),
                # except shared
            #    uiOutput("mr_dsem2_r_outputmainsummary_1"),
            #    uiOutput("mr_dsem2_r_outputmainsummary_2"),
            #    uiOutput("mr_dsem2_r_outputmainsummary_3"),
            #    uiOutput("mr_dsem2_r_outputmainsummary_4"),
            #    uiOutput("mr_dsem2_r_outputmainsummary_5"),
            #    uiOutput("mr_dsem2_r_outputmainsummary_6"),
                uiOutput("mr_dsem2_r_cor"),
                uiOutput("mr_dsem2_r_outputlikelihood"),
                
                # shared
                uiOutput("mr_dsem2_r_sharedparameter"),
                uiOutput("mr_dsem2_r_slikelihood")
                #               paste0("Response ", ii),
                
            )
    })))
})

output$mr_dsem2_var1_select<-renderUI({
  res3<-c("Response 1","Response 2","Response 3","Response 4","Response 5","Response 6")
  counts=mr_dsem2_results()$option$RespCount+1
  res3<-res3[1:counts]
  selectInput("mr_dsem2_var1_select","",choices =as.list(res3),selected = "Response 1", multiple = FALSE)
})

output$mr_dsem2_var2_select<-renderUI({
  res3<-c("Overall","Response 1","Response 2","Response 3","Response 4","Response 5","Response 6")
  counts=mr_dsem2_results()$option$RespCount+1
  res3<-res3[1:counts]
  selectInput("mr_dsem2_var2_select","",choices =as.list(res3),selected = "Overall", multiple = FALSE)
})


output$mr_dsem2_g_mcplots <- renderUI({
    input$mr_dsem2_var2_select
    res3<-c("Response 1","Response 2","Response 3","Response 4","Response 5","Response 6")
    res3<-res3[1:mr_dsem2_results()$option$RespCount]
    for (j in 1:mr_dsem2_results()$option$RespCount) {
        if (input$mr_dsem2_var2_select == res3[j]) ii=j
        if (input$mr_dsem2_var2_select == "Overall") ii=2
    }
    do.call(tabsetPanel, c(id = 'plot', lapply(1:1, function(i) {
        tabPanel(
            #paste0("Response ", i),
            div(style = "padding:10px"),
            uiOutput(paste0('mr_dsem2_m_showplot1_', ii)),
            uiOutput(paste0('mr_dsem2_p_showplot1_', ii)),
            uiOutput(paste0('mr_dsem2_l_showplot1_', ii))
        )
    })))
})



output$mr_dsem2_r_plot01 <- renderPlot({
  input$file
  input$resetData
  input$di_option_run
  input$mr_dsem2_g_run
  
  if (g_resetresult == FALSE || is.null(mr_dsem2_results()$Plot01))
    return()
  
  local({
    output$mr_dsem2_r_downloadplot01 <- downloadPlot(
      mr_dsem2_results()$Plot01
    )
  })
  mr_dsem2_results()$Plot01
})


output$mr_dsem2_r_showplot01 <- renderUI({
  input$file
  input$resetData
  input$di_option_run
  input$mr_dsem2_g_run
  
  if (g_resetresult == FALSE || is.null(mr_dsem2_results()$Plot01))
    return()
  div(
    style = "position: relative; border: 1px solid #D3D3D3;",
    withSpinner(
      plotOutput("mr_dsem2_r_plot01", height="600px"),
      type = 1,
      color = "#2c3e50",
      size = 1.2
    ),
    div(
      style = "position: absolute; left:0.5em; bottom: 0.5em;",
      dropdown(
        downloadButton(outputId = "mr_dsem2_r_downloadplot01", label = "Download Plot"),
        size = "xs",
        icon = icon("download", class = "opt"),
        up = TRUE
      )
    )
  )    
  
})



output$mr_dsem2_r_plot02 <- renderPlot({
  input$file
  input$resetData
  input$di_option_run
  input$mr_dsem2_g_run
  
  if (g_resetresult == FALSE || is.null(mr_dsem2_results()$Plot02))
    return()
  
  local({
    output$mr_dsem2_r_downloadplot02 <- downloadPlot(
      mr_dsem2_results()$Plot02
    )
  })
  mr_dsem2_results()$Plot02
})


output$mr_dsem2_r_showplot02 <- renderUI({
  input$file
  input$resetData
  input$di_option_run
  input$mr_dsem2_g_run
  
  if (g_resetresult == FALSE || is.null(mr_dsem2_results()$Plot02))
    return()
  div(
    style = "position: relative; border: 1px solid #D3D3D3;",
    withSpinner(
      plotOutput("mr_dsem2_r_plot02", height="600px"),
      type = 1,
      color = "#2c3e50",
      size = 1.2
    ),
    div(
      style = "position: absolute; left:0.5em; bottom: 0.5em;",
      dropdown(
        downloadButton(outputId = "mr_dsem2_r_downloadplot02", label = "Download Plot"),
        size = "xs",
        icon = icon("download", class = "opt"),
        up = TRUE
      )
    )
  )    
  
})


output$mr_dsem2_r_showplot03 <- renderUI({
  input$file
  input$resetData
  input$di_option_run
  input$mr_dsem2_g_run
  
  div(
    style = "position: relative; border: 1px solid #D3D3D3;",
    withSpinner(
      plotOutput("mr_dsem2_r_plot03", height="600px"),
      type = 1,
      color = "#2c3e50",
      size = 1.2
    ),
    div(
      style = "position: absolute; left:0.5em; bottom: 0.5em;",
      dropdown(
        downloadButton(outputId = "mr_dsem2_r_downloadplot03", label = "Download Plot"),
        size = "xs",
        icon = icon("download", class = "opt"),
        up = TRUE
      )
    )
  )    
  
})


output$mr_dsem2_r_showplot04 <- renderUI({
  input$file
  input$resetData
  input$di_option_run
  input$mr_dsem2_g_run
  
  div(
    style = "position: relative; border: 1px solid #D3D3D3;",
    withSpinner(
      plotOutput("mr_dsem2_r_plot04", height="600px"),
      type = 1,
      color = "#2c3e50",
      size = 1.2
    ),
    div(
      style = "position: absolute; left:0.5em; bottom: 0.5em;",
      dropdown(
        downloadButton(outputId = "mr_dsem2_r_downloadplot04", label = "Download Plot"),
        size = "xs",
        icon = icon("download", class = "opt"),
        up = TRUE
      )
    )
  )    
  
})



output$mr_dsem2_r_showplot05 <- renderUI({
  input$file
  input$resetData
  input$di_option_run
  input$mr_dsem2_g_run
  
  div(
    style = "position: relative; border: 1px solid #D3D3D3;",
    withSpinner(
      plotOutput("mr_dsem2_r_plot05", height="600px"),
      type = 1,
      color = "#2c3e50",
      size = 1.2
    ),
    div(
      style = "position: absolute; left:0.5em; bottom: 0.5em;",
      dropdown(
        downloadButton(outputId = "mr_dsem2_r_downloadplot05", label = "Download Plot"),
        size = "xs",
        icon = icon("download", class = "opt"),
        up = TRUE
      )
    )
  )    
  
})


output$mr_dsem2_r_plot03 <- renderPlot({
  if (mr_dsem2_g_resp_value()==2) yy_name=c("y1","y2")
  if (mr_dsem2_g_resp_value()==3) yy_name=c("y1","y2","y3")
  if (mr_dsem2_g_resp_value()==4) yy_name=c("y1","y2","y3","y4")
  if (mr_dsem2_g_resp_value()==5) yy_name=c("y1","y2","y3","y4","y5")
  if (mr_dsem2_g_resp_value()==6) yy_name=c("y1","y2","y3","y4","y5","y6")
  res <- mr_dsem2_results()
  par(mfrow=c(1,1))
  path_diagram1(res,1,yy_name)
})

output$mr_dsem2_r_plot04 <- renderPlot({
  if (mr_dsem2_g_resp_value()==2) yy_name=c("y1","y2")
  if (mr_dsem2_g_resp_value()==3) yy_name=c("y1","y2","y3")
  if (mr_dsem2_g_resp_value()==4) yy_name=c("y1","y2","y3","y4")
  if (mr_dsem2_g_resp_value()==5) yy_name=c("y1","y2","y3","y4","y5")
  if (mr_dsem2_g_resp_value()==6) yy_name=c("y1","y2","y3","y4","y5","y6")
  res <- mr_dsem2_results()
  par(mfrow=c(1,1))
  path_diagram1(res,2,yy_name)
})

output$mr_dsem2_r_plot05 <- renderPlot({
  if (mr_dsem2_g_resp_value()==2) yy_name=c("y1","y2")
  if (mr_dsem2_g_resp_value()==3) yy_name=c("y1","y2","y3")
  if (mr_dsem2_g_resp_value()==4) yy_name=c("y1","y2","y3","y4")
  if (mr_dsem2_g_resp_value()==5) yy_name=c("y1","y2","y3","y4","y5")
  if (mr_dsem2_g_resp_value()==6) yy_name=c("y1","y2","y3","y4","y5","y6")
  res <- mr_dsem2_results()
  par(mfrow=c(1,1))
  path_diagram2(res,yy_name)
})
  