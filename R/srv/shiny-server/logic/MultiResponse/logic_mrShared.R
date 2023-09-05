# Shared Components ####

output$mr_shared_g_respslider <- renderUI({
    sliderInput(
        inputId = sprintf("mr_shared_g_resp_count"),
        label = "Number of Response Variables",
        value = 2,
        min = 2,
        max = 6
    )
})

output$mr_shared_g_corr_structure<-renderUI({
    selectInput(
        "mr_shared_g_corr_structure",
        "Correlation Structure",
        choices = c("shared"),
        selected = "shared",
        multiple = FALSE
    )
})

mr_shared_g_resp_value <- reactive({
    if (is.null(input$mr_shared_g_resp_count))
        return(2)
    else
        return(input$mr_shared_g_resp_count)
})


# ¦§ Mean, Phi, Lambda ####
source('logic/MultiResponse/Shared/logic_mrShared_m.R', local = T)
source('logic/MultiResponse/Shared/logic_mrShared_p.R', local = T)
source('logic/MultiResponse/Shared/logic_mrShared_l.R', local = T)

# ¦¦ Additional Settings ####
output$mr_shared_a_mord<-renderUI({
    selectInput(
        "mr_shared_a_mord",
        "Order for Mean Model",
        choices=c("0"="0","1"="1"),
        multiple=FALSE
    )
})

output$mr_shared_a_dord<-renderUI({
    selectInput(
        "mr_shared_a_dord",
        "Order for Dispersion Model",
        choices=c("1"="1","2"="2"),
        multiple=FALSE
    )
})

output$mr_shared_a_order<-renderUI({
    selectInput(
        "mr_shared_a_order",
        "Order for Laplace's Approximation",
        choices=c("1"="1","2"="2"),
        multiple=FALSE
    )
})

output$mr_shared_a_reml<-renderUI({
    selectInput(
        "mr_shared_a_reml",
        "REML or ML",
        choice = c("REML" = TRUE,"ML" = FALSE),
        multiple=FALSE
    )
})

# Shared Run ####

observeEvent(input$mr_shared_g_run, {
    g_resetresult <<- TRUE
})

mr_shared_results <- eventReactive(input$mr_shared_g_run, {
    # Start Progressbar  
    withProgress(message = 'Shared', style = "notification", value = 0.1, {
          
    # ¦§ Variable Declaration 1 ####
    MM_list <- vector(mode = "list", length = input$mr_shared_g_resp_count)
    DM_list <- vector(mode = "list", length = input$mr_shared_g_resp_count)
    f <- data2
    
    for(k in 1:input$mr_shared_g_resp_count) {
        m_form <- input[[paste0("mr_shared_m_model_", k)]]
        p_form <- NULL
        l_form <- NULL

        if (!is.null(input[[paste0("mr_shared_p_check_phi_", k)]]) && input[[paste0("mr_shared_p_check_phi_", k)]])
            p_form <- input[[paste0("mr_shared_p_model_", k)]]
        else
            p_form <- phi~1

        if (!is.null(input[[paste0("mr_shared_l_check_lambda_", k)]]) && input[[paste0("mr_shared_l_check_lambda_", k)]])
            l_form <- input[[paste0("mr_shared_l_model_", k)]]
        else
            l_form <- lambda~1

        res <- m_form
        m_RandDistM <- NULL
        p_RandDistM <- NULL
        l_RandDistM <- NULL
        m_nRand <- length(input[[paste0("mr_shared_m_rand_", k)]])
        p_nRand <- length(input[[paste0("mr_shared_p_rand_", k)]])
        l_nRand <- length(input[[paste0("mr_shared_l_rand_", k)]])

        if (m_nRand > 0) {
            for (i in 1:m_nRand) {
                m_RandDistM <- c(m_RandDistM,input[[paste0("mr_shared_m_rf_",i,"_", k)]])
            }
        }
        if (!is.null(input[[paste0("mr_shared_p_check_phi_", k)]]) && input[[paste0("mr_shared_p_check_phi_", k)]]) {
            if (p_nRand > 0) {
                for (i in 1:p_nRand) {
                    p_RandDistM <- c(p_RandDistM,input[[paste0("mr_shared_p_rf_",i,"_", k)]])
                }
            }
        }
        LinPredRandVariance=NULL
        if (!is.null(input[[paste0("mr_shared_l_check_lambda_", k)]]) && input[[paste0("mr_shared_l_check_lambda_", k)]]) {
            if (l_nRand > 0) {
                for (i in 1:l_nRand) {
                    l_RandDistM <- c(l_RandDistM,input[[paste0("mr_shared_l_rf_",i,"_", k)]])
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
                Link = isolate(input[[paste0("mr_shared_m_link_", k)]]),
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
                  Link = isolate(input[[paste0("mr_shared_m_link_", k)]]),
                  LinPred = as.formula(m_form),
                  RandDist = m_RandDistM
                )
        }
        
        MM_list[[k]] <- MM
        DM_list[[k]] <- DM
    }
    
    # Increase Progressbar
    incProgress(0.4, detail = paste("Loading..."))
    input_datamain <- vector(mode = "list", length = input$mr_shared_g_resp_count)
    input_respdist = NULL
    input_factor = NULL

    for(k in 1:input$mr_shared_g_resp_count) {
        input_respdist <- c(input_respdist, input[[paste0("mr_shared_m_dist_", k)]])
        input_datamain[[k]] <- f
    }
    
    # ¦§ Variable Declaration 2 ####
    
    mord = NULL
    dord = NULL
    orderLaplace = NULL
    reml = NULL
    if (is.null(input$mr_shared_a_mord))
        mord = 0
    else
        mord = as.numeric(input$re_dhglm_a_mord)
    
    if (is.null(input$mr_shared_a_dord))
        dord = 1
    else
        dord = as.numeric(input$mr_shared_a_dord)
    
    if (is.null(input$mr_shared_a_order))
        orderLaplace = 1
    else
        orderLaplace = as.numeric(input$mr_shared_a_order)
    
    if (is.null(input$mr_shared_a_reml))
        reml = TRUE
    else 
        reml = FALSE
    
    corrStructure = input$mr_shared_g_corr_structure
    
    # ¦§ Fitted Model ####
    
    fittedModel <<- jointfit(
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
    
    # ¦§ Beta_Coefficient ####
    betaCoeff <- cbind(fittedModel$Beta_h, fittedModel$SE_Beta, fittedModel$t_Beta, fittedModel$p_Beta)
    colnames(betaCoeff) <- c("Estimate", "Std. Error", "t_value", "p.value")
    
    # ¦§ Phi_Coefficient ####
    phiCoeff <- cbind(fittedModel$Log_Phi, fittedModel$SE_Log_Phi)
    colnames(phiCoeff) <- c("Estimate", "Std. Error")
    
    # ¦§ Lambda_Coefficient ####
    lambdaCoeff <- cbind(fittedModel$Log_Lambda, fittedModel$SE_Log_Lambda)
    colnames(lambdaCoeff) <- c("Estimate", "Std. Error")
    
    # ¦§ Shared Parameter ####
    sharedParameter <- cbind(fittedModel$Shared, fittedModel$SE_Shared)
    colnames(sharedParameter) <- c("Estimate", "Std. Error")
    
    # ¦§ Likelihood ####
    
    likelihood = matrix(fittedModel$CAIC)
    colnames(likelihood) <- c("cAIC")
    
    # ¦¦ Option ####
    fittedModel$option<<-list(
        betaCoeff = betaCoeff,
        phiCoeff = phiCoeff,
        lambdaCoeff = lambdaCoeff,
        sharedParameter = sharedParameter,
        Likelihood = likelihood
    )
    
    # Set Progressbar
    setProgress(1, detail = "Finish")
    return(fittedModel)
    
    }) # End Progressbar
})

# Shared tabpanel ####
output$mr_shared_g_tabpanel <- renderUI({
    do.call(tabsetPanel, c(id = 'mr_shared_g_tabp', lapply(1:mr_shared_g_resp_value(), function(i) {
        tabPanel(
            title = paste0('Response', i),
            br(),
            uiOutput(paste0("mr_shared_m_model_", i)),
            uiOutput(paste0("mr_shared_p_model_", i)),
            uiOutput(paste0("mr_shared_l_model_", i)),
            uiOutput(paste0("mr_shared_r_accordion_", i))
        )
    })))
})

# Shared Accordion ####

lapply(1:6, function(k) {  
    output[[paste0('mr_shared_r_accordion_',k)]]<-renderUI({
        input$file
        input$resetData
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        sharedAccordion <- bs_accordion_sidebar(
            id = paste0("mr_shared_r_accordion_", k),
            spec_side = c(width = 3, offset = 0),
            spec_main = c(width = 9, offset = 0)    
        )
        
        sharedAccordion <- sharedAccordion %>%
            bs_append(
                title_side = "Mean",
                content_side = NULL,
                content_main = div(
                    h3(strong("Model for Mean")),
                    uiOutput(paste0("mr_shared_m_resp_", k)),
                    uiOutput(paste0("mr_shared_m_variable_", k)),
                    fluidRow(
                        column(
                            8,
                            uiOutput(paste0("mr_shared_m_interaction_", k))
                        ),
                        column(
                            4,
                            uiOutput(paste0("mr_shared_m_interactionappend_", k)),
                            style = "text-align:right; padding:15px"
                        )
                    ),
                    uiOutput(paste0("mr_shared_m_rand_", k)),
                    uiOutput(paste0("mr_shared_m_check_slope_", k)),
                    uiOutput(paste0("mr_shared_m_check_withoutslope_", k)),
                    uiOutput(paste0("mr_shared_m_slope_", k)),
                    uiOutput(paste0("mr_shared_m_check_randinteraction_", k)),
                    fluidRow(
                        column(
                            8,
                            uiOutput(paste0("mr_shared_m_randinteraction_", k))
                        ),
                        column(
                            4,
                            uiOutput(paste0("mr_shared_m_randinteractionappend_", k)),
                            style = "text-align:right; padding:15px"
                        )
                    ), #new components
                    uiOutput(paste0("mr_shared_m_randfamily_", k)),
                    hr(style = "border-color: #2C3E50;"),
                    uiOutput(paste0("mr_shared_m_dist_", k)),
                    # uiOutput(paste0("mr_shared_m_check_binomd_", k)),
                    # uiOutput(paste0("mr_shared_m_binomd_", k)),
                    uiOutput(paste0("mr_shared_m_link_", k)),
                    uiOutput(paste0("mr_shared_m_check_nointercept_", k)),
                    uiOutput(paste0("mr_shared_m_check_offset_", k)),
                    uiOutput(paste0("mr_shared_m_offset_", k)),
                    uiOutput(paste0("mr_shared_m_factor_", k)),
                    hr(style = "border-color: #2C3E50;"),
                    uiOutput(paste0("mr_shared_m_check_exp_", k))
                )
            )
        
        sharedAccordion <- sharedAccordion %>%
        bs_append(
            title_side = "Phi",
            content_side = uiOutput(paste0("mr_shared_p_check_phi_", k)),
            content_main = div(
                h3(strong("Model for phi")),
                
                uiOutput(paste0("mr_shared_p_resp_", k)),
                uiOutput(paste0("mr_shared_p_variable_", k)),
                fluidRow(
                    column(
                        8,
                        uiOutput(paste0("mr_shared_p_interaction_", k))
                    ),
                    column(
                        4,
                        uiOutput(paste0("mr_shared_p_interactionappend_", k)),
                        style = "text-align:right; padding:15px"
                    )
                ),
                uiOutput(paste0("mr_shared_p_rand_", k)),
                uiOutput(paste0("mr_shared_p_check_slope_", k)),
                uiOutput(paste0("mr_shared_p_check_withoutslope_", k)),
                uiOutput(paste0("mr_shared_p_slope_", k)),
                uiOutput(paste0("mr_shared_p_check_randinteraction_", k)),
                fluidRow(
                    column(
                        8,
                        uiOutput(paste0("mr_shared_p_randinteraction_", k))
                    ),
                    column(
                        4,
                        uiOutput(paste0("mr_shared_p_randinteractionappend_", k)),
                        style = "text-align:right; padding:15px"
                    )
                ), #new components
                uiOutput(paste0("mr_shared_p_randfamily_", k)),
                uiOutput(paste0("mr_shared_p_link_", k)),
                uiOutput(paste0("mr_shared_p_check_offset_", k)),
                uiOutput(paste0("mr_shared_p_offset_", k))
            )
        )
        
        sharedAccordion <- sharedAccordion %>%
        bs_append(
            title_side = "Lambda",
            content_side = uiOutput(paste0("mr_shared_l_check_lambda_", k)),
            content_main = div(
                h3(strong("Model for Lambda")),
                
                uiOutput(paste0("mr_shared_l_resp_", k)),
                uiOutput(paste0("mr_shared_l_variable_", k)),
                fluidRow(
                    column(
                        8,
                        uiOutput(paste0("mr_shared_l_interaction_", k))
                    ),
                    column(
                        4,
                        uiOutput(paste0("mr_shared_l_interactionappend_", k)),
                        style = "text-align:right; padding:15px"
                    )
                ),
                uiOutput(paste0("mr_shared_l_rand_", k)),
                uiOutput(paste0("mr_shared_l_check_randinteraction_", k)),
                fluidRow(
                    column(
                        8,
                        uiOutput(paste0("mr_shared_l_randinteraction_", k))
                    ),
                    column(
                        4,
                        uiOutput(paste0("mr_shared_l_randinteractionappend_", k)),
                        style = "text-align:right; padding:15px"
                    )
                ), #new components
                uiOutput(paste0("mr_shared_l_randfamily_", k)),
                uiOutput(paste0("mr_shared_l_link_", k))
            )
        )
        
        if (k == 1) {
            sharedAccordion <- sharedAccordion %>%
            bs_append(
                title_side = "Setting",
                content_side = NULL,
                content_main = div(
                    h3(strong("Additional Settings")),
                    uiOutput("mr_shared_a_mord"),
                    uiOutput("mr_shared_a_dord"),
                    uiOutput("mr_shared_a_order"),
                    uiOutput("mr_shared_a_reml")
                )
            )
        }
        
        div(
            sharedAccordion,
            use_bs_tooltip(),
            use_bs_accordion_sidebar() # needs to be at end, for some reason
        )
    })
})

# Shared Results ####

output$mr_shared_r_betacoeff <- renderTable({
    input$file
    input$resetData
    input$mr_shared_g_run
    
    if (g_resetresult == FALSE || is.null(mr_shared_results()$option$betaCoeff))
        return()
  
    mr_shared_results()$option$betaCoeff
}, rownames = TRUE, bordered = TRUE, caption = "Beta Coefficient", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_shared_r_phicoeff <- renderTable({
    input$file
    input$resetData
    input$mr_shared_g_run
    
    if (g_resetresult == FALSE || is.null(mr_shared_results()$option$phiCoeff))
        return()
  
    mr_shared_results()$option$phiCoeff
}, rownames = TRUE, bordered = TRUE, caption = "Phi Coefficient", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_shared_r_lambdacoeff <- renderTable({
    input$file
    input$resetData
    input$mr_shared_g_run
    
    if (g_resetresult == FALSE || is.null(mr_shared_results()$option$lambdaCoeff))
        return()
  
    mr_shared_results()$option$lambdaCoeff
}, rownames = TRUE, bordered = TRUE, caption = "Lambda Coefficient", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_shared_r_sharedparameter <- renderTable({
    input$file
    input$resetData
    input$mr_shared_g_run
    
    if (g_resetresult == FALSE || is.null(mr_shared_results()$option$sharedParameter))
        return()
  
    mr_shared_results()$option$sharedParameter
}, rownames = TRUE, bordered = TRUE, caption = "Shared Parameter", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_shared_r_likelihood <- renderTable({
    input$file
    input$resetData
    input$mr_shared_g_run
    
    if (g_resetresult == FALSE || is.null(mr_shared_results()$option$Likelihood))
        return()
  
    mr_shared_results()$option$Likelihood
}, rownames = TRUE, bordered = TRUE, caption = "Likelihood", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)




    

