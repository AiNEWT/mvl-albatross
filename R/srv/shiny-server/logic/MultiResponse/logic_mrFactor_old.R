# Factor Slider ####

output$mr_factor_g_respslider <- renderUI({
    sliderInput(
        inputId = sprintf("mr_factor_g_resp_count"),
        label = "Number of Factors",
        value = 2,
        min = 2,
        max = 6
    )
})

mr_factor_g_resp_value <- reactive({
    if (is.null(input$mr_factor_g_resp_count))
        return(2)
    else
        return(input$mr_factor_g_resp_count)
})

# Factor Correlation Structure ####

output$mr_factor_g_corr_structure<-renderUI({
    selectInput(
        "mr_factor_g_corr_structure",
        "Correlation Structure",
        choices = c("factor"),
        selected = "factor",
        multiple = FALSE
    )
})

# Factor Run ####

observeEvent(input$mr_factor_g_run, {
    g_resetresult <<- TRUE
})

mr_factor_results <- eventReactive(input$mr_factor_g_run, {
    withProgress(message = 'Factor Model', style = "notification", value = 0, {
    MM_list <- vector(mode = "list", length = input$mr_factor_g_resp_count)
    f <- data2
    
    # ¦§ MeanModel (MM) ####
    for(k in 1:input$mr_factor_g_resp_count) {
        m_form <- input[[paste0("mr_factor_m_model_", k)]]

        res <- m_form
        m_RandDistM <- NULL
        m_nRand <- length(input[[paste0("mr_factor_m_rand_", k)]])
        comp <- NULL
        if (m_nRand > 0) {
            for (i in 1:m_nRand) {
                m_RandDistM <- c(m_RandDistM,input[[paste0("mr_factor_m_rf_",i,"_", k)]])
            }
        }
        MM <- DHGLMMODELING(
            Model = "mean", 
            Link = input[[paste0("mr_factor_m_link_", k)]], 
            LinPred = as.formula(m_form), 
            RandDist = m_RandDistM
        )
        MM_list[[k]] <- MM
    }
    incProgress(0.5, detail = paste("Loading..."))
    input_datamain <- vector(mode = "list", length = input$mr_factor_g_resp_count)
    input_respdist = NULL
    input_factor = NULL

    # ¦§ Factor, Data Declaration ####
    
    for(k in 1:input$mr_factor_g_resp_count) {
        if (input[['mr_factor_g_corr_structure']] == "factor") {
            input_factor <- as.numeric(c(input_factor, input[[paste0("mr_factor_m_factor_", k)]]))
        }
        input_respdist <- c(input_respdist, input[[paste0("mr_factor_m_dist_", k)]])
        input_datamain[[k]] <- f
    }
    
    corrStructure = input$mr_factor_g_corr_structure
    
    # ¦§ Fitted Model #### 
    
    fittedModel <<- jointfit(
        RespDist = input_respdist,
        DataMain = input_datamain,
        structure = corrStructure,
        MeanModel = MM_list,
        order = as.numeric(input$mr_factor_a_order),
        REML = input$mr_factor_a_reml,
        factor = input_factor
    )

    factorLevels <- levels(as.factor(input_factor))
    lengthOfFactorLevels <- length(factorLevels)    
    
    # You Need to Remove This Code
    Sys.sleep(5)
    
    # ¦§ Lambda ###
    
    lambdaCount <- input$mr_factor_g_resp_count - length(levels(as.factor(input_factor)))
    lambdaResults <- matrix(c(fittedModel$lambda, fittedModel$se_lambda), ncol = lambdaCount, byrow = TRUE)
    rownames(lambdaResults) <- c("Coefficient", "SE")
    lambdaColnames <- NULL
    # Find First Factor Index [ex) c(1, 1, 1, 2, 2, 2) -> 1, 4]
    lambdaIndex <- NULL
    for (i in 1:lengthOfFactorLevels)
        lambdaIndex <- c(lambdaIndex, which(factorLevels[i] == input_factor)[1])
    for (i in 1:input$mr_factor_g_resp_count) {
        if (any(i == lambdaIndex))
            next
        lambdaColnames <- c(lambdaColnames, paste0("\\(\\lambda_{", i, "}\\)"))
    }
    colnames(lambdaResults) <- lambdaColnames
    
    # ¦§ Gamma ####
    
    gammaCount <- lengthOfFactorLevels * (lengthOfFactorLevels + 1) / 2
    gammaResults <- matrix(c(fittedModel$gamma, fittedModel$se_gamma), ncol = gammaCount, byrow = TRUE)
    rownames(gammaResults) <- c("Coefficient", "SE")
    gammaColnames <- NULL
    for (i in 1:lengthOfFactorLevels)
        for (j in i:lengthOfFactorLevels)
            gammaColnames <- c(gammaColnames, paste0("\\(\\gamma_{", i, j, "}\\)"))
    colnames(gammaResults) <- gammaColnames
    
    # ¦§ Beta ####
    
    betaCount <- input$mr_factor_g_resp_count
    betaResults <- matrix(c(fittedModel$beta, fittedModel$se_beta), ncol = betaCount, byrow = TRUE)
    rownames(betaResults) <- c("Coefficient", "SE")
    betaColnames <- NULL
    for (i in 1:betaCount)
        betaColnames <- c(betaColnames, paste0("\\(\\beta_{0", i,"}\\)"))
    colnames(betaResults) <- betaColnames
    
    # ¦¦ Likelihood ####
    
    likelihoodResults <- NULL
    if (is.null(fittedModel$caic)) {
        likelihoodResults <- matrix(c(fittedModel$deviance, fittedModel$df), nrow = 1, ncol = 2)
        colnames(likelihoodResults) <- c("deviance", "df")
    }
    else {
        likelihoodResults <- matrix(c(fittedModel$deviance, fittedModel$df, fittedModel$caic), nrow = 1, ncol = 3)
        colnames(likelihoodResults) <- c("deviance", "df", "cAIC")
    }
    
    fittedModel$option<<-list(
        lambda = lambdaResults,
        gamma = gammaResults,
        beta = betaResults,
        likelihood = likelihoodResults
    )
    
    setProgress(1, detail = "Finish")
    return(fittedModel)
    })
})

# Factor Tab Panel ####
output$mr_factor_g_tabpanel <- renderUI({
    do.call(tabsetPanel, c(id = 'factortabp', lapply(1:mr_factor_g_resp_value(), function(i) {
        tabPanel(
            title = paste0('Factor', i),
            br(),
            uiOutput(paste0('mr_factor_r_accordion_', i))
        )
    })))
})

# Factor Components ####
lapply(1:6, function(k) {  
    output[[paste0('mr_factor_r_accordion_',k)]]<-renderUI({
        input$file
        input$resetData
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        factorAccordion <- bs_accordion_sidebar(
            id = paste0("mr_factor_r_accordion_", k),
            spec_side = c(width = 3, offset = 0),
            spec_main = c(width = 9, offset = 0)    
        )
        
        factorAccordion <- factorAccordion %>%
            bs_append(
                title_side = "Factor",
                content_side = NULL,
                content_main = div(
                    h3(strong("Factor")),
                    uiOutput(paste0("mr_factor_m_model_", k)),
                    uiOutput(paste0("mr_factor_m_resp_", k)),
                    uiOutput(paste0("mr_factor_m_variable_", k)),
                    fluidRow(
                        column(
                            8,
                            uiOutput(paste0("mr_factor_m_interaction_", k))
                        ),
                        column(
                            4,
                            uiOutput(paste0("mr_factor_m_interactionappend_", k)),
                            style = "text-align:right; padding:15px"
                        )
                    ),
                    uiOutput(paste0("mr_factor_m_rand_", k)),
                    uiOutput(paste0("mr_factor_m_randfamily_", k)),
                    uiOutput(paste0("mr_factor_m_dist_", k)),
                    # uiOutput(paste0("mr_factor_m_check_binomd_", k)),
                    # uiOutput(paste0("mr_factor_m_binomd_", k)),
                    uiOutput(paste0("mr_factor_m_link_", k)),
                    uiOutput(paste0("mr_factor_m_factor_", k)),
                    uiOutput(paste0("mr_factor_m_check_nointercept_", k)),
                    uiOutput(paste0("mr_factor_m_check_offset_", k)),
                    uiOutput(paste0("mr_factor_m_offset_", k))
                )
            )
        
        div(
            factorAccordion,
            use_bs_tooltip(),
            use_bs_accordion_sidebar() # needs to be at end, for some reason
        )
    })
  
    output[[paste0('mr_factor_m_model_',k)]] <- renderUI({
        input[[paste0('mr_factor_m_resp_',k)]]
        input[[paste0('mr_factor_m_variable_', k)]]$right
        input[[paste0('mr_factor_m_rand_', k)]]
        input[[paste0('mr_factor_m_check_nointercept_', k)]]
        input[[paste0('mr_factor_m_binomd_', k)]]
        input[[paste0('mr_factor_m_slope_', k)]]
        
        responsePart = input[[paste0('mr_factor_m_resp_',k)]]
        variablePart = paste0(input[[paste0('mr_factor_m_variable_', k)]]$right, collapse="+")
        randomPart = NULL
        meanModel = NULL
        
        # binomd not working
        # if (all(!is.null(input[[paste0('mr_factor_m_check_binomd_', k)]]), 
        #         input[[paste0('mr_factor_m_check_binomd_', k)]], 
        #         !is.null(input[[paste0('mr_factor_m_binomd_', k)]]), 
        #         input[[paste0('mr_factor_m_binomd_', k)]] != ""))
        #     responsePart = paste0("cbind(", responsePart, ",", input[[paste0('mr_factor_m_binomd_', k)]], "-", responsePart, ")")
        
        if(!is.null(input[[paste0('mr_factor_m_rand_', k)]]) && length(input[[paste0('mr_factor_m_rand_', k)]]) > 0) {
            randomPart = paste0(input[[paste0('mr_factor_m_rand_', k)]], collapse = ")+(1|")
            randomPart = paste0("+(1|", randomPart, ")")
            
            #random slope (a|b) form
            if(all(!is.null(input[[paste0('mr_factor_m_check_slope_', k)]]), input[[paste0('mr_factor_m_check_slope_', k)]], length(input[[paste0('mr_factor_m_slope_', k)]]) > 0)) 
                for(i in 1:length(input[[paste0('mr_factor_m_rand_', k)]])) 
                    for(j in 1:length(input[[paste0('mr_factor_m_slope_', k)]])) 
                        randomPart = paste0(randomPart,"+(", input[[paste0('mr_factor_m_slope_', k)]][j], "|", input[[paste0('mr_factor_m_rand_', k)]][i], ")")
        }
        
        if(!is.null(input[[paste0('mr_factor_m_check_nointercept_', k)]]) && !input[[paste0('mr_factor_m_check_nointercept_', k)]]) {
            if(variablePart == "")
                meanModel = paste0(responsePart, "~", "1", randomPart)
            else
                meanModel = paste0(responsePart, "~", variablePart, randomPart)
        }
        
        if(!is.null(input[[paste0('mr_factor_m_check_nointercept_', k)]]) && input[[paste0('mr_factor_m_check_nointercept_', k)]]) {
            if(variablePart == "")
                meanModel = paste0(responsePart, "~", "-1", randomPart)
            else
                meanModel = paste0(responsePart, "~", "-1+", variablePart, randomPart)
        }
        textAreaInput(sprintf("mr_factor_m_model_%s", k), "Model",value = meanModel, height = "80px")
    })
    
    output[[paste0('mr_factor_m_resp_',k)]] <- renderUI({
        input$file
        input$resetData
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        nameValue = names(select_if(data2, is.numeric))
        
        selectInput(
            sprintf("mr_factor_m_resp_%s", k),
            "Response Variable", 
            choices = as.list(c("", nameValue)),
            selected = NULL,
            multiple = FALSE
        )
    })
    
    output[[paste0('mr_factor_m_variable_',k)]] <- renderUI({
        input$file
        input$resetData
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        input[[paste0('mr_factor_m_interactionappend_',k)]]
    
        nameValue <- c(names(data2), g_mr_factor_m_interaction[[k]])
    
        chooserInput(
            sprintf("mr_factor_m_variable_%s", k), 
            "Variable", 
            "Selected", 
            nameValue, 
            c(), 
            size = 15, 
            multiple = TRUE, 
            width = 150
        )
    })
    
    output[[paste0('mr_factor_m_interaction_',k)]] <- renderUI({
        input$file
        input$resetData
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
    
        nameValue <- names(data2)
        
        selectInput(
            sprintf("mr_factor_m_interaction_%s", k),
            "Make Interaction Variable", 
            choices = as.list(c("", nameValue)),
            selected = NULL,
            multiple = TRUE
        )
    })
    
    output[[paste0('mr_factor_m_interactionappend_',k)]] <- renderUI({
        input$file
        input$resetData
    
        actionButton(sprintf('mr_factor_m_interactionappend_%s',k), "Append")
    })
    
    observeEvent(input[[paste0('mr_factor_m_interactionappend_',k)]], {
        g_mr_factor_m_interaction[[k]] <<- c(g_mr_factor_m_interaction[[k]], paste(input[[paste0('mr_factor_m_interaction_',k)]], collapse=":"))
    })

    output[[paste0('mr_factor_m_rand_',k)]] <- renderUI({
        input$file
        input$resetData
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        input[[paste0('mr_factor_m_randinteractionappend_',k)]]
        
        nameValue <- c(names(data2), g_mr_factor_m_randinteraction[[k]])
        
        selectInput(
            sprintf("mr_factor_m_rand_%s", k), 
            "Random Effects", 
            choices = as.list(nameValue),
            multiple = TRUE
        )
    })
    
    output[[paste0('mr_factor_m_check_slope_',k)]] <- renderUI({
        input$file
        input$resetData
        
        checkboxInput(
            sprintf("mr_factor_m_check_slope_%s", k), 
            "Random slope model", 
            value = FALSE
        )
    })
    
    output[[paste0('mr_factor_m_slope_',k)]] <- renderUI({
        input$file
        input$resetData
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        if(is.null(input[[paste0('mr_factor_m_check_slope_',k)]]) || !input[[paste0('mr_factor_m_check_slope_',k)]])
            return()
        
        nameValue <- names(data2)
        
        selectInput(
            sprintf("mr_factor_m_slope_%s", k),
            "Random slope", 
            choices = as.list(c("", nameValue)),
            multiple = TRUE
        )
    })
    
    output[[paste0('mr_factor_m_check_randinteraction_',k)]] <- renderUI({
        input$file
        input$resetData
    
        checkboxInput(
            sprintf("mr_factor_m_check_randinteraction_%s", k), 
            "Interaction in the Random effect", 
            value = FALSE
        )
    })
    
    output[[paste0('mr_factor_m_randinteraction_',k)]] <- renderUI({
        input$file
        input$resetData
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        if(is.null(input[[paste0('mr_factor_m_check_randinteraction_',k)]]) || !input[[paste0('mr_factor_m_check_randinteraction_',k)]])
            return()
        
        nameValue <- names(data2)
    
        selectInput(
            sprintf("mr_factor_m_randinteraction_%s", k), 
            "Make Interaction Variable",
            choices = as.list(c("", nameValue)),
            selected = NULL, 
            multiple = TRUE
        )
    })
    
    output[[paste0('mr_factor_m_randinteractionappend_',k)]] <- renderUI({
        input$file
        input$resetData
    
        if(is.null(input[[paste0('mr_factor_m_check_randinteraction_',k)]]) || !input[[paste0('mr_factor_m_check_randinteraction_',k)]])
            return()
        
        actionButton(sprintf("mr_factor_m_randinteractionappend_%s", k), "Append")
    })
    
    observeEvent(input[[paste0('mr_factor_m_randinteractionappend_', k)]], {
        if(length(input[[paste0('mr_factor_m_randinteraction_', k)]]) > 1)
            g_mr_factor_m_randinteraction[[k]] <<- c(g_mr_factor_m_randinteraction[[k]], paste(input[[paste0('mr_factor_m_randinteraction_', k)]], collapse=":"))
    })
    
    output[[paste0('mr_factor_m_dist_', k)]] <- renderUI({
        input$file
        input$resetData
        
        dist = c(
            "binomial" = "binomial",
            "gaussian" = "gaussian",
            "poisson" = "poisson",
            "gamma" = "gamma"
        )
        
        selectInput(
            sprintf("mr_factor_m_dist_%s", k),
            "Distribution", 
            choices = dist, 
            multiple = FALSE
        )
    })
    
    output[[paste0('mr_factor_m_check_binomd_', k)]] <- renderUI({
        input$file
        input$resetData
        
        if (is.null(input[[paste0('mr_factor_m_dist_', k)]]) || input[[paste0('mr_factor_m_dist_', k)]] != "binomial")
            return()
        
        checkboxInput(sprintf("mr_factor_m_check_binomd_%s", k), "Binomial denominator", value = FALSE)
    })

    output[[paste0('mr_factor_m_binomd_', k)]] <- renderUI({
        input$file
        input$resetData
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        if (is.null(input[[paste0('mr_factor_m_check_binomd_', k)]]) || !input[[paste0('mr_factor_m_check_binomd_', k)]])
            return()
        
        nameValue = names(select_if(data2, is.numeric))
        
        selectInput(
            sprintf("mr_factor_m_binomd_%s", k),
            "Binomial denominator",
            choices = as.list(c("", nameValue)),
            multiple = FALSE
        )
    })
    
    output[[paste0('mr_factor_m_link_', k)]] <- renderUI({
        input$file
        input$resetData
        
        if(is.null(input[[paste0('mr_factor_m_dist_', k)]]))
            return()
        
        link = c(
            "identity" = "identity",
            "log" = "log",
            "logit" = "logit",
            "probit" = "probit",
            "cloglog" = "cloglog",
            "inverse" = "inverse"
        )
        
        if (input[[paste0('mr_factor_m_dist_', k)]] == "binomial")
            selection = "logit"
        else if (input[[paste0('mr_factor_m_dist_', k)]] == "poisson")
            selection = "log"
        else if (input[[paste0('mr_factor_m_dist_', k)]] == "gamma")
            selection = "log"
        else
            selection = "identity"
        
        selectInput(
            sprintf("mr_factor_m_link_%s", k),
            "Link Function",
            choices = link,
            selected = selection,
            multiple = FALSE
        )
    })
    
    output[[paste0('mr_factor_m_check_nointercept_', k)]] <- renderUI({
        input$file
        input$resetData
        
        checkboxInput(sprintf("mr_factor_m_check_nointercept_%s", k), "No intercept model", value = FALSE)
    })

    output[[paste0('mr_factor_m_check_offset_', k)]] <- renderUI({
        input$file
        input$resetData
        
        checkboxInput(sprintf("mr_factor_m_check_offset_%s", k), "Offset Variable", value = FALSE)
    })
    
    output[[paste0('mr_factor_m_offset_', k)]] <- renderUI({
        input$file
        input$resetData
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        if (is.null(input[[paste0('mr_factor_m_check_offset_', k)]]) || !input[[paste0('mr_factor_m_check_offset_', k)]])
            return()
        
        nameValue = names(select_if(data2, is.numeric))
    
        selectInput(
            sprintf("mr_factor_m_offset_%s", k),
            "Offset Variable (Numeric Only)",
            choices = as.list(c("", nameValue)),
            selected = NULL,
            multiple = FALSE
        )
    })
    
    output[[paste0('mr_factor_m_factor_', k)]] <- renderUI({
        if(is.null(input$mr_factor_g_corr_structure) || input$mr_factor_g_corr_structure != "factor")
            return()
        
        factor=c(
            "1st factor" = "1",
            "2nd factor" = "2",
            "3rd factor" = "3",
            "4th factor" = "4", 
            "5th factor" = "5"
        )
        
        selectInput(
            sprintf("mr_factor_m_factor_%s", k), 
            'Factor', 
            choices = factor,
            multiple = FALSE
        )
    })
})

observe({
    input$file
    input$resetData
    lapply(1:6, function(k) {
        if(is.null(input[[paste0('mr_factor_m_model_', k)]]) 
            || is.null(input[[paste0('mr_factor_m_rand_',k)]]) 
            || length(input[[paste0('mr_factor_m_rand_',k)]]) == 0) {
            output[[paste0('mr_factor_m_randfamily_',k)]] <- renderUI({
                return()
            })
        }
        
        input[[paste0('mr_factor_m_rand_',k)]]
        
        if(length(input[[paste0('mr_factor_m_rand_',k)]]) > 0) {
            output[[paste0('mr_factor_m_randfamily_',k)]] <- renderUI({
                nRand <- length(strsplit(input[[paste0('mr_factor_m_model_',k)]], "|", fixed = TRUE)[[1]]) - 1
                if (nRand >= 1) {
                    dynamic_selection_list <- lapply(1:nRand, function(i) {
                        selectInput(
                            paste0("mr_factor_m_rf_", i, "_", k),
                            paste("Distribution for Random effects",i),
                            choices=c("gaussian","gamma","beta")
                        )
                    })
                    do.call(tagList, dynamic_selection_list)
                }
            })
        } else {
            output[[paste0('mr_factor_m_randfamily_',k)]] <- renderUI({
                return()
            })
        }
    })
})

output$mr_factor_a_order<-renderUI({
    selectInput("mr_factor_a_order","Order for Laplace's Approximation",choices=c("1"="1","2"="2"),multiple=FALSE)
})

output$mr_factor_a_reml<-renderUI({
    selectInput("mr_factor_a_reml","REML method or not",choice = c(TRUE, FALSE),multiple=FALSE)
})

# Factor Results ####

output$mr_factor_r_lambda <- renderUI({
    input$file
    input$resetData
    input$mr_factor_g_run
    
    tagList(
        withMathJax(),
        withMathJax(tableOutput("mr_factor_r_lambdatable"))
    )
})

output$mr_factor_r_gamma <- renderUI({
    input$file
    input$resetData
    input$mr_factor_g_run
    
    tagList(
        withMathJax(),
        withMathJax(tableOutput("mr_factor_r_gammatable"))
    )
})

output$mr_factor_r_beta <- renderUI({
    input$file
    input$resetData
    input$mr_factor_g_run
    
    tagList(
        withMathJax(),
        withMathJax(tableOutput("mr_factor_r_betatable"))
    )
})

output$mr_factor_r_lambdatable <- renderTable({
    input$file
    input$resetData
    input$mr_factor_g_run
    
    if (g_resetresult == FALSE || is.null(mr_factor_results()$option$lambda))
        return()
  
    mr_factor_results()$option$lambda
}, rownames = TRUE, bordered = TRUE, caption = "Lambda", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_factor_r_gammatable <- renderTable({
    input$file
    input$resetData
    input$mr_factor_g_run
    
    if (g_resetresult == FALSE || is.null(mr_factor_results()$option$gamma))
        return()
  
    mr_factor_results()$option$gamma
}, rownames = TRUE, bordered = TRUE, caption = "Gamma", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_factor_r_betatable <- renderTable({
    input$file
    input$resetData
    input$mr_factor_g_run
    
    if (g_resetresult == FALSE || is.null(mr_factor_results()$option$beta))
        return()
  
    mr_factor_results()$option$beta
}, rownames = TRUE, bordered = TRUE, caption = "Beta", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_factor_r_likelihood <- renderTable({
    input$file
    input$resetData
    input$mr_factor_g_run
    
    if (g_resetresult == FALSE || is.null(mr_factor_results()$option$likelihood))
        return()
  
    mr_factor_results()$option$likelihood
}, rownames = FALSE, bordered = TRUE, caption = "Likelihood", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)