lapply(1:6, function(k) {    
    output[[paste0('mr_mdhglm_m_model_',k)]] <- renderUI({
        input[[paste0('mr_mdhglm_m_resp_',k)]]
        input[[paste0('mr_mdhglm_m_variable_', k)]]$right
        input[[paste0('mr_mdhglm_m_rand_', k)]]
        input[[paste0('mr_mdhglm_m_check_nointercept_', k)]]
        input[[paste0('mr_mdhglm_m_binomd_', k)]]
        input[[paste0('mr_mdhglm_m_slope_', k)]]
        
        responsePart = input[[paste0('mr_mdhglm_m_resp_',k)]]
        variablePart = paste0(input[[paste0('mr_mdhglm_m_variable_', k)]]$right, collapse=" + ")
        randomPart = NULL
        meanModel = NULL
        
        # binomd not working
        # if (all(!is.null(input[[paste0('mr_mdhglm_m_check_binomd_', k)]]), 
        #         input[[paste0('mr_mdhglm_m_check_binomd_', k)]], 
        #         !is.null(input[[paste0('mr_mdhglm_m_binomd_', k)]]), 
        #         input[[paste0('mr_mdhglm_m_binomd_', k)]] != ""))
        #     responsePart = paste0("cbind(", responsePart, ",", input[[paste0('mr_mdhglm_m_binomd_', k)]], "-", responsePart, ")")
        
        if(!is.null(input[[paste0('mr_mdhglm_m_rand_', k)]]) && length(input[[paste0('mr_mdhglm_m_rand_', k)]]) > 0) {
            if (is.null(input[[paste0('mr_mdhglm_m_check_withoutslope_', k)]]) || !input[[paste0('mr_mdhglm_m_check_withoutslope_', k)]]) {
                randomPart = paste0(input[[paste0('mr_mdhglm_m_rand_', k)]], collapse = ") + (1|")
                randomPart = paste0(" + (1|", randomPart, ")")
            }
            
            #random slope (a|b) form
            if(all(!is.null(input[[paste0('mr_mdhglm_m_check_slope_', k)]]), input[[paste0('mr_mdhglm_m_check_slope_', k)]], length(input[[paste0('mr_mdhglm_m_slope_', k)]]) > 0)) 
                for(i in 1:length(input[[paste0('mr_mdhglm_m_rand_', k)]])) 
                    for(j in 1:length(input[[paste0('mr_mdhglm_m_slope_', k)]])) 
                        randomPart = paste0(randomPart, " + (", input[[paste0('mr_mdhglm_m_slope_', k)]][j], "|", input[[paste0('mr_mdhglm_m_rand_', k)]][i], ")")
        }
        
        if(!is.null(input[[paste0('mr_mdhglm_m_check_nointercept_', k)]]) && !input[[paste0('mr_mdhglm_m_check_nointercept_', k)]]) {
            if(variablePart == "")
                meanModel = paste(responsePart, "~", "1", randomPart)
            else
                meanModel = paste(responsePart, "~", variablePart, randomPart)
        }
        
        if(!is.null(input[[paste0('mr_mdhglm_m_check_nointercept_', k)]]) && input[[paste0('mr_mdhglm_m_check_nointercept_', k)]]) {
            if(variablePart == "")
                meanModel = paste(responsePart, "~", "-1", randomPart)
            else
                meanModel = paste(responsePart, "~", "-1 +", variablePart, randomPart)
        }
        textAreaInput(sprintf("mr_mdhglm_m_model_%s", k), "Model for Mean", value = meanModel, height = "60px")
    })
    
    output[[paste0('mr_mdhglm_m_resp_',k)]] <- renderUI({
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
                sprintf("mr_mdhglm_m_resp_%s", k),
                "Response Variable", 
                choices = as.list(c("", nameValue)),
                selected = NULL,
                multiple = FALSE
            ),
            bsTooltip(
                sprintf("mr_mdhglm_m_resp_%s", k),
                "Numeric Only",
                "right", 
                options = list(container = "body")
            )
        )
    })
    
    output[[paste0('mr_mdhglm_m_variable_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        input[[paste0('mr_mdhglm_m_interactionappend_',k)]]
    
        nameValue <- c(names(data2), g_mr_mdhglm_m_interaction[[k]])
    
        chooserInput(
            sprintf("mr_mdhglm_m_variable_%s", k), 
            "Variable", 
            "Selected", 
            nameValue, 
            c(), 
            size = 15, 
            multiple = TRUE, 
            width = 150
        )
    })
    
    output[[paste0('mr_mdhglm_m_interaction_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
    
        nameValue <- names(data2)
        
        selectInput(
            sprintf("mr_mdhglm_m_interaction_%s", k),
            "Make Interaction Variable", 
            choices = as.list(c("", nameValue)),
            selected = NULL,
            multiple = TRUE
        )
    })
    
    output[[paste0('mr_mdhglm_m_interactionappend_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
    
        actionButton(sprintf('mr_mdhglm_m_interactionappend_%s',k), "Append")
    })
    
    observeEvent(input[[paste0('mr_mdhglm_m_interactionappend_',k)]], {
        g_mr_mdhglm_m_interaction[[k]] <<- c(g_mr_mdhglm_m_interaction[[k]], paste(input[[paste0('mr_mdhglm_m_interaction_',k)]], collapse=":"))
    })

    output[[paste0('mr_mdhglm_m_rand_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        input[[paste0('mr_mdhglm_m_randinteractionappend_',k)]]
        
        nameValue <- c(names(data2), g_mr_mdhglm_m_randinteraction[[k]])
        
        selectInput(
            sprintf("mr_mdhglm_m_rand_%s", k), 
            "Random Effects", 
            choices = as.list(nameValue),
            multiple = TRUE
        )
    })
    
    output[[paste0('mr_mdhglm_m_check_slope_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        
        checkboxInput(
            sprintf("mr_mdhglm_m_check_slope_%s", k), 
            "Random slope model", 
            value = FALSE
        )
    })
    
    output[[paste0('mr_mdhglm_m_check_withoutslope_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        
        if(is.null(input[[paste0('mr_mdhglm_m_check_slope_',k)]]) || !input[[paste0('mr_mdhglm_m_check_slope_',k)]])
            return()
        
        checkboxInput(
            sprintf("mr_mdhglm_m_check_withoutslope_%s", k), 
            "Without Random Intercept", 
            value = FALSE
        )
    })
    
    output[[paste0('mr_mdhglm_m_slope_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        if(is.null(input[[paste0('mr_mdhglm_m_check_slope_',k)]]) || !input[[paste0('mr_mdhglm_m_check_slope_',k)]])
            return()
        
        nameValue <- names(data2)
        
        selectInput(
            sprintf("mr_mdhglm_m_slope_%s", k),
            "Random slope", 
            choices = as.list(c("", nameValue)),
            multiple = TRUE
        )
    })
    
    output[[paste0('mr_mdhglm_m_check_randinteraction_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
    
        checkboxInput(
            sprintf("mr_mdhglm_m_check_randinteraction_%s", k), 
            "Interaction in the Random effect", 
            value = FALSE
        )
    })
    
    output[[paste0('mr_mdhglm_m_randinteraction_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        if(is.null(input[[paste0('mr_mdhglm_m_check_randinteraction_',k)]]) || !input[[paste0('mr_mdhglm_m_check_randinteraction_',k)]])
            return()
        
        nameValue <- names(data2)
    
        selectInput(
            sprintf("mr_mdhglm_m_randinteraction_%s", k), 
            "Make Interaction Variable",
            choices = as.list(c("", nameValue)),
            selected = NULL, 
            multiple = TRUE
        )
    })
    
    output[[paste0('mr_mdhglm_m_randinteractionappend_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
    
        if(is.null(input[[paste0('mr_mdhglm_m_check_randinteraction_',k)]]) || !input[[paste0('mr_mdhglm_m_check_randinteraction_',k)]])
            return()
        
        actionButton(sprintf("mr_mdhglm_m_randinteractionappend_%s", k), "Append")
    })
    
    observeEvent(input[[paste0('mr_mdhglm_m_randinteractionappend_', k)]], {
        if(length(input[[paste0('mr_mdhglm_m_randinteraction_', k)]]) > 1)
            g_mr_mdhglm_m_randinteraction[[k]] <<- c(g_mr_mdhglm_m_randinteraction[[k]], paste(input[[paste0('mr_mdhglm_m_randinteraction_', k)]], collapse=":"))
    })
    
    output[[paste0('mr_mdhglm_m_dist_', k)]] <- renderUI({
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
            sprintf("mr_mdhglm_m_dist_%s", k),
            "Distribution", 
            choices = dist, 
            multiple = FALSE
        )
    })
    
    output[[paste0('mr_mdhglm_m_check_binomd_', k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        
        if (is.null(input[[paste0('mr_mdhglm_m_dist_', k)]]) || input[[paste0('mr_mdhglm_m_dist_', k)]] != "binomial")
            return()
        
        checkboxInput(sprintf("mr_mdhglm_m_check_binomd_%s", k), "Binomial denominator", value = FALSE)
    })

    output[[paste0('mr_mdhglm_m_binomd_', k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        if (is.null(input[[paste0('mr_mdhglm_m_check_binomd_', k)]]) || !input[[paste0('mr_mdhglm_m_check_binomd_', k)]])
            return()
        
        nameValue = names(select_if(data2, is.numeric))
        
        selectInput(
            sprintf("mr_mdhglm_m_binomd_%s", k),
            "Binomial denominator",
            choices = as.list(c("", nameValue)),
            multiple = FALSE
        )
    })
    
    output[[paste0('mr_mdhglm_m_link_', k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        
        if(is.null(input[[paste0('mr_mdhglm_m_dist_', k)]]))
            return()
        
        link = c(
            "identity" = "identity",
            "log" = "log",
            "logit" = "logit",
            "probit" = "probit",
            "cloglog" = "cloglog",
            "inverse" = "inverse"
        )
        
        if (input[[paste0('mr_mdhglm_m_dist_', k)]] == "binomial")
            selection = "logit"
        else if (input[[paste0('mr_mdhglm_m_dist_', k)]] == "poisson")
            selection = "log"
        else if (input[[paste0('mr_mdhglm_m_dist_', k)]] == "gamma")
            selection = "log"
        else
            selection = "identity"
        
        selectInput(
            sprintf("mr_mdhglm_m_link_%s", k),
            "Link Function",
            choices = link,
            selected = selection,
            multiple = FALSE
        )
    })
    
    output[[paste0('mr_mdhglm_m_check_nointercept_', k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        
        checkboxInput(sprintf("mr_mdhglm_m_check_nointercept_%s", k), "No intercept model", value = FALSE)
    })

    output[[paste0('mr_mdhglm_m_check_offset_', k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        
        checkboxInput(sprintf("mr_mdhglm_m_check_offset_%s", k), "Offset Variable", value = FALSE)
    })
    
    output[[paste0('mr_mdhglm_m_offset_', k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        if (is.null(input[[paste0('mr_mdhglm_m_check_offset_', k)]]) || !input[[paste0('mr_mdhglm_m_check_offset_', k)]])
            return()
        
        nameValue = names(select_if(data2, is.numeric))
    
        div(
            selectInput(
                sprintf("mr_mdhglm_m_offset_%s", k),
                "Offset Variable",
                choices = as.list(c("", nameValue)),
                selected = NULL,
                multiple = FALSE
            ),
            bsTooltip(
                sprintf("mr_mdhglm_m_offset_%s", k),
                "Numeric Only",
                "right", 
                options = list(container = "body")
            )
        )
    })
    
    output[[paste0('mr_mdhglm_m_factor_', k)]] <- renderUI({
        if(is.null(input$mr_mdhglm_g_corr_structure) || input$mr_mdhglm_g_corr_structure != "factor")
            return()
        
        factor=c(
            "1st factor" = "1",
            "2nd factor" = "2",
            "3rd factor" = "3",
            "4th factor" = "4", 
            "5th factor" = "5"
        )
        
        selectInput(
            sprintf("mr_mdhglm_m_factor_%s", k), 
            'factor', 
            choices = factor,
            multiple = FALSE
        )
    })
    
    output[[paste0('mr_mdhglm_m_check_vif_', k)]] <-renderUI({
        input$file
        input$resetData
        input$di_option_run
        
        checkboxInput(sprintf("mr_mdhglm_m_check_vif_%s", k), "VIF", value = FALSE)
    })
    
    output[[paste0('mr_mdhglm_m_check_robustse_', k)]] <-renderUI({
        input$file
        input$resetData
        input$di_option_run
        
        checkboxInput(sprintf("mr_mdhglm_m_check_robustse_%s", k), "Robust Standard Errors", value = FALSE)
    })
    
    output[[paste0('mr_mdhglm_m_check_confint_', k)]] <-renderUI({
        input$file
        input$resetData
        input$di_option_run
        
        checkboxInput(sprintf("mr_mdhglm_m_check_confint_%s", k), "Confidence Intervals for Coefficients", value = FALSE)
    })
    
    output[[paste0('mr_mdhglm_m_check_exp_', k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        
        checkboxInput(sprintf("mr_mdhglm_m_check_exp_%s", k), "Exponential scale", value = FALSE)
    })
    
    output[[paste0('mr_mdhglm_m_check_rcodes_', k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        
        checkboxInput(sprintf("mr_mdhglm_m_check_rcodes_%s", k), "R Codes", value = FALSE)
    })
})

output$mr_mdhglm_m_check_comparison_1<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("mr_mdhglm_m_check_comparison_1", "Comparison With Other Model", value = FALSE)
})

observe({
    input$file
    input$resetData
    input$di_option_run
    lapply(1:6, function(k) {
        if(is.null(input[[paste0('mr_mdhglm_m_model_', k)]]) 
            || is.null(input[[paste0('mr_mdhglm_m_rand_',k)]]) 
            || length(input[[paste0('mr_mdhglm_m_rand_',k)]]) == 0) {
            output[[paste0('mr_mdhglm_m_randfamily_',k)]] <- renderUI({
                return()
            })
        }
        
        input[[paste0('mr_mdhglm_m_rand_',k)]]
        
        if(length(input[[paste0('mr_mdhglm_m_rand_',k)]]) > 0) {
            output[[paste0('mr_mdhglm_m_randfamily_',k)]] <- renderUI({
                nRand <- length(strsplit(input[[paste0('mr_mdhglm_m_model_',k)]], "|", fixed = TRUE)[[1]]) - 1
                if (nRand >= 1) {
                    dynamic_selection_list <- lapply(1:nRand, function(i) {
                        selectInput(
                            paste0("mr_mdhglm_m_rf_", i, "_", k),
                            paste("Distribution for Random effects",i),
                            choices=c("gaussian","gamma","beta")
                        )
                    })
                    do.call(tagList, dynamic_selection_list)
                }
            })
        } else {
            output[[paste0('mr_mdhglm_m_randfamily_',k)]] <- renderUI({
                return()
            })
        }
    })
})