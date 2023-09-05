lapply(1:6, function(k) {
    
    output[[paste0('mr_joint_p_check_phi_',k)]]<-renderUI({
        input$file
        input$resetData
        input$di_option_run
    
        checkboxInput(sprintf("mr_joint_p_check_phi_%s", k), strong("Use"), value = FALSE)
    })
    
    output[[paste0('mr_joint_p_model_',k)]] <- renderUI({
        input[[paste0('mr_joint_p_resp_',k)]]
        input[[paste0('mr_joint_p_variable_', k)]]$right
        input[[paste0('mr_joint_p_rand_', k)]]
        input[[paste0('mr_joint_p_slope_', k)]]
        
        if(is.null(input[[paste0('mr_joint_p_check_phi_',k)]]) || !input[[paste0('mr_joint_p_check_phi_',k)]])
            return()
        
        if(any(is.null(input[[paste0('mr_joint_p_resp_',k)]]), is.null(input[[paste0('mr_joint_p_variable_', k)]])))
            return()
        
        responsePart = input[[paste0('mr_joint_p_resp_',k)]]
        variablePart = paste0(input[[paste0('mr_joint_p_variable_', k)]]$right, collapse=" + ")
        randomPart = NULL
        phiModel = NULL
        
        if(!is.null(input[[paste0('mr_joint_p_rand_', k)]]) && length(input[[paste0('mr_joint_p_rand_', k)]]) > 0) {
            randomPart = paste0(input[[paste0('mr_joint_p_rand_', k)]], collapse = ") + (1|")
            randomPart = paste0(" + (1|", randomPart, ")")
        }
        
        if(variablePart == "")
            phiModel = paste(responsePart, "~", "1", randomPart)
        else
            phiModel = paste(responsePart, "~", variablePart, randomPart)
        
        textAreaInput(sprintf("mr_joint_p_model_%s", k), "Model for Phi", value = phiModel, height = "60px")
    })
    
    output[[paste0('mr_joint_p_resp_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        if(is.null(input[[paste0('mr_joint_p_check_phi_',k)]]) || !input[[paste0('mr_joint_p_check_phi_',k)]])
            return()
        
        nameValue = names(select_if(data2, is.numeric))
        
        div(
            selectInput(
                sprintf("mr_joint_p_resp_%s", k),
                "Response Variable", 
                choices = as.list(c("", nameValue)),
                selected = NULL,
                multiple = FALSE
            ),
            bsTooltip(
                sprintf("mr_joint_p_resp_%s", k),
                "Numeric Only",
                "right", 
                options = list(container = "body")
            )
        )
    })
    
    output[[paste0('mr_joint_p_variable_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        input[[paste0('mr_joint_p_interactionappend_',k)]]
    
        if(is.null(input[[paste0('mr_joint_p_check_phi_',k)]]) || !input[[paste0('mr_joint_p_check_phi_',k)]])
            return()
        
        nameValue <- c(names(data2), g_mr_joint_p_interaction[[k]])
    
        chooserInput(
            sprintf("mr_joint_p_variable_%s", k), 
            "Variable", 
            "Selected", 
            nameValue, 
            c(), 
            size = 15, 
            multiple = TRUE, 
            width = 150
        )
    })
    
    output[[paste0('mr_joint_p_interaction_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
    
        if(is.null(input[[paste0('mr_joint_p_check_phi_',k)]]) || !input[[paste0('mr_joint_p_check_phi_',k)]])
            return()
        
        nameValue <- names(data2)
        
        selectInput(
            sprintf("mr_joint_p_interaction_%s", k),
            "Make Interaction Variable", 
            choices = as.list(c("", nameValue)),
            selected = NULL,
            multiple = TRUE
        )
    })
    
    output[[paste0('mr_joint_p_interactionappend_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
    
        if(is.null(input[[paste0('mr_joint_p_check_phi_',k)]]) || !input[[paste0('mr_joint_p_check_phi_',k)]])
            return()
        
        actionButton(sprintf('mr_joint_p_interactionappend_%s',k), "Append")
    })
    
    observeEvent(input[[paste0('mr_joint_p_interactionappend_',k)]], {
        g_mr_joint_p_interaction[[k]] <<- c(g_mr_joint_p_interaction[[k]], paste(input[[paste0('mr_joint_p_interaction_',k)]], collapse=":"))
    })

    output[[paste0('mr_joint_p_rand_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        input[[paste0('mr_joint_p_randinteractionappend_',k)]]
        
        if(is.null(input[[paste0('mr_joint_p_check_phi_',k)]]) || !input[[paste0('mr_joint_p_check_phi_',k)]])
            return()
        
        nameValue <- c(names(data2), g_mr_joint_p_randinteraction[[k]])
        
        selectInput(
            sprintf("mr_joint_p_rand_%s", k), 
            "Random Effects", 
            choices = as.list(nameValue),
            multiple = TRUE
        )
    })
    
    output[[paste0('mr_joint_p_check_randinteraction_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
    
        if(is.null(input[[paste0('mr_joint_p_check_phi_',k)]]) || !input[[paste0('mr_joint_p_check_phi_',k)]])
            return()
        
        checkboxInput(
            sprintf("mr_joint_p_check_randinteraction_%s", k), 
            "Interaction in the Random effect", 
            value = FALSE
        )
    })
    
    output[[paste0('mr_joint_p_randinteraction_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        if(is.null(input[[paste0('mr_joint_p_check_phi_',k)]]) || !input[[paste0('mr_joint_p_check_phi_',k)]])
            return()
        
        if(is.null(input[[paste0('mr_joint_p_check_randinteraction_',k)]]) || !input[[paste0('mr_joint_p_check_randinteraction_',k)]])
            return()
        
        nameValue <- names(data2)
    
        selectInput(
            sprintf("mr_joint_p_randinteraction_%s", k), 
            "Make Interaction Variable",
            choices = as.list(c("", nameValue)),
            selected = NULL, 
            multiple = TRUE
        )
    })
    
    output[[paste0('mr_joint_p_randinteractionappend_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
    
        if(is.null(input[[paste0('mr_joint_p_check_phi_',k)]]) || !input[[paste0('mr_joint_p_check_phi_',k)]])
            return()
        
        if(is.null(input[[paste0('mr_joint_p_check_randinteraction_',k)]]) || !input[[paste0('mr_joint_p_check_randinteraction_',k)]])
            return()
        
        actionButton(sprintf("mr_joint_p_randinteractionappend_%s", k), "Append")
    })
    
    observeEvent(input[[paste0('mr_joint_p_randinteractionappend_', k)]], {
        if(length(input[[paste0('mr_joint_p_randinteraction_', k)]]) > 1)
            g_mr_joint_p_randinteraction[[k]] <<- c(g_mr_joint_p_randinteraction[[k]], paste(input[[paste0('mr_joint_p_randinteraction_', k)]], collapse=":"))
    })
   
    output[[paste0('mr_joint_p_link_', k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        
        if(is.null(input[[paste0('mr_joint_p_check_phi_',k)]]) || !input[[paste0('mr_joint_p_check_phi_',k)]])
            return()
        
        link = c("log" = "log","inverse" = "inverse")
        
        selectInput(
            sprintf("mr_joint_p_link_%s", k),
            "Link Function",
            choices = link,
            multiple = FALSE
        )
    })
    
    output[[paste0('mr_joint_p_check_offset_', k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        
        if(is.null(input[[paste0('mr_joint_p_check_phi_',k)]]) || !input[[paste0('mr_joint_p_check_phi_',k)]])
            return()
        
        checkboxInput(sprintf("mr_joint_p_check_offset_%s", k), "Offset Variable", value = FALSE)
    })
    
    output[[paste0('mr_joint_p_offset_', k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        if(is.null(input[[paste0('mr_joint_p_check_phi_',k)]]) || !input[[paste0('mr_joint_p_check_phi_',k)]])
            return()
        
        if (is.null(input[[paste0('mr_joint_p_check_offset_', k)]]) || !input[[paste0('mr_joint_p_check_offset_', k)]])
            return()
        
        nameValue = names(select_if(data2, is.numeric))
    
        div(
            selectInput(
                sprintf("mr_joint_p_offset_%s", k),
                "Offset Variable",
                choices = as.list(c("", nameValue)),
                selected = NULL,
                multiple = FALSE
            ),
            bsTooltip(
                sprintf("mr_joint_p_offset_%s", k),
                "Numeric Only",
                "right", 
                options = list(container = "body")
            )
        )
    })
})

observe({
    input$file
    input$resetData
    input$di_option_run
    lapply(1:6, function(k) {
        if(is.null(input[[paste0('mr_joint_p_model_', k)]]) 
            || is.null(input[[paste0('mr_joint_p_rand_',k)]]) 
            || length(input[[paste0('mr_joint_p_rand_',k)]]) == 0) {
            output[[paste0('mr_joint_p_randfamily_',k)]] <- renderUI({
                return()
            })
        }
        
        input[[paste0('mr_joint_p_rand_',k)]]
        
        if(length(input[[paste0('mr_joint_p_rand_',k)]]) > 0) {
            output[[paste0('mr_joint_p_randfamily_',k)]] <- renderUI({
                nRand <- length(strsplit(input[[paste0('mr_joint_p_model_',k)]], "|", fixed = TRUE)[[1]]) - 1
                if (nRand >= 1) {
                    dynamic_selection_list <- lapply(1:nRand, function(i) {
                        selectInput(
                            paste0("mr_joint_p_rf_", i, "_", k),
                            paste("Distribution for Random effects",i),
                            choices=c("gaussian","gamma","beta")
                        )
                    })
                    do.call(tagList, dynamic_selection_list)
                }
            })
        } else {
            output[[paste0('mr_joint_p_randfamily_',k)]] <- renderUI({
                return()
            })
        }
    })
})