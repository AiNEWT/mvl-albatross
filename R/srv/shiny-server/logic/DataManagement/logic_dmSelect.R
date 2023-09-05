# Select Rows Preview ####

observeEvent(input$dm_sr_preview, {
    if (is.null(input$dm_sr_expression)) {
        return()
    }
    
    g_sr_previewdata <<- data2
    lapply(1:length(g_sr_previewdata), function(k) {
        if (g_attributes[[k]] == "factor") {
            g_sr_previewdata[[k]] <<- as.character(g_sr_previewdata[[k]])
        }
    })
    tryCatch(expr=filterExpr <- rlang::parse_expr(input$dm_sr_expression),
            error = function(e) {
                showNotification("Wrong Expression.", type="warning")
                g_error <<- TRUE
            }
    )
    tryCatch(expr=g_sr_previewdata <<- dplyr::filter(g_sr_previewdata,!!filterExpr), 
            error = function(e) {
                showNotification("Wrong Expression.", type="warning")
                g_error <<- TRUE
            }
    )
    
    if (g_error == TRUE) {
        g_error <<- FALSE
        return()
    }
})

# Select Rows Run ####

observeEvent(input$dm_sr_run, {
    if (is.null(input$dm_sr_expression)) {
        return()
    }
    
    lapply(1:length(data2), function(k) {
        if (g_attributes[[k]] == "factor") {
            data2[[k]] <<- as.character(data2[[k]])
        }
    })
    tryCatch(expr=filterExpr <- rlang::parse_expr(input$dm_sr_expression),
            error = function(e) {
                showNotification("Wrong Expression.", type="warning")
                g_error <<- TRUE
            }
    )
    tryCatch(expr=data2 <<- dplyr::filter(data2,!!filterExpr), 
            error = function(e) {
                showNotification("Wrong Expression.", type="warning")
                g_error <<- TRUE
            }
    )
    
    if (g_error == TRUE) {
        g_error <<- FALSE
        return()
    }
    
    lapply(1:length(data2), function(k) {
        if (g_attributes[[k]] == "factor") {
            data2[[k]] <<- as.factor(data2[[k]])
        }
    })
})

# Select Rows Components ####

output$dm_sr_expression<-renderUI({
    input$dm_sr_varname1
    input$dm_sr_varname2
    input$dm_sr_arith
    
    expressionValue=NULL
    expressionValue=paste0(c(input$dm_sr_varname1, input$dm_sr_varname2), collapse=paste0(" ", input$dm_sr_arith, " "))
        
    textAreaInput("dm_sr_expression","Expression",value=expressionValue)
})

output$dm_sr_check_tools <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("dm_sr_check_tools", "Tools for Expression", value = FALSE)
})

output$dm_sr_selectvarname1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_mvi_run
    input$dm_ct_run
    input$dm_sr_run  
    input$dm_rn_run
    input$dm_del_run
    input$dm_md_run
    
    if(is.null(input$dm_sr_check_tools) || input$dm_sr_check_tools == FALSE)
        return()
    
    nameValue = names(data2)
    
    selectInput(
        "dm_sr_varname1",
        "Variable 1",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = FALSE
    )
})

output$dm_sr_selectvarname2 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_mvi_run
    input$dm_ct_run
    input$dm_sr_run  
    input$dm_rn_run
    input$dm_del_run
    input$dm_md_run
    
    if(is.null(input$dm_sr_check_tools) || input$dm_sr_check_tools == FALSE)
        return()
    
    nameValue = names(data2)
    if (is.null(input$dm_sr_varname1)) {
        NULL
    } else if (class(data2[[input$dm_sr_varname1]]) == "numeric" || class(data2[[input$dm_sr_varname1]]) == "integer") {
        nameValue = names(select_if(data2, is.numeric))
    } else if (class(data2[[input$dm_sr_varname1]]) == "factor") {
        nameValue = names(select_if(data2, is.factor))
    }
    
    selectInput(
        "dm_sr_varname2",
        "Variable 2",
        choices = as.list(c("", nameValue)), 
        selected = NULL,
        multiple = FALSE
    )
})

output$dm_sr_selectarith <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_mvi_run
    input$dm_ct_run
    input$dm_sr_run  
    input$dm_rn_run
    input$dm_del_run
    input$dm_md_run
    
    if(is.null(input$dm_sr_check_tools) || input$dm_sr_check_tools == FALSE)
        return()
    
    Operators=c("", "==", "!=", ">", ">=", "<", "<=")
    selectInput("dm_sr_arith","Operator", choices=Operators, selected=NULL, multiple = FALSE)
})