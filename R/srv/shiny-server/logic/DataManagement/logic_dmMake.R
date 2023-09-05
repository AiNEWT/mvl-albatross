# Make Variable Run ####

observeEvent(input$dm_mv_run, {
    if (input$dm_mv_expression == "" || input$dm_mv_newvariablename == "") {
        return()
    }
    
    newVariableName <- input$dm_mv_newvariablename
    
    nameAccepted = which(names(data2)==newVariableName)
    if (length(nameAccepted) != 0) {
        showNotification("Type another variable name.", type="warning")
        return()
    }
    
    tryCatch(expr = data2<<-mutate(data2, newVariable=!!rlang::parse_expr(input$dm_mv_expression)),
        error = function(e) {
            showNotification("Wrong Expression.", type="warning")
            g_error <<- TRUE
        }
    )
    
    if (g_error == TRUE) {
        g_error <<- FALSE
        return()
    }
    
    names(data2)[length(data2)]<<-newVariableName
    
    if (class(data2[[length(data2)]]) == "integer") {
        g_attributes <<- c(g_attributes, "integer")
        g_resetattributes <<- c(g_resetattributes, "integer")
    } else {
        g_attributes <<- c(g_attributes, "numeric")
        g_resetattributes <<- c(g_resetattributes, "numeric")
    }
    
    g_colattr$cols[[newVariableName]] <<- col_double()
    g_resetcolattr$cols[[newVariableName]] <<- col_double()
    
    # run mutate -> attribute removed
    # if attribute == NULL editData is not working
    # make attribute
    attr(data2, "spec") <<- g_colattr
})

# Make Variable Components ####

output$dm_mv_newvariablename <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    textAreaInput("dm_mv_newvariablename","New Variable Name",value=NULL)
})

output$dm_mv_expression<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_varname
    input$dm_mv_arith
    
    expressionValue=NULL
    if (is.null(input$dm_mv_arith)) {
        NULL
    }
    else if (input$dm_mv_arith == "log") {
        expressionValue=paste0(input$dm_mv_varname, collapse=" + ")
        expressionValue=paste0('log(', expressionValue, ')') 
    } else {
        expressionValue=paste0(input$dm_mv_varname, collapse=paste0(" ", input$dm_mv_arith, " "))
    }
        
    textAreaInput("dm_mv_expression","Expression",value=expressionValue)
})

output$dm_mv_check_tools <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("dm_mv_check_tools", "Tools for Expression", value = FALSE)
})

output$dm_mv_selectvarname <- renderUI({
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
    
    if(is.null(input$dm_mv_check_tools) || input$dm_mv_check_tools == FALSE)
        return()
    
    div(
        selectInput(
            "dm_mv_varname",
            "Variable",
            choices = as.list(c("", names(select_if(data2, is.numeric)))), 
            selected = NULL, 
            multiple = TRUE
        ),
        bsTooltip(
            "dm_mv_varname", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$dm_mv_selectarith <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$dm_mv_check_tools) || input$dm_mv_check_tools == FALSE)
        return()
    
    Operators=c("","+","-","*","/","^","%/%","%%","log")
    selectInput("dm_mv_arith","Operator", choices=Operators, selected=NULL, multiple = FALSE)
})