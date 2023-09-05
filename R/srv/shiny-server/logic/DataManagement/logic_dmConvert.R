# Convert Type Run ####

observeEvent(input$dm_ct_run, {
    
    if (is.null(input$dm_ct_variable) || input$dm_ct_type == "")
        return()
    
    cfVariable <- input$dm_ct_variable
    cfType <- input$dm_ct_type
    for(k in 1:length(cfVariable)){
        cfIndex <- which(names(data2) == cfVariable[k])
        if (cfType == "Factor") {
            data2[[cfIndex]] <<- as.factor(data2[[cfIndex]])
            g_attributes[cfIndex] <<- "factor"
            g_colattr$cols[[cfIndex]] <<- col_character()
        } else if (cfType == "Integer") {
            data2[[cfIndex]] <<- as.integer(data2[[cfIndex]])
            g_attributes[cfIndex] <<- "integer"
            g_colattr$cols[[cfIndex]] <<- col_double()
        } else if (cfType == "Numeric") {
            data2[[cfIndex]] <<- as.numeric(data2[[cfIndex]])
            g_attributes[cfIndex] <<- "numeric"
            g_colattr$cols[[cfIndex]] <<- col_double()
        }
    }
    attr(data2, "spec") <<- g_colattr
})

# Convert Type Components ####

output$dm_ct_selecttype <- renderUI({
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
    
    Types=c("", "Factor", "Integer", "Numeric")
    selectInput("dm_ct_type","Data Type", choices = Types, selected = NULL, multiple = FALSE)
})

output$dm_ct_variables <- renderUI({
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
    
    if (is.null(input$dm_ct_type))
        return()
    
    nameValue = NULL
    if (input$dm_ct_type == "Factor") {
        nameValue = names(select_if(data2, function(x) {return(is.numeric(x) || is.character(x))}))
    } else {
        nameValue = names(data2)
    }
    
    selectInput(
        "dm_ct_variable",
        "Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})