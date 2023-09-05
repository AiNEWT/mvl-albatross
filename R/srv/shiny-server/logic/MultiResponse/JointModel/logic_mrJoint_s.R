output$mr_joint_s1_model<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$mr_joint_s1_survivaltime
    input$mr_joint_s1_variable$right
    input$mr_joint_s1_indicator
    input$mr_joint_s1_initialtime
    
    if(any(is.null(input$mr_joint_s1_survivaltime), is.null(input$mr_joint_s1_variable), is.null(input$mr_joint_s1_indicator)))
        return()
    
    responsePart = NULL
    variablePart = paste0(input$mr_joint_s1_variable$right, collapse=" + ")
    randomPart = NULL
    modelEquation = NULL
    
    if (!is.null(input$mr_joint_s1_initialtime) && input$mr_joint_s1_initialtime != "") {
        responsePart = paste0("Surv(", input$mr_joint_s1_survivaltime," - ",input$mr_joint_s1_initialtime,", ",input$mr_joint_s1_indicator," == ", input$mr_joint_s1_status, ")")
    } else {
        responsePart = paste0("Surv(", input$mr_joint_s1_survivaltime,", ",input$mr_joint_s1_indicator," == ", input$mr_joint_s1_status, ")")
    }
    
    if(length(input$mr_joint_s1_variable$right) == 0) {
        variablePart = "1"
    }
    
    if(!is.null(input$mr_joint_s1_rand) && length(input$mr_joint_s1_rand) > 0) {
        randomPart = paste0(input$mr_joint_s1_rand, collapse = ") + (1|")
        randomPart = paste0(" + (1|", randomPart, ")")
    }
    
    modelEquation = paste0(responsePart, " ~ ", variablePart, randomPart)
    
    textAreaInput("mr_joint_s1_model","Model for Interesting Event",value = modelEquation, height = "60px")
})

output$mr_joint_s1_survivaltime<-renderUI({
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
            "mr_joint_s1_survivaltime",
            "Survival Time", 
            choices = as.list(c("",nameValue)),
            selected = NULL,
            multiple = FALSE
        ),
        bsTooltip(
            "mr_joint_s1_survivaltime", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$mr_joint_s1_check_initial<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("mr_joint_s1_check_initial", "Initial Time", value = FALSE)
})

output$mr_joint_s1_initialtime<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$mr_joint_s1_check_initial) || !input$mr_joint_s1_check_initial) 
        return()
    
    nameValue = names(select_if(data2, is.numeric))

    div(
        selectInput(
            "mr_joint_s1_initialtime",
            "Initial Time",
            choices = as.list(c("",nameValue)),
            multiple = FALSE
        ),
        bsTooltip(
            "mr_joint_s1_initialtime", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$mr_joint_s1_variable<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$mr_joint_s1_interactionappend
    
    nameValue <- c(names(data2), g_mr_joint_s1_interaction)
    
    chooserInput("mr_joint_s1_variable", "Variable", "Selected", nameValue, c(), size = 15, multiple = TRUE, width=150)
})


output$mr_joint_s1_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue <- names(data2)
    
    selectInput(
        "mr_joint_s1_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$mr_joint_s1_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    actionButton("mr_joint_s1_interactionappend", "Append")
})

observeEvent(input$mr_joint_s1_interactionappend, {
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    g_mr_joint_s1_interaction <<- c(g_mr_joint_s1_interaction, paste(input$mr_joint_s1_interaction, collapse=":"))
})

output$mr_joint_s1_rand <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$mr_joint_s1_randinteractionappend
    
    nameValue <- c(names(data2), g_sa_fm_randinteraction)
    
    selectInput(
        "mr_joint_s1_rand",
        "Random Effects",
        choices = as.list(nameValue), 
        multiple = TRUE
    )
})

output$mr_joint_s1_indicator<-renderUI({
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
            "mr_joint_s1_indicator",
            "Censoring Indicator",
            choices = as.list(c("",nameValue)),
            multiple = FALSE
        ),
        bsTooltip(
            "mr_joint_s1_indicator", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$mr_joint_s1_status <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    selectInput(
        "mr_joint_s1_status",
        "Censoring Value",
        choices = as.list(c("0" = "0", "1" = "1", "2" = "2")),
        selected = "1",
        multiple = FALSE
    )
})

output$mr_joint_s1_dist <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    selectInput(
        "mr_joint_s1_dist",
        "Model for Baseline Hazard",
        choices = as.list(c("Semi-parametric" = "FM", "Weibull model" = "Weibull")),
        selected = "FM",
        multiple = FALSE
    )
})

output$mr_joint_s2_check_competing <- renderUI({
    input$file
    input$resetData
    input$di_option_run

    checkboxInput("mr_joint_s2_check_competing", strong("Use"), value = FALSE)
})

output$mr_joint_s2_model<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$mr_joint_s2_survivaltime
    input$mr_joint_s2_variable$right
    input$mr_joint_s2_indicator
    input$mr_joint_s2_initialtime
    
    if(is.null(input$mr_joint_s2_check_competing) || !input$mr_joint_s2_check_competing)
        return()
    
    if(any(is.null(input$mr_joint_s2_survivaltime), is.null(input$mr_joint_s2_variable), is.null(input$mr_joint_s2_indicator)))
        return()
    
    responsePart = NULL
    variablePart = paste0(input$mr_joint_s2_variable$right, collapse=" + ")
    randomPart = NULL
    modelEquation = NULL
    
    if (!is.null(input$mr_joint_s2_initialtime) && input$mr_joint_s2_initialtime != "") {
        responsePart = paste0("Surv(", input$mr_joint_s2_survivaltime," - ",input$mr_joint_s2_initialtime,", ",input$mr_joint_s2_indicator," == ", input$mr_joint_s2_status, ")")
    } else {
        responsePart = paste0("Surv(", input$mr_joint_s2_survivaltime,", ",input$mr_joint_s2_indicator," == ", input$mr_joint_s2_status, ")")
    }
    
    if(length(input$mr_joint_s2_variable$right) == 0) {
        variablePart = "1"
    }
    
    if(!is.null(input$mr_joint_s2_rand) && length(input$mr_joint_s2_rand) > 0) {
        randomPart = paste0(input$mr_joint_s2_rand, collapse = ") + (1|")
        randomPart = paste0(" + (1|", randomPart, ")")
    }
    
    modelEquation = paste0(responsePart, " ~ ", variablePart, randomPart)
    
    textAreaInput("mr_joint_s2_model", "Model for Competing Event", value = modelEquation, height = "60px")
})

output$mr_joint_s2_survivaltime<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
        
    if(is.null(input$mr_joint_s2_check_competing) || !input$mr_joint_s2_check_competing)
        return()
    
    nameValue = names(select_if(data2, is.numeric))
    
    div(
        selectInput(
            "mr_joint_s2_survivaltime",
            "Survival Time", 
            choices = as.list(c("",nameValue)),
            selected = NULL,
            multiple = FALSE
        ),
        bsTooltip(
            "mr_joint_s2_survivaltime", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$mr_joint_s2_check_initial<-renderUI({
    input$file
    input$resetData
    input$di_option_run
        
    if(is.null(input$mr_joint_s2_check_competing) || !input$mr_joint_s2_check_competing)
        return()
    
    checkboxInput("mr_joint_s2_check_initial", "Initial Time", value = FALSE)
})

output$mr_joint_s2_initialtime<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
        
    if(is.null(input$mr_joint_s2_check_competing) || !input$mr_joint_s2_check_competing)
        return()
    
    if(is.null(input$mr_joint_s2_check_initial) || !input$mr_joint_s2_check_initial) 
        return()
    
    nameValue = names(select_if(data2, is.numeric))
    
    div(
        selectInput(
            "mr_joint_s2_initialtime",
            "Initial Time",
            choices = as.list(c("",nameValue)),
            multiple = FALSE
        ),
        bsTooltip(
            "mr_joint_s2_initialtime", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$mr_joint_s2_variable<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$mr_joint_s2_interactionappend
        
    if(is.null(input$mr_joint_s2_check_competing) || !input$mr_joint_s2_check_competing)
        return()
    
    nameValue <- c(names(data2), g_mr_joint_s2_interaction)
    
    chooserInput("mr_joint_s2_variable", "Variable", "Selected", nameValue, c(), size = 15, multiple = TRUE, width=150)
})


output$mr_joint_s2_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
        
    if(is.null(input$mr_joint_s2_check_competing) || !input$mr_joint_s2_check_competing)
        return()
    
    nameValue <- names(data2)
    
    selectInput(
        "mr_joint_s2_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$mr_joint_s2_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
        
    if(is.null(input$mr_joint_s2_check_competing) || !input$mr_joint_s2_check_competing)
        return()
    
    actionButton("mr_joint_s2_interactionappend", "Append")
})

observeEvent(input$mr_joint_s2_interactionappend, {
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
        
    if(is.null(input$mr_joint_s2_check_competing) || !input$mr_joint_s2_check_competing)
        return()
    
    g_mr_joint_s2_interaction <<- c(g_mr_joint_s2_interaction, paste(input$mr_joint_s2_interaction, collapse=":"))
})

output$mr_joint_s2_rand <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$mr_joint_s2_randinteractionappend
        
    if(is.null(input$mr_joint_s2_check_competing) || !input$mr_joint_s2_check_competing)
        return()
    
    nameValue <- c(names(data2), g_sa_fm_randinteraction)
    
    selectInput(
        "mr_joint_s2_rand",
        "Random Effects",
        choices = as.list(nameValue), 
        multiple = TRUE
    )
})

output$mr_joint_s2_indicator<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
        
    if(is.null(input$mr_joint_s2_check_competing) || !input$mr_joint_s2_check_competing)
        return()
    
    nameValue = names(select_if(data2, is.numeric))
    
    div(
        selectInput(
            "mr_joint_s2_indicator",
            "Censoring Indicator",
            choices = as.list(c("",nameValue)),
            multiple = FALSE
        ),
        bsTooltip(
            "mr_joint_s2_indicator", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$mr_joint_s2_status <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
        
    if(is.null(input$mr_joint_s2_check_competing) || !input$mr_joint_s2_check_competing)
        return()
    
    selectInput(
        "mr_joint_s2_status",
        "Censoring Value",
        choices = as.list(c("0" = "0", "1" = "1", "2" = "2")),
        selected = "2",
        multiple = FALSE
    )
})

output$mr_joint_s2_dist <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
        
    if(is.null(input$mr_joint_s2_check_competing) || !input$mr_joint_s2_check_competing)
        return()
    
    selectInput(
        "mr_joint_s2_dist",
        "Model for Baseline Hazard",
        choices = as.list(c("Semi-parametric" = "FM", "Weibull model" = "Weibull")),
        selected = "FM",
        multiple = FALSE
    )
})