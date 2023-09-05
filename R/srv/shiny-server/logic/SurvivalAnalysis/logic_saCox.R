# Cox Run ####

observeEvent(input$sa_cox_run, {
    g_resetresult <<- TRUE
})

sa_cox_results <- eventReactive(input$sa_cox_run, {
    if (input$sa_cox_survivaltime == "" || input$sa_cox_indicator == "") {
        showNotification("Please select survival time and indicator", type="warning")
        return()
    }

    if (input$sa_cox_check_initial == TRUE && input$sa_cox_initialtime == "") {
        showNotification("Please choose initial time", type="warning")
        return()
    }
    
    # Start Progressbar
    withProgress(message = 'Cox Model', style = "notification", value = 0.1, {
    
    modelFormula <- input$sa_cox_model

    coxresults <- coxph(
        formula = as.formula(modelFormula), 
        data = data2, 
        robust = input$sa_cox_check_robustse,
        ties = input$sa_cox_ties
    )
    coxresults$call$formula <- as.formula(modelFormula)
    summaryresults<-summary(coxresults)
    
    incProgress(0.5, detail = paste("Loading..."))
    
    zphresults = NULL
    if(length(input$sa_cox_variable$right)>0)
        zphresults <- cox.zph(coxresults)
    
    incProgress(0.5, detail = paste("Loading..."))
    
    return(list(
        Cox = coxresults, 
        Summary = summaryresults,
        Zph = zphresults,
        Formula = modelFormula
    ))
    
    })
})

# Cox Components ####

output$sa_cox_model<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$sa_cox_survivaltime
    input$sa_cox_variable$right
    input$sa_cox_indicator
    input$sa_cox_initialtime
    
    if(any(is.null(input$sa_cox_survivaltime), is.null(input$sa_cox_variable), is.null(input$sa_cox_indicator)))
        return()

    responsePart = NULL
    variablePart = paste0(input$sa_cox_variable$right, collapse="+")
    coxModel = NULL

    
    if (!is.null(input$sa_cox_initialtime) && input$sa_cox_initialtime != "") {
        responsePart = paste0("Surv(", input$sa_cox_survivaltime," - ",input$sa_cox_initialtime,", ",input$sa_cox_indicator," == 1)")
    } else {
        responsePart = paste0("Surv(", input$sa_cox_survivaltime,", ",input$sa_cox_indicator," == 1)")
    }
    
    if(length(input$sa_cox_variable$right) == 0) {
        variablePart = "1"
    }
    coxModel = paste0(responsePart, " ~ ", variablePart)
    
    textAreaInput("sa_cox_model","Model",value = coxModel)
})

output$sa_cox_survivaltime<-renderUI({
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
            "sa_cox_survivaltime",
            "Survival Time", 
            choices = as.list(c("",nameValue)),
            selected = NULL,
            multiple = FALSE
        ),
        bsTooltip(
            "sa_cox_survivaltime", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$sa_cox_check_initial<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("sa_cox_check_initial", "Initial Time", value = FALSE)
})

output$sa_cox_initialtime<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$sa_cox_check_initial) || !input$sa_cox_check_initial) 
        return()
    
    nameValue = names(select_if(data2, is.numeric))

    div(
        selectInput(
            "sa_cox_initialtime",
            "Initial Time",
            choices = as.list(c("",nameValue)),
            multiple = FALSE
        ),
        bsTooltip(
            "sa_cox_initialtime", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$sa_cox_variable<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$sa_cox_interactionappend
    
    nameValue <- c(names(data2), g_sa_cox_interaction)
    
    chooserInput("sa_cox_variable", "Variable", "Selected", nameValue, c(), size = 15, multiple = TRUE, width=150)
})


output$sa_cox_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue <- names(data2)
    
    selectInput(
        "sa_cox_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$sa_cox_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    actionButton("sa_cox_interactionappend", "Append")
})

observeEvent(input$sa_cox_interactionappend, {
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    g_sa_cox_interaction <<- c(g_sa_cox_interaction, paste(input$sa_cox_interaction, collapse=":"))
})

output$sa_cox_indicator<-renderUI({
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
            "sa_cox_indicator",
            "Censoring Indicator",
            choices = as.list(c("",nameValue)),
            multiple = FALSE
        ),
        bsTooltip(
            "sa_cox_initialtime", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

#if we choose "exact", it will not work in most situations
output$sa_cox_ties <-renderUI({
    input$file
    input$resetData
    input$di_option_run

    tieList <- c(
        "efron" = "efron",
        "breslow" = "breslow"
        # "exact" = "exact"
    )

    selectInput(
        "sa_cox_ties",
        "Method for tie handling",
        choices = tieList,
        selected = "efron",
        multiple = FALSE
    )
})

output$sa_cox_check_robustse <- renderUI({
    input$file
    input$resetData
    input$di_option_run

    checkboxInput("sa_cox_check_robustse", "Robust Standard Error", value = FALSE)
})

# Cox Results ####

output$sa_cox_title <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$sa_cox_run
    
    titleOutput = NULL
    if (g_resetresult == FALSE || is.null(sa_cox_results()$Cox))
        titleOutput = ""
    else
        titleOutput = sa_cox_results()$Formula
    
    paste0("Model : ", titleOutput)
})

output$sa_cox_nullmodelNotification <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$sa_cox_run
    
    if(g_resetresult == FALSE || is.null(sa_cox_results()$Cox) || !is.null(sa_cox_results()$Zph))
        return()
        
    nullNotification = "Null model"
    
    nullNotification
})

output$sa_cox_numberresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_cox_run
    
    if (g_resetresult == FALSE || is.null(sa_cox_results()$Cox))
        return()
    
    numbermatrix <- matrix(data=c(sa_cox_results()$Cox$n, sa_cox_results()$Cox$nevent), ncol=2)
    colnames(numbermatrix) <- c("n", "number of events")
    
    numbermatrix
}, bordered = TRUE, caption = "", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 0)

output$sa_cox_concordanceresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_cox_run
    
    if (g_resetresult == FALSE || is.null(sa_cox_results()$Cox)  
        || is.null(sa_cox_results()$Summary))
        return()
    
    if (is.null(sa_cox_results()$Zph)) {
        concordancelength <- length(sa_cox_results()$Summary$concordance)
        concordancematrix <-matrix(
            data=sa_cox_results()$Summary$concordance[c(concordancelength-1, concordancelength)], ncol=2)
        colnames(concordancematrix) <- c("Concordance", "se")
        return(concordancematrix)
    }
    
    concordancematrix <-matrix(
        data=sa_cox_results()$Summary$concordance, ncol=2)
    colnames(concordancematrix) <- c("Concordance", "se")
    
    concordancematrix
}, bordered = TRUE, caption = "", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 3)

output$sa_cox_robustresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_cox_run
    
    if (g_resetresult == FALSE || is.null(sa_cox_results()$Cox)  
        || is.null(sa_cox_results()$Summary) || is.null(sa_cox_results()$Zph) 
        || !sa_cox_results()$Summary$used.robust)
        return()
    
    robustmatrix <-matrix(
        data=c(
            sa_cox_results()$Summary$robscore[1],
            sa_cox_results()$Summary$robscore[length(sa_cox_results()$Summary$robscore)]
        ), ncol=2)
    colnames(robustmatrix) <- c("Robust", "p-value")
    
    robustmatrix
}, bordered = TRUE, caption = "", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 3)

output$sa_cox_testresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_cox_run
    
    if (g_resetresult == FALSE || is.null(sa_cox_results()$Cox)  
        || is.null(sa_cox_results()$Summary) || is.null(sa_cox_results()$Zph))
        return()
    
    testmatrix <-matrix(
        data=c(
            sa_cox_results()$Summary$logtest, 
            sa_cox_results()$Summary$waldtest, 
            sa_cox_results()$Summary$sctest
        ), 
        ncol=3, 
        byrow = TRUE
    )
    colnames(testmatrix) <- c("Statistic", "df", "p-value")
    rownames(testmatrix) <- c("Likelihood ratio test", "Wald test", "Score(logrank) test")
    
    testmatrix
}, bordered = TRUE, rownames = TRUE, caption = "Test Results", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"),digits = 3)

output$sa_cox_coefresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_cox_run
    
    if (g_resetresult == FALSE || is.null(sa_cox_results()$Cox)  
        || is.null(sa_cox_results()$Summary) || is.null(sa_cox_results()$Zph))
        return()
    
    sa_cox_results()$Summary$coefficients
}, bordered = TRUE, rownames = TRUE, caption = "Coefficients", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"),digits = 6)

output$sa_cox_confintresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_cox_run
    
    if (g_resetresult == FALSE || is.null(sa_cox_results()$Cox)  
        || is.null(sa_cox_results()$Summary) || is.null(sa_cox_results()$Zph))
        return()
    
    sa_cox_results()$Summary$conf.int
}, bordered = TRUE, rownames = TRUE, caption = "Confidence Interval", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"),digits = 4)

output$sa_cox_zphresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_cox_run
    
    if (g_resetresult == FALSE || is.null(sa_cox_results()$Cox)  
        || is.null(sa_cox_results()$Zph))
        return()
    
    zphmatrix <- sa_cox_results()$Zph$table
    colnames(zphmatrix) <- c("Statistic", "df", "p-value")
    
    
    zphmatrix
}, bordered = TRUE, rownames = TRUE, caption = "Proportional Hazard Assumption Check", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 3)
