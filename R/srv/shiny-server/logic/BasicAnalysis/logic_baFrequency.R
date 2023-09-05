# Frequency Run ####

observeEvent(input$ba_fa_run, {
    g_resetresult <<- TRUE
})

ba_fa_results <- eventReactive(input$ba_fa_run, {
    faTitle = NULL
    faResults = NULL
    faRpresults = NULL
    faCpresults = NULL
    faPresults = NULL
    faCsresults = NULL
    faFeresults = NULL
    faMnresults = NULL
    
    if (input$ba_fa_varname1 == "" || input$ba_fa_checkcolumn == TRUE && is.null(input$ba_fa_varname2)) {
        showNotification("Please choose variable", type="warning")
        return()
    }
    
    # Start Progressbar
    withProgress(message = 'Freqeuncy Analysis', style = "notification", value = 0, {
        
    # ?? Single Variable ####
    if (input$ba_fa_checkcolumn == FALSE) {
        faRowVariableLevels = levels(as.factor(data2[[input$ba_fa_varname1]]))
        faTitle <- paste0("Row Variable : ", input$ba_fa_varname1)
        faResults <- data.frame(table(data2[[input$ba_fa_varname1]]))["Freq"]
        faFreqpercent <- data.frame(apply(faResults, 1, function(x) return(x / sum(faResults) * 100)))
        faResults <- cbind(faResults, faFreqpercent)
        rownames(faResults) <- c(faRowVariableLevels)
        colnames(faResults) <- c("Freqeuncy", "Freq Percent")
        
        faPlot1 <- ggplot(data2[input$ba_fa_varname1], 
            aes_string(x=input$ba_fa_varname1)) + 
            geom_bar() + 
            ggtitle("Frequency Plot") + 
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5))
        
        # Increase Progressbar
        incProgress(1.0, detail = paste("Loading..."))
        
        return(list(
            Title = faTitle,
            Favalue = faResults,
            Plot1 = faPlot1
        ))
    }

    # Increase Progressbar
    incProgress(0.5, detail = paste("Loading..."))
    
    # ?? Title ####        
    faTitle <- paste0("Row Variable : ", input$ba_fa_varname1, " , Column Variable : ", input$ba_fa_varname2)

    # ?? Frequency Table ####        
    faRowVariableLevels = levels(as.factor(data2[[input$ba_fa_varname1]]))
    faColumnVarialbeLevels = levels(as.factor(data2[[input$ba_fa_varname2]]))
    if (length(faRowVariableLevels) < 2 ||
        length(faColumnVarialbeLevels) < 2) {
        showNotification("Variable must have at least 2 levels", type="warning")
        return()
    }

    faResults <- matrix(
        data = table(data2[[input$ba_fa_varname1]], data2[[input$ba_fa_varname2]]), 
        nrow = length(faRowVariableLevels),
        ncol = length(faColumnVarialbeLevels)
    )
    faCalResults <- data.frame(faResults)
    faResults <- cbind(faResults, rowSums(faResults))
    faResults <- rbind(faResults, colSums(faResults))
    
    rownames(faResults) <- c(faRowVariableLevels, "Total")
    colnames(faResults) <- c(faColumnVarialbeLevels, "Total")
    
    # ?? Plot1 ####
    faPlotdata <- data.frame(table(data2[[input$ba_fa_varname1]], data2[[input$ba_fa_varname2]]))
    colnames(faPlotdata) <- c(input$ba_fa_varname1, input$ba_fa_varname2, "Freq")
    faPlot1 <- ggplot(faPlotdata, aes_string(fill=input$ba_fa_varname2, y="Freq", x=input$ba_fa_varname1)) + 
        geom_bar(position="dodge", stat="identity") + 
        ggtitle("Frequency Plot") + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    
    # ?? Row Percent ####
    if (input$ba_fa_rpcheck == TRUE) {
        faRpresults <- apply(faCalResults, 1, function(x) return(x/sum(x) * 100))
        faRpresults <- t(data.frame(faRpresults))
        
        rownames(faRpresults) <- faRowVariableLevels
        colnames(faRpresults) <- faColumnVarialbeLevels
    }
    
    # ?? Column Percent ####
    if (input$ba_fa_cpcheck == TRUE) {
        faCpresults <- sapply(faCalResults, function(x) return(x/sum(x) * 100))
        faCpresults <- data.frame(faCpresults)
    
        rownames(faCpresults) <- faRowVariableLevels
        colnames(faCpresults) <- faColumnVarialbeLevels
    }
    
    # ?? Percent ####
    if (input$ba_fa_pcheck == TRUE) {
        faPresults <- sapply(faCalResults, function(x) return(x/sum(faCalResults) * 100))
        faPresults <- data.frame(faPresults)
    
        rownames(faPresults) <- faRowVariableLevels
        colnames(faPresults) <- faColumnVarialbeLevels
    }
    
    # ?? Chi-squared Test ####
    if (input$ba_fa_cscheck == TRUE) {
        faCstest <- chisq.test(data2[[input$ba_fa_varname1]], data2[[input$ba_fa_varname2]])
        faCsresults <- matrix(data = c(faCstest$statistic, faCstest$parameter, faCstest$p.value)
                                      ,nrow = 1, ncol = 3)
        colnames(faCsresults) <- c("X-squared", "df", "p-value")
        rownames(faCsresults) <- "Chi-squared"
    }
    
    # ?? Fisher Exact Test ####
    if (input$ba_fa_fecheck == TRUE) {
        al=input$ba_fa_fecheck_option
        faFetest <- fisher.test(data2[[input$ba_fa_varname1]], data2[[input$ba_fa_varname2]],alternative=al)
        faFeresults <- matrix(data = c(round(faFetest$p.value, digits = 4), faFetest$alternative), nrow = 1, ncol = 2)
        colnames(faFeresults) <- c("p-value", "Alternative")
        rownames(faFeresults) <- "Fisher"
    }
    
    # ?? McNemar's test ####
    if(input$ba_fa_mncheck == TRUE) {
        faMntest <- mcnemar.test(data2[[input$ba_fa_varname1]], data2[[input$ba_fa_varname2]], correct = TRUE)
        faMnresults <- matrix(data = c(round(faMntest$statistic, digits = 4), 
                                       faMntest$parameter, 
                                       faMntest$p.value
                                       ), 
                              nrow = 1, ncol = 3)
        colnames(faMnresults) <- c("McNemar's chi-squared", "df", "p-value")
        rownames(faMnresults) <- "McNemar"
    }
    
    # Increase Progressbar
    incProgress(0.5, detail = paste("Loading..."))
    
    return(list(
        Title = faTitle,
        Favalue = faResults,
        Rowpercent = faRpresults,
        Columnpercent = faCpresults,
        Allpercent = faPresults,
        Chisquared = faCsresults,
        Fisher = faFeresults,
        McNemar = faMnresults,
        Plot1 = faPlot1
    ))
    
    }) # End Progressbar
})

# Frequency Components ####

output$ba_fa_selectvarname1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue <- names(select_if(data2, is.factor))
    nameValue <- c(nameValue, names(select_if(data2, is.integer)))
    
    selectInput(
        "ba_fa_varname1",
        "Row Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = FALSE
    )
})

output$ba_fa_checkcolumn <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput(
        inputId = 'ba_fa_checkcolumn',
        label = 'Use Column Variable',
        value = FALSE
    )
})

output$ba_fa_selectvarname2 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if (is.null(input$ba_fa_checkcolumn) || input$ba_fa_checkcolumn == FALSE)
        return()
    
    nameValue <- names(select_if(data2, is.factor))
    nameValue <- c(nameValue, names(select_if(data2, is.integer)))
    nameValue <- nameValue[-which(nameValue==input$ba_fa_varname1)]
    
    selectInput(
        "ba_fa_varname2",
        "Column Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = FALSE
    )
})

output$ba_fa_rpcheck <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if (is.null(input$ba_fa_checkcolumn) || input$ba_fa_checkcolumn == FALSE)
        return()
    
    checkboxInput("ba_fa_rpcheck", "Row Percent")
})

output$ba_fa_cpcheck <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if (is.null(input$ba_fa_checkcolumn) || input$ba_fa_checkcolumn == FALSE)
        return()
    
    checkboxInput("ba_fa_cpcheck", "Column Percent")
})

output$ba_fa_pcheck <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if (is.null(input$ba_fa_checkcolumn) || input$ba_fa_checkcolumn == FALSE)
        return()
    
    checkboxInput("ba_fa_pcheck", "Percent")
})

output$ba_fa_cscheck <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if (is.null(input$ba_fa_checkcolumn) || input$ba_fa_checkcolumn == FALSE)
        return()
    
    checkboxInput("ba_fa_cscheck", "Chi-squared Test")
})

output$ba_fa_fecheck <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if (is.null(input$ba_fa_checkcolumn) || input$ba_fa_checkcolumn == FALSE)
        return()
    
    checkboxInput("ba_fa_fecheck", "Fisher Exact Test")
})

output$ba_fa_fecheck_option <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$ba_fa_fecheck
    
    if (is.null(input$ba_fa_checkcolumn) || input$ba_fa_checkcolumn == FALSE)
        return()
    if (is.null(input$ba_fa_fecheck) || input$ba_fa_fecheck == FALSE)
      return()
    ah=c("two.sided","greater","less")
    selection="two.sided"
    selectInput("ba_fa_fecheck_option", "Alternative of Fisher Exact Test", choices = ah, selected = selection, multiple = FALSE)
})

output$ba_fa_mncheck <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if (is.null(input$ba_fa_checkcolumn) || input$ba_fa_checkcolumn == FALSE)
        return()
    
    checkboxInput("ba_fa_mncheck", "McNemar's Test")
})

# Frequency Results ####

output$ba_fa_title <- renderUI({
    input$ba_fa_run
    
    titleOutput = NULL
    if (g_resetresult == FALSE || is.null(ba_fa_results()$Title))
        titleOutput = "Row Variable : "
    else
        titleOutput = ba_fa_results()$Title
    
    paste(titleOutput)
})

output$ba_fa_faresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_fa_run
    
    if (g_resetresult == FALSE || is.null(ba_fa_results()$Favalue))
        return()
    
    ba_fa_results()$Favalue
}, rownames = TRUE, colnames = TRUE, bordered = TRUE, caption = "Frequency Table", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 0)

output$ba_fa_rpresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_fa_run
    
    if (g_resetresult == FALSE || is.null(ba_fa_results()$Rowpercent))
        return()
    
    ba_fa_results()$Rowpercent
}, rownames = TRUE, colnames = TRUE, bordered = TRUE, caption = "Row Percent Table", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

output$ba_fa_cpresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_fa_run
    
    if (g_resetresult == FALSE || is.null(ba_fa_results()$Columnpercent))
        return()
    
    ba_fa_results()$Columnpercent
}, rownames = TRUE, colnames = TRUE, bordered = TRUE, caption = "Column Percent Table", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

output$ba_fa_presults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_fa_run
    
    if (g_resetresult == FALSE || is.null(ba_fa_results()$Allpercent))
        return()
    
    ba_fa_results()$Allpercent
}, rownames = TRUE, colnames = TRUE, bordered = TRUE, caption = "Percent Table", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

output$ba_fa_csresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_fa_run
    
    if (g_resetresult == FALSE || is.null(ba_fa_results()$Chisquared))
        return()
    
    ba_fa_results()$Chisquared
}, bordered = TRUE, caption = "Pearson's Chi-squared", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

output$ba_fa_feresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_fa_run
    
    if (g_resetresult == FALSE || is.null(ba_fa_results()$Fisher))
        return()
    
    ba_fa_results()$Fisher
}, bordered = TRUE, caption = "Fisher Exact Test", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

output$ba_fa_mnresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_fa_run
    
    if (g_resetresult == FALSE || is.null(ba_fa_results()$McNemar))
        return()
    
    ba_fa_results()$McNemar
}, bordered = TRUE, caption = "McNemar's Chi-squared test with continuity correction", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

# Frequency Plots ####

output$ba_fa_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$ba_fa_run

    if (g_resetresult == FALSE || is.null(ba_fa_results()$Plot1))
        return()
    
    local({
        output$ba_fa_downloadplot1 <- downloadPlot(
            ba_fa_results()$Plot1
        )
    })
    
    ba_fa_results()$Plot1
})

output$ba_fa_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$ba_fa_run

    if (g_resetresult == FALSE || is.null(ba_fa_results()$Plot1))
        return()
    
    div(
        style = "position: relative; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("ba_fa_plot1"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left:0.5em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "ba_fa_downloadplot1", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})
