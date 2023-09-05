# ANOVA Run ####

observeEvent(input$ba_an_run, {
    g_resetresult <<- TRUE
})

ba_an_results <- eventReactive(input$ba_an_run, {
    if (is.null(input$ba_an_varname1) || is.null(input$ba_an_varname2)) {
        showNotification("Please choose variable", type="warning")
        return()
    }
    
    # Start Progressbar
    withProgress(message = 'ANOVA', style = "notification", value = 0, {
        
    anSpecificresults = NULL
    anShapiroresults = NULL
    anBpresults = NULL
    anKwresults = NULL
    
    temp=c(data2[input$ba_an_varname1])
    data2[input$ba_an_varname1]=ifelse(is.na(temp),"NA",temp)
    length1=length(input$ba_an_varname2)
    for (i in 1:length1) {
           temp=c(data2[input$ba_an_varname2[i]])
           data2[input$ba_an_varname2[i]]=ifelse(is.na(temp),"NA",temp)
    }
    anFormula <- paste(isolate(paste(input$ba_an_varname1,"~")), paste0(input$ba_an_varname2,collapse="+"))
    anLength <- length(input$ba_an_varname2)
    anLm <- lm(anFormula, data = data2, na.rm = TRUE)
    anLm$call$formula <- as.formula(anFormula)
    
    # ¦§ Specific Results ####
    anSpecificresults <- data.frame(unclass(summary.aov(aov(as.formula(anFormula), data = data2))))
    anSpecificresults[nrow(anSpecificresults) + 1,] <- 
        c(sum(anSpecificresults[[1]]), sum(anSpecificresults[[2]]), NA, NA, NA)
    rownames(anSpecificresults)[[nrow(anSpecificresults)]] <- "Total"
    anSpecificresults[[1]] <- as.integer(anSpecificresults[[1]])
    
    # Increase Progressbar
    incProgress(0.5, detail = paste("Loading..."))
    
    # ¦§ Shapiro-Wilk Test ####
    if (input$ba_an_shapirocheck == FALSE) {
        anShapiroresults <- NULL
    } else if (all(anLm$residuals == 0) || length(anLm$residuals) >= 5000) {
        showNotification("Too many data in selected groups.", type="warning")
    } else {
        anShapiro <- shapiro.test(stdres(anLm))
        anShapiroresults <- matrix(data=c(anShapiro$statistic, anShapiro$p.value),nrow = 1, ncol = 2)
        colnames(anShapiroresults) <- c("W", "p-value")
        rownames(anShapiroresults) <- "Shapiro-Wilks"
    }

    # ¦§ Breusch-Pagan Test ####
    if (input$ba_an_bpcheck == TRUE) {
        anBptest <- bptest(anLm, data = data2)
        anBpresults <- matrix(data=c(anBptest$statistic, anBptest$parameter, anBptest$p.value),nrow = 1, ncol = 3)
        colnames(anBpresults) <- c("BP", "df", "p-value")
        rownames(anBpresults) <- "Breusch-Pagan"
    }
    
    # ¦§ Kruskal-Wallis Test ####
    if (input$ba_an_kwcheck == TRUE) {
        if (length(input$ba_an_varname2) >= 2) {
            showNotification("Choose one group parameter (Kruskal-Wallis)", type="warning")
        } else if (length(grep(":", anFormula)) != 0 && 
                   any(sapply(strsplit(input$ba_an_varname2, ":")[[1]], function(x) return(is.factor(data2[[x]]))))){
            # if interation variable include factor variable kruskal-wallis doesn't support
            showNotification("Kruskal-Wallis test doesn't support this formula.", type="warning")
        } else {
            anKwtest <- kruskal.test(as.formula(anFormula), data = data2)
            anKwresults <- matrix(data=c(anKwtest$statistic, anKwtest$parameter, anKwtest$p.value),nrow = 1, ncol = 3)
            colnames(anKwresults) <- c("Statistic", "df", "p-value")
            rownames(anKwresults) <- "Kruskal-Wallis"
        }
    }

    # Increase Progressbar
    incProgress(0.5, detail = paste("Loading..."))
    
    # ¦§ Anova Results ####
    anAnovaresults <- anSpecificresults[anLength + 1,]
    anAnovaresults[2, 1] <- sum(anSpecificresults[1:anLength, 1])
    anAnovaresults[2, 2] <- sum(anSpecificresults[1:anLength, 2])
    anAnovaresults[3, 1] <- sum(anAnovaresults[1, 1], anAnovaresults[2, 1])
    anAnovaresults[3, 2] <- sum(anAnovaresults[1, 2], anAnovaresults[2, 2])
    row.names(anAnovaresults) <- c("Residuals", "Model", "Total")
    anAnovaresults <- anAnovaresults[c("Model", "Residuals", "Total"),]
    anAnovaresults[1, 3] <- anAnovaresults[1, 2] / anAnovaresults[1, 1]
    anAnovaresults[1, 4] <- anAnovaresults[1, 3] / anAnovaresults[2, 3]
    anAnovaresults[1, 5] <- 1-pf(anAnovaresults[1, 4], anAnovaresults[1, 1], anAnovaresults[2, 1])
    
    # ¦¦ Plot ####
    anPlot1 <- NULL
    anPlot2 <- NULL
    if (all(anLm$residuals == 0) || length(anLm$residuals) >= 5000) {
        anPlot1 <- NULL
        anPlot2 <- NULL
    } else {
        anPlot1 <- autoplot(anLm)[[2]] + 
            ggtitle("Normal Probability Plot") +
            theme_bw() + 
            theme(plot.title = element_text(hjust = 0.5))
        anPlot2 <- autoplot(anLm)[[3]] + 
            ggtitle("Scale-Location") +
            theme_bw() + 
            theme(plot.title = element_text(hjust = 0.5))
    }
    
    return(list(
        Formula = anFormula,
        Anova = anAnovaresults,
        Specific = anSpecificresults,
        Shapiro = anShapiroresults,
        Bp = anBpresults,
        Kw = anKwresults,
        Plot1 = anPlot1,
        Plot2 = anPlot2
    ))
    
    }) # End Progressbar
})

# ANOVA Components ####

output$ba_an_selectvarname1 <- renderUI({
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
            "ba_an_varname1",
            "Variable",
            choices = as.list(c("", nameValue)), 
            selected = NULL, 
            multiple = FALSE
        ),
        bsTooltip(
            "ba_an_varname1", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$ba_an_selectvarname2 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$ba_an_interactionappend
    
    nameValue <- c(names(data2), g_ba_an_interaction)
    nameValue <- nameValue[-which(nameValue==input$ba_an_varname1)]
    
    selectInput(
        "ba_an_varname2",
        "Group Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$ba_an_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue <- c(names(data2), g_ba_an_interaction)
    
    selectInput(
        "ba_an_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$ba_an_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$ba_an_run
    
    actionButton("ba_an_interactionappend", "Append")
})

observeEvent(input$ba_an_interactionappend, {
    g_ba_an_interaction <<- c(g_ba_an_interaction, paste(input$ba_an_interaction, collapse=":"))
})

output$ba_an_shapirocheck <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("ba_an_shapirocheck", "Shapiro-Wilk Test (Normality)")
})

output$ba_an_bpcheck <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("ba_an_bpcheck", "Breusch-Pagan Test (Homoscedasticity)")
})

output$ba_an_kwcheck <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("ba_an_kwcheck", "Kruskal-Wallis Test")
})

# ANOVA Results ####

output$ba_an_title <- renderText({
    input$file
    input$resetData
    input$di_option_run
    input$ba_an_run
    
    titleOutput = NULL
    if (g_resetresult == FALSE || is.null(ba_an_results()$Formula))
        titleOutput = ""
    else
        titleOutput = ba_an_results()$Formula
    
    paste0("Model : ", titleOutput)
})

output$ba_an_anovaresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_an_run
    
    if (g_resetresult == FALSE || is.null(ba_an_results()$Anova))
        return()
    
    ba_an_results()$Anova
}, rownames = TRUE, bordered = TRUE, caption = "ANOVA", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

output$ba_an_specificresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_an_run
    
    if (g_resetresult == FALSE || is.null(ba_an_results()$Specific))
        return()
    
    ba_an_results()$Specific
}, rownames = TRUE, bordered = TRUE, caption = "Specific Results", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

output$ba_an_shapiroresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_an_run
    
    if (g_resetresult == FALSE || is.null(ba_an_results()$Shapiro))
        return()
    
    ba_an_results()$Shapiro
}, rownames = TRUE, bordered = TRUE, caption = "Shapiro-Wilks Results", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

output$ba_an_bpresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_an_run
    
    if (g_resetresult == FALSE || is.null(ba_an_results()$Bp))
        return()
    
    ba_an_results()$Bp
}, rownames = TRUE, bordered = TRUE, caption = "Breusch-Pagan Results", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

output$ba_an_kwresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_an_run
    
    if (g_resetresult == FALSE || is.null(ba_an_results()$Kw))
        return()
    
    ba_an_results()$Kw
}, bordered = TRUE, caption = "Kruskal-Wallis Results", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

# ANOVA Plot ####

output$ba_an_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$ba_an_run
    
    if (g_resetresult == FALSE || is.null(ba_an_results()$Plot1))
        return()

    local({
        output$ba_an_downloadplot1 <- downloadPlot(
            ba_an_results()$Plot1
        )
    })
    
    ba_an_results()$Plot1
})

output$ba_an_plot2 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$ba_an_run
    
    if (g_resetresult == FALSE || is.null(ba_an_results()$Plot2))
        return()
    
    local({
        output$ba_an_downloadplot2 <- downloadPlot(
            ba_an_results()$Plot2
        )
    })
    
    ba_an_results()$Plot2
})

output$ba_an_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$ba_an_run

    if (g_resetresult == FALSE || is.null(ba_an_results()$Plot1))
        return()
    
    div(
        style = "position: relative; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("ba_an_plot1"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left:0.5em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "ba_an_downloadplot1", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})

output$ba_an_showplot2 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$ba_an_run

    if (g_resetresult == FALSE || is.null(ba_an_results()$Plot2))
        return()
    
    div(
        style = "position: relative; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("ba_an_plot2"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left:0.5em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "ba_an_downloadplot2", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})