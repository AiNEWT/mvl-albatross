# T-Test Run ####

observeEvent(input$ba_tt_run, {
    g_resetresult <<- TRUE
})

ba_tt_results <- eventReactive(input$ba_tt_run, {
    ttX = NULL
    ttY = NULL
    ttMu <- input$ba_tt_mu
    ttPaired = FALSE
    ttEqual = FALSE
    ttType = input$ba_tt_type
    ttLmx = NULL
    ttLmy = NULL
    ttShapiroresults = NULL
    ttLeveneresults = NULL
    ttWilcoxresults = NULL
    ttDesc = NULL
    ttDescnames <- c("n", "mean", "median", "sd", "se")
    
    ttPlot1 = NULL
    ttPlot2 = NULL
    ttPlot3 = NULL
    
    # Checking Blank Inputs
    if (ttType == "1S" && input$ba_tt_varname1 == "") {
        showNotification("Please choose variable", type="warning")
        return()
    }
    else if (ttType == "2P" && (input$ba_tt_varname1 == "" || input$ba_tt_varname2 == "")) {
        showNotification("Please choose variable", type="warning")
        return()
    }
    else if (ttType == "2Up" && (input$ba_tt_varname1 == "" || input$ba_tt_varname2 == "" || 
             input$ba_tt_groupvarname1 == "" || input$ba_tt_groupvarname2 == "")) {
        showNotification("Please choose variable", type="warning")
        return()
    }
    
    # Start Progressbar
    withProgress(message = 't-test', style = "notification", value = 0, {
        
    # Confidence level between 0 and 1
    if (input$ba_tt_level > 1 || input$ba_tt_level < 0) {
        showNotification("Significance level must be a single number between 0 and 1.", type="warning")
        return()
    }
      
    # ¦§ Plot1 ####
    if (ttType == "1S") {
        ttX = data2[[input$ba_tt_varname1]]
        
        if (length(ttX) <= 2) {
            showNotification("There is not enough data in selected groups.", type="warning")
            return()
        }
        
        if (var(ttX) == 0) {
            showNotification("Data are essentially constant", type="warning")
            return()
        }
        
        ttLmx = lm(ttX ~ 1, data = data2)
        ttPlot1 <- autoplot(ttLmx)[[2]] + 
            ggtitle("Normal Probability Plot") + 
            theme_bw() + 
            theme(plot.title = element_text(hjust = 0.5))
        ttPlot2 <- NULL
    }
    else if (ttType == "2P") {
        ttX = data2[[input$ba_tt_varname1]]
        ttY = data2[[input$ba_tt_varname2]]
        
        if (length(ttX) <= 2 || length(ttY) <= 2) {
            showNotification("There is not enough data in selected groups.", type="warning")
            return()
        }
        
        if (var(ttX - ttY) == 0) {
            showNotification("Data are essentially constant", type="warning")
            return()
        }
        
        ttLmx = lm((ttX - ttY) ~ 1, data = data2)
        ttPlot1 <- autoplot(ttLmx)[[2]] + 
            ggtitle("Normal Probability Plot") + 
            theme_bw() + 
            theme(plot.title = element_text(hjust = 0.5))
        ttPlot2 <- NULL
    } else {
        ttX = data2[[input$ba_tt_varname1]][data2[[input$ba_tt_varname2]] == input$ba_tt_groupvarname1]
        ttY = data2[[input$ba_tt_varname1]][data2[[input$ba_tt_varname2]] == input$ba_tt_groupvarname2]
        
        if (length(ttX) <= 2 || length(ttY) <= 2) {
            showNotification("There is not enough data in selected groups.", type="warning")
            return()
        }
        
        if (var(ttX, na.rm = TRUE) == 0 || var(ttY, na.rm = TRUE) == 0) {
            showNotification("Data are essentially constant", type="warning")
            return()
        }
        
        ttLmx = lm(ttX ~ 1, data = data2)
        ttLmy = lm(ttY ~ 1, data = data2)
        ttPlot1 <- autoplot(ttLmx)[[2]] + 
            ggtitle(paste0("Normal Probability Plot of ", input$ba_tt_groupvarname1)) + 
            theme_bw() + 
            theme(plot.title = element_text(hjust = 0.5))
        ttPlot2 <- autoplot(ttLmy)[[2]] + 
            ggtitle(paste0("Normal Probability Plot of ", input$ba_tt_groupvarname2)) + 
            theme_bw() + 
            theme(plot.title = element_text(hjust = 0.5))
    }

    # Increase Progressbar
    incProgress(0.5, detail = paste("Loading..."))
        
    if (input$ba_tt_type == "2P")
        ttPaired = TRUE
    
    # ¦§ Shapiro-Wilk Test ####
    if (input$ba_tt_shapirocheck == FALSE) {
        ttShapiroresults <- NULL
    }
    else if (length(ttX) >= 5000 || length(ttY) >= 5000) {
        showNotification("Too many data in selected groups.", type="warning")
    }
    else if (ttType == "1S" || ttType == "2P") {
        ttShapiro <- shapiro.test(stdres(ttLmx))
        ttShapiroresults <-
            matrix(
                data = c(ttShapiro$statistic, ttShapiro$p.value),
                nrow = 1,
                ncol = 2
            )
        colnames(ttShapiroresults) <- c("W", "p-value")
        rownames(ttShapiroresults) <- "Shapiro-Wilks"
    } else {
        ttShapiro1 <- shapiro.test(stdres(ttLmx))
        ttShapiro2 <- shapiro.test(stdres(ttLmy))
        ttShapiroresults <-
            matrix(
                data = c(
                    ttShapiro1$statistic,
                    ttShapiro1$p.value,
                    ttShapiro2$statistic,
                    ttShapiro2$p.value
                ),
                nrow = 2,
                ncol = 2
            )
        colnames(ttShapiroresults) <- c("W", "p-value")
        rownames(ttShapiroresults) <- c(input$ba_tt_groupvarname1, input$ba_tt_groupvarname2)
    }
    
    # ¦§ Check Equal Variance ####
    if (input$ba_tt_type == "2Up" && input$ba_tt_equal == "TRUE") {
        ttEqual = TRUE
    }
    else if (input$ba_tt_type == "2Up" && input$ba_tt_equal == "Levene Result") {
        ttLevenematrix <- list(ttX,ttY)
        names(ttLevenematrix) <- c(input$ba_tt_groupvarname1, input$ba_tt_groupvarname2)
        ttLevenematrix <- reshape2::melt(ttLevenematrix)

        ttLeveneresults <- leveneTest(ttLevenematrix$value, as.factor(ttLevenematrix$L1))
        if (ttLeveneresults[["Pr(>F)"]][[1]] < 0.05)
            ttEqual = FALSE
        else
            ttEqual = TRUE
    }
    
    # ¦§ Levene Test ####
    if (!is.null(input$ba_tt_levenecheck) && input$ba_tt_levenecheck == TRUE && input$ba_tt_type == "2Up") {
        ttLevenematrix <- list(ttX,ttY)
        names(ttLevenematrix) <- c(input$ba_tt_groupvarname1, input$ba_tt_groupvarname2)
        ttLevenematrix <- reshape2::melt(ttLevenematrix)

        ttLeveneresults <- leveneTest(ttLevenematrix$value, as.factor(ttLevenematrix$L1))
    }
    
    # ¦§ Wilcoxon Test ####
    if (!is.null(input$ba_tt_wilcoxcheck) && input$ba_tt_wilcoxcheck == TRUE) {
        ttWilcox <- wilcox.test(
            x = ttX,
            y = ttY,
            alternative = input$ba_tt_alternative,
            mu = ttMu,
            paired = ttPaired,
            conf.level = 1 - input$ba_tt_level
        )
        ttWilcoxresults <-
            matrix(
                data = c(ttWilcox$statistic, ttWilcox$p.value),
                nrow = 1,
                ncol = 2
            )
        if (ttType == "2P")
            colnames(ttWilcoxresults) <- c("V", "p-value")
        else
            colnames(ttWilcoxresults) <- c("W", "p-value")
        rownames(ttWilcoxresults) <- "Wilcox"
    }
    
    # ¦§ Descriptive Statistics ####
    if (ttType == "1S") {
        ttDesc <- psych::describe(ttX)
    } else {
        ttXdesc <- psych::describe(ttX)
        ttYdesc <- psych::describe(ttY)
        ttDesc <- rbind(ttXdesc, ttYdesc)
    }
    
    if (ttType == "1S")
        rownames(ttDesc) <- input$ba_tt_varname1
    else if (ttType == "2P")
        rownames(ttDesc) <- c(input$ba_tt_varname1, input$ba_tt_varname2)
    else
        rownames(ttDesc) <- c(input$ba_tt_groupvarname1, input$ba_tt_groupvarname2)
    
    ttDesc <- data.frame(ttDesc)
    ttDesc <- subset(ttDesc, select=-c(vars, min, max, trimmed, mad, range, skew, kurtosis))
    ttDesc[["n"]] <- as.integer(ttDesc[["n"]])
    
    ttDescresults <- ttDesc[ttDescnames]
    
    # ¦§ Plot2 ####
    if (ttType == "1S") {
        ttEstmatrix <- reshape2::melt(data2[[input$ba_tt_varname1]])

        ttPlot3 <- ggplot(ttEstmatrix, aes(x=value)) + 
            geom_density(adjust = 0.1, alpha = .3) + 
            geom_vline(xintercept = ttMu, color = "red", size = 1.5) + 
            ggtitle("Density Plot") +
            theme_bw() + 
            theme(plot.title = element_text(hjust = 0.5))
    } else if (ttType == "2P") {
        ttEstmatrix <- reshape2::melt(data2[[input$ba_tt_varname1]]-data2[[input$ba_tt_varname2]])

        ttPlot3 <- ggplot(ttEstmatrix, aes(x=value)) + 
            geom_density(adjust = 0.1, alpha = .3) + 
            geom_vline(xintercept = ttMu, color = "red",size = 1.5) + 
            ggtitle("Density Plot") +
            theme_bw() + 
            theme(plot.title = element_text(hjust = 0.5))
    } else {
        ttEstmatrix <- list(ttX,ttY)
        names(ttEstmatrix) <- c(input$ba_tt_groupvarname1, input$ba_tt_groupvarname2)
        ttEstmatrix <- reshape2::melt(ttEstmatrix)

        ttPlot3 <- ggplot(ttEstmatrix, aes(x=value, fill=L1)) + 
            geom_density(adjust = 0.1, alpha = .3) + 
            ggtitle("Density Plot") +
            theme_bw() + 
            theme(plot.title = element_text(hjust = 0.5))
    }
    
    # ¦§ T-Test ####
    ttResults <- t.test(
        x = ttX,
        y = ttY,
        alternative = input$ba_tt_alternative,
        mu = ttMu,
        paired = ttPaired,
        var.equal = ttEqual,
        conf.level = 1 - input$ba_tt_level,
        na.rm = TRUE
    )
    
    ttConf.int <- ttResults$conf.int
    names(ttConf.int) <- c("lowerCI", "upperCI")
    ttValue <- c(ttResults$statistic, ttResults$parameter, p.value = ttResults$p.value, ttResults$estimate, se = ttResults$stderr, ttConf.int)
    
    # ¦§ Hypotheses ####
    hypoValue1 = NULL
    hypoValue2 = NULL
    if (input$ba_tt_type == "1S") {
        hypoValue1 = paste0("\\(\\mu_{", input$ba_tt_varname1, "}\\) = ", input$ba_tt_mu)
        if (input$ba_tt_alternative == "two.sided") {
            hypoValue2 = paste0("\\(\\mu_{", input$ba_tt_varname1, "}\\) \\(\\neq\\) ", input$ba_tt_mu)
        } else if (input$ba_tt_alternative == "less") {
            hypoValue2 = paste0("\\(\\mu_{", input$ba_tt_varname1, "}\\) < ", input$ba_tt_mu)
        } else {
            hypoValue2 = paste0("\\(\\mu_{", input$ba_tt_varname1, "}\\) > ", input$ba_tt_mu)
        }
    } else if (input$ba_tt_type == "2P") {
        hypoValue1 = paste0("\\(\\mu_{", input$ba_tt_varname1, "}\\) - \\(\\mu_{", input$ba_tt_varname2, "}\\) = ", input$ba_tt_mu)
        if (input$ba_tt_alternative == "two.sided") {
            hypoValue2 = paste0("\\(\\mu_{", input$ba_tt_varname1, "}\\) - \\(\\mu_{", input$ba_tt_varname2, "}\\) \\(\\neq\\) ", input$ba_tt_mu)
        } else if (input$ba_tt_alternative == "less") {
            hypoValue2 = paste0("\\(\\mu_{", input$ba_tt_varname1, "}\\) - \\(\\mu_{", input$ba_tt_varname2, "}\\) < ", input$ba_tt_mu)
        } else {
            hypoValue2 = paste0("\\(\\mu_{", input$ba_tt_varname1, "}\\) - \\(\\mu_{", input$ba_tt_varname2, "}\\) > ", input$ba_tt_mu)
        }
    } else {
        hypoValue1 = paste0("\\(\\mu_{", input$ba_tt_groupvarname1, "}\\) - \\(\\mu_{", input$ba_tt_groupvarname2, "}\\) = ", input$ba_tt_mu)
        if (input$ba_tt_alternative == "two.sided") {
            hypoValue2 = paste0("\\(\\mu_{", input$ba_tt_groupvarname1, "}\\) - \\(\\mu_{", input$ba_tt_groupvarname2, "}\\) \\(\\neq\\) ", input$ba_tt_mu)
        } else if (input$ba_tt_alternative == "less") {
            hypoValue2 = paste0("\\(\\mu_{", input$ba_tt_groupvarname1, "}\\) - \\(\\mu_{", input$ba_tt_groupvarname2, "}\\) < ", input$ba_tt_mu)
        } else {
            hypoValue2 = paste0("\\(\\mu_{", input$ba_tt_groupvarname1, "}\\) - \\(\\mu_{", input$ba_tt_groupvarname2, "}\\) > ", input$ba_tt_mu)
        }
    }
    
    ttHypothesis <- matrix(data = c("\\(\\ H_0 \\)", "\\(\\ H_1 \\)", hypoValue1, hypoValue2), nrow = 2, ncol = 2)
    
    # ¦§ Title ####
    titleOutput = NULL
    if (input$ba_tt_type == "1S") {
        titleOutput = "One Sample t-test"
    } else if (input$ba_tt_type == "2P") {
        titleOutput = "Paired t-test"
    } else {
        titleOutput = "Unpaired t-test"
    }
    
    # ¦¦ Remove & Append Tab ####
    removeTab(inputId = "ba_tt_resulttabset", target = "Normal Q-Q")
    if (input$ba_tt_shapirocheck == TRUE) {
        appendTab(inputId = "ba_tt_resulttabset",
            tabPanel("Normal Q-Q",
                br(),
                uiOutput("ba_tt_showplot1"),
                uiOutput("ba_tt_showplot2")
            )
        )
    }
    
    # Increase Progressbar
    incProgress(0.5, detail = paste("Loading..."))
        
    return(list(
        Title = titleOutput,
        Hypotheses = ttHypothesis,
        Desc = ttDescresults,
        Ttest = data.frame(t(ttValue)),
        Shapiro = ttShapiroresults,
        Levene = ttLeveneresults,
        Wilcox = ttWilcoxresults,
        Plot1 = ttPlot1,
        Plot2 = ttPlot2,
        Plot3 = ttPlot3
    ))
    
    }) # End Progressbar
})

# T-Test Components ####

output$ba_tt_type <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    shiny::radioButtons(
        "ba_tt_type",
        "Type of t-test",
        choices = list("One sample" = "1S", "Paired" = "2P", "Unpaired" = "2Up"),
        selected = "1S",
        inline = TRUE
    )
})

output$ba_tt_selectvarname1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if (is.null(input$ba_tt_type))
        return()
    
    label = NULL
    nameValue = names(select_if(data2, is.numeric))
    
    if (input$ba_tt_type == "2P")
        label = "Variable 1"
    else
        label = "Variable"
    

    div(
        selectInput(
            "ba_tt_varname1",
            label,
            choices = as.list(c("", nameValue)), 
            selected = NULL, 
            multiple = FALSE
        ),
        bsTooltip(
            "ba_tt_varname1", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$ba_tt_selectvarname2 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    label = NULL
    nameValue = NULL
    
    if (is.null(input$ba_tt_type) || input$ba_tt_type == "1S")
        return()
    else if (input$ba_tt_type == "2P") {
        label <- "Variable 2"
        nameValue = names(select_if(data2, is.numeric))
        nameValue <- nameValue[-which(nameValue==input$ba_tt_varname1)]
    }
    else {
        label <- "Group Variable"
        nameValue <- names(select_if(data2, is.factor))
        nameValue <- c(nameValue, names(select_if(data2, is.integer)))
        
        # is.null to ignore Error
        if (!is.null(input$ba_tt_varname1) && is.integer(data2[[input$ba_tt_varname1]])) {
            nameValue <- nameValue[-which(nameValue==input$ba_tt_varname1)]
        }
    }
    
    if (input$ba_tt_type == "2P") {
        div(
            selectInput(
                "ba_tt_varname2",
                label,
                choices = as.list(c("", nameValue)), 
                selected = NULL, 
                multiple = FALSE
            ),
            bsTooltip(
                "ba_tt_varname2", 
                "Numeric Only",
                "right", 
                options = list(container = "body")
            )
        )
    } else {
        selectInput(
            "ba_tt_varname2",
            label,
            choices = as.list(c("", nameValue)), 
            selected = NULL, 
            multiple = FALSE
        )
    }
    
})

output$ba_tt_groupvarname1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if (is.null(input$ba_tt_type) || input$ba_tt_type != "2Up")
        return()
    
    nameValue = NULL
    
    if (!is.null(input$ba_tt_varname2))
        nameValue = levels(as.factor(data2[[input$ba_tt_varname2]]))
    
    selectInput(
        "ba_tt_groupvarname1",
        "Group 1",
        choices = as.list(c("", nameValue)),
        selected = nameValue[1],
        multiple = FALSE
    )
})

output$ba_tt_groupvarname2 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if (is.null(input$ba_tt_type) || input$ba_tt_type != "2Up")
        return()
    
    nameValue = NULL
    
    if (!is.null(input$ba_tt_varname2)) {
        nameValue = levels(as.factor(data2[[input$ba_tt_varname2]]))
        nameValue <- nameValue[-which(nameValue==input$ba_tt_groupvarname1)]
    }

    selectInput(
        "ba_tt_groupvarname2",
        "Group 2",
        choices = as.list(c("", nameValue)),
        selected = nameValue[1],
        multiple = FALSE
    )
})

output$ba_tt_mu <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if (is.null(input$ba_tt_type))
        return()
    
    label = NULL
    
    if (input$ba_tt_type == "1S" || input$ba_tt_type == "2Up")
        label = "Null Value"
    else if (input$ba_tt_type == "2P")
        label = "Null Difference"
    
    numericInput("ba_tt_mu", label, value = 0)
})

output$ba_tt_equal <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if (is.null(input$ba_tt_type) || input$ba_tt_type != "2Up")
        return()
    
    selectInput(
        "ba_tt_equal",
        "Equal Variances",
        choices = list("TRUE", "FALSE", "Levene Result"),
        selected = NULL, 
        multiple = FALSE
    )
})

output$ba_tt_alternative <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    selectInput(
        "ba_tt_alternative",
        "Alternative Hypothesis",
        choices = list("two.sided", "less", "greater"),
        selected = NULL, 
        multiple = FALSE
    )
})

output$ba_tt_level <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    numericInput("ba_tt_level", "Significance Level", value = 0.05, min = 0, max = 1, step = 0.05)
})

output$ba_tt_shapirocheck <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("ba_tt_shapirocheck", "Shapiro-Wilk Test (Normality)")
})

output$ba_tt_levenecheck <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if (is.null(input$ba_tt_type) || input$ba_tt_type != "2Up")
        return()
    
    checkboxInput("ba_tt_levenecheck", "Levene Test (Variance Equality)")
})

output$ba_tt_wilcoxcheck <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if (is.null(input$ba_tt_type))
        return()
    
    label = NULL
    
    if (input$ba_tt_type == "1S")
        label = "Wilcoxon Test"
    else if (input$ba_tt_type == "2P")
        label = "Wilcoxon Signed Rank Test"
    else if (input$ba_tt_type == "2Up")
        label = "Wilcoxon Rank Sum Test"
    
    checkboxInput("ba_tt_wilcoxcheck", label)
})

# T-Test Results ####

output$ba_tt_title <- renderText({
    input$file
    input$resetData
    input$di_option_run
    input$ba_tt_run
    
    if (g_resetresult == FALSE || is.null(ba_tt_results()$Title)) {
        return("t-test")
    }

    ba_tt_results()$Title
})

output$ba_tt_hypotheses <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$ba_tt_run
    
    tagList(
        withMathJax(),
        withMathJax(tableOutput("ba_tt_hypothesestable"))
    )
})

output$ba_tt_hypothesestable <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_tt_run
    
    if (g_resetresult == FALSE || is.null(ba_tt_results()$Hypotheses))
        return()
    
    ba_tt_results()$Hypotheses
}, rownames = FALSE, colnames = FALSE, bordered = TRUE, caption = "Hypotheses", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

output$ba_tt_Desc <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_tt_run
    
    if (g_resetresult == FALSE || is.null(ba_tt_results()$Desc))
        return()
    
    ba_tt_results()$Desc
}, rownames = TRUE, bordered = TRUE, caption = "Descriptive Statistics", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

output$ba_tt_ttestresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_tt_run
    
    if (g_resetresult == FALSE || is.null(ba_tt_results()$Ttest))
        return()
    
    ba_tt_results()$Ttest
}, bordered = TRUE, caption = "t-test", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

output$ba_tt_shapiroresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_tt_run
    
    if (g_resetresult == FALSE || is.null(ba_tt_results()$Shapiro))
        return()
    
    ba_tt_results()$Shapiro
}, rownames = TRUE, bordered = TRUE, caption = "Shapiro-Wilks Results", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

output$ba_tt_leveneresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_tt_run
    
    if (g_resetresult == FALSE || is.null(ba_tt_results()$Levene))
        return()
    
    ba_tt_results()$Levene
}, bordered = TRUE, caption = "Levene Results", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

output$ba_tt_wilcoxresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_tt_run
    
    if (g_resetresult == FALSE || is.null(ba_tt_results()$Wilcox))
        return()
    
    ba_tt_results()$Wilcox
}, bordered = TRUE, caption = "Wilcoxon Test", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

# T-Test Plots ####

output$ba_tt_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$ba_tt_run

    if (g_resetresult == FALSE || is.null(ba_tt_results()$Plot1))
        return()

    local({
        output$ba_tt_downloadplot1 <- downloadPlot(
            ba_tt_results()$Plot1
        )
    })
    
    ba_tt_results()$Plot1
})

output$ba_tt_plot2 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$ba_tt_run

    if (g_resetresult == FALSE || is.null(ba_tt_results()$Plot2))
        return()

    local({
        output$ba_tt_downloadplot2 <- downloadPlot(
            ba_tt_results()$Plot2
        )
    })
    
    ba_tt_results()$Plot2
})

output$ba_tt_plot3 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$ba_tt_run

    if (g_resetresult == FALSE || is.null(ba_tt_results()$Plot3))
        return()

    local({
        output$ba_tt_downloadplot3 <- downloadPlot(
            ba_tt_results()$Plot3
        )
    })
    
    ba_tt_results()$Plot3
})

output$ba_tt_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$ba_tt_run

    if (g_resetresult == FALSE || is.null(ba_tt_results()$Plot1))
        return()
    
    div(
        style = "position: relative; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("ba_tt_plot1"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left:0.5em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "ba_tt_downloadplot1", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})

output$ba_tt_showplot2 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$ba_tt_run

    if (g_resetresult == FALSE || is.null(ba_tt_results()$Plot2))
        return()
    
    div(
        style = "position: relative; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("ba_tt_plot2"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left:0.5em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "ba_tt_downloadplot2", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})

output$ba_tt_showplot3 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$ba_tt_run

    if (g_resetresult == FALSE || is.null(ba_tt_results()$Plot3))
        return()
    
    div(
        style = "position: relative; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("ba_tt_plot3"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left:0.5em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "ba_tt_downloadplot3", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})