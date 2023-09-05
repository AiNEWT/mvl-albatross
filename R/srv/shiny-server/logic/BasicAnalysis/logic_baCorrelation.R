# Correlation Run ####

observeEvent(input$ba_co_run, {
    g_resetresult <<- TRUE
})

ba_co_results <- eventReactive(input$ba_co_run, {
    if (length(input$ba_co_varname1$right) < 2) {
        showNotification("Please select at least 2 data", type = "warning")
        return()
    }
    
    # Start Progressbar
    withProgress(message = 'Correlation Analysis', style = "notification", value = 0, {
        
    coDataframe <- data2[input$ba_co_varname1$right]

    coCofftype <- NULL
    coCorrtest <- NULL
    coPlottype <- NULL
    coPlot <- NULL
    
    coTitle <- NULL
    coRresults <- NULL
    coNresults <- NULL
    coPresults <- NULL
    
    # ¦§ Title ####
    coTitle <- paste(input$ba_co_varname1$right, collapse = ', ')
    
    # ¦§ Correlation Analysis ####
    if (input$ba_co_cofftype == "Pearson")
        coCofftype = "pearson"
    else if (input$ba_co_cofftype == "Kendall")
        coCofftype = "kendall"
    else if (input$ba_co_cofftype == "Spearman")
        coCofftype = "spearman"
        
    coCoresults <- corr.test(coDataframe, method = coCofftype)
    coRresults <- coCoresults[["r"]]
    coNresults <- coCoresults[["n"]]
    if (length(coNresults) == 1) {
        coNresults <- data.frame(coNresults)
        rownames(coNresults) <- "Total"
        colnames(coNresults) <- "Count"
    }
    coPresults <- abs(coCoresults[["p"]])
    
    # Increase Progressbar
    incProgress(0.5, detail = paste("Loading..."))
    
    # ¦¦ Plot ####
    if (input$ba_co_plottype == "Scatter") {
        lm_eqn <- function(df){
            m <- lm(y ~ x, df)
            eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                             list(a = format(unname(coef(m)[1]), digits = 2),
                                  b = format(unname(coef(m)[2]), digits = 2),
                                  r2 = format(summary(m)$r.squared, digits = 3)))
            as.character(as.expression(eq))
        }
        
        lowerFn <- function(data, mapping, method = "lm", showFormula = FALSE, ...) {
            p <- ggplot(data = data, mapping = mapping) +
                geom_point() +
                geom_smooth(method = method, color = "blue", ...)
            if(showFormula == TRUE)
                p = p + geom_text(x = max(data[[1]]), y = max(data[[2]]), 
                                  label = lm_eqn(data), parse = TRUE)
            p
        }
        bins = 30
        check_formula = FALSE
        
        if(!is.null(input$ba_co_bins)) bins = input$ba_co_bins
        if(!is.null(input$ba_co_check_formula)) check_formula = input$ba_co_check_formula
        
        coPlot1 = ggpairs(coDataframe, 
                          lower = list(continuous = wrap(lowerFn, 
                                                         method = "lm", 
                                                         showFormula = check_formula)),
                          diag = list(continuous = wrap("barDiag", bins = bins))) + 
            theme_bw() +
            ggtitle("Correlation Plot") + 
            theme(plot.title = element_text(hjust = 0.5))
    }
    else if (input$ba_co_plottype == "Heat Map") {
        if (sum(is.na(coRresults)) >= 1) {
            showNotification("Standard deviation is zero", type = "warning")
            return()
        }
        coPlot1 <- ggcorrplot(coRresults, hc.order = TRUE, outline.col = "white") +
            ggtitle("Correlation Plot") + 
            theme(plot.title = element_text(hjust = 0.5))
    }
    else if (input$ba_co_plottype == "Heat Map (with Value)") {
        if (sum(is.na(coRresults)) >= 1) {
            showNotification("Standard deviation is zero", type = "warning")
            return()
        }
        coPlot1 <- ggcorrplot(coRresults, hc.order = TRUE, outline.col = "white", lab = TRUE) + 
            ggtitle("Correlation Plot") + 
            theme(plot.title = element_text(hjust = 0.5))
    }
    
    # Increase Progressbar
    incProgress(0.5, detail = paste("Loading..."))
    
    return(list(
        Title = coTitle,
        Covalue = coRresults,
        Nvalue = coNresults,
        Pvalue = coPresults,
        Plot1 = coPlot1
    ))
    
    }) # End Progressbar
})

# Correlation Components ####
output$ba_co_selectvarname1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue = names(select_if(data2, is.numeric))
    
    div(
        chooserInput(
            "ba_co_varname1", 
            "Variable", 
            "Selected", 
            nameValue, 
            c(), 
            size = 15, 
            multiple = TRUE, 
            width = 150
        ),
        bsTooltip(
            "ba_co_varname1", 
            "Numeric Only",
            "bottom", 
            options = list(container = "body")
        )
    )
})

output$ba_co_cofftype <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    selectInput(
        "ba_co_cofftype",
        "Coefficient Type",
        choices = list("Pearson", "Kendall", "Spearman"),
        selected = NULL,
        multiple = FALSE
    )
})

output$ba_co_plottype <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    selectInput(
        "ba_co_plottype",
        "Plot Type",
        choices = list("Scatter", "Heat Map", "Heat Map (with Value)"),
        selected = NULL,
        multiple = FALSE
    )
})

output$ba_co_bins <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$ba_co_plottype
    
    if (is.null(input$ba_co_plottype) || input$ba_co_plottype != "Scatter")
        return()
    
    sliderInput(
        inputId = sprintf("ba_co_bins"),
        label = "Number of Groups for Histogram",
        value = 30,
        min = 1,
        max = 100
    )
})

# output$ba_co_check_formula <- renderUI({
#     input$file
#     input$resetData
#     input$di_option_run
#     input$ba_co_plottype
#     
#     if (is.null(input$ba_co_plottype) || input$ba_co_plottype != "Scatter")
#         return()
#     
#     checkboxInput(
#         inputId = 'ba_co_check_formula',
#         label = 'Report Regression line formula',
#         value = FALSE
#     )
# })

# Correlation Results ####

output$ba_co_title <- renderText({
    input$file
    input$resetData
    input$di_option_run
    input$ba_co_run
    
    titleOutput = NULL
    if (g_resetresult == FALSE || is.null(ba_co_results()$Title))
        titleOutput = ""
    else
        titleOutput = ba_co_results()$Title
    
    paste("Variable : ", titleOutput)
})

output$ba_co_coresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_co_run
    
    if (g_resetresult == FALSE || is.null(ba_co_results()$Covalue))
        return()
    
    ba_co_results()$Covalue
}, rownames = TRUE, bordered = TRUE, caption = "Correlation Matrix", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

output$ba_co_nresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_co_run
    
    if (g_resetresult == FALSE || is.null(ba_co_results()$Nvalue))
        return()
    
    ba_co_results()$Nvalue
}, rownames = TRUE, bordered = TRUE, caption = "Sample Size", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 0)

output$ba_co_presults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_co_run
    
    if (g_resetresult == FALSE || is.null(ba_co_results()$Pvalue))
        return()
    
    ba_co_results()$Pvalue
}, rownames = TRUE, bordered = TRUE, caption = "p-values", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

# Correlation Plots ####



output$ba_co_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$ba_co_run

    if (g_resetresult == FALSE || is.null(ba_co_results()$Plot1))
        return()
    
    local({
        output$ba_co_downloadplot1 <- downloadPlot(
            ba_co_results()$Plot1
        )
    })
    
    ba_co_results()$Plot1
})

output$ba_co_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$ba_co_run

    if (g_resetresult == FALSE || is.null(ba_co_results()$Plot1))
        return()
    
    div(
        style = "position: relative; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("ba_co_plot1"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left:0.5em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "ba_co_downloadplot1", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})


