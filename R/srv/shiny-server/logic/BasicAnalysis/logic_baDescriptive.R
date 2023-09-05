observeEvent(input$ba_ds_help, {
    rintrojs::introjs(session, options = list(
        steps = data.frame(
            element = c(
                NA, 
                "#ba_ds_selectvarname1", 
                "#ba_ds_check_quantile",
                "#ba_ds_select_quantile",
                "#ba_ds_checkgroup", 
                "#ba_ds_selectvarname2", 
                "#ba_ds_run"
            ),
            intro = c(
                "Hello there this is Descriptive Statistics.",
                "Choose variable.",
                "If you want choose group variable. you need to check this combobox.",
                "Choose group variable.",
                "Press run button!"
            ),
            position = c("right")
        )
    ))
})

# Descriptive Run ####

observeEvent(input$ba_ds_run, {
    g_resetresult <<- TRUE
})

ba_ds_results <- eventReactive(input$ba_ds_run, {
    if (input$ba_ds_varname1 == "") {
        showNotification("Please choose variable", type="warning")
        return()
    }
    
    if (input$ba_ds_checkgroup == TRUE && is.null(input$ba_ds_varname2)) {
        showNotification("Please choose group variable", type="warning")
        return()
    }
    
    # Start Progressbar
    withProgress(message = 'Descriptive Statistics', style = "notification", value = 0, {
        
    removeTab(inputId = "ba_ds_resulttabset", target = "Group Histogram")
    removeTab(inputId = "ba_ds_resulttabset", target = "Group Box-plot")

    # ¦§ Variable Results ####
    dsVariableresults <- psych::describe(data2[[input$ba_ds_varname1]])
    dsVariableresults <- data.frame(dsVariableresults)[c("n", "mean", "min", "median", "max", "sd", "se")]
    dsVariableresults <- cbind(length(data2[[input$ba_ds_varname1]]),
                               sum(is.na(data2[[input$ba_ds_varname1]])),
                               dsVariableresults[1:3],
                               as.vector(summary(data2[[input$ba_ds_varname1]]))[2],
                               dsVariableresults[4],
                               as.vector(summary(data2[[input$ba_ds_varname1]]))[5],
                               dsVariableresults[5:7]
                               )
    colnames(dsVariableresults)[c(1,2,6,8)] = c("N", "missing data", "1st Qu.", "3rd Qu.")
    rownames(dsVariableresults) <- input$ba_ds_varname1
    
    # ¦§ Percentile Results ####
    
    
    dsPercentileresults <- NULL
    if(!is.null(input$ba_ds_check_percentile) && input$ba_ds_check_percentile && !is.null(input$ba_ds_select_percentile)) {
        
        try(eval(parse(text = paste0("dsPercentileresults = t(as.data.frame(quantile(
                                     x = data2[[input$ba_ds_varname1]], probs = ", 
                                     input$ba_ds_select_percentile, 
                                     ", na.rm = TRUE)))")
        )))
        if(!is.null(dsPercentileresults)) rownames(dsPercentileresults) <- input$ba_ds_varname1
    }
    
    # ¦§ Variable Plot ####
    var=data2[input$ba_ds_varname1]
    dsPlot <- ggplot(var, aes_string(x = input$ba_ds_varname1)) + 
        geom_histogram(bins = input$ba_ds_bins) + 
        ggtitle(paste0("Histogram of ", input$ba_ds_varname1)) + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    
    # ¦§ Variable Box-Plot ####
    
    var=data2[input$ba_ds_varname1]
    dsBoxplot <- ggplot(var, aes_string(y = input$ba_ds_varname1)) + 
        geom_boxplot() + 
        ggtitle(paste0("Box-plot of ", input$ba_ds_varname1)) + 
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    

    if (input$ba_ds_checkgroup == FALSE) {
        
        # Increase Progressbar
        incProgress(1.0, detail = paste("Loading..."))
        
        return(list(
            Variable = input$ba_ds_varname1,
            Variabletable = dsVariableresults,
            Percentile = dsPercentileresults,
            Plot1 = dsPlot,
            Boxplot = dsBoxplot
        ))
    }
    
        
    # ¦§ Group Results 1 ####
    dsDataframe <- data.frame(data2)
    dsNames <- c(input$ba_ds_varname2, "n", "mean", "min", "median", "max", "sd", "se")
    
    dsResults <- psych::describeBy(dsDataframe[input$ba_ds_varname1], group = dsDataframe[input$ba_ds_varname2])
    
    # dsResults2
    eval(parse(text = paste0(
        "dsResults2 = data2 %>% dplyr::group_by(", 
        paste0(input$ba_ds_varname2, collapse = ", "), 
        ") %>% dplyr::summarize(missing = sum(is.na(", 
        input$ba_ds_varname1, 
        ")), Q1 = quantile(", 
        input$ba_ds_varname1, 
        ", 0.25, na.rm = TRUE), Q3 = quantile(", 
        input$ba_ds_varname1, 
        ", 0.75, na.rm = TRUE))"
    )))
    
    dsResults2 = na.omit(dsResults2)

    dims <- attr(dsResults, "dimnames")
    dsMatrix <- NULL
    dsNA <- c(NA, 0, rep(NA, 11))
    names(dsNA) <- c("vars", "n", "mean", "sd", "median", "trimmed", "mad", "min", "max", "range", "skew", "kurtosis", "se")

    for (i in 1:length(dsResults)) {
        if (is.null(dsResults[[i]])) {
            dsMatrix <- rbind(dsMatrix, dsNA)
        } else {
            dsMatrix <- rbind(dsMatrix, dsResults[[i]])
        }
    }
    
    # Increase Progressbar
    incProgress(0.5, detail = paste("Loading..."))
    
    # ¦§ Group Plot ####
    
    dimLength <- length(dims)
    groupNamesMatrix <- matrix(c(rep("", length(dsResults))), nrow=length(dsResults), ncol=dimLength)
    maxChangecount <- matrix(c(rep(0, dimLength)), nrow=dimLength, ncol=1)
    maxChangecount[1] <- 1
    if (dimLength > 1) {
        for (i in 2:dimLength) {
            maxChangecount[i] <- length(c(dims[[i - 1]])) * maxChangecount[i - 1]
        } 
    }
    
    g_ds_grouphistselect <<- matrix(c(rep("", length(dsResults))), nrow=length(dsResults), ncol=1)
    ds_histselectnames <- c(input$ba_ds_varname2)
    for (i in 1:dimLength) {
        maxDimcount <- length(dims[[i]])
        dimCount <- 1
        changeCount <- 0
        
        for (j in 1:length(dsResults)) {
            if (g_ds_grouphistselect[j] == "") {
                g_ds_grouphistselect[j] <<- paste0(ds_histselectnames[i], " : ", dims[[i]][dimCount])
            } else {
                g_ds_grouphistselect[j] <<- paste0(g_ds_grouphistselect[j], " / ", ds_histselectnames[i], " : ", dims[[i]][dimCount])
            }
            groupNamesMatrix[j, i] <- dims[[i]][dimCount]
            
            changeCount <- changeCount + 1
            if (changeCount == maxChangecount[i]) {
                dimCount <- dimCount + 1
                changeCount <- 0
            }
            
            if (dimCount > maxDimcount) {
                dimCount <- 1
            }
        }
    }
    
    # ¦§ Group Results 2 ####
    
    colnames(groupNamesMatrix) <- c(input$ba_ds_varname2)
    rownames(groupNamesMatrix) <- NULL
    rownames(dsMatrix) <- NULL
    dsResults <- cbind(dsMatrix, groupNamesMatrix)
    
    dsResults[["n"]] <- as.integer(dsResults[["n"]])
    dsResults <- subset(dsResults, select=-c(vars, trimmed, mad, range, skew, kurtosis))
    dsResults <- dsResults[dsNames]
    
    dsResults2$N = dsResults$n + dsResults2$missing
    
    dsResults = cbind(dsResults[,1:dimLength], 
                      dsResults2$N, 
                      dsResults2$missing, 
                      dsResults[,dimLength + 1:3],
                      dsResults2$Q1, 
                      dsResults[,dimLength + 4],
                      dsResults2$Q3, 
                      dsResults[,dimLength + 5:7])
    
    colnames(dsResults) = c(input$ba_ds_varname2, "N", "missing data", "n", "mean", 
                            "min", "1st Qu.", "median", "3rd Qu.", "max", "sd", "se")

    # ¦§ Group Percentile Results ####
    
    dsGroupPercentileresults = NULL
    
    if(!is.null(input$ba_ds_check_percentile) && 
       input$ba_ds_check_percentile && 
       !is.null(input$ba_ds_select_percentile)) {
        
        temp_text = NULL
        try(eval(parse(text = paste0(
            "temp_text = ", 
            input$ba_ds_select_percentile
        ))))
        
        if(!is.null(temp_text)) {
            try(eval(parse(text = paste0(
                "dsGroupPercentileresults = data2 %>% group_by(", 
                paste0(input$ba_ds_varname2, collapse = ", "), 
                ") %>% dplyr::summarize(quantile(", 
                input$ba_ds_varname1, 
                ", probs = ", 
                paste0(temp_text, 
                       collapse = paste0(
                           ", na.rm = TRUE), quantile(", 
                           input$ba_ds_varname1, 
                           ", probs = "
                       )), 
                ", na.rm = TRUE))" 
            ))))
        }
        
        if(!is.null(dsGroupPercentileresults)) {
            colnames(dsGroupPercentileresults) <- c(input$ba_ds_varname2, colnames(dsPercentileresults))
            #dsGroupPercentileresults = na.omit(dsGroupPercentileresults)
        }
    }
    
    # ¦§ Group Box-Plot ####
    
    dsGroupBoxplot = NULL
    if(!is.null(input$ba_ds_check_groupboxplot) &&
       input$ba_ds_check_groupboxplot && dimLength == 1) {
        dsGroupBoxplot <- ggplot(data = data2, 
                                 mapping = aes_string(
                                     x = input$ba_ds_varname2,
                                     y = input$ba_ds_varname1, 
                                     group = input$ba_ds_varname2)) + 
            geom_boxplot() + 
            ggtitle(paste0("Box-plot of ", input$ba_ds_varname1)) + 
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5))
    }
    
    # Increase Progressbar
    incProgress(0.5, detail = paste("Loading..."))
    
    # ¦¦ Group Append Tab ####
    appendTab(inputId = "ba_ds_resulttabset",
        tabPanel(
            "Group Histogram",
            br(),
            fluidRow(
                column(
                    9,
                    uiOutput("ba_ds_grouphistselect")
                ),
                column(
                    3,
                    uiOutput("ba_ds_grouphistappend"),
                    style = "text-align:right; padding:15px"
                )
            ),
            uiOutput("ba_ds_showgrouphist")
        )
    )
    
    if(!is.null(dsGroupBoxplot)) {
        appendTab(inputId = "ba_ds_resulttabset",
                  tabPanel(
                      "Group Box-plot",
                      br(),
                      uiOutput("ba_ds_showgroupboxplot")
                  )
        )
    }
    
    return(list(
        Variable = input$ba_ds_varname1,
        Variabletable = dsVariableresults,
        Percentile = dsPercentileresults,
        Groupnames = groupNamesMatrix,
        Grouptable = dsResults,
        GroupPercentile = dsGroupPercentileresults,
        Plot1 = dsPlot,
        Boxplot = dsBoxplot, 
        GroupBoxplot = dsGroupBoxplot
    ))
    
    }) # End Progressbar
})

# Descriptive Components ####

output$ba_ds_selectvarname1 <- renderUI({
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
            "ba_ds_varname1",
            "Variable",
            choices = as.list(c("", nameValue)), 
            selected = NULL, 
            multiple = FALSE
        ),
        bsTooltip(
            "ba_ds_varname1", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$ba_ds_check_percentile <- renderUI({
    input$file
    input$resetData
    input$di_option_run

    checkboxInput(
        inputId = 'ba_ds_check_percentile',
        label = 'Report Percentiles',
        value = FALSE
    )
})

output$ba_ds_select_percentile <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$ba_ds_check_percentile

     if (is.null(input$ba_ds_check_percentile) || input$ba_ds_check_percentile == FALSE)  return()

     values="c(0,0.25,0.5,0.75,1)"
     textAreaInput("ba_ds_select_percentile", "Percentiles", value = values, height = "80px")
})


output$ba_ds_checkgroup <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput(
        inputId = 'ba_ds_checkgroup',
        label = 'Use Group',
        value = FALSE
    )
})

output$ba_ds_selectvarname2 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if (is.null(input$ba_ds_checkgroup) || input$ba_ds_checkgroup == FALSE)
        return()
    
    nameValue <- names(select_if(data2, is.factor))
    nameValue <- c(nameValue, names(select_if(data2, is.integer)))

    # is.null to ignore Error
    if (!is.null(input$ba_ds_varname1) && is.integer(data2[[input$ba_ds_varname1]])) {
        nameValue <- nameValue[-which(nameValue==input$ba_ds_varname1)]
    }
    
    selectInput(
        "ba_ds_varname2",
        "Groups",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$ba_ds_check_groupboxplot <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if (is.null(input$ba_ds_checkgroup) || input$ba_ds_checkgroup == FALSE)
        return()
    
    checkboxInput(
        inputId = 'ba_ds_check_groupboxplot',
        label = 'Show Group Box-Plot (only for 1 group)',
        value = FALSE
    )
})

# Descriptive Results ####

output$ba_ds_title <- renderUI({
    input$ba_ds_run
    
    titleOutput = NULL
    if (g_resetresult == FALSE || is.null(ba_ds_results()$Variable))
        titleOutput = ""
    else
        titleOutput = ba_ds_results()$Variable
    
    paste("Variable : ", titleOutput)
})

output$ba_ds_variableresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_ds_run
    
    if (g_resetresult == FALSE || is.null(ba_ds_results()$Variabletable))
        return()
    
    ba_ds_results()$Variabletable
}, rownames = TRUE, bordered = TRUE, caption = "Variable Results", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)


output$ba_ds_percentileresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_ds_run
    
    if (g_resetresult == FALSE || is.null(ba_ds_results()$Percentile))
        return()
    
    ba_ds_results()$Percentile
}, rownames = TRUE, bordered = TRUE, caption = "Percentile Results", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)



output$ba_ds_groupresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_ds_run
    
    if (g_resetresult == FALSE || is.null(ba_ds_results()$Grouptable))
        return()
    
    ba_ds_results()$Grouptable
}, bordered = TRUE, caption = "Group Results", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

output$ba_ds_grouppercentileresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$ba_ds_run
    
    if (g_resetresult == FALSE || is.null(ba_ds_results()$GroupPercentile))
        return()
    
    ba_ds_results()$GroupPercentile
}, rownames = FALSE, bordered = TRUE, caption = "Group Percentile Results", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

# Descriptive Plot ####

output$ba_ds_variablehist <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$ba_ds_run
    
    if (g_resetresult == FALSE || is.null(ba_ds_results()$Variable))
        return()
    
    local({
        output$ba_ds_downloadvariablehist <- downloadPlot(
            ba_ds_results()$Plot1
        )
    })
    
    ba_ds_results()$Plot1
})

output$ba_ds_showvariablehist <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$ba_ds_run

    if (g_resetresult == FALSE || is.null(ba_ds_results()$Variable))
        return()
    
    div(
        style = "position: relative; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("ba_ds_variablehist"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left:0.5em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "ba_ds_downloadvariablehist", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})

output$ba_ds_boxplot <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$ba_ds_run
    
    if (g_resetresult == FALSE || is.null(ba_ds_results()$Variable))
        return()
    
    local({
        output$ba_ds_downloadboxplot <- downloadPlot(
            ba_ds_results()$Boxplot
        )
    })
    
    ba_ds_results()$Boxplot
})

output$ba_ds_showboxplot <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$ba_ds_run
    
    if (g_resetresult == FALSE || is.null(ba_ds_results()$Variable))
        return()
    
    div(
        style = "position: relative; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("ba_ds_boxplot"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left:0.5em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "ba_ds_downloadboxplot", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})



output$ba_ds_grouphistselect <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$ba_ds_run
    
    if (g_resetresult == FALSE || is.null(ba_ds_results()$Grouptable))
        return()
    
    choicehist <- g_ds_grouphistselect
    for (i in 1:nrow(ba_ds_results()$Grouptable))
        if (ba_ds_results()$Grouptable[["n"]][i] <= 1)
            choicehist[i] = ""
        
    choicehist = subset(choicehist, choicehist != "")
    
    selectInput("ba_ds_grouphistselect","Select Histogram", choices=as.list(choicehist), width="100%", multiple = TRUE)
})

output$ba_ds_grouphistappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$ba_ds_run
    
    if (g_resetresult == FALSE || is.null(ba_ds_results()$Grouptable))
        return()
    
    actionButton("ba_ds_grouphistappend", "Append Histogram")
})

ba_ds_grouphistresult <- eventReactive(input$ba_ds_grouphistappend, {
    g_ds_grouphist <<- input$ba_ds_grouphistselect
    g_ds_grouphistindex <<- NULL
    
    for (i in 1:length(g_ds_grouphist))
        g_ds_grouphistindex <<- c(g_ds_grouphistindex, which(g_ds_grouphistselect==g_ds_grouphist[i]))
    
    return(list(
        g_ds_grouphist = g_ds_grouphist
    ))
})

output$ba_ds_bins <- renderUI({
#  if(is.null(input$dm_mvi_check_tools) || input$dm_mvi_check_tools == FALSE)
#      return()
  sliderInput(
    inputId = sprintf("ba_ds_bins"),
    label = "Number of Groups",
    value = 30,
    min = 1,
    max = 100
  )
})

output$ba_ds_min<-renderUI({
    textAreaInput("ba_ds_min","Minimum",value=NULL)
})

output$ba_ds_max<-renderUI({
    textAreaInput("ba_ds_max","Maxmum",value=NULL)
})

ba_ds_grouphist <- reactive({
    input$ba_ds_run
    input$ba_ds_bins
    input$ba_ds_min
    input$ba_ds_max
    
    if (g_resetresult == FALSE || is.null(ba_ds_results()$Grouptable))
        return()
    
    if (is.null(ba_ds_grouphistresult()$g_ds_grouphist))
        return()
    
    grouphist <- vector(mode = "list", length = length(ba_ds_grouphistresult()$g_ds_grouphist))
    
    for (k in 1:length(ba_ds_grouphistresult()$g_ds_grouphist)) {
        selectedGroupindex = g_ds_grouphistindex[k]
        selectedGroupvalue = ba_ds_results()$Groupnames[selectedGroupindex,]
        selectedGroupname = c(ba_ds_results()$Variable, names(selectedGroupvalue))
        selectedGroupmatrix <- data2[selectedGroupname]
        selectedGrouptitle = paste("Group : ", g_ds_grouphistselect[selectedGroupindex])
        for(i in 1:length(selectedGroupvalue)) {
            number = 1 + i
            selectedGroupmatrix <- subset(selectedGroupmatrix, selectedGroupmatrix[number]==selectedGroupvalue[[i]])
        }
    
        # Histogram
        # hist(selectedGroupmatrix[[1]], main=g_ds_grouphistselect[selectedGroupindex], xlab=selectedGroupname[1])
        var=ba_ds_results()$Variable
        if (is.null(input$ba_ds_min) == FALSE)  var=var[var>=is.numeric(input$ba_ds_min)]
        if (is.null(input$ba_ds_max) == FALSE)  var=var[var<=is.numeric(input$ba_ds_max)]
        grouphist[[k]] <- ggplot(selectedGroupmatrix, aes_string(x=var)) + 
            geom_histogram(bins = input$ba_ds_bins) + 
            ggtitle(paste0("Histogram of ", selectedGrouptitle)) + 
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5))
    }

    return(grouphist)
})

lapply(1:100, function(k) {
    output[[paste0('ba_ds_grouphist_', k)]] <- renderPlot({
        input$ba_ds_run
        
        local({
            output[[paste0('ba_ds_downloadgrouphist_', k)]] <- downloadPlot(
                ba_ds_grouphist()[[k]]
            )
        })
        
        ba_ds_grouphist()[[k]]
    })
})

output$ba_ds_showgrouphist <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$ba_ds_run
    
    if (g_resetresult == FALSE || is.null(ba_ds_results()$Grouptable))
        return()
    
    if (is.null(ba_ds_grouphistresult()$g_ds_grouphist))
        return()
    
    do.call(tagList, c(lapply(1:length(ba_ds_grouphistresult()$g_ds_grouphist), function(k) {
        div(
            style = "position: relative; border: 1px solid #D3D3D3;",
            withSpinner(
                plotOutput(paste0("ba_ds_grouphist_", k)),
                type = 1,
                color = "#2c3e50",
                size = 1.2
            ),
            div(
                style = "position: absolute; left:0.5em; bottom: 0.5em;",
                dropdown(
                    downloadButton(outputId = paste0('ba_ds_downloadgrouphist_', k), label = "Download Plot"),
                    size = "xs",
                    icon = icon("download", class = "opt"),
                    up = TRUE
                )
            )
        )
    })))
})


# Group Box-plot

output$ba_ds_groupboxplot <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$ba_ds_run
    
    if (g_resetresult == FALSE || is.null(ba_ds_results()$Variable) || is.null(ba_ds_results()$GroupBoxplot))
        return()
    
    local({
        output$ba_ds_downloadgroupboxplot <- downloadPlot(
            ba_ds_results()$GroupBoxplot
        )
    })
    
    ba_ds_results()$GroupBoxplot
})

output$ba_ds_showgroupboxplot <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$ba_ds_run
    
    if (g_resetresult == FALSE || is.null(ba_ds_results()$Variable) || is.null(ba_ds_results()$GroupBoxplot))
        return()
    
    div(
        style = "position: relative; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("ba_ds_groupboxplot"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left:0.5em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "ba_ds_downloadgroupboxplot", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})