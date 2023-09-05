# Frailty Run ####

observeEvent(input$sa_fm_g_run, {
    g_resetresult <<- TRUE
})

sa_fm_results <- eventReactive(input$sa_fm_g_run, {
    # ¦§ Warning & Notify ####
    if (input$sa_fm_s_survivaltime == "" || input$sa_fm_s_indicator == "") {
        showNotification("Please select survival time and indicator", type="warning")
        return()
    }
    
    if (input$sa_fm_s_check_initial == TRUE && input$sa_fm_s_initialtime == "") {
        showNotification("Please choose initial time", type="warning")
        return()
    }
        
    withProgress(message = 'Frailty Model', style = "notification", value = 0.1, {
    
    # ¦§ Variable Declaration ####
    # in FrailtyHL, they don't support non-random effect model. 
    # So in revised frailtyHL(aa), we made constant value named "idid". 
    # if we set "idid" as random effect, frailtyHL will automatically estimate idid as 0 during maximizing likelihood. 
    tempData = data2
    dataLength = nrow(tempData)
    tempData$idid <- 1:dataLength
        
    modelFormula <- input$sa_fm_s_model
    tempModelFormula <- modelFormula
    
    modelRandFamily <- "Normal"
    modelNumberRandom <- length(input$sa_fm_s_rand)
    modelMord <- input$sa_fm_a_mord
    modelDord <- input$sa_fm_a_dord
    modelVarfixed <- TRUE
    modelVarinit <- 0
    
    if(modelNumberRandom > 0) {
        modelRandFamily <- input$sa_fm_s_randfamily
        modelVarfixed <- FALSE
        modelVarinit <- 0.163
    } else {
        tempModelFormula <- paste0(modelFormula, " + (1|idid)")
        # in FrailtyHL, they don't support non-random effect model. 
        # So in revised frailtyHL(aa), they made constant value named "idid". 
        # if we set "idid" as random effect, frailtyHL will automatically estimate idid as 0 during maximizing likelihood. 
    }
    
    # ¦§ Fitted Model ####
    
    incProgress(0.4, detail = paste("Loading..."))
    if(input$sa_fm_s_check_group) 
        fittedModel <<- frailtyHL_grouped(
            formula = as.formula(tempModelFormula), 
            data = tempData, 
            RandDist = modelRandFamily, 
            mord = modelMord, 
            dord = modelDord, 
            varfixed = modelVarfixed, 
            varinit = modelVarinit,
            grouped = TRUE
        )
    else 
        fittedModel <<- aa(
            formula = as.formula(tempModelFormula), 
            data = tempData, 
            RandDist = modelRandFamily, 
            mord = modelMord, 
            dord = modelDord, 
            varfixed = modelVarfixed, 
            varinit = modelVarinit
        )
    fittedModel$formula <- as.formula(modelFormula)
    incProgress(0.3, detail = paste("Loading..."))
    
    # ¦§ Description ####
    
    dataLength = length(data2[[input$sa_fm_s_survivaltime]])
    eventNumber = sum(data2[[input$sa_fm_s_indicator]] == 1)
    modelDesc <- matrix(data = c(fittedModel$Model, dataLength, eventNumber, fittedModel$Method), ncol=4)
    colnames(modelDesc) <- c("Model", "Number of data", "Number of events", "Method")
    
    comparisonLikelihood <- NULL
    # ¦§ Likelihood ####
    if(modelMord == 0 && modelDord == 1) {
        colnames(fittedModel$likelihood) <- c("\\(\\ -2h_0 \\)", "\\(\\ -2h_p \\)", "\\(\\ -2p_{\\beta, \\upsilon}(h_p) \\)")
        comparisonLikelihood <- fittedModel$likelihood
    }
    else if(modelMord == 0 && modelDord == 2) {
        colnames(fittedModel$likelihood) <- c("\\(\\ -2h_0 \\)", "\\(\\ -2h_p \\)", 
                                              "\\(\\ -2p_{\\beta, \\upsilon}(h_p) \\)", "\\(\\ -2s_{\\beta, \\upsilon}(h_p) \\)")
        comparisonLikelihood <- fittedModel$likelihood[,1:3]
    }
    else if(modelMord == 1 && modelDord == 1) {
        colnames(fittedModel$likelihood) <- c("\\(\\ -2h_0 \\)", "\\(\\ -2h_p \\)", 
                                              "\\(\\ -2p_{\\upsilon}(h_p) \\)", "\\(\\ -2p_{\\beta, \\upsilon}(h_p) \\)")
        comparisonLikelihood <- fittedModel$likelihood[,c(1,2,4)]
    }
    else {      # modelMord == 1 && modelDord == 2
        colnames(fittedModel$likelihood) <- c("\\(\\ -2h_0 \\)", "\\(\\ -2h_p \\)", 
                                              "\\(\\ -2p_{\\upsilon}(h_p) \\)", "\\(\\ -2s_{\\upsilon}(h_p) \\)", 
                                              "\\(\\ -2p_{\\beta, \\upsilon}(h_p) \\)", "\\(\\ -2s_{\\beta, \\upsilon}(h_p) \\)")
        comparisonLikelihood <- fittedModel$likelihood[,c(1,2,5)]
    }
    
    # ¦§ Comparison Model ####
    
    fmComparisonmodel <- NULL
    if(input$sa_fm_s_check_comparison) {
        fmComparisonmodel <- matrix(
            data = c(modelFormula, NA, fittedModel$Method, round(comparisonLikelihood, digits = 3), round(fittedModel$aic, digits = 3)), 
            ncol = 9
        )
        if(modelNumberRandom > 0)
            fmComparisonmodel[1,2] <- modelRandFamily
        fmComparisonmodel <- rbind(g_sa_fm_comparisonmodel, fmComparisonmodel)
        colnames(fmComparisonmodel) <- c("Model", "Rand", "Method", 
                                         "\\(\\ -2h_0 \\)", "\\(\\ -2h_p \\)", "\\(\\ -2p_{\\beta, \\upsilon}(h_p) \\)", 
                                         colnames(fittedModel$aic))
        g_sa_fm_comparisonmodel <<- fmComparisonmodel
    }
    
    # ¦§ R Codes ####
    
    Rcodes <- NULL
    if(input$sa_fm_s_check_rcodes) {
        Rcodes <- paste0(
            "frailtyHL(Formula = ", 
            modelFormula, 
            ", data = dataName, RandDist = \"", 
            modelRandFamily, 
            "\", mord = ", 
            modelMord, 
            ", dord = ", 
            modelDord, 
            ")"
        )
        Rcodes <- matrix(Rcodes)
        colnames(Rcodes) <- "Call"
    }

    # ¦¦ Options ####
    
    fittedModel$option <- list(
        Formula = modelFormula, 
        RandFamily = modelRandFamily, 
        Mord = modelMord, 
        Dord = modelDord, 
        Grouped = input$sa_fm_s_check_group,
        nRand = modelNumberRandom,
        ModelDesc = modelDesc,
        ComparisonModel = fmComparisonmodel,
        Rcodes = Rcodes
    )
    
    removeTab(inputId = "sa_fm_resulttabset", target = "Random Effect Inferences")
    if (modelNumberRandom > 0) {
        appendTab(inputId = "sa_fm_resulttabset",
            tabPanel("Random Effect Inferences",
                br(),
                uiOutput("sa_fm_r_showplot1")
            )
        )
    }
        
    setProgress(1, detail = "Finish")
    return(fittedModel)
    })
})

# Frailty Components ####

output$sa_fm_s_model<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$sa_fm_s_survivaltime
    input$sa_fm_s_variable$right
    input$sa_fm_s_indicator
    input$sa_fm_s_initialtime
    
    if(any(is.null(input$sa_fm_s_survivaltime), is.null(input$sa_fm_s_variable), is.null(input$sa_fm_s_indicator)))
        return()
    
    responsePart = NULL
    variablePart = paste0(input$sa_fm_s_variable$right, collapse="+")
    randomPart = NULL
    modelEquation = NULL
    
    if (!is.null(input$sa_fm_s_initialtime) && input$sa_fm_s_initialtime != "" && input$sa_fm_s_check_initial) {
        responsePart = paste0("Surv(", input$sa_fm_s_survivaltime," - ",input$sa_fm_s_initialtime,", ",input$sa_fm_s_indicator," == 1)")
    } else {
        responsePart = paste0("Surv(", input$sa_fm_s_survivaltime,", ",input$sa_fm_s_indicator," == 1)")
    }
    
    if(length(input$sa_fm_s_variable$right) == 0) {
        variablePart = "1"
    }
    
    if(!is.null(input$sa_fm_s_rand) && length(input$sa_fm_s_rand) > 0) {
        randomPart = paste0(input$sa_fm_s_rand, collapse = ") + (1|")
        randomPart = paste0(" + (1|", randomPart, ")")
    }
    
    modelEquation = paste0(responsePart, " ~ ", variablePart, randomPart)
    
    textAreaInput("sa_fm_s_model","Model",value = modelEquation)
})

output$sa_fm_s_survivaltime<-renderUI({
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
            "sa_fm_s_survivaltime",
            "Survival Time", 
            choices = as.list(c("",nameValue)),
            selected = NULL,
            multiple = FALSE
        ),
        bsTooltip(
            "sa_fm_s_survivaltime", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$sa_fm_s_check_initial<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("sa_fm_s_check_initial", "Initial Time", value = FALSE)
})

output$sa_fm_s_initialtime<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$sa_fm_s_check_initial) || !input$sa_fm_s_check_initial)
        return()
      
    nameValue = names(select_if(data2, is.numeric))
    
    div(
        selectInput(
            "sa_fm_s_initialtime",
            "Initial Time",
            choices = as.list(c("",nameValue)),
            multiple = FALSE
        ),
        bsTooltip(
            "sa_fm_s_initialtime", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$sa_fm_s_variable<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$sa_fm_s_interactionappend
    
    nameValue <- c(names(data2), g_sa_fm_interaction)
    
    chooserInput("sa_fm_s_variable", "Variable", "Selected", nameValue, c(), size = 15, multiple = TRUE, width=150)
})


output$sa_fm_s_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue <- names(data2)
    
    selectInput(
        "sa_fm_s_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$sa_fm_s_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    actionButton("sa_fm_s_interactionappend", "Append")
})

observeEvent(input$sa_fm_s_interactionappend, {
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    g_sa_fm_interaction <<- c(g_sa_fm_interaction, paste(input$sa_fm_s_interaction, collapse=":"))
})

output$sa_fm_s_rand <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$sa_fm_s_randinteractionappend
    
    nameValue <- c(names(data2), g_sa_fm_randinteraction)
    
    selectInput(
        "sa_fm_s_rand",
        "Random Effects",
        choices = as.list(nameValue), 
        multiple = TRUE
    )
})

output$sa_fm_s_randinteraction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue <- names(data2)
    
    selectInput(
        "sa_fm_s_randinteraction",
        "Interaction in Random Effects",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$sa_fm_s_randinteractionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    actionButton("sa_fm_s_randinteractionappend", "Append")
})

observeEvent(input$sa_fm_s_randinteractionappend, {
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    g_sa_fm_randinteraction <<- c(g_sa_fm_randinteraction, paste(input$sa_fm_s_randinteraction, collapse=":"))
})

# in frailtyHL(aa also), they don't support specific randfamily. all gaussian or all gamma. 
# RandDist : Distribution for random effect ("Normal" or "Gamma").
output$sa_fm_s_randfamily <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$sa_fm_s_model) || is.null(input$sa_fm_s_rand) || length(input$sa_fm_s_rand) == 0)
        return()
    
    selectInput(
        "sa_fm_s_randfamily",
        "Distribution for Random effects",
        choices = c(
            "gaussian"="Normal",
            "gamma"="Gamma"
        )
    )
})

output$sa_fm_s_indicator<-renderUI({
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
            "sa_fm_s_indicator",
            "Censoring Indicator",
            choices = as.list(c("",nameValue)),
            multiple = FALSE
        ),
        bsTooltip(
            "sa_fm_s_indicator", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$sa_fm_s_check_group <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("sa_fm_s_check_group", "Grouped Duration", value = FALSE)
})

output$sa_fm_s_check_robustse <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("sa_fm_s_check_robustse", "Robust Standard Error", value = FALSE)
})

output$sa_fm_s_check_confint <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("sa_fm_s_check_confint", "Confidence Intervals for Coefficients", value = FALSE)
})

output$sa_fm_s_check_exp <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("sa_fm_s_check_exp", "Exponential Scale", value = FALSE)
})

output$sa_fm_s_check_comparison <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("sa_fm_s_check_comparison", "model comparison", value = FALSE)
})

output$sa_fm_s_check_rcodes <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("sa_fm_s_check_rcodes", "R code", value = FALSE)
})

output$sa_fm_a_mord<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    selectInput("sa_fm_a_mord","Order for Mean",choices=c("0"="0","1"="1"),multiple=FALSE)
})

output$sa_fm_a_dord<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    selectInput("sa_fm_a_dord","Order for Dispersion",choices=c("1"="1","2"="2"),multiple=FALSE)
})

# Frailty Results ####

output$sa_fm_r_title <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$sa_fm_g_run
    
    titleOutput = NULL
    if (g_resetresult == FALSE || is.null(sa_fm_results()))
        titleOutput = ""
    else
        titleOutput = sa_fm_results()$option$Formula
    
    return(paste0("Model : ", titleOutput))
})

output$sa_fm_r_mainsummary <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_fm_g_run
    
    if (g_resetresult == FALSE || is.null(sa_fm_results()$option$ModelDesc))
        return()
    
    sa_fm_results()$option$ModelDesc
}, rownames = TRUE, bordered = TRUE, caption = "Model Description", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"))

output$sa_fm_fix_coef <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_fm_g_run
    
    if (g_resetresult == FALSE || is.null(sa_fm_results()$FixCoef))
        return()
    
    sa_fm_results()$FixCoef
}, rownames = TRUE, bordered = TRUE, caption = "Estimates from the mean model", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$sa_fm_rand_coef <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_fm_g_run
    
    if (g_resetresult == FALSE || is.null(sa_fm_results()$RandCoef))
        return()
    
    # In non-random effect model, randCoef has 'idid' which is dummay value. 
    # So we will return NULL in that case
    if(length(strsplit(sa_fm_results()$option$Formula, "|", fixed = TRUE)[[1]]) == 1)
        return()
    
    
    sa_fm_results()$RandCoef
}, rownames = TRUE, bordered = TRUE, caption = "Estimates from the dispersion model", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$sa_fm_r_likelihood <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$sa_fm_g_run
    
    tagList(
        withMathJax(),
        withMathJax(tableOutput("sa_fm_r_likelihoodtable"))
    )
})

output$sa_fm_r_likelihoodtable <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_fm_g_run
    
    if (g_resetresult == FALSE || is.null(sa_fm_results()$likelihood))
        return()
    
    sa_fm_results()$likelihood
}, colnames = TRUE, bordered = TRUE, caption = "Likelihood", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$sa_fm_r_aic <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_fm_g_run
    
    if (g_resetresult == FALSE || is.null(sa_fm_results()$aic))
        return()
    
    sa_fm_results()$aic
}, colnames = TRUE, bordered = TRUE, caption = "AIC", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$sa_fm_r_comparisonmodel <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$sa_fm_g_run
    
    tagList(
        withMathJax(),
        withMathJax(tableOutput("sa_fm_r_comparisonmodeltable"))
    )
})

output$sa_fm_r_comparisonmodeltable <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_fm_g_run
    
    if (g_resetresult == FALSE || is.null(sa_fm_results()$option$ComparisonModel))
        return()
    
    sa_fm_results()$option$ComparisonModel
}, colnames = TRUE, bordered = TRUE, caption = "Comparison Model", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"))

output$sa_fm_r_rcodes<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_fm_g_run
    
    if (g_resetresult == FALSE || is.null(sa_fm_results()$option$Rcodes))
        return()
  
    sa_fm_results()$option$Rcodes
}, rownames = FALSE, bordered = TRUE, caption = "R Codes", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

# Frailty Plots ####

sa_fm_r_reactiveplot1 <- reactive({
    if(g_resetresult == FALSE || is.null(sa_fm_results()))
        return()
    
    if(is.null(input$sa_fm_r_selectplot1))
        plotType = 1
    else
        plotType = as.integer(input$sa_fm_r_selectplot1)
    
    res <- sa_fm_results()
    p <- res$p
    q <- res$q
    qLength <- length(q)
    qcum <- cumsum(q)
    qcum <- c(0, qcum)
    
    var <- diag(res$Hinv)[(p+qcum[plotType]+1):(p+qcum[plotType+1])]
    SE <- sqrt(var)
    v_h <- res$v_h[(qcum[plotType]+1):(qcum[plotType+1])]
    lb <- v_h - 1.96*SE
    ub <- v_h + 1.96*SE
    v_h <- sort(v_h)
    
    # Plot 1
    plot1 <- ggplot(data = data.frame(v_h), aes(sample=v_h)) + 
      stat_qq() + 
      stat_qq_line(linetype = 2) + 
      labs(x = "Theoretical Quantiles", y = "Sample Quantiles") + 
      ggtitle("Normal Probability Plot") + 
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5))
              
    # Plot 2
    plot2 <- ggplot(data.frame(v_h), aes(v_h)) + geom_histogram() +
      labs(x = "Estimated Frailty Effects", y = "Frequency") +
      ggtitle(paste0("Histogram of Frailty Effects", plotType)) +
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5))
    
    lb <- sort(lb)
    ub <- sort(ub)
    randCount <- 1:q[plotType]
        
    # Plot 3
    ciTable <- cbind(lb, ub, randCount)
    plot3 <- ggplot(data.frame(cbind(v_h, randCount)), aes(randCount, v_h)) + 
        geom_point() + 
        geom_line() +
        geom_segment(data = ciTable, aes(x = randCount, y = lb, xend = randCount, yend = ub)) +
        labs(x = "Frailty Number", y = "Estimated Frailty Effects") + 
        ggtitle(paste0("95% CI for Frailty Effect", plotType)) +
        theme_bw() + 
        theme(plot.title = element_text(hjust = 0.5))
    
    return(ggarrange(ggarrange(plot1, plot2, ncol = 2), plot3, nrow = 2) )
})

output$sa_fm_r_selectplot1<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$sa_fm_g_run
    
    if (g_resetresult == FALSE || is.null(sa_fm_results()))
        return()
    
    numberRand <- sa_fm_results()$option$nRand
    
    if (numberRand == 0)
        return()
    
    choiceName = NULL
    choiceValue = NULL
    
    for (i in 1:numberRand) {
        choiceName = c(choiceName, paste0("Random Effect", i))
        choiceValue = c(choiceValue, i)
    }

    radioGroupButtons(
        inputId = "sa_fm_r_selectplot1",
        label = "Select Plot", 
        choiceNames = choiceName,
        choiceValues = choiceValue, 
        selected = "r1",
        direction = "vertical"
    )
})

output$sa_fm_r_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$sa_fm_g_run
    
    local({
        output$sa_fm_r_downloadplot1 <- downloadPlot({
            sa_fm_r_reactiveplot1()
        })      
    })
    
    sa_fm_r_reactiveplot1()
})

output$sa_fm_r_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$sa_fm_g_run

    if(g_resetresult == FALSE || is.null(sa_fm_results()))
        return()
    
    div(
        style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("sa_fm_r_plot1", height="600px"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left: 0.5em; bottom: 0.5em;",
            dropdown(
                uiOutput("sa_fm_r_selectplot1"),
                size = "xs",
                icon = icon("chart-bar", class = "opt"), 
                up = TRUE
            )
        ),
        div(
            style = "position: absolute; left:4em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "sa_fm_r_downloadplot1", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})


        