# LM Run ####

observeEvent(input$rg_lm_g_run, {
    g_resetresult <<- TRUE
})

rg_lm_results <- eventReactive(input$rg_lm_g_run, {
    # ¦§ Warning & Notify ####
    if (input$rg_lm_m_resp == "") {
      showNotification("Please choose response", type="warning")
      return()
    }
    
    if (input$rg_lm_m_check_vif == TRUE && length(input$rg_lm_m_variable$right) < 2) {
        showNotification("Model contains fewer than 2 terms.", type="warning")
        return()
    }
    
    # Start Progressbar
    withProgress(message = 'LM', style = "notification", value = 0, {
      
    # ¦§ LM ####
      
    lmResults = NULL
    PredictionlmResults = NULL
    
    if(input$rg_lm_m_check_offset && input$rg_lm_m_offset == ""){
        showNotification("Choose offset variable.", type="warning")
        return()
    }
    else if(input$rg_lm_m_check_offset){
        lmResults <-lm(as.formula(input$rg_lm_m_model), data = data2, offset = data2[[input$rg_lm_m_offset]])
 
        offsetFormula = paste(input$rg_lm_m_model, "+offset(", input$rg_lm_m_offset, ")")
        lmResultsExcludeNA <- lm(as.formula(offsetFormula), data = data2, na.action = na.exclude)
    } else {
        lmResults <-lm(as.formula(input$rg_lm_m_model),data = data2)
        lmResultsExcludeNA <- lm(as.formula(input$rg_lm_m_model), data = data2, na.action = na.exclude)
    }
    
    # ¦§ Summary ####
    
    summaryResults <- summary(lmResults)
    
    summaryFstats <- summaryResults$fstatistic["value"]
    summaryNumdf <- summaryResults$fstatistic["numdf"]
    summaryDendf <- summaryResults$fstatistic["dendf"]
    
    if (is.null(summaryFstats) || is.null(summaryNumdf) || is.null(summaryDendf))
        pfResults <- NA
    else pfResults <- pf(summaryFstats, summaryNumdf, summaryDendf, lower.tail = F)
    
    attributes(pfResults) <- NULL
    
    M.summary <- c(
      "Num of obs" = nobs(lmResults),
      "Fstats" = summaryFstats,
      "Prob > F" = pfResults,
      "R-squared" = summaryResults$r.squared,
      "Adj R-squared" = summaryResults$adj.r.squared,
      "Residual Std.Error" = summaryResults$sigma,
      "AIC" = AIC(lmResults)
    )
    
    if (is.null(summaryNumdf) || is.null(summaryDendf))
        names(M.summary)[2] <- "F"
    else names(M.summary)[2] <- sprintf("F(%g, %g)", summaryNumdf, summaryDendf)
    
    if(input$rg_lm_m_check_bptest){
        bpResults <- bptest(lmResults)
        bpResults <- c(bpResults[1]$statistic, bpResults[2]$parameter, bpResults[4]$p.value)
        names(bpResults) <- c("BP.stats", "BP.df", "BP.p.value")
        M.summary <- c(M.summary, bpResults)
    }
    
    # ¦§ Coefficients ####
    
    coefResults <- summaryResults$coefficients
    if(input$rg_lm_m_check_robustse) {
        coefResults <- unclass(coeftest(lmResults, vcov = vcovHC(lmResults, "HC1")))
    }
    
    if(input$rg_lm_m_check_confint) {
        coefResults <- cbind(coefResults, confint(lmResults))
    }
    
    if(input$rg_lm_m_check_vif) {
        coefResults <- cbind(coefResults, VIF = c(NA, vif(lmResults)))
    }
    
    # ¦§ Margins ####
    
    lmMargins <- NULL
    if (all(input$rg_lm_m_check_margins, !all(is.null(input$rg_lm_m_fmargins), is.null(input$rg_lm_m_cmargins)))) {
        lmMarginlist <- list()
        if(input$rg_lm_m_check_offset == FALSE) {
            if(!is.null(input$rg_lm_m_fmargins)){
                lmFmargins <- input$rg_lm_m_fmargins
                lmFmarginlist <- lapply(lmFmargins, function(x) levels(as.factor(data2[[x]])))
                names(lmFmarginlist) <- lmFmargins
                lmMarginlist <- append(lmMarginlist, lmFmarginlist)
            }
            if(!is.null(input$rg_lm_m_cmargins)){
                lmCmargins <- input$rg_lm_m_cmargins
                lmCmarginlist <- lapply(1:length(lmCmargins), function(x){
                    eval(parse(text = paste("input$rg_lm_m_marginvalue", x, sep = "")))
                })
                names(lmCmarginlist) <- lmCmargins
                lmMarginlist <- append(lmMarginlist, lmCmarginlist)
            }
        lmMargins <- append(lmMargins, summary(margins(lmResults, at = lmMarginlist)))
        }
    }
    
    # ¦§ Comparision Model ####
    
    lmComparisonmodel <- NULL
    lmAnova <- anova(lmResults)
    if (input$rg_lm_m_check_comparison) {
        if (input$rg_lm_m_check_offset == TRUE) {
            lmComparisonmodel <- matrix(c(
                input$rg_lm_m_model,
                input$rg_lm_m_offset,
                lmResults$df.residual,
                round(lmAnova[["Sum Sq"]][[nrow(lmAnova)]], digits = 5),
                round(summaryResults$r.squared, digits = 5),
                round(summaryResults$adj.r.squared, digits = 5),
                round(AIC(lmResults), digits = 5),
                round(BIC(lmResults), digits = 5)),
                ncol = 8
            )
        } else {
            lmComparisonmodel <- matrix(c(
                input$rg_lm_m_model,
                NA,
                lmResults$df.residual,
                round(lmAnova[["Sum Sq"]][[nrow(lmAnova)]], digits = 5),
                round(summaryResults$r.squared, digits = 5),
                round(summaryResults$adj.r.squared, digits = 5),
                round(AIC(lmResults), digits = 5),
                round(BIC(lmResults), digits = 5)),
                ncol = 8
            )
        }
        lmComparisonmodel <- rbind(g_rg_lm_comparisonmodel, lmComparisonmodel)
        colnames(lmComparisonmodel) <- c("Model", "Offset", "Res.df", "Sum.sq", "R-squared", "Adj R-squared", "AIC", "BIC")
        g_rg_lm_comparisonmodel <<- lmComparisonmodel
    }
        
    # ¦§ R Codes ####
    
    Rcodes <- NULL
    if (input$rg_lm_m_check_rcodes) {
        option <- list(form = input$rg_lm_m_model)
        rform <- as.character(as.formula(option$form))
        Rcodes <- matrix(paste0("lm(formula = ", paste(rform[2], rform[1], rform[3]), ", data = data)"))
        colnames(Rcodes) <- "Call"
    }
    
    # Increase Progressbar
    incProgress(0.5, detail = paste("Loading..."))
    
    # ¦¦ Plots ####
    lmPlots <- ggplotglm(lmResults)
    
    # Increase Progressbar
    incProgress(0.5, detail = paste("Loading..."))
    
    return(list(
        Title = input$rg_lm_m_model,
        Model = lmResults,
        ANOVA = anova(lmResults),
        M.summary = M.summary,
        Coefficients = coefResults,
        Margins = lmMargins,
        ModelExcludeNA = lmResultsExcludeNA,
        Comparisonmodel = lmComparisonmodel,
        Rcodes = Rcodes,
        Plot1 = lmPlots
    ))
    
    }) # End Progressbar
})

# LM Components ####

output$rg_lm_m_model<-renderUI({
    input$rg_lm_m_resp
    input$rg_lm_m_variable$right
    
    responsePart = input$rg_lm_m_resp
    variablePart = paste0(input$rg_lm_m_variable$right, collapse="+")
    modelEquation = NULL
    
    if(!is.null(input$rg_lm_m_check_nointercept) && !input$rg_lm_m_check_nointercept) {
        if(variablePart == "")
            modelEquation = paste(responsePart, "~", "1")
        else
            modelEquation = paste(responsePart, "~", variablePart)
    }
    
    if(!is.null(input$rg_lm_m_check_nointercept) && input$rg_lm_m_check_nointercept) {
        if(variablePart == "")
            modelEquation = paste(responsePart, "~", "-1")
        else
            modelEquation = paste(responsePart, "~", "-1 +", variablePart)
    }
    
    textAreaInput("rg_lm_m_model", "Model", value = modelEquation, height = "80px")
})

output$rg_lm_m_resp <- renderUI({
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
            "rg_lm_m_resp",
            "Response Variable",
            choices = as.list(c("", nameValue)), 
            selected = NULL, 
            multiple = FALSE
        ),
        bsTooltip(
            "rg_lm_m_resp", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$rg_lm_m_variable <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$rg_lm_m_interactionappend
    
    nameValue <- c(names(data2), g_rg_lm_interaction)
    
    chooserInput("rg_lm_m_variable", "Variable", "Selected", nameValue, c(), size = 15, multiple = TRUE, width = 150)
})

output$rg_lm_m_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue <- names(data2)
    
    selectInput(
        "rg_lm_g_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$rg_lm_m_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    actionButton("rg_lm_m_interactionappend", "Append")
})

observeEvent(input$rg_lm_m_interactionappend, {
    g_rg_lm_interaction <<- c(g_rg_lm_interaction, paste(input$rg_lm_g_interaction, collapse=":"))
})

output$rg_lm_m_check_nointercept<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_lm_m_check_nointercept", "No Intercept Model", value = FALSE)
})

output$rg_lm_m_check_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_lm_m_check_offset", "Offset Variable", value = FALSE)
})
    
output$rg_lm_m_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    nameValue = names(select_if(data2, is.numeric))
    
    if (is.null(input$rg_lm_m_check_offset) || !input$rg_lm_m_check_offset)
        return()
    
    div(
        selectInput(
            "rg_lm_m_offset",
            "Offset Variable",
            choices = as.list(c("", nameValue)),
            selected = NULL,
            multiple = FALSE
        ),
        bsTooltip(
            "rg_lm_m_offset", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$rg_lm_m_check_margins<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_lm_m_check_margins", "Margins", value = FALSE)
})

output$rg_lm_m_margins1<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if (is.null(input$rg_lm_m_check_margins) || !input$rg_lm_m_check_margins || is.null(input$rg_lm_m_variable$right))
        return()
    
    Factornames = names(select_if(data2, is.factor))
    Contiousnames = names(select_if(data2, is.numeric))
    Factornames = as.list(intersect(Factornames, input$rg_lm_m_variable$right))
    Contiousnames = as.list(intersect(Contiousnames, input$rg_lm_m_variable$right))
    
    tagList(
        selectInput(
            "rg_lm_m_fmargins",
            "Factor Variables",
            choices = Factornames,
            multiple = TRUE
        ),
        bsTooltip(
            "rg_lm_m_fmargins", 
            "Factor Only",
            "right", 
            options = list(container = "body")
        ),
        selectInput(
            "rg_lm_m_cmargins",
            "Contiuous Variables", 
            choices = Contiousnames,
            multiple = TRUE
        ),
        bsTooltip(
            "rg_lm_m_cmargins", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$rg_lm_m_margins2<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if (is.null(input$rg_lm_m_check_margins) || !input$rg_lm_m_check_margins || is.null(input$rg_lm_m_cmargins))
        return()
    
    lapply(1:length(input$rg_lm_m_cmargins), function(i) {
        numericInput(
            inputId = paste0("rg_lm_m_marginvalue", i), 
            label = input$rg_lm_m_cmargins[i],
            value = mean(data2[[input$rg_lm_m_cmargins[i]]]),
            max = max(data2[[input$rg_lm_m_cmargins[i]]]), 
            min = min(data2[[input$rg_lm_m_cmargins[i]]])
        )
    })
})

output$rg_lm_m_check_vif<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_lm_m_check_vif", "VIF", value = FALSE)
})

output$rg_lm_m_check_robustse<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_lm_m_check_robustse", "Robust Standard Errors", value = FALSE)
})

output$rg_lm_m_check_confint<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_lm_m_check_confint", "Confidence Intervals for Coefficients", value = FALSE)
})

output$rg_lm_m_check_bptest<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_lm_m_check_bptest", "Bruesch-Pagan Test", value = FALSE)
})

output$rg_lm_m_check_comparison<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_lm_m_check_comparison", "model comparison", value = FALSE)
})

output$rg_lm_m_check_rcodes<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_lm_m_check_rcodes", "R Codes", value = FALSE)
})

# LM Results ####

output$rg_lm_r_model <- renderText({
    input$file
    input$resetData
    input$di_option_run
    input$rg_lm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_lm_results()$Title))
        titleOutput = ""
    else
        titleOutput = rg_lm_results()$Title
    
    paste("Model : ", titleOutput)
})

output$rg_lm_r_anova<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_lm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_lm_results()$ANOVA))
        return()
  
    rg_lm_results()$ANOVA
}, rownames = TRUE, bordered = TRUE, caption = "ANOVA", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_lm_r_modelsummary<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_lm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_lm_results()$M.summary))
        return()
  
    rg_lm_results()$M.summary 
}, colnames = FALSE, rownames = TRUE, bordered = TRUE, caption = "Model Summary", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_lm_r_coefficients<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_lm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_lm_results()$Coefficients))
        return()
  
    rg_lm_results()$Coefficients 
}, rownames = TRUE, bordered = TRUE, caption = "Coefficients", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_lm_r_margins<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_lm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_lm_results()$Margins))
        return()
  
    rg_lm_results()$Margins 
}, rownames = TRUE, bordered = TRUE, caption = "Margins", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_lm_r_comparisonmodel<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_lm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_lm_results()$Comparisonmodel))
        return()
  
    rg_lm_results()$Comparisonmodel
}, rownames = TRUE, bordered = TRUE, caption = "Comparison Model", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_lm_r_rcodes<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_lm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_lm_results()$Rcodes))
        return()
  
    rg_lm_results()$Rcodes
}, rownames = FALSE, bordered = TRUE, caption = "R Codes", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

# LM Plots ####

output$rg_lm_r_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$rg_lm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_lm_results()$Plot1))
        return()
      
    local({
        output$rg_lm_r_downloadplot1 <- downloadPlot(
            rg_lm_results()$Plot1
        )
    })
    rg_lm_results()$Plot1
})

output$rg_lm_r_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$rg_lm_g_run

    if (g_resetresult == FALSE || is.null(rg_lm_results()$Plot1))
        return()
    
    div(
        style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("rg_lm_r_plot1", height="600px"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left:0.5em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "rg_lm_r_downloadplot1", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})

# LM Prediction ####

output$rg_lm_r_prediction <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_lm_g_run
  
    if (g_resetresult == FALSE || is.null(rg_lm_results()$ModelExcludeNA)) {
        return()
    }
    regression <<- rg_lm_results()$ModelExcludeNA
    predictionOutput <- data2
    
    mu <- predict(regression, predictionOutput, interval="confidence")
    StudentResidual <- rstandard(regression)
    predictionOutput$pred <- mu[,1]
    predictionOutput$predLL <- mu[,2]
    predictionOutput$predUL <- mu[,3]
    predictionOutput$StudentResidual <- StudentResidual

    return(predictionOutput)
}, rownames = TRUE, bordered = TRUE, spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)