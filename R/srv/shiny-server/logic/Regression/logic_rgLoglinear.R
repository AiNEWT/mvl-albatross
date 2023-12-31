# Log-linear Run ####

observeEvent(input$rg_loglinear_g_run, {
    g_resetresult <<- TRUE
})

rg_loglinear_results <- eventReactive(input$rg_loglinear_g_run, {
    # �� Warning & Notify ####
    if (input$rg_loglinear_m_resp == "") {
      showNotification("Please choose response", type="warning")
      return()
    }
    
    if (input$rg_loglinear_m_dist == "binomial" && input$rg_loglinear_m_check_binomd == TRUE && input$rg_loglinear_m_binomd == "") {
        showNotification("Please choose binomial denominator.", type="warning")
        return()
    }
  
    if (input$rg_loglinear_m_check_vif == TRUE && length(input$rg_loglinear_m_variable$right) < 2) {
        showNotification("Model contains fewer than 2 terms.", type="warning")
        return()
    }
    
    # Start Progressbar
    withProgress(message = 'Log-linear', style = "notification", value = 0, {
      
    # �� GLM ####
        
    loglinearResults = NULL
    
    if(input$rg_loglinear_m_check_offset){
        loglinearResults <- glm(
            as.formula(input$rg_loglinear_m_model),
            family = eval(parse(text = paste(input$rg_loglinear_m_dist,"(link = ",input$rg_loglinear_m_link,")", sep = ""))),
            data = data2,
            offset = data2[[input$rg_loglinear_m_offset]]
        )
    } else {
        loglinearResults <- glm(
            as.formula(input$rg_loglinear_m_model),
            family = eval(parse(text = paste(input$rg_loglinear_m_dist,"(link = ",input$rg_loglinear_m_link,")", sep = ""))),
            data = data2
        )
    }
    
    # �� Summary ####
    
    summaryResults <- summary(loglinearResults)

    M.summary1 <- c(
        "Family" = input$rg_loglinear_m_dist,
        "Link" = input$rg_loglinear_m_link,
        "Optimization" = "IWLS",
        "Num of iteration" = loglinearResults$iter
    )
    
    M.summary2 <- c(
        "Num of obs" = nobs(loglinearResults),
        "Res. deviance" = loglinearResults$deviance,
        "Null deviance" = loglinearResults$null.deviance,
        "Log likelihood" = logLik(loglinearResults)[[1]],
        "AIC" = AIC(loglinearResults),
        "BIC" = BIC(loglinearResults)
    )
    names(M.summary2)[2] <- sprintf("Res. deviance(df=%g)", loglinearResults$df.residual)
    names(M.summary2)[3] <- sprintf("Null deviance(df=%g)", loglinearResults$df.null)
    
    # �� Coefficients ####
    
    coefResults = NULL
    if(input$rg_loglinear_m_check_exp){
        lengthCoef <- length(coef(loglinearResults))
        CoefMatrix1 <- matrix(exp(coef(loglinearResults)), lengthCoef, 1)
        CoefMatrix2 <- matrix(summaryResults$coefficients[,c(-1,-2)], lengthCoef, 2)
        
        if (loglinearResults$df.residual==0 && (input$rg_loglinear_m_dist=="gaussian" | input$rg_loglinear_m_dist=="gamma")) {
            CoefMatrix3 <- matrix(NA, lengthCoef, 2)
        } else {
            CoefMatrix3 <- matrix(exp(confint(loglinearResults)), lengthCoef, 2)
        }
        coefResults<-cbind(CoefMatrix1, CoefMatrix2)
        colnames(coefResults) <- cbind("Estimate(Exp)","t-value","Pr(>|t|)")
        if(input$rg_loglinear_m_check_confint) {
            coefResults <- cbind(coefResults, CoefMatrix3)
            colnames(coefResults) <- cbind("Estimate(Exp)","t-value","Pr(>|t|)","2.5%","97.5%")
        }
        rownames(coefResults) <- names(coef(loglinearResults))
    } else {
        coefResults <- summaryResults$coefficients
        lengthCoef <- nrow(coefResults)
        CoefMatrix1 <- matrix(coefResults, lengthCoef, 4)
        if(input$rg_loglinear_m_check_confint) {
            if (loglinearResults$df.residual==0 && (input$rg_loglinear_m_dist=="gaussian" || input$rg_loglinear_m_dist=="gamma")) {
                CoefMatrix3 <- matrix(NA, lengthCoef, 2)
            } else CoefMatrix3<-matrix(confint(loglinearResults), lengthCoef, 2)
            if(input$rg_loglinear_m_check_robustse) {
                coefResults <- unclass(coeftest(loglinearResults, vcov=vcovHC(loglinearResults,"HC1")))
            }
            coefResults <- cbind(CoefMatrix1 ,CoefMatrix3)
            colnames(coefResults) <- cbind("Estimate","Std.Error","t-value","Pr(>|t|)","2.5%","97.5%")
            rownames(coefResults) <- rownames(summaryResults$coefficients)
        }
    }
    
    if(input$rg_loglinear_m_check_vif) {
        coefResults <- cbind(coefResults, VIF = c(NA, vif(loglinearResults)))
    }
    
    # �� Margins ####
    
    loglinearMargins <- NULL
    if (all(input$rg_loglinear_m_check_margins, !all(is.null(input$rg_loglinear_m_fmargins), is.null(input$rg_loglinear_m_cmargins)))) {
        loglinearMarginlist <- list()
        if(input$rg_loglinear_m_check_offset == FALSE) {
            if(!is.null(input$rg_loglinear_m_fmargins)){
                loglinearFmargins <- input$rg_loglinear_m_fmargins
                loglinearFmarginlist <- lapply(loglinearFmargins, function(x) levels(as.factor(data2[[x]])))
                names(loglinearFmarginlist) <- loglinearFmargins
                loglinearMarginlist <- append(loglinearMarginlist, loglinearFmarginlist)
            }
            if(!is.null(input$rg_loglinear_m_cmargins)){
                loglinearCmargins <- input$rg_loglinear_m_cmargins
                loglinearCmarginlist <- lapply(1:length(loglinearCmargins), function(x){
                    eval(parse(text = paste("input$rg_loglinear_m_marginvalue", x, sep = "")))
                })
                names(loglinearCmarginlist) <- loglinearCmargins
                loglinearMarginlist <- append(loglinearMarginlist, loglinearCmarginlist)
            }
        loglinearMargins <- append(loglinearMargins, summary(margins(loglinearResults, at = loglinearMarginlist)))
        }
    }

    # �� Comparision Model ####
    
    loglinearComparisonmodel <- NULL
    if (input$rg_loglinear_m_check_comparison) {
        if (input$rg_loglinear_m_check_offset == TRUE) {
            loglinearComparisonmodel <- matrix(c(
                input$rg_loglinear_m_model, 
                input$rg_loglinear_m_dist, 
                input$rg_loglinear_m_link, 
                input$rg_loglinear_m_offset,
                loglinearResults$df.residual, 
                round(loglinearResults$deviance, digits=5), 
                round(AIC(loglinearResults), digits=5), 
                round(BIC(loglinearResults), digits=5)),
                ncol = 8
            )
        } else {
            loglinearComparisonmodel <- matrix(c(
                input$rg_loglinear_m_model, 
                input$rg_loglinear_m_dist, 
                input$rg_loglinear_m_link, 
                NA,
                loglinearResults$df.residual, 
                round(loglinearResults$deviance, digits=5),
                round(AIC(loglinearResults), digits=5), 
                round(BIC(loglinearResults), digits=5)),
                ncol = 8
            )
        }
        loglinearComparisonmodel <- rbind(g_rg_glm_comparisonmodel, loglinearComparisonmodel)
        colnames(loglinearComparisonmodel) <- c("Model", "Distribution", "Link", "Offset", "Res.df", "Res.dev", "AIC", "BIC")
        g_rg_glm_comparisonmodel <<- loglinearComparisonmodel
    }

    # �� R Codes ####
    
    Rcodes <- NULL
    option <- list(
        form = input$rg_loglinear_m_model,
        dist = input$rg_loglinear_m_dist,
        link = input$rg_loglinear_m_link,
        checkbinomd = input$rg_loglinear_m_check_binomd,
        binomd = input$rg_loglinear_m_binomd
    )
    if (input$rg_loglinear_m_check_rcodes) {
        rform <- as.character(as.formula(option$form))
        Rcodes <- paste0("glm(formula = ", paste(rform[2], rform[1], rform[3]), ", ")
        family <- paste0("family = ", option$dist, "(", option$link, ")")
        Rcodes <- paste0(Rcodes, family,", data = data)")
        Rcodes <- matrix(Rcodes)
        colnames(Rcodes) <- "Call"
    }
        
    # Increase Progressbar
    incProgress(0.5, detail = paste("Loading..."))
    
    # �� Plots ####
    loglinearPlots <- ggplotglm(loglinearResults)

    # Increase Progressbar
    incProgress(0.5, detail = paste("Loading..."))
    
    return(list(
        Title = input$rg_loglinear_m_model,
        Model = loglinearResults,
        M.summary1 = M.summary1,
        M.summary2 = M.summary2,
        Coefficients = coefResults,
        Margins = loglinearMargins,
        Rcodes = Rcodes,
        Comparisonmodel = loglinearComparisonmodel,
        Option = option,
        Plot1 = loglinearPlots
    ))
    
    }) # End Progressbar
})

# Log-linear Components ####

output$rg_loglinear_m_model<-renderUI({
    input$rg_loglinear_m_resp
    input$rg_loglinear_m_variable$right
    
    if (all(!is.null(input$rg_loglinear_m_check_binomd), input$rg_loglinear_m_check_binomd)) {
        modelEquation1 = paste("cbind(", input$rg_loglinear_m_resp, ",", input$rg_loglinear_m_binomd, "-", input$rg_loglinear_m_resp, ")")
    } else {
        modelEquation1 = input$rg_loglinear_m_resp
    }
    
    if (length(input$rg_loglinear_m_variable$right) >= 1) {
        if (!is.null(input$rg_loglinear_m_check_nointercept) && input$rg_loglinear_m_check_nointercept == TRUE)
            modelEquation2 = paste0(c("-1", input$rg_loglinear_m_variable$right), collapse="+")
        else
            modelEquation2 = paste0(input$rg_loglinear_m_variable$right, collapse="+")
    } else {
        if (!is.null(input$rg_loglinear_m_check_nointercept) && input$rg_loglinear_m_check_nointercept == TRUE)
            modelEquation2 = "-1"
        else
            modelEquation2 = "1"
    }
    
    modelEquation = paste(modelEquation1, "~", modelEquation2)
    
    textAreaInput("rg_loglinear_m_model", "Model", value = modelEquation, height = "80px")
})

output$rg_loglinear_m_resp <- renderUI({
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
            "rg_loglinear_m_resp",
            "Response Variable",
            choices = as.list(c("", nameValue)), 
            selected = NULL, 
            multiple = FALSE
        ),
        bsTooltip(
            "rg_loglinear_m_resp", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$rg_loglinear_m_variable <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$rg_loglinear_m_interactionappend
    
    nameValue <- c(names(data2), g_rg_loglinear_interaction)
    
    chooserInput("rg_loglinear_m_variable", "Variable", "Selected", nameValue, c(), size = 15, multiple = TRUE, width = 150)
})

output$rg_loglinear_m_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue <- names(data2)
    
    selectInput(
        "rg_loglinear_g_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$rg_loglinear_m_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    actionButton("rg_loglinear_m_interactionappend", "Append")
})

observeEvent(input$rg_loglinear_m_interactionappend, {
    g_rg_loglinear_interaction <<- c(g_rg_loglinear_interaction, paste(input$rg_loglinear_g_interaction, collapse=":"))
})

output$rg_loglinear_m_check_binomd<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$rg_loglinear_m_dist) || input$rg_loglinear_m_dist != "binomial")
        return()
    
    checkboxInput("rg_loglinear_m_check_binomd", "Binomial Denominator", value = FALSE)
})
    
output$rg_loglinear_m_binomd<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
        
    if (is.null(input$rg_loglinear_m_check_binomd) || !input$rg_loglinear_m_check_binomd)
        return()
    
    nameValue = names(data2)

    selectInput(
        "rg_loglinear_m_binomd",
        "Binomial Denominator",
        choices = as.list(c("", nameValue)),
        multiple = FALSE
    )
})
    
output$rg_loglinear_m_dist<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    dist = c(
        "poisson"="poisson"
    )
    
    selectInput(
        "rg_loglinear_m_dist",
        "Distribution", 
        choices = dist, 
        multiple = FALSE
    )
})

output$rg_loglinear_m_link<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    link = NULL
    selection = NULL
    
    if(is.null(input$rg_loglinear_m_dist))
        return()

    selection="log"
    link = c("log"="log")
    
    selectInput("rg_loglinear_m_link", "Link Function", choices = link, selected = selection, multiple = FALSE)
})

output$rg_loglinear_m_check_nointercept<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_loglinear_m_check_nointercept", "No Intercept Model", value = FALSE)
})

output$rg_loglinear_m_check_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_loglinear_m_check_offset", "Offset Variable", value = FALSE)
})
    
output$rg_loglinear_m_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
        
    if (is.null(input$rg_loglinear_m_check_offset) || !input$rg_loglinear_m_check_offset)
        return()
    
    nameValue = names(select_if(data2, is.numeric))
    
    div(
        selectInput(
            "rg_loglinear_m_offset",
            "Offset Variable",
            choices = as.list(c("", nameValue)),
            selected = NULL,
            multiple = FALSE
        ),
        bsTooltip(
            "rg_loglinear_m_offset", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$rg_loglinear_m_check_margins<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_loglinear_m_check_margins", "Margins", value = FALSE)
})

output$rg_loglinear_m_margins1<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$rg_loglinear_m_check_margins) || !input$rg_loglinear_m_check_margins || is.null(input$rg_loglinear_m_variable$right))
        return()
    
    Factornames = names(select_if(data2, is.factor))
    Contiousnames = names(select_if(data2, is.numeric))
    Factornames = as.list(intersect(Factornames, input$rg_loglinear_m_variable$right))
    Contiousnames = as.list(intersect(Contiousnames, input$rg_loglinear_m_variable$right))
  
    tagList(
        selectInput(
            "rg_loglinear_m_fmargins",
            "Factor Variables", 
            choices = Factornames, 
            multiple = TRUE
        ),
        bsTooltip(
            "rg_loglinear_m_fmargins", 
            "Factor Only",
            "right", 
            options = list(container = "body")
        ),
        selectInput(
            "rg_loglinear_m_cmargins",
            "Contiuous Variables", 
            choices = Contiousnames, 
            multiple = TRUE
        ),
        bsTooltip(
            "rg_loglinear_m_cmargins", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$rg_loglinear_m_margins2<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$rg_loglinear_m_check_margins) || !input$rg_loglinear_m_check_margins || is.null(input$rg_loglinear_m_cmargins))
        return()

    lapply(1:length(input$rg_loglinear_m_cmargins), function(i) {
        numericInput(
            inputId = paste0("rg_loglinear_m_marginvalue", i), 
            label = input$rg_loglinear_m_cmargins[i],
            value = mean(data2[[input$rg_loglinear_m_cmargins[i]]]),
            max = max(data2[[input$rg_loglinear_m_cmargins[i]]]), 
            min = min(data2[[input$rg_loglinear_m_cmargins[i]]])
        )
    })
})

output$rg_loglinear_m_check_vif<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_loglinear_m_check_vif", "VIF", value = FALSE)
})

output$rg_loglinear_m_check_robustse<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_loglinear_m_check_robustse", "Robust Standard Errors", value = FALSE)
})

output$rg_loglinear_m_check_confint<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_loglinear_m_check_confint", "Confidence Intervals for Coefficients", value = FALSE)
})

output$rg_loglinear_m_check_bptest<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_loglinear_m_check_exp", "Exponential Scale", value = FALSE)
})

output$rg_loglinear_m_check_comparison<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_loglinear_m_check_comparison", "model comparison", value = FALSE)
})
    
output$rg_loglinear_m_check_rcodes<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_loglinear_m_check_rcodes", "R Codes", value = FALSE)
})

# Log-linear Results ####

output$rg_loglinear_r_model <- renderText({
    input$file
    input$resetData
    input$di_option_run
    input$rg_loglinear_g_run
    
    if (g_resetresult == FALSE || is.null(rg_loglinear_results()$Title))
        titleOutput = ""
    else
        titleOutput = rg_loglinear_results()$Title
    
    paste("Model : ", titleOutput)
})

output$rg_loglinear_r_modelsummary1<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_loglinear_g_run
    
    if (g_resetresult == FALSE || is.null(rg_loglinear_results()$M.summary1))
        return()
  
    rg_loglinear_results()$M.summary1
}, rownames = TRUE, bordered = TRUE, caption = "Model Summary 1", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_loglinear_r_modelsummary2<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_loglinear_g_run
    
    if (g_resetresult == FALSE || is.null(rg_loglinear_results()$M.summary2))
        return()
  
    rg_loglinear_results()$M.summary2
}, colnames = FALSE, rownames = TRUE, bordered = TRUE, caption = "Model Summary 2", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_loglinear_r_coefficients<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_loglinear_g_run
    
    if (g_resetresult == FALSE || is.null(rg_loglinear_results()$Coefficients))
        return()
  
    rg_loglinear_results()$Coefficients 
}, rownames = TRUE, bordered = TRUE, caption = "Coefficients", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_loglinear_r_margins<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_loglinear_g_run
    
    if (g_resetresult == FALSE || is.null(rg_loglinear_results()$Margins))
        return()
  
    rg_loglinear_results()$Margins 
}, rownames = TRUE, bordered = TRUE, caption = "Margins", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_loglinear_r_comparisonmodel<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_loglinear_g_run
    
    if (g_resetresult == FALSE || is.null(rg_loglinear_results()$Comparisonmodel))
        return()
  
    rg_loglinear_results()$Comparisonmodel
}, rownames = TRUE, bordered = TRUE, caption = "Comparison Model", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_loglinear_r_rcodes<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_loglinear_g_run
    
    if (g_resetresult == FALSE || is.null(rg_loglinear_results()$Rcodes))
        return()
  
    rg_loglinear_results()$Rcodes
}, rownames = FALSE, bordered = TRUE, caption = "R Codes", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

# Log-linear Plots ####

output$rg_loglinear_r_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$rg_loglinear_g_run
    
    if (g_resetresult == FALSE || is.null(rg_loglinear_results()$Plot1))
        return()
    
    local({
        output$rg_loglinear_r_downloadplot1 <- downloadPlot(
            rg_loglinear_results()$Plot1
        )
    })
    
    rg_loglinear_results()$Plot1
})

output$rg_loglinear_r_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$rg_loglinear_g_run

    if (g_resetresult == FALSE || is.null(rg_loglinear_results()$Plot1))
        return()
    
    div(
        style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("rg_loglinear_r_plot1", height="600px"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left:0.5em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "rg_loglinear_r_downloadplot1", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})

# Log-linear Prediction ####

output$rg_loglinear_r_prediction <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_loglinear_g_run

    if (g_resetresult == FALSE || is.null(rg_loglinear_results()$Model)) {
        return()
    }
    regression <- rg_loglinear_results()$Model
    predictionOutput <- data2
    
    pred <- predict(regression, se.fit=TRUE)
    mu1<-pred$fit
    se1<-pred$se.fit
    LL1<-mu1-1.96*se1
    UL1<-mu1+1.96*se1
    
    if (rg_loglinear_results()$Option$link == "log") {
        mu1<-exp(mu1)
        LL1<-exp(LL1)
        UL1<-exp(UL1)
    }
      
    StudentResidual <- rstandard(regression)
    predictionOutput$pred <- mu1
    predictionOutput$pred95LL <- LL1
    predictionOutput$pred95UL <- UL1
    predictionOutput$StudentResidual <- StudentResidual
    predictionOutput$leverage <- glm.diag(regression)$h
    predictionOutput$cook <- glm.diag(regression)$cook

    return(predictionOutput)
}, rownames = TRUE, bordered = TRUE, spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)