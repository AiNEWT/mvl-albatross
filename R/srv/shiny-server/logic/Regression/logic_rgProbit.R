# Probit Run ####

observeEvent(input$rg_probit_g_run, {
    g_resetresult <<- TRUE
})

rg_probit_results <- eventReactive(input$rg_probit_g_run, {
    # ¦§ Warning & Notify ####
    if (input$rg_probit_m_resp == "") {
      showNotification("Please choose response", type="warning")
      return()
    }
    
    if (input$rg_probit_m_dist == "binomial" && input$rg_probit_m_check_binomd == TRUE && input$rg_probit_m_binomd == "") {
        showNotification("Please choose binomial denominator.", type="warning")
        return()
    }
  
    if (input$rg_probit_m_check_vif == TRUE && length(input$rg_probit_m_variable$right) < 2) {
        showNotification("Model contains fewer than 2 terms.", type="warning")
        return()
    }
    
    # Start Progressbar
    withProgress(message = 'Probit', style = "notification", value = 0, {
      
    # ¦§ GLM ####
        
    probitResults = NULL
    
    if(input$rg_probit_m_check_offset){
        probitResults <- glm(
            as.formula(input$rg_probit_m_model),
            family = eval(parse(text = paste(input$rg_probit_m_dist,"(link = ",input$rg_probit_m_link,")", sep = ""))),
            data = data2,
            offset = data2[[input$rg_probit_m_offset]]
        )
    } else {
        probitResults <- glm(
            as.formula(input$rg_probit_m_model),
            family = eval(parse(text = paste(input$rg_probit_m_dist,"(link = ",input$rg_probit_m_link,")", sep = ""))),
            data = data2
        )
    }
    
    # ¦§ Summary ####
    
    summaryResults <- summary(probitResults)

    M.summary1 <- c(
        "Family" = input$rg_probit_m_dist,
        "Link" = input$rg_probit_m_link,
        "Optimization" = "IWLS",
        "Num of iteration" = probitResults$iter
    )
    
    M.summary2 <- c(
        "Num of obs" = nobs(probitResults),
        "Res. deviance" = probitResults$deviance,
        "Null deviance" = probitResults$null.deviance,
        "Log likelihood" = logLik(probitResults)[[1]],
        "AIC" = AIC(probitResults),
        "BIC" = BIC(probitResults)
    )
    names(M.summary2)[2] <- sprintf("Res. deviance(df=%g)", probitResults$df.residual)
    names(M.summary2)[3] <- sprintf("Null deviance(df=%g)", probitResults$df.null)
    
    # ¦§ Coefficients ####
    
    coefResults = NULL
    if(input$rg_probit_m_check_exp){
        lengthCoef <- length(coef(probitResults))
        CoefMatrix1 <- matrix(exp(coef(probitResults)), lengthCoef, 1)
        CoefMatrix2 <- matrix(summaryResults$coefficients[,c(-1,-2)], lengthCoef, 2)
        
        if (probitResults$df.residual==0 && (input$rg_probit_m_dist=="gaussian" | input$rg_probit_m_dist=="gamma")) {
            CoefMatrix3 <- matrix(NA, lengthCoef, 2)
        } else {
            CoefMatrix3 <- matrix(exp(confint(probitResults)), lengthCoef, 2)
        }
        coefResults<-cbind(CoefMatrix1, CoefMatrix2)
        colnames(coefResults) <- cbind("Estimate(Exp)","t-value","Pr(>|t|)")
        if(input$rg_probit_m_check_confint) {
            coefResults <- cbind(coefResults, CoefMatrix3)
            colnames(coefResults) <- cbind("Estimate(Exp)","t-value","Pr(>|t|)","2.5%","97.5%")
        }
        rownames(coefResults) <- names(coef(probitResults))
    } else {
        coefResults <- summaryResults$coefficients
        lengthCoef <- nrow(coefResults)
        CoefMatrix1 <- matrix(coefResults, lengthCoef, 4)
        if(input$rg_probit_m_check_confint) {
            if (probitResults$df.residual==0 && (input$rg_probit_m_dist=="gaussian" || input$rg_probit_m_dist=="gamma")) {
                CoefMatrix3 <- matrix(NA, lengthCoef, 2)
            } else CoefMatrix3<-matrix(confint(probitResults), lengthCoef, 2)
            if(input$rg_probit_m_check_robustse) {
                coefResults <- unclass(coeftest(probitResults, vcov=vcovHC(probitResults,"HC1")))
            }
            coefResults <- cbind(CoefMatrix1 ,CoefMatrix3)
            colnames(coefResults) <- cbind("Estimate","Std.Error","t-value","Pr(>|t|)","2.5%","97.5%")
            rownames(coefResults) <- rownames(summaryResults$coefficients)
        }
    }
    
    if(input$rg_probit_m_check_vif) {
        coefResults <- cbind(coefResults, VIF = c(NA, vif(probitResults)))
    }
    
    # ¦§ Margins ####
    
    probitMargins <- NULL
    if (all(input$rg_probit_m_check_margins, !all(is.null(input$rg_probit_m_fmargins), is.null(input$rg_probit_m_cmargins)))) {
        probitMarginlist <- list()
        if(input$rg_probit_m_check_offset == FALSE) {
            if(!is.null(input$rg_probit_m_fmargins)){
                probitFmargins <- input$rg_probit_m_fmargins
                probitFmarginlist <- lapply(probitFmargins, function(x) levels(as.factor(data2[[x]])))
                names(probitFmarginlist) <- probitFmargins
                probitMarginlist <- append(probitMarginlist, probitFmarginlist)
            }
            if(!is.null(input$rg_probit_m_cmargins)){
                probitCmargins <- input$rg_probit_m_cmargins
                probitCmarginlist <- lapply(1:length(probitCmargins), function(x){
                    eval(parse(text = paste("input$rg_probit_m_marginvalue", x, sep = "")))
                })
                names(probitCmarginlist) <- probitCmargins
                probitMarginlist <- append(probitMarginlist, probitCmarginlist)
            }
        probitMargins <- append(probitMargins, summary(margins(probitResults, at = probitMarginlist)))
        }
    }

    # ¦§ Comparision Model ####
    
    probitComparisonmodel <- NULL
    if (input$rg_probit_m_check_comparison) {
        if (input$rg_probit_m_check_offset == TRUE) {
            probitComparisonmodel <- matrix(c(
                input$rg_probit_m_model, 
                input$rg_probit_m_dist, 
                input$rg_probit_m_link, 
                input$rg_probit_m_offset,
                probitResults$df.residual, 
                round(probitResults$deviance, digits=5), 
                round(AIC(probitResults), digits=5), 
                round(BIC(probitResults), digits=5)),
                ncol = 8
            )
        } else {
            probitComparisonmodel <- matrix(c(
                input$rg_probit_m_model, 
                input$rg_probit_m_dist, 
                input$rg_probit_m_link, 
                NA,
                probitResults$df.residual, 
                round(probitResults$deviance, digits=5),
                round(AIC(probitResults), digits=5), 
                round(BIC(probitResults), digits=5)),
                ncol = 8
            )
        }
        probitComparisonmodel <- rbind(g_rg_glm_comparisonmodel, probitComparisonmodel)
        colnames(probitComparisonmodel) <- c("Model", "Distribution", "Link", "Offset", "Res.df", "Res.dev", "AIC", "BIC")
        g_rg_glm_comparisonmodel <<- probitComparisonmodel
    }

    # ¦§ R Codes ####
    
    Rcodes <- NULL
    option <- list(
        form = input$rg_probit_m_model,
        dist = input$rg_probit_m_dist,
        link = input$rg_probit_m_link,
        checkbinomd = input$rg_probit_m_check_binomd,
        binomd = input$rg_probit_m_binomd
    )
    if (input$rg_probit_m_check_rcodes) {
        rform <- as.character(as.formula(option$form))
        Rcodes <- paste0("glm(formula = ", paste(rform[2], rform[1], rform[3]), ", ")
        family <- paste0("family = ", option$dist, "(", option$link, ")")
        Rcodes <- paste0(Rcodes, family,", data = data)")
        Rcodes <- matrix(Rcodes)
        colnames(Rcodes) <- "Call"
    }
    
    # Increase Progressbar
    incProgress(0.5, detail = paste("Loading..."))
    
    # ¦¦ Plots ####
    probitPlots <- ggplotglm(probitResults)
    
    # Increase Progressbar
    incProgress(0.5, detail = paste("Loading..."))    
    
    return(list(
        Title = input$rg_probit_m_model,
        Model = probitResults,
        M.summary1 = M.summary1,
        M.summary2 = M.summary2,
        Coefficients = coefResults,
        Margins = probitMargins,
        Rcodes = Rcodes,
        Comparisonmodel = probitComparisonmodel,
        Option = option,
        Plot1 = probitPlots
    ))
    
    }) # End Progressbar
})

# Probit Components ####

output$rg_probit_m_model<-renderUI({
    input$rg_probit_m_resp
    input$rg_probit_m_variable$right
    
    if (all(!is.null(input$rg_probit_m_check_binomd), input$rg_probit_m_check_binomd)) {
        modelEquation1 = paste("cbind(", input$rg_probit_m_resp, ",", input$rg_probit_m_binomd, "-", input$rg_probit_m_resp, ")")
    } else {
        modelEquation1 = input$rg_probit_m_resp
    }
    
    if (length(input$rg_probit_m_variable$right) >= 1) {
        if (!is.null(input$rg_probit_m_check_nointercept) && input$rg_probit_m_check_nointercept == TRUE)
            modelEquation2 = paste0(c("-1", input$rg_probit_m_variable$right), collapse="+")
        else
            modelEquation2 = paste0(input$rg_probit_m_variable$right, collapse="+")
    } else {
        if (!is.null(input$rg_probit_m_check_nointercept) && input$rg_probit_m_check_nointercept == TRUE)
            modelEquation2 = "-1"
        else
            modelEquation2 = "1"
    }
    
    modelEquation = paste(modelEquation1, "~", modelEquation2)
    
    textAreaInput("rg_probit_m_model", "Model", value = modelEquation, height = "80px")
})

output$rg_probit_m_resp <- renderUI({
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
            "rg_probit_m_resp",
            "Response Variable",
            choices = as.list(c("", nameValue)), 
            selected = NULL, 
            multiple = FALSE
        ),
        bsTooltip(
            "rg_probit_m_resp", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$rg_probit_m_variable <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$rg_probit_m_interactionappend
    
    nameValue <- c(names(data2), g_rg_probit_interaction)
    
    chooserInput("rg_probit_m_variable", "Variable", "Selected", nameValue, c(), size = 15, multiple = TRUE, width = 150)
})

output$rg_probit_m_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue <- names(data2)
    
    selectInput(
        "rg_probit_g_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$rg_probit_m_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    actionButton("rg_probit_m_interactionappend", "Append")
})

observeEvent(input$rg_probit_m_interactionappend, {
    g_rg_probit_interaction <<- c(g_rg_probit_interaction, paste(input$rg_probit_g_interaction, collapse=":"))
})

output$rg_probit_m_check_binomd<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$rg_probit_m_dist) || input$rg_probit_m_dist != "binomial")
        return()
    
    checkboxInput("rg_probit_m_check_binomd", "Binomial Denominator", value = FALSE)
})
    
output$rg_probit_m_binomd<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if (is.null(input$rg_probit_m_check_binomd) || !input$rg_probit_m_check_binomd)
        return()
    
    nameValue = names(data2)
    
    selectInput(
        "rg_probit_m_binomd",
        "Binomial Denominator",
        choices = as.list(c("", nameValue)),
        multiple = FALSE
    )
})
    
output$rg_probit_m_dist<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    dist = c(
        "binomial"="binomial"
    )
    
    selectInput(
        "rg_probit_m_dist",
        "Distribution", 
        choices = dist, 
        multiple = FALSE
    )
})

output$rg_probit_m_link<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    link = NULL
    selection = NULL
    
    if(is.null(input$rg_probit_m_dist))
        return()

    selection="probit"
    link = c("probit"="probit")
    
    selectInput("rg_probit_m_link", "Link Function", choices = link, selected = selection, multiple = FALSE)
})

output$rg_probit_m_check_nointercept<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_probit_m_check_nointercept", "No Intercept Model", value = FALSE)
})

output$rg_probit_m_check_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_probit_m_check_offset", "Offset Variable", value = FALSE)
})
    
output$rg_probit_m_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if (is.null(input$rg_probit_m_check_offset) || !input$rg_probit_m_check_offset)
        return()
    
    nameValue = names(select_if(data2, is.numeric))
    
    div(
        selectInput(
            "rg_probit_m_offset",
            "Offset Variable",
            choices = as.list(c("", nameValue)),
            selected = NULL,
            multiple = FALSE
        ),
        bsTooltip(
            "rg_probit_m_offset", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$rg_probit_m_check_margins<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_probit_m_check_margins", "Margins", value = FALSE)
})

output$rg_probit_m_margins1<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$rg_probit_m_check_margins) || !input$rg_probit_m_check_margins || is.null(input$rg_probit_m_variable$right))
        return()
    
    Factornames = names(select_if(data2, is.factor))
    Contiousnames = names(select_if(data2, is.numeric))
    Factornames = as.list(intersect(Factornames, input$rg_probit_m_variable$right))
    Contiousnames = as.list(intersect(Contiousnames, input$rg_probit_m_variable$right))
    
    tagList(
        selectInput(
            "rg_probit_m_fmargins",
            "Factor Variables", 
            choices = Factornames, 
            multiple = TRUE
        ),
        bsTooltip(
            "rg_probit_m_fmargins", 
            "Factor Only",
            "right", 
            options = list(container = "body")
        ),
        selectInput(
            "rg_probit_m_cmargins",
            "Contiuous Variables", 
            choices = Contiousnames, 
            multiple = TRUE
        ),
        bsTooltip(
            "rg_probit_m_cmargins", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$rg_probit_m_margins2<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$rg_probit_m_check_margins) || !input$rg_probit_m_check_margins || is.null(input$rg_probit_m_cmargins))
        return()
    
    lapply(1:length(input$rg_probit_m_cmargins), function(i) {
        numericInput(
            inputId = paste0("rg_probit_m_marginvalue", i), 
            label = input$rg_probit_m_cmargins[i],
            value = mean(data2[[input$rg_probit_m_cmargins[i]]]),
            max = max(data2[[input$rg_probit_m_cmargins[i]]]), 
            min = min(data2[[input$rg_probit_m_cmargins[i]]])
        )
    })
})

output$rg_probit_m_check_vif<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_probit_m_check_vif", "VIF", value = FALSE)
})

output$rg_probit_m_check_robustse<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_probit_m_check_robustse", "Robust Standard Errors", value = FALSE)
})

output$rg_probit_m_check_confint<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_probit_m_check_confint", "Confidence Intervals for Coefficients", value = FALSE)
})

output$rg_probit_m_check_bptest<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_probit_m_check_exp", "Exponential Scale", value = FALSE)
})

output$rg_probit_m_check_comparison<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_probit_m_check_comparison", "model comparison", value = FALSE)
})
    
output$rg_probit_m_check_rcodes<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_probit_m_check_rcodes", "R Codes", value = FALSE)
})

# Probit Results ####

output$rg_probit_r_model <- renderText({
    input$file
    input$resetData
    input$di_option_run
    input$rg_probit_g_run
    
    if (g_resetresult == FALSE || is.null(rg_probit_results()$Title))
        titleOutput = ""
    else
        titleOutput = rg_probit_results()$Title
    
    paste("Model : ", titleOutput)
})

output$rg_probit_r_modelsummary1<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_probit_g_run
    
    if (g_resetresult == FALSE || is.null(rg_probit_results()$M.summary1))
        return()
  
    rg_probit_results()$M.summary1
}, rownames = TRUE, bordered = TRUE, caption = "Model Summary 1", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_probit_r_modelsummary2<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_probit_g_run
    
    if (g_resetresult == FALSE || is.null(rg_probit_results()$M.summary2))
        return()
  
    rg_probit_results()$M.summary2
}, colnames = FALSE, rownames = TRUE, bordered = TRUE, caption = "Model Summary 2", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_probit_r_coefficients<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_probit_g_run
    
    if (g_resetresult == FALSE || is.null(rg_probit_results()$Coefficients))
        return()
  
    rg_probit_results()$Coefficients 
}, rownames = TRUE, bordered = TRUE, caption = "Coefficients", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_probit_r_margins<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_probit_g_run
    
    if (g_resetresult == FALSE || is.null(rg_probit_results()$Margins))
        return()
  
    rg_probit_results()$Margins 
}, rownames = TRUE, bordered = TRUE, caption = "Margins", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_probit_r_comparisonmodel<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_probit_g_run
    
    if (g_resetresult == FALSE || is.null(rg_probit_results()$Comparisonmodel))
        return()
  
    rg_probit_results()$Comparisonmodel
}, rownames = TRUE, bordered = TRUE, caption = "Comparison Model", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_probit_r_rcodes<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_probit_g_run
    
    if (g_resetresult == FALSE || is.null(rg_probit_results()$Rcodes))
        return()
  
    rg_probit_results()$Rcodes
}, rownames = FALSE, bordered = TRUE, caption = "R Codes", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

# Probit Plots ####

output$rg_probit_r_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$rg_probit_g_run
    
    if (g_resetresult == FALSE || is.null(rg_probit_results()$Plot1))
        return()

    local({
        output$rg_probit_r_downloadplot1 <- downloadPlot(
            rg_probit_results()$Plot1
        )
    })
    
    rg_probit_results()$Plot1
})

output$rg_probit_r_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$rg_probit_g_run

    if (g_resetresult == FALSE || is.null(rg_probit_results()$Plot1))
        return()
    
    div(
        style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("rg_probit_r_plot1", height="600px"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left:0.5em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "rg_probit_r_downloadplot1", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})

# Probit Prediction ####

output$rg_probit_r_prediction <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_probit_g_run

    if (g_resetresult == FALSE || is.null(rg_probit_results()$Model)) {
        return()
    }
    regression <- rg_probit_results()$Model
    predictionOutput <- data2
    
    pred <- predict(regression, se.fit=TRUE)
    mu1<-pred$fit
    se1<-pred$se.fit
    LL1<-mu1-1.96*se1
    UL1<-mu1+1.96*se1
    
    if (rg_probit_results()$Option$link == "probit") {
        mu1<-pnorm(mu1)
        LL1<-pnorm(LL1)
        UL1<-pnorm(UL1)
    } 
    
    if (!is.null(rg_probit_results()$Option$checkbinomd) && rg_probit_results()$Option$checkbinomd) {
        binomD<-rg_probit_results()$Option$binomd
        mu1<-mu1*data2[[binomD]]
        LL1<-LL1*data2[[binomD]]
        UL1<-UL1*data2[[binomD]]
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