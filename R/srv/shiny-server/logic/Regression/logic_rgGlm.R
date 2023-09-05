# GLM Run ####

observeEvent(input$rg_glm_g_run, {
    g_resetresult <<- TRUE
})

rg_glm_results <- eventReactive(input$rg_glm_g_run, {
    # ¦§ Warning & Notify ####
    if (input$rg_glm_m_resp == "") { 
      showNotification("Please choose response", type="warning")
      return()
    }
    
    if (input$rg_glm_m_dist == "binomial" && input$rg_glm_m_check_binomd == TRUE && input$rg_glm_m_binomd == "") {
        showNotification("Please choose binomial denominator.", type="warning")
        return()
    }
    
    if (input$rg_glm_m_check_vif == TRUE && length(input$rg_glm_m_variable$right) < 2) {
        showNotification("Model contains fewer than 2 terms.", type="warning")
        return()
    }
    
    # Start Progressbar
    withProgress(message = 'GLM', style = "notification", value = 0, {
      
    # ¦§ GLM ####
      
    glmResults = NULL
    
    if(input$rg_glm_m_check_offset){
        glmResults <- glm(
            as.formula(input$rg_glm_m_model),
            family = eval(parse(text = paste(input$rg_glm_m_dist,"(link = ",input$rg_glm_m_link,")", sep = ""))),
            data = data2,
            offset = data2[[input$rg_glm_m_offset]]
        )
    } else {
        glmResults <- glm(
            as.formula(input$rg_glm_m_model),
            family = eval(parse(text = paste(input$rg_glm_m_dist,"(link = ",input$rg_glm_m_link,")", sep = ""))),
            data = data2
        )
    }
    
    # ¦§ Summary ####
    
    summaryResults <- summary(glmResults)

    M.summary1 <- c(
        "Family" = input$rg_glm_m_dist,
        "Link" = input$rg_glm_m_link,
        "Optimization" = "IWLS",
        "Num of iteration" = glmResults$iter
    )
    
    M.summary2 <- c(
        "Num of obs" = nobs(glmResults),
        "Res. deviance" = glmResults$deviance,
        "Null deviance" = glmResults$null.deviance,
        "Log likelihood" = logLik(glmResults)[[1]],
        "AIC" = AIC(glmResults),
        "BIC" = BIC(glmResults)
    )
    names(M.summary2)[2] <- sprintf("Res. deviance(df=%g)", glmResults$df.residual)
    names(M.summary2)[3] <- sprintf("Null deviance(df=%g)", glmResults$df.null)
    
    # ¦§ Coefficients ####
    
    coefResults = NULL
    if(input$rg_glm_m_check_exp){
        lengthCoef <- length(coef(glmResults))
        CoefMatrix1 <- matrix(exp(coef(glmResults)), lengthCoef, 1)
        CoefMatrix2 <- matrix(summaryResults$coefficients[,c(-1,-2)], lengthCoef, 2)
        
        if (glmResults$df.residual==0 && (input$rg_glm_m_dist=="gaussian" | input$rg_glm_m_dist=="gamma")) {
            CoefMatrix3 <- matrix(NA, lengthCoef, 2)
        } else {
            CoefMatrix3 <- matrix(exp(confint(glmResults)), lengthCoef, 2)
        }
        coefResults<-cbind(CoefMatrix1, CoefMatrix2)
        colnames(coefResults) <- cbind("Estimate(Exp)","t-value","Pr(>|t|)")
        if(input$rg_glm_m_check_confint) {
            coefResults <- cbind(coefResults, CoefMatrix3)
            colnames(coefResults) <- cbind("Estimate(Exp)","t-value","Pr(>|t|)","2.5%","97.5%")
        
        }
        rownames(coefResults) <- names(coef(glmResults))
    } else {
        coefResults <- summaryResults$coefficients
        lengthCoef <- nrow(coefResults)
        CoefMatrix1 <- matrix(coefResults, lengthCoef, 4)
        if(input$rg_glm_m_check_confint) {
            if (glmResults$df.residual==0 && (input$rg_glm_m_dist=="gaussian" || input$rg_glm_m_dist=="gamma")) {
                CoefMatrix3 <- matrix(NA, lengthCoef, 2)
            } else CoefMatrix3<-matrix(confint(glmResults), lengthCoef, 2)
            if(input$rg_glm_m_check_robustse) {
                coefResults <- unclass(coeftest(glmResults, vcov=vcovHC(glmResults,"HC1")))
            }
            coefResults <- cbind(CoefMatrix1 ,CoefMatrix3)
            colnames(coefResults) <- cbind("Estimate","Std.Error","t-value","Pr(>|t|)","2.5%","97.5%")
            rownames(coefResults) <- rownames(summaryResults$coefficients)
        }
    }
    
    if(input$rg_glm_m_check_vif) {
        coefResults <- cbind(coefResults, VIF = c(NA, vif(glmResults)))
    }
    
    # ¦§ Margins ####
    
    glmMargins <- NULL
    if (all(input$rg_glm_m_check_margins, !all(is.null(input$rg_glm_m_fmargins), is.null(input$rg_glm_m_cmargins)))) {
        glmMarginlist <- list()
        if(input$rg_glm_m_check_offset == FALSE) {
            if(!is.null(input$rg_glm_m_fmargins)){
                glmFmargins <- input$rg_glm_m_fmargins
                glmFmarginlist <- lapply(glmFmargins, function(x) levels(as.factor(data2[[x]])))
                names(glmFmarginlist) <- glmFmargins
                glmMarginlist <- append(glmMarginlist, glmFmarginlist)
            }
            if(!is.null(input$rg_glm_m_cmargins)){
                glmCmargins <- input$rg_glm_m_cmargins
                glmCmarginlist <- lapply(1:length(glmCmargins), function(x){
                    eval(parse(text = paste("input$rg_glm_m_marginvalue", x, sep = "")))
                })
                names(glmCmarginlist) <- glmCmargins
                glmMarginlist <- append(glmMarginlist, glmCmarginlist)
            }
        glmMargins <- append(glmMargins, summary(margins(glmResults, at = glmMarginlist)))
        }
    }

    # ¦§ Comparision Model ####
    
    glmComparisonmodel <- NULL
    if (input$rg_glm_m_check_comparison) {
      if (input$rg_glm_m_dist=="gaussian" | input$rg_glm_m_dist=="Gamma" | input$rg_glm_m_dist=="inverse.gaussian") {
        if (input$rg_glm_m_check_offset == TRUE) {
            glmComparisonmodel <- matrix(c(
                input$rg_glm_m_model, 
                input$rg_glm_m_dist, 
                input$rg_glm_m_link, 
                input$rg_glm_m_offset,
                round(glmResults$df.residual, digits=5),
                glmResults$df.residual,
                round(AIC(glmResults), digits=5), 
                round(BIC(glmResults), digits=5)),
                ncol = 8
            )
        } else {
            glmComparisonmodel <- matrix(c(
                input$rg_glm_m_model, 
                input$rg_glm_m_dist, 
                input$rg_glm_m_link, 
                NA,
                round(glmResults$df.residual, digits=5),
                glmResults$df.residual, 
                round(AIC(glmResults), digits=5), 
                round(BIC(glmResults), digits=5)),
                ncol = 8
            )
        }
      } else {
        if (input$rg_glm_m_check_offset == TRUE) {
          glmComparisonmodel <- matrix(c(
            input$rg_glm_m_model, 
            input$rg_glm_m_dist, 
            input$rg_glm_m_link, 
            input$rg_glm_m_offset,
            round(glmResults$deviance, digits=5), 
            glmResults$df.residual,
            round(AIC(glmResults), digits=5), 
            round(BIC(glmResults), digits=5)),
            ncol = 8
          )
        } else {
          glmComparisonmodel <- matrix(c(
            input$rg_glm_m_model, 
            input$rg_glm_m_dist, 
            input$rg_glm_m_link, 
            NA,
            round(glmResults$deviance, digits=5),
            glmResults$df.residual, 
            round(AIC(glmResults), digits=5), 
            round(BIC(glmResults), digits=5)),
            ncol = 8
          )
        }        
      }
        glmComparisonmodel <- rbind(g_rg_glm_comparisonmodel, glmComparisonmodel)
        colnames(glmComparisonmodel) <- c("Model", "Distribution", "Link", "Offset","Scaled Deviance","df", "AIC", "BIC")
        g_rg_glm_comparisonmodel <<- glmComparisonmodel
    }
    
    # ¦§ R Codes ####
    
    Rcodes <- NULL
    option <- list(
        form = input$rg_glm_m_model,
        dist = input$rg_glm_m_dist,
        link = input$rg_glm_m_link,
        checkbinomd = input$rg_glm_m_check_binomd,
        binomd = input$rg_glm_m_binomd
    )
    if (input$rg_glm_m_check_rcodes) {
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
    glmPlots <- ggplotglm(glmResults)
    
    # Increase Progressbar
    incProgress(0.5, detail = paste("Loading..."))
     
    return(list(
        Title = input$rg_glm_m_model,
        Model = glmResults,
        M.summary1 = M.summary1,
        M.summary2 = M.summary2,
        Coefficients = coefResults,
        Margins = glmMargins,
        Rcodes = Rcodes,
        Comparisonmodel = glmComparisonmodel,
        Option = option,
        Plot1 = glmPlots
    ))
    
    }) # End Progressbar
})

# GLM Components ####

output$rg_glm_m_model<-renderUI({
    input$rg_glm_m_resp
    input$rg_glm_m_variable$right
    
    if (all(!is.null(input$rg_glm_m_check_binomd), input$rg_glm_m_check_binomd)) {
        modelEquation1 = paste("cbind(", input$rg_glm_m_resp, ",", input$rg_glm_m_binomd, "-", input$rg_glm_m_resp, ")")
    } else {
        modelEquation1 = input$rg_glm_m_resp
    }
    
    if (length(input$rg_glm_m_variable$right) >= 1) {
        if (!is.null(input$rg_glm_m_check_nointercept) && input$rg_glm_m_check_nointercept == TRUE)
            modelEquation2 = paste0(c("-1", input$rg_glm_m_variable$right), collapse="+")
        else
            modelEquation2 = paste0(input$rg_glm_m_variable$right, collapse="+")
    } else {
        if (!is.null(input$rg_glm_m_check_nointercept) && input$rg_glm_m_check_nointercept == TRUE)
            modelEquation2 = "-1"
        else
            modelEquation2 = "1"
    }
    
    modelEquation = paste(modelEquation1, "~", modelEquation2)
    
    textAreaInput("rg_glm_m_model", "Model", value = modelEquation, height = "80px")
})

output$rg_glm_m_resp <- renderUI({
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
            "rg_glm_m_resp",
            "Response Variable",
            choices = as.list(c("", nameValue)), 
            selected = NULL, 
            multiple = FALSE
        ),
        bsTooltip(
            "rg_glm_m_resp", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$rg_glm_m_variable <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$rg_glm_m_interactionappend
    
    nameValue <- c(names(data2), g_rg_glm_interaction)
    
    chooserInput("rg_glm_m_variable", "Variable", "Selected", nameValue, c(), size = 15, multiple = TRUE, width = 125)
})

output$rg_glm_m_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue <- names(data2)
    
    selectInput(
        "rg_glm_g_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$rg_glm_m_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    actionButton("rg_glm_m_interactionappend", "Append")
})

observeEvent(input$rg_glm_m_interactionappend, {
    g_rg_glm_interaction <<- c(g_rg_glm_interaction, paste(input$rg_glm_g_interaction, collapse=":"))
})

output$rg_glm_m_check_binomd<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$rg_glm_m_dist) || input$rg_glm_m_dist != "binomial")
        return()
    
    checkboxInput("rg_glm_m_check_binomd", "Binomial Denominator", value = FALSE)
})
    
output$rg_glm_m_binomd<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if (is.null(input$rg_glm_m_check_binomd) || !input$rg_glm_m_check_binomd)
        return()
    
    nameValue = names(data2)
    
    selectInput(
        "rg_glm_m_binomd",
        "Binomial Denominator",
        choices = as.list(c("", nameValue)),
        multiple = FALSE
    )
})
    
output$rg_glm_m_dist<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    dist = c(
        "gaussian"="gaussian",
        "binomial"="binomial",
        "poisson"="poisson",
        "gamma"="Gamma",
        "inverse.gaussian"="inverse.gaussian"
    )
    
    selectInput(
        "rg_glm_m_dist",
        "Distribution", 
        choices = dist, 
        multiple = FALSE
    )
})

output$rg_glm_m_link<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    link = NULL
    selection = NULL
    
    if(is.null(input$rg_glm_m_dist))
        return()

    if (input$rg_glm_m_dist == "gaussian"){
        selection ="identity"
        link = c("identity"="identity", "log"="log", "inverse"="inverse")
    } else if (input$rg_glm_m_dist == "binomial"){
        selection="logit"
        link = c("logit"="logit", "probit"="probit", "cauchit"="cauchit", "log"="log", "cloglog"="cloglog", "identity"="identity")
    } else if (input$rg_glm_m_dist == "poisson"){
        selection="log"
        link = c("log"="log", "identity"="identity", "sqrt"="sqrt")
    } else if (input$rg_glm_m_dist == "Gamma"){
        selection="log"
        link = c("inverse"="inverse", "identity"="identity", "log"="log")
    } else if (input$rg_glm_m_dist == "inverse.gaussian"){
        selection="1/mu^2"
        link = c("1/mu^2"="1/mu^2", "inverse"="inverse", "log"="log", "identity"="identity")
    } 
    
    selectInput("rg_glm_m_link", "Link Function", choices = link, selected = selection, multiple = FALSE)
})

output$rg_glm_m_check_nointercept<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_glm_m_check_nointercept", "No Intercept Model", value = FALSE)
})

output$rg_glm_m_check_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_glm_m_check_offset", "Offset Variable", value = FALSE)
})
    
output$rg_glm_m_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if (is.null(input$rg_glm_m_check_offset) || !input$rg_glm_m_check_offset)
        return()
    
    nameValue = names(select_if(data2, is.numeric))
    
    div(
        selectInput(
            "rg_glm_m_offset",
            "Offset Variable",
            choices = as.list(c("", nameValue)),
            selected = NULL,
            multiple = FALSE
        ),
        bsTooltip(
            "rg_glm_m_offset", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$rg_glm_m_check_margins<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_glm_m_check_margins", "Margins", value = FALSE)
})

output$rg_glm_m_margins1<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$rg_glm_m_check_margins) || !input$rg_glm_m_check_margins || is.null(input$rg_glm_m_variable$right))
        return()
    
    Factornames = names(select_if(data2, is.factor))
    Contiousnames = names(select_if(data2, is.numeric))
    Factornames = as.list(intersect(Factornames, input$rg_glm_m_variable$right))
    Contiousnames = as.list(intersect(Contiousnames, input$rg_glm_m_variable$right))
    
    tagList(
        selectInput(
            "rg_glm_m_fmargins",
            "Factor Variables", 
            choices = Factornames, 
            multiple = TRUE
        ),
        bsTooltip(
            "rg_glm_m_fmargins", 
            "Factor Only",
            "right", 
            options = list(container = "body")
        ),
        selectInput(
            "rg_glm_m_cmargins",
            "Contiuous Variables", 
            choices = Contiousnames, 
            multiple = TRUE
        ),
        bsTooltip(
            "rg_glm_m_cmargins", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$rg_glm_m_margins2<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$rg_glm_m_check_margins) || !input$rg_glm_m_check_margins || is.null(input$rg_glm_m_cmargins))
        return()
    
    lapply(1:length(input$rg_glm_m_cmargins), function(i) {
        numericInput(
            inputId = paste0("rg_glm_m_marginvalue", i), 
            label = input$rg_glm_m_cmargins[i],
            value = mean(data2[[input$rg_glm_m_cmargins[i]]]),
            max = max(data2[[input$rg_glm_m_cmargins[i]]]), 
            min = min(data2[[input$rg_glm_m_cmargins[i]]])
        )
    })
})

output$rg_glm_m_check_vif<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_glm_m_check_vif", "VIF", value = FALSE)
})

output$rg_glm_m_check_robustse<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_glm_m_check_robustse", "Robust Standard Errors", value = FALSE)
})

output$rg_glm_m_check_confint<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_glm_m_check_confint", "Confidence Intervals for Coefficients", value = FALSE)
})

output$rg_glm_m_check_exp<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_glm_m_check_exp", "Exponential Scale", value = FALSE)
})

output$rg_glm_m_check_comparison<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_glm_m_check_comparison", "model comparison", value = FALSE)
})
    
output$rg_glm_m_check_rcodes<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_glm_m_check_rcodes", "R Codes", value = FALSE)
})

# GLM Results ####

output$rg_glm_r_model <- renderText({
    input$file
    input$resetData
    input$di_option_run
    input$rg_glm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_glm_results()$Title))
        titleOutput = ""
    else
        titleOutput = rg_glm_results()$Title
    
    paste("Model : ", titleOutput)
})

output$rg_glm_r_modelsummary1<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_glm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_glm_results()$M.summary1))
        return()
  
    rg_glm_results()$M.summary1
}, rownames = TRUE, bordered = TRUE, caption = "Model Summary 1", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_glm_r_modelsummary2<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_glm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_glm_results()$M.summary2))
        return()
  
    rg_glm_results()$M.summary2
}, colnames = FALSE, rownames = TRUE, bordered = TRUE, caption = "Model Summary 2", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_glm_r_coefficients<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_glm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_glm_results()$Coefficients))
        return()
  
    rg_glm_results()$Coefficients 
}, rownames = TRUE, bordered = TRUE, caption = "Coefficients", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_glm_r_margins<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_glm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_glm_results()$Margins))
        return()
  
    rg_glm_results()$Margins 
}, rownames = TRUE, bordered = TRUE, caption = "Margins", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_glm_r_comparisonmodel<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_glm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_glm_results()$Comparisonmodel))
        return()
  
    rg_glm_results()$Comparisonmodel
}, rownames = TRUE, bordered = TRUE, caption = "Comparison Model", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_glm_r_rcodes<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_glm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_glm_results()$Rcodes))
        return()
  
    rg_glm_results()$Rcodes
}, rownames = FALSE, bordered = TRUE, caption = "R Codes", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

# GLM Plots ####

output$rg_glm_r_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$rg_glm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_glm_results()$Plot1))
        return()

    local({
        output$rg_glm_r_downloadplot1 <- downloadPlot(
            rg_glm_results()$Plot1
        )
    })
    rg_glm_results()$Plot1
})

output$rg_glm_r_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$rg_glm_g_run

    if (g_resetresult == FALSE || is.null(rg_glm_results()$Plot1))
        return()
    
    div(
        style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("rg_glm_r_plot1", height="600px"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left:0.5em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "rg_glm_r_downloadplot1", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})

# GLM Prediction ####

output$rg_glm_r_prediction <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_glm_g_run

    if (g_resetresult == FALSE || is.null(rg_glm_results()$Model)) {
        return()
    }
    regression <- rg_glm_results()$Model
    predictionOutput <- data2
    
    pred <- predict(regression, se.fit=TRUE)
    mu1<-pred$fit
    se1<-pred$se.fit
    LL1<-mu1-1.96*se1
    UL1<-mu1+1.96*se1
    
    if (rg_glm_results()$Option$link == "identity") {
        mu1<-mu1
        LL1<-LL1
        UL1<-UL1
    } else if (rg_glm_results()$Option$link == "log") {
        mu1<-exp(mu1)
        LL1<-exp(LL1)
        UL1<-exp(UL1)
    } else if (rg_glm_results()$Option$link == "logit") {
        mu1<-1/(1+exp(-mu1))
        LL1<-1/(1+exp(-LL1))
        UL1<-1/(1+exp(-UL1))
    } else if (rg_glm_results()$Option$link == "probit") {
        mu1<-pnorm(mu1)
        LL1<-pnorm(LL1)
        UL1<-pnorm(UL1)
    } else if (rg_glm_results()$Option$link == "cloglog") {
        mu1<-1-exp(-exp(mu1))
        LL1<-1-exp(-exp(LL1))
        UL1<-1-exp(-exp(UL1))
    } else if (rg_glm_results()$Option$link == "inverse") {
        mu1<-1/mu1
        LL1<-1/LL1
        UL1<-1/UL1
    } else if (rg_glm_results()$Option$link == "sqrt") {
        mu1<-mu1^2
        LL1<-LL1^2
        UL1<-UL1^2
    } else if (rg_glm_results()$Option$link == "cauchit") {
        mu1<-cauchitlink(mu1,inverse=TRUE)
        LL1<-cauchitlink(LL1,inverse=TRUE)
        UL1<-cauchitlink(UL1,inverse=TRUE)
    }
      
    if (!is.null(rg_glm_results()$Option$checkbinomd) && rg_glm_results()$Option$checkbinomd) {
        binomD<-rg_glm_results()$Option$binomd
        mu1<-mu1*data2[[binomD]]
        LL1<-LL1*data2[[binomD]]
        UL1<-UL1*data2[[binomD]]
    }
      
    StudentResidual <- rstandard(regression)
    predictionOutput$pred <- mu1
    predictionOutput$pred95LL <- LL1
    predictionOutput$pred95UL <- UL1
    tryCatch(expr = predictionOutput$StudentResidual <- StudentResidual,
        error = function(e) {
            showNotification("Wrong Expression.", type="warning")
        }
    )
    predictionOutput$leverage <- glm.diag(regression)$h
    predictionOutput$cook <- glm.diag(regression)$cook

    return(predictionOutput)
}, rownames = TRUE, bordered = TRUE, spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)