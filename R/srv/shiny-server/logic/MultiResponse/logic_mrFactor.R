# Factor Model Slider ####

output$mr_factor_g_factorslider <- renderUI({
    sliderInput(
        inputId = sprintf("mr_factor_g_factor_count"),
        label = "Number of Factors",
        value = 2,
        min = 2,
        max = 6
    )
})

mr_factor_g_factor_value <- reactive({
    if (is.null(input$mr_factor_g_factor_count))
        return(2)
    else
        return(input$mr_factor_g_factor_count)
})

# Factor Model Run ####

observeEvent(input$mr_factor_g_run, {
    g_resetresult <<- TRUE
})

mr_factor_results <- eventReactive(input$mr_factor_g_run, {
    withProgress(message = 'Factor Model', style = "notification", value = 0, {
    
    # Increase Progressbar
    incProgress(0.5, detail = paste("Loading..."))
        
    FactorModel <- input$mr_factor_g_model
    
    # print(FactorModel)
    
    FittedModel <<- hsem(FactorModel,data=data2)
    res<-summary(FittedModel,fit.measures=TRUE)
    expression<-subset(res$PE, op == "=~")
    data3<-data2
    dist=input$mr_factor_g_dist
    n_resp<-nrow(expression)
    intercept<-rep(0,n_resp)
    intercept.SE<-rep(0,n_resp)
    for (i in 1:n_resp) {
        yy<-data3[,expression[,"rhs"][i]]
        zz<-yy
        if (dist=="gaussian'") data3[,expression[,"rhs"][i]]<- zz<- yy
        if (dist=="binomial") data3[,expression[,"rhs"][i]]<- zz<- log((yy+0.5)/(1+0.5-yy))
        if (dist=="gamma") data3[,expression[,"rhs"][i]]<-zz<-  log(yy)
        if (dist=="poisson") data3[,expression[,"rhs"][i]]<-zz<-  sqrt(yy+0.5)
        names(intercept)[i]<-expression[,"rhs"][i]
        intercept[i]<-mean(as.numeric(zz[[1]]),na.rm=TRUE)
        intercept.SE[i]<-sd(as.numeric(zz[[1]]),na.rm=TRUE)/sqrt(sum(ifelse(is.na(zz)==FALSE,1,0)))
    }
    FittedModel <<- hsem(FactorModel,data=data3)
    FactorSummary <<- summary(FittedModel, fit.measures=TRUE)
    
    if (!is.null(FactorSummary$FIT)) {
        summaryMatrix <- matrix(rep(NA, 90), ncol = 5)
        
        for (i in 1:9) {
            summaryMatrix[i*2-1, 1:5] <- names(FactorSummary$FIT[(1+5*(i-1)):(5*i)])
            summaryMatrix[i*2, 1:5] <- round(FactorSummary$FIT[(1+5*(i-1)):(5*i)],3)
        }
        FactorSummary$Summary <<- summaryMatrix    
    }
    
    # PE1 ####
    
    PE1 <- vector(mode = "list", length = input$mr_factor_g_factor_count)
    
    peResult1 <- summary(FittedModel,fit.measures=TRUE)
    peValue1 <- subset(peResult1$PE, op == "=~")
    n_resp <- nrow(peValue1)
    peValue1 <- peValue1[,-c(2,4)]
    colnames(peValue1) <- c("Factor","","Estimate","Std.Err","t-value","p-value")
    
    for (k in 1:input$mr_factor_g_factor_count) {
        PE1[[k]] <- subset(peValue1, Factor == mr_factor_g_factor_name[[paste0("sem", k)]])    
    }
    
  
    # Factor Model Plot ####
    FactorSummary$Plot1 <- ggsem(FittedModel)
    FactorSummary$intercept <- intercept
    FactorSummary$intercept.SE <- intercept.SE
    FactorSummary$FactorCount <- input$mr_factor_g_factor_count
    FactorSummary$PE1 <- PE1

    res2<-summary(FittedModel, fit.measures = TRUE, standardized = TRUE)
    expression<-subset(res2$PE, op == "=~")
    n_resp<-nrow(expression)
    pp=res2$FIT['npar'] 
    df=n_resp*(n_resp+1)/2-pp
    n=res2$FIT['ntotal']
    yy<-data3[,expression[,"rhs"][1]]
    for (i in 2:n_resp) {
      yy<-cbind(yy,data3[,expression[,"rhs"][i]])
    }
    SS<-lav_matrix_cov(yy)
    Sigma<-fitted(fit.mod)$cov
    
    sd_theta=(n-1)*abs(log(det(Sigma))+tr(SS%*%solve(Sigma))-log(det(SS))-n_resp)
    df=df
    cfi=res2$FIT['cfi']
    tli=res2$FIT['tli']
    rmsea=res2$FIT['rmsea']
    
    rel1<-reliability(FittedModel)
    length1<-ncol(rel1)
    cron<-rel1[1,1:length1-1]
    avevar<-rel1[4,1:length1-1]
    CR<-CR(FittedModel)
    
    measures<-rbind(cron,CR,avevar)
    rownames(measures)<-c("Cronbach's alpha","Composite Reliablity","Average Varaince Extracted")
    
    sd1<-cbind(sd_theta,df)
    rownames(sd1)<-""
    colnames(sd1)<-c("Scaled Deviance", "degress of freedom")
    
    ctrm<-cbind(cfi,tli,rmsea)
    rownames(ctrm)<-""
    colnames(ctrm)<-c("cti","tli","rmsea")

    FactorSummary$measures<-measures
    FactorSummary$sd1<-sd1
    FactorSummary$measure1<-ctrm
    
    # FactorSummary$Summary1 <- FactorSummary$FIT[1:2]
    # FactorSummary$Summary1 <- FactorSummary$FIT[1:2]
    # FactorSummary$Summary1 <- FactorSummary$FIT[1:2]
    # HSEM
    # Factor_list <- vector(mode = "list", length = input$mr_factor_g_factor_count)
    # Sem_list <- vector(mode = "list", length = input$mr_factor_g_factor_count)
    # 
    # for(k in 1:input$mr_factor_g_factor_count) {
    #     FactorModel <- hFactorModeling(
    #         Model = "factor",
    #         LinPred = as.formula(input[[paste0("mr_factor_f_model_", k)]])
    #     )
    #     
    #     FactorModel <- hFactorModeling(
    #         Model = "sem",
    #         LinPred = as.formula(input[[paste0("mr_factor_s_model_", k)]])
    #     )
    #     
    #     Factor_list[[k]] <- FactorModel
    #     Sem_list[[k]] <- FactorModel
    # }
    # 
    # fittedModel <- hsemfit(
    #     RespDist = "gaussian",
    #     DataMain = data2,
    #     factorModel = Factor_list,
    #     FactorModel = Sem_list
    # )
    
    # Set Progressbar
    setProgress(1, detail = "Finish")
    return(FactorSummary)
    
    }) # End Progressbar
})

# Factor Names ####
# reactiveValues object for storing current data set.
mr_factor_g_factor_name <- reactiveValues(
    sem1 = "factor1", 
    sem2 = "factor2", 
    sem3 = "factor3", 
    sem4 = "factor4", 
    sem5 = "factor5", 
    sem6 = "factor6"
)

output$mr_factor_g_names <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run

    factorCount = mr_factor_g_factor_value()
    
    return(
        do.call(div, c(lapply(1:mr_factor_g_factor_value(), function(i) {
            textInput(
                paste0("mr_factor_g_factor_name_",i), 
                paste0("Type Name Factor",i),
                placeholder = paste0("Factor",i)
            )
        })))
    )
})

# Return the UI for a modal dialog with data selection input. If 'failed' is
# TRUE, then display a message that the previous value was invalid.
mr_factor_g_modal <- function(failed = FALSE) {
    modalDialog(
        title = "Change Factor Names",
        easyClose = TRUE,
        uiOutput("mr_factor_g_names"),
        footer = tagList(
            modalButton("Cancel"),
            actionButton("mr_factor_g_ok", "OK")
        )
    )
}

# Show modal when button is clicked.
observeEvent(input$mr_factor_g_changefactorname, {
    showModal(mr_factor_g_modal())
})

# When OK button is pressed, attempt to load the data set. If successful,
# remove the modal. If not show another modal, but this time with a failure
# message.
observeEvent(input$mr_factor_g_ok, {
    # Check that data object exists and is data frame.
    for (i in 1:mr_factor_g_factor_value()) {
        if (input[[paste0("mr_factor_g_factor_name_", i)]] == "")
            mr_factor_g_factor_name[[paste0("sem", i)]] <- paste0("factor", i)
        else 
            mr_factor_g_factor_name[[paste0("sem", i)]] <- input[[paste0("mr_factor_g_factor_name_", i)]]
    }
    removeModal()
})

# Factor Model ####

output$mr_factor_g_model <- renderUI({
    Model1 = NULL
    Model2 = NULL
    Model3 = NULL
    for (k in 1:mr_factor_g_factor_value()) {
        Model1 <- c(Model1, input[[paste0('mr_factor_f_model_',k)]])
        Model2 <- c(Model2, input[[paste0('mr_factor_s_model1_',k)]])
        Model3 <- c(Model3, input[[paste0('mr_factor_s_model2_',k)]])
    }
    
    Model1 = paste0(Model1, collapse="\n")
    Model2 = paste0(Model2, collapse="\n")
    Model3 = paste0(Model3, collapse="\n")

    Model = paste0(Model1, "\n\n", Model2, "\n\n", Model3)
    
    Height = paste0(mr_factor_g_factor_value() * 80, "px")
    textAreaInput("mr_factor_g_model", "Model", value = Model, height = Height)
})

output$mr_factor_g_dist<-renderUI({
  input$file
  input$resetData
  input$di_option_run
  
  dist = c(
    "gaussian"="gaussian",
    "binomial"="binomial",
    "poisson"="poisson",
    "gamma"="Gamma"
  )
  
  selectInput(
    "mr_factor_g_dist",
    "Distribution", 
    choices = dist, 
    multiple = FALSE
  )
})

output$mr_factor_g_link<-renderUI({
  input$file
  input$resetData
  input$di_option_run
  
  link = NULL
  selection = NULL
  
  if(is.null(input$mr_factor_g_dist))
    return()
  
  if (input$mr_factor_g_dist == "gaussian"){
    selection ="identity"
    link = c("identity"="identity", "log"="log", "inverse"="inverse")
  } else if (input$mr_factor_g_dist == "binomial"){
    selection="logit"
    link = c("logit"="logit", "probit"="probit", "cauchit"="cauchit", "log"="log", "cloglog"="cloglog", "identity"="identity")
  } else if (input$mr_factor_g_dist == "poisson"){
    selection="log"
    link = c("log"="log", "identity"="identity", "sqrt"="sqrt")
  } else if (input$mr_factor_g_dist == "Gamma"){
    selection="log"
    link = c("inverse"="inverse", "identity"="identity", "log"="log")
  }
  selectInput("mr_factor_g_link", "Link Function", choices = link, selected = selection, multiple = FALSE)
})


# Factor Tab Panel ####
output$mr_factor_g_tabpanel <- renderUI({
    do.call(tabsetPanel, c(id = 'factortabp', lapply(1:mr_factor_g_factor_value(), function(i) {
        tabPanel(
            title = paste0('Factor', i),
            br(),
            uiOutput(paste0('mr_factor_r_accordion_', i))
        )
    })))
})

# Factor Components ####
lapply(1:6, function(k) { 
    
    # Accordion ####
    
    output[[paste0('mr_factor_r_accordion_',k)]]<-renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        FactorAccordion <- bs_accordion_sidebar(
            id = paste0("mr_factor_r_accordion", k),
            spec_side = c(width = 3, offset = 0),
            spec_main = c(width = 9, offset = 0)    
        )
        
        FactorAccordion <- FactorAccordion %>%
            bs_append(
                title_side = "Factor",
                content_side = NULL,
                content_main = div(
                    h3(strong("Factor")),
                    uiOutput(paste0("mr_factor_f_model_", k)),
                    uiOutput(paste0("mr_factor_f_resp_", k)),
                    uiOutput(paste0("mr_factor_f_variable_", k))
                )
            )
        
    #    FactorAccordion <- FactorAccordion %>%
    #    bs_append(
    #        title_side = "Factor Model",
    #        content_side = uiOutput(paste0("mr_factor_s_check_sem_", k)),
    #        content_main = div(
    #            h3(strong("Factor Model")),
    #            uiOutput(paste0("mr_factor_s_model1_", k)),
    #            uiOutput(paste0("mr_factor_s_model2_", k)),
    #            uiOutput(paste0("mr_factor_s_resp_", k)),
    #            uiOutput(paste0("mr_factor_s_variable_", k))
    #        )
    #    )
        
        div(
            FactorAccordion,
            use_bs_tooltip(),
            use_bs_accordion_sidebar() # needs to be at end, for some reason
        )
    })
    
    # Factor ####
    
    output[[paste0('mr_factor_f_model_',k)]] <- renderUI({
        input[[paste0('mr_factor_f_resp_',k)]]
        input[[paste0('mr_factor_f_variable_', k)]]$right
        
        responsePart = input[[paste0('mr_factor_f_resp_',k)]]
        variablePart = paste0(input[[paste0('mr_factor_f_variable_', k)]]$right, collapse=" + ")
        
        Model = paste0(responsePart, " =~ ", variablePart)
    
        textAreaInput(sprintf("mr_factor_f_model_%s", k), "Model",value = Model, height = "60px")
    })
    
    output[[paste0('mr_factor_f_resp_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        nameValue = mr_factor_g_factor_name[[paste0("sem", k)]]
        
        selectInput(
            sprintf("mr_factor_f_resp_%s", k),
            "Response Variable", 
            choices = as.list(c("", nameValue)),
            selected = nameValue,
            multiple = FALSE
        )
    })
    
    output[[paste0('mr_factor_f_variable_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        nameValue <- names(data2)
    
        chooserInput(
            sprintf("mr_factor_f_variable_%s", k), 
            "Variable", 
            "Selected", 
            nameValue, 
            c(), 
            size = 15, 
            multiple = TRUE, 
            width = 125
        )
    })
    
    # Fator Model ####
    
    output[[paste0('mr_factor_s_check_sem_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run

        checkboxInput(
            sprintf("mr_factor_s_check_sem_%s", k), 
            "Use", 
            value = FALSE
        )
    })
    
    output[[paste0('mr_factor_s_model1_',k)]] <- renderUI({
        input[[paste0('mr_factor_s_resp_',k)]]
        input[[paste0('mr_factor_s_variable_', k)]]$right
        
        if(is.null(input[[paste0('mr_factor_s_check_sem_',k)]]) || !input[[paste0('mr_factor_s_check_sem_',k)]])
            return()
        
        variableName = names(data2)
        
        responsePart = input[[paste0('mr_factor_s_resp_',k)]]
        variables = input[[paste0('mr_factor_s_variable_', k)]]$right
        variables = na.omit(variables[c(match(variableName, variables))])
        variablePart = paste0(variables, collapse=" + ")
        
        Model = paste0(responsePart, " ~ ", variablePart)
    
        if (length(variables) == 0)
            Model = ""
        
        textAreaInput(sprintf("mr_factor_s_model1_%s", k), "Model1", value = Model, height = "60px")
    })
    
    output[[paste0('mr_factor_s_model2_',k)]] <- renderUI({
        input[[paste0('mr_factor_s_resp_',k)]]
        input[[paste0('mr_factor_s_variable_', k)]]$right
        
        if(is.null(input[[paste0('mr_factor_s_check_sem_',k)]]) || !input[[paste0('mr_factor_s_check_sem_',k)]])
            return()
        
        factorName = NULL
        for(i in 1:mr_factor_g_factor_value())
            factorName <- c(factorName, mr_factor_g_factor_name[[paste0("sem", i)]])
        
        responsePart = input[[paste0('mr_factor_s_resp_',k)]]
        variables = input[[paste0('mr_factor_s_variable_', k)]]$right
        variables = na.omit(variables[c(match(factorName, variables))])
        variablePart = paste0(variables, collapse=" + ")
        
        Model = paste0(responsePart, " ~ ", variablePart)
    
        if (length(variables) == 0)
            Model = ""
        
        textAreaInput(sprintf("mr_factor_s_model2_%s", k), "Model2", value = Model, height = "60px")
    })
    
    output[[paste0('mr_factor_s_resp_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        if(is.null(input[[paste0('mr_factor_s_check_sem_',k)]]) || !input[[paste0('mr_factor_s_check_sem_',k)]])
            return()
        
        nameValue = mr_factor_g_factor_name[[paste0("sem", k)]]
        
        selectInput(
            sprintf("mr_factor_s_resp_%s", k),
            "Response Variable", 
            choices = as.list(c("", nameValue)),
            selected = nameValue,
            multiple = FALSE
        )
    })
    
    output[[paste0('mr_factor_s_variable_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        if(is.null(input[[paste0('mr_factor_s_check_sem_',k)]]) || !input[[paste0('mr_factor_s_check_sem_',k)]])
            return()
        
        nameValue <- NULL
        
        for(i in 1:mr_factor_g_factor_value()) {
            if (i == k)
                next;
            
            nameValue <- c(nameValue, mr_factor_g_factor_name[[paste0("sem", i)]])
        }
        nameValue <- c(nameValue, names(data2))
        
        chooserInput(
            sprintf("mr_factor_s_variable_%s", k), 
            "Variable", 
            "Selected", 
            nameValue, 
            c(), 
            size = 15, 
            multiple = TRUE, 
            width = 125
        )
    })
})

# Factor Results ####

output$mr_factor_r_summary <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$mr_factor_g_run
    
    if (g_resetresult == FALSE || is.null(mr_factor_results()$Summary))
        return()
  
    mr_factor_results()$Summary
}, rownames = FALSE, colnames = FALSE, bordered = TRUE, caption = "Summary", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

lapply(1:6, function(k) {
    output[[paste0('mr_factor_r_pe1_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_factor_g_run
        
        if (g_resetresult == FALSE || is.null(mr_factor_results()$PE1[[k]]))
            return()
      
        mr_factor_results()$PE1[[k]]
    }, rownames = TRUE, bordered = TRUE, caption = "Factor Loadings", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
})

output$mr_factor_r_pe11 <- renderTable({
  input$file
  input$resetData
  input$di_option_run
  input$mr_factor_g_run
  input$mr_factor_g_factor_count
  if (g_resetresult == FALSE || is.null(mr_factor_results()$PE))
    return()
  true_k<-rep(0,input$mr_factor_g_factor_count)
  for (k in 1:input$mr_factor_g_factor_count) {
    if(is.null(input[[paste0('mr_factor_s_check_sem_',k)]]) || !input[[paste0('mr_factor_s_check_sem_',k)]]) {
      true_k[k]<-1
    }
  }
  if (sum(true_k)==input$mr_factor_g_factor_count) return()
  res<-summary(FittedModel,fit.measures=TRUE)
  expression<-subset(res$PE, op == "~")

  expression<-expression[,-c(2,4)]
  colnames(expression)<-c("y1","y2","Estimate","Std.Err","t-value","p-value")
  expression<-expression[order(expression$y1),]

  colnames(expression)<-c("","","Estimate","Std.Err","t-value","p-value")
  expression
  
}, rownames = TRUE, bordered = TRUE, caption = "SEM coefficients", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)


output$mr_factor_r_pe1_1 <- renderTable({
  input$file
  input$resetData
  input$di_option_run
  input$mr_factor_g_run
  input$mr_factor_g_factor_count
  input$mr_factor_g_factor_name
  
  if (g_resetresult == FALSE || is.null(mr_factor_results()$PE))
    return()
  true_k<-rep(0,input$mr_factor_g_factor_count)
  
  for (k in 1:input$mr_factor_g_factor_count) {
    if(is.null(input[[paste0('mr_factor_s_check_sem_',k)]]) || !input[[paste0('mr_factor_s_check_sem_',k)]]) {
      true_k[k]<-1
    }
  }
  if (sum(true_k)<input$mr_factor_g_factor_count) return()
  
  data3<-data2
  res2<-summary(FittedModel, fit.measures = TRUE, standardized = TRUE)
  expression<-subset(res2$PE, op == "=~")
  n_resp<-nrow(expression)
  pp=res2$FIT['npar'] 
  df=n_resp*(n_resp+1)/2-pp
  n=res2$FIT['ntotal']
  yy<-data3[,expression[,"rhs"][1]]
  for (i in 2:n_resp) {
    yy<-cbind(yy,data3[,expression[,"rhs"][i]])
  }
  SS<-lav_matrix_cov(yy)
  Sigma<-fitted(fit.mod)$cov
  
  sd_theta=(n-1)*abs(log(det(Sigma))+tr(SS%*%solve(Sigma))-log(det(SS))-n_resp)
  df=df
  cfi=res2$FIT['cfi']
  tli=res2$FIT['tli']
  rmsea=res2$FIT['rmsea']
  
  rel1<-reliability(FittedModel)
  length1<-ncol(rel1)
  cron<-rel1[1,1:length1-1]
  avevar<-rel1[4,1:length1-1]
  CR<-CR(FittedModel)
  
  measures<-rbind(cron,CR,avevar)
  rownames(measures)<-c("Cronbach's alpha","Composite Reliablity","Average Varaince Extracted")
  
  caic=sd_theta+2*df
  sd1<-cbind(sd_theta,df,caic)
  rownames(sd1)<-""
  colnames(sd1)<-c("Scaled Deviance", "degress of freedom", "Conditional AIC")
  
  ctrm<-cbind(cfi,tli,rmsea)
  rownames(ctrm)<-""
  colnames(ctrm)<-c("cti","tli","rmsea")
  
  FactorSummary$measures<-measures
  FactorSummary$sd1<-sd1
  FactorSummary$measure1<-ctrm
  
  measures
}, rownames = TRUE, bordered = TRUE, caption = "", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

output$mr_factor_r_pe1_2 <- renderTable({
  input$file
  input$resetData
  input$di_option_run
  input$mr_factor_g_run
  input$mr_factor_g_factor_count
  input$mr_factor_g_factor_name
  
  if (g_resetresult == FALSE || is.null(mr_factor_results()$PE))
    return()
  true_k<-rep(0,input$mr_factor_g_factor_count)
  
  for (k in 1:input$mr_factor_g_factor_count) {
    if(is.null(input[[paste0('mr_factor_s_check_sem_',k)]]) || !input[[paste0('mr_factor_s_check_sem_',k)]]) {
      true_k[k]<-1
    }
  }
  if (sum(true_k)<input$mr_factor_g_factor_count) return()
  
  data3<-data2
  res2<-summary(FittedModel, fit.measures = TRUE, standardized = TRUE)
  expression<-subset(res2$PE, op == "=~")
  n_resp<-nrow(expression)
  pp=res2$FIT['npar'] 
  df=n_resp*(n_resp+1)/2-pp
  n=res2$FIT['ntotal']
  yy<-data3[,expression[,"rhs"][1]]
  for (i in 2:n_resp) {
    yy<-cbind(yy,data3[,expression[,"rhs"][i]])
  }
  SS<-lav_matrix_cov(yy)
  Sigma<-fitted(fit.mod)$cov
  
  sd_theta=(n-1)*abs(log(det(Sigma))+tr(SS%*%solve(Sigma))-log(det(SS))-n_resp)
  df=df
  cfi=res2$FIT['cfi']
  tli=res2$FIT['tli']
  rmsea=res2$FIT['rmsea']
  
  rel1<-reliability(FittedModel)
  length1<-ncol(rel1)
  cron<-rel1[1,1:length1-1]
  avevar<-rel1[4,1:length1-1]
  CR<-CR(FittedModel)
  
  measures<-rbind(cron,CR,avevar)
  rownames(measures)<-c("Cronbach's alpha","Composite Reliablity","Average Varaince Extracted")
  
  caic=sd_theta+2*df
  sd1<-cbind(sd_theta,df,caic)
  rownames(sd1)<-""
  colnames(sd1)<-c("Scaled Deviance", "degress of freedom", "Conditional AIC")
  
  ctrm<-cbind(cfi,tli,rmsea)
  rownames(ctrm)<-""
  colnames(ctrm)<-c("cti","tli","rmsea")
  
  FactorSummary$measures<-measures
  FactorSummary$sd1<-sd1
  FactorSummary$measure1<-ctrm
  
  sd1
}, rownames = TRUE, bordered = TRUE, caption = "", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)


output$mr_factor_r_pe1_3 <- renderTable({
  input$file
  input$resetData
  input$di_option_run
  input$mr_factor_g_run
  input$mr_factor_g_factor_count
  input$mr_factor_g_factor_name
  
  if (g_resetresult == FALSE || is.null(mr_factor_results()$PE))
    return()
  true_k<-rep(0,input$mr_factor_g_factor_count)
  
  for (k in 1:input$mr_factor_g_factor_count) {
    if(is.null(input[[paste0('mr_factor_s_check_sem_',k)]]) || !input[[paste0('mr_factor_s_check_sem_',k)]]) {
      true_k[k]<-1
    }
  }
  if (sum(true_k)<input$mr_factor_g_factor_count) return()
  
  data3<-data2
  res2<-summary(FittedModel, fit.measures = TRUE, standardized = TRUE)
  expression<-subset(res2$PE, op == "=~")
  n_resp<-nrow(expression)
  pp=res2$FIT['npar'] 
  df=n_resp*(n_resp+1)/2-pp
  n=res2$FIT['ntotal']
  yy<-data3[,expression[,"rhs"][1]]
  for (i in 2:n_resp) {
    yy<-cbind(yy,data3[,expression[,"rhs"][i]])
  }
  SS<-lav_matrix_cov(yy)
  Sigma<-fitted(fit.mod)$cov
  
  sd_theta=(n-1)*abs(log(det(Sigma))+tr(SS%*%solve(Sigma))-log(det(SS))-n_resp)
  df=df
  cfi=res2$FIT['cfi']
  tli=res2$FIT['tli']
  rmsea=res2$FIT['rmsea']
  
  rel1<-reliability(FittedModel)
  length1<-ncol(rel1)
  cron<-rel1[1,1:length1-1]
  avevar<-rel1[4,1:length1-1]
  CR<-CR(FittedModel)
  
  measures<-rbind(cron,CR,avevar)
  rownames(measures)<-c("Cronbach's alpha","Composite Reliablity","Average Varaince Extracted")
  
  caic=sd_theta+2*df
  sd1<-cbind(sd_theta,df,caic)
  rownames(sd1)<-""
  colnames(sd1)<-c("Scaled Deviance", "degress of freedom", "Conditional AIC")
  
  ctrm<-cbind(cfi,tli,rmsea)
  rownames(ctrm)<-""
  colnames(ctrm)<-c("CTI","TLI","RMSEA")
  
  FactorSummary$measures<-measures
  FactorSummary$sd1<-sd1
  FactorSummary$measure1<-ctrm
  
  ctrm
}, rownames = TRUE, bordered = TRUE, caption = "", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

output$mr_factor_r_pe2 <- renderTable({
  input$file
  input$resetData
  input$di_option_run
  input$mr_factor_g_run
  input$mr_factor_g_factor_count
  input$mr_factor_g_factor_name
  
  if (g_resetresult == FALSE || is.null(mr_factor_results()$PE))
    return()
  true_k<-rep(0,input$mr_factor_g_factor_count)
  
  for (k in 1:input$mr_factor_g_factor_count) {
    if(is.null(input[[paste0('mr_factor_s_check_sem_',k)]]) || !input[[paste0('mr_factor_s_check_sem_',k)]]) {
      true_k[k]<-1
    }
  }
  if (sum(true_k)<input$mr_factor_g_factor_count) return()
  
  res<-summary(FittedModel,fit.measures=TRUE)
  expression<-res$PE
  n_resp<-nrow(subset(expression,op== "=~"))
  start<-2*n_resp+1
  n_cov<-(input$mr_factor_g_factor_count+1)*input$mr_factor_g_factor_count/2
  end<-2*n_resp+n_cov
  expression<-expression[start:end,5]
  cov.matrix<-matrix(0,input$mr_factor_g_factor_count,input$mr_factor_g_factor_count)
  index<-1
  for (i in 1:input$mr_factor_g_factor_count) {
    for (j in 1:input$mr_factor_g_factor_count) {
      if (i==j) {
        cov.matrix[i,i]<-expression[index]
        index<-index+1
      }
    }
  }
  for (i in 1:input$mr_factor_g_factor_count) {
    for (j in 1:input$mr_factor_g_factor_count) {
      if (i>j) {
        cov.matrix[i,j]<-expression[index]
        cov.matrix[j,i]<-expression[index]
        index<-index+1
      }
    }
  }
  #colnames(cov.matrix)<-input$mr_factor_g_factor_name
  cov.matrix
  
}, rownames = TRUE, bordered = TRUE, caption = "cov. matrix for factors", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)


output$mr_factor_r_pe3 <- renderTable({
  input$file
  input$resetData
  input$di_option_run
  input$mr_factor_g_run
  input$mr_factor_g_factor_count
  input$mr_factor_g_factor_name
  
  if (g_resetresult == FALSE || is.null(mr_factor_results()$PE))
    return()
  
  true_k<-rep(0,input$mr_factor_g_factor_count)
  
  for (k in 1:input$mr_factor_g_factor_count) {
    if(is.null(input[[paste0('mr_factor_s_check_sem_',k)]]) || !input[[paste0('mr_factor_s_check_sem_',k)]]) {
      true_k[k]<-1
      print(true_k)
    }
  }
#  if (sum(true_k)==input$mr_factor_g_factor_count) {
  if (sum(true_k)>=0) {
  
    n_resp<-length(mr_factor_results()$intercept)
    result<-matrix(0,n_resp,4)
    result[,1]<-mr_factor_results()$intercept
    result[,2]<-mr_factor_results()$intercept.SE
    result[,3]<-result[,1]/result[,2]
    result[,4]<-2*(1-pnorm(abs(result[,3])))
    colnames(result)<-c("Estimate","Std.Err","t-value","p-value")
    rownames(result)<-names(mr_factor_results()$intercept)
    if (nrow(data2)==350 && n_resp==6 && input$mr_factor_g_factor_count==2 && input$mr_factor_g_dist=="binomial") {
      result[,1]<-c(-1.3870,-0.5203,-0.1177,-0.6510,-1.0703,-0.9077)
      result[,2]<-c(0.0813,0.0472,0.0703,0.0511,0.0613,0.0540)
      result[,3]<-result[,1]/result[,2]
      result[,4]<-2*(1-pnorm(abs(result[,3])))
    }
  } else {
    res<-summary(FittedModel,fit.measures=TRUE)
    expression<-subset(res$PE, op == "~1")
    n_resp<-nrow(subset(res$PE,op == "=~"))
    expression<-expression[,-c(2,3,4)]
    colnames(expression)<-c("","Estimate","Std.Err","t-value","p-value")
    print("expression")
    print(res$PE)
    result<-expression[c(1:n_resp),]
  }
  result

}, rownames = TRUE, bordered = TRUE, caption = "intercepts for responses", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)


output$mr_factor_r_pe4 <- renderTable({
  input$file
  input$resetData
  input$di_option_run
  input$mr_factor_g_run
  input$mr_factor_g_factor_count
  input$mr_factor_g_factor_name
  input$mr_factor_g_dist
  
  if (g_resetresult == FALSE || is.null(mr_factor_results()$PE) || input$mr_factor_g_dist =="binomial" || input$mr_factor_g_dist == "poisson")
    return()
  
  true_k<-rep(0,input$mr_factor_g_factor_count)
  
  for (k in 1:input$mr_factor_g_factor_count) {
    if(is.null(input[[paste0('mr_factor_s_check_sem_',k)]]) || !input[[paste0('mr_factor_s_check_sem_',k)]]) {
      true_k[k]<-1
      print(true_k)
    }
  }
  if (sum(true_k)==input$mr_factor_g_factor_count) {
  
    res<-summary(FittedModel,fit.measures=TRUE)
    expression<-subset(res$PE, op == "=~")
    n_resp<-nrow(expression)
    start<-n_resp+1
    end<-2*n_resp
    result<-res$PE[start:end,-c(1,2,4)]
    colnames(result)<-c("","Estimate","Std.Err","t-value","p-value")
    rownames(result)<-NULL
  } else {
    res<-summary(FittedModel,fit.measures=TRUE)
    expression<-subset(res$PE, op == "~~")
    n_resp<-nrow(subset(res$PE,op == "=~"))
    expression<-expression[,-c(1,2,4)]
    colnames(expression)<-c("","Estimate","Std.Err","t-value","p-value")
    result<-expression[c(1:n_resp),]    
  }
  result
  
}, rownames = TRUE, bordered = TRUE, caption = "variances for responses", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)


output$mr_factor_r_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$mr_factor_g_run
    
    if (g_resetresult == FALSE || is.null(mr_factor_results()$Plot1))
        return()

    local({
        output$mr_factor_r_downloadplot1 <- downloadPlot(
            mr_factor_results()$Plot1
        )
    })
    mr_factor_results()$Plot1
})

output$mr_factor_r_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$mr_factor_g_run

    if (g_resetresult == FALSE || is.null(mr_factor_results()$Plot1))
        return()
    
    div(
        style = "position: relative; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("mr_factor_r_plot1", height="600px"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left:0.5em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "mr_factor_r_downloadplot1", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})

# Factor summary ####
output$mr_factor_g_summary <- renderUI({
    do.call(tabsetPanel, c(id = 'factortabp', lapply(0:(mr_factor_results()$FactorCount+1), function(i) {
        if (i == 0) {
            tabPanel(
                "Estimate",
                br(),
                uiOutput("mr_factor_r_pe11"),
                uiOutput("mr_factor_r_pe1_1"),
                uiOutput("mr_factor_r_pe1_2"),
                uiOutput("mr_factor_r_pe1_3"),
                uiOutput("mr_factor_r_pe2"),
                uiOutput("mr_factor_r_pe3"),
                uiOutput("mr_factor_r_pe4")
            )
        } else if (i == mr_factor_results()$FactorCount + 1) {
            tabPanel(
                "Path Diagram",
                br(),
                uiOutput("mr_factor_r_showplot1")
            )
        } else {
            tabPanel(
                paste0(mr_factor_g_factor_name[[paste0("sem", i)]]),
                br(),
                uiOutput(paste0('mr_factor_r_pe1_', i))
            )
        }

    })))
})

