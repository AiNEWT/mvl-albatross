# SEM Slider ####

output$mr_sem_g_factorslider <- renderUI({
    sliderInput(
        inputId = sprintf("mr_sem_g_factor_count"),
        label = "Number of Factors",
        value = 2,
        min = 2,
        max = 6
    )
})

mr_sem_g_factor_value <- reactive({
    if (is.null(input$mr_sem_g_factor_count))
        return(2)
    else
        return(input$mr_sem_g_factor_count)
})

# SEM Run ####

observeEvent(input$mr_sem_g_run, {
    g_resetresult <<- TRUE
})

mr_sem_results <- eventReactive(input$mr_sem_g_run, {
    withProgress(message = 'SEM', style = "notification", value = 0, {
    removeTab(inputId = "mr_sem_resulttabset", target = "Model Comparison")
    # Increase Progressbar
    incProgress(0.5, detail = paste("Loading..."))
        
    SemModel <- paste0(input$mr_sem_g_model,"\n",input$mr_sem_g_model_2)
    
    # print(SemModel)
    
#    FittedModel <<- sem(
#        SemModel,
#        data = data2, 
#        std.lv=TRUE, 
#        missing='direct',
#        estimator="ML",
#        se='robust.huber.white',
#        test='Yuan.Bentler' 
#    )
    nfactor=input$mr_sem_g_factor_count
    FittedModel <<- hsem(SemModel,data=data2)
    res<-summary(FittedModel,fit.measures=TRUE)
    expression<-subset(res$PE, op == "=~")
    data3<-data2
    dist=input$mr_sem_g_dist
    n_resp<-nrow(expression)
    intercept<-rep(0,n_resp)
    intercept.SE<-rep(0,n_resp)
    muhat3<-matrix(0,nrow(data3),n_resp)
    yyy<-matrix(0,nrow(data3),n_resp)
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
        muhat3[,i]<-zz[[1]]
        yyy[,i]<-yy[[1]]
    }
    FittedModel <<- hsem(SemModel,data=data3)
    SemSummary <<- summary(FittedModel, fit.measures=TRUE)
    nresp<-n_resp
    yy<-data3[,expression[,"rhs"][c(1:nresp)]]
    if (!is.null(SemSummary$FIT)) {
        summaryMatrix <- matrix(rep(NA, 90), ncol = 5)
        
        for (i in 1:9) {
            summaryMatrix[i*2-1, 1:5] <- names(SemSummary$FIT[(1+5*(i-1)):(5*i)])
            summaryMatrix[i*2, 1:5] <- round(SemSummary$FIT[(1+5*(i-1)):(5*i)],3)
        }
        SemSummary$Summary <<- summaryMatrix    
    }
    
    SemSummary$res<-res
    # SEM Plot ####
    SemSummary$Plot1 <- ggsem(FittedModel)
    SemSummary$Plot01 <- ggcorrplot::ggcorrplot(cor(yy))
    SemSummary$intercept <- intercept
    SemSummary$intercept.SE <- intercept.SE
    SemSummary$dist<-dist
    
    vv<-predict(FittedModel)
    samplesize<-nrow(vv)
    SemSummary$Plot02 <- ggcorrplot::ggcorrplot(cor(vv))
    expression<-subset(res$PE, op == "=~")
    nana<-ifelse(is.na(yy),0,1)
    yy1<-matrix(0,samplesize,nresp)
    ee<-matrix(0,samplesize,nresp)
    mu1<-matrix(0,samplesize,nresp)
    x1<-rep(0,samplesize)
    
    row1<-1
    for (i in 1:nrow(data3)) {
      nana1<-sum(nana[i,])
      if (nana1==nresp) {
        yy1[row1,1:nresp]<-as.numeric(yy[i,1:nresp])
        if (nrow(data3)==714) x1[row1]<-data3$sex[i]
        row1<-row1+1
      } 
    }

    mutotal<-rep(0,samplesize*nresp)
    eetotal<-rep(0,samplesize*nresp)
    SemSummary$residual <- eetotal
    SemSummary$muhat<- mutotal
    SemSummary$residual<- c(vv)
    SemSummary$muhat<-muhat3
    expression1=expression[1,1]
    j=1
    expression_var<-subset(res$PE, op == "~~")   
      for (i in 1:nresp) {
        expression_comp=expression[i,1]
        if (expression_comp==expression1) {
          j=j
        } else j=j+1
        expression1=expression_comp
        end<-samplesize*i
        start<-samplesize*(i-1)+1
#        ee[,i]<-(yyy[,i]-mean(yyy[,i])-res$PE[i,5]*vv[,j])/sqrt(var(yyy[,i]-mean(yyy[,i])-res$PE[i,5]*vv[,j]))
        ee[,i]<-(yyy[,i]-mean(yyy[,i])-res$PE[i,5]*vv[,j])/sqrt(expression_var[i,5])
        mu1[,i]<-res$PE[i,5]*vv[,j]
        if (dist=="binomial") {
          mu1[,i]<-1/(1+exp(-1*mu1[,i]))
          ee[,i]<-(yyy[,i]-mu1[,i])/sqrt(var(yyy[,i]-mu1[,i]))
          absee<-abs(ee[,i])
          ee[,i]<-(absee-mean(absee))/sd(absee)
          mu1[,i]<-log(mu1[1:samplesize,i]/(1-mu1[1:samplesize,i]))
        }
        mutotal[start:end]<-mu1[1:samplesize,i]
        eetotal[start:end]<-ee[1:samplesize,i]
      }
    SemSummary$residual <- eetotal
    SemSummary$muhat<- mutotal    
    n_sem=0
    mu2total=0
    ee2total=0
    ee2=0
    mu2=0
    expression_sem<-subset(res$PE, op == "~")
    expression_sem<-expression_sem[order(expression_sem$'lhs'),]
    if (nrow(expression_sem)==0) n_sem=0
    if (nrow(expression_sem)!=0) {
      init=expression_sem[1,1]
      n_sem=1
      for (i in 1:nrow(expression_sem)) {
        if (expression_sem[i,1]!=init) {
          n_sem=n_sem+1
          init=expression_sem[i,1]
        }
      }
      mu2total<-rep(0,samplesize*n_sem)
      ee2total<-rep(0,samplesize*n_sem)
      ee2<-matrix(0,samplesize,n_sem)
      mu2<-matrix(0,samplesize,n_sem)
    }
    #for (i in 1:5) {
    #  end<-samplesize*i
    #  start<-samplesize*(i-1)+1
    #  mu2total[start:end]<-mu2[1:samplesize,i]
    #  ee2total[start:end]<-ee2[1:samplesize,i]
    #}
    if (n_sem>=1) {
        i<-1
        end<-samplesize*i
        start<-samplesize*(i-1)+1
        xx<-rep(0,samplesize)
        init=expression_sem[1,1]
        for (j in 1:nrow(expression_sem)) {
          if (expression_sem[j,1]==init) {
            if (sum(regexpr(expression_sem[j,3],colnames(vv))) == -1*ncol(vv)) {
                row1=1
                for (iii in 1:nrow(data3)) {
                  nana1<-sum(nana[iii,])
                  if (nana1==nresp) {
                    xx[row1]=data3[iii,expression_sem[j,3]]
                    row1=row1+1
                  }
                }
                xx1<-rep(0,samplesize)
                for (iii in 1:samplesize) {
                  xx1[iii]=xx[[iii]][1]
                }
              mu2[,i]=mu2[,i]+expression_sem[j,5]*xx1[1:samplesize]
            } else {
              mu2[,i]=mu2[,i]+expression_sem[j,5]*vv[,expression_sem[j,3]]
            }
            ee2[,i]=(vv[,i]-mu2[,i])/sqrt(var(vv[,i]-mu2[,i]))
            mu2total[start:end]<-mu2[1:samplesize,i]
            ee2total[start:end]<-ee2[1:samplesize,i]
          } else {
            init=expression_sem[j,1]
            i<-i+1
            if (sum(regexpr(expression_sem[j,3],colnames(vv))) == -1*ncol(vv)) {
              row1=1
              for (iii in 1:nrow(data3)) {
                nana1<-sum(nana[iii,])
                if (nana1==nresp) {
                  xx[row1]<-data3[iii,expression_sem[j,3]][1,1]
                  row1<-row1+1
                }
              } 
              xx1<-rep(0,samplesize)
              for (iii in 1:samplesize) {
                xx1[iii]=xx[[iii]][1]
              }
              mu2[,i]=mu2[,i]+expression_sem[j,5]*xx1
            } else {
              mu2[,i]=mu2[,i]+expression_sem[j,5]*vv[,expression_sem[j,3]]
            }
            ee2[,i]=(vv[,i]-mu2[,i])/sqrt(var(vv[,i]-mu2[,i]))
            mu2total[start:end]<-mu2[1:samplesize,i]
            ee2total[start:end]<-ee2[1:samplesize,i]
          }
        }
    }
    
    n_sem<-nfactor
    mu2total<-rep(0,samplesize*n_sem)
    ee2total<-rep(0,samplesize*n_sem)
    ee2<-matrix(0,samplesize,n_sem)
    mu2<-matrix(0,samplesize,n_sem)    
    for (i in 1:n_sem) {
      end<-samplesize*i
      start<-samplesize*(i-1)+1
      ee2total[start:end]<-vv[,i]
    }
    SemSummary$samplesize<-samplesize
    SemSummary$nresp<-nresp
    SemSummary$nfactor<-nfactor
    SemSummary$n_sem<-n_sem

    SemSummary$residual2 <- ee2total
    SemSummary$muhat2 <- mu2total
    
    
    # Comparison Model ####
    #  semComparisonmodel1 = NULL
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
    Sigma<-fitted(FittedModel)$cov[1:n_resp,1:n_resp]
    
    sd_theta=(n-1)*abs(log(det(Sigma))+tr(SS%*%solve(Sigma))-log(det(SS))-n_resp)
    caic=sd_theta+2*df
    cfi=res2$FIT['cfi']
    tli=res2$FIT['tli']
    rmsea=res2$FIT['rmsea']
    
#    rel1<-reliability(FittedModel)
#    length1<-ncol(rel1)
#    cron<-rel1[1,1:length1-1]
#    avevar<-rel1[4,1:length1-1]
#    CR<-CR(FittedModel)
    
#    measures<-rbind(cron,CR,avevar)
#    rownames(measures)<-c("Cronbach's alpha","Composite Reliablity","Average Varaince Extracted")
    
    sd1<-cbind(sd_theta,df,caic)
    rownames(sd1)<-""
    colnames(sd1)<-c("Scaled Deviance", "df", "cAIC")
    
    ctrm<-cbind(cfi,tli,rmsea)
    rownames(ctrm)<-""
    colnames(ctrm)<-c("CTI","TLI","RMSEA")
    
    selection<-cbind(sd_theta,df,caic,cfi,tli,rmsea)
    colnames(selection)<-c("Scaled Deviance", "df", "cAIC","CTI","TLI","RMSEA")
    rownames(selection)<-""
    
    semComparisonmodel2 = NULL
    semComparisonmodel1 <- matrix(SemModel,ncol=1)
    ml=cfi*100*n
    rl=tli*100*n
    semComparisonmodel2 <- matrix(
      c(
        ml,rl,sd_theta,df,caic,cfi,tli,rmsea
      ),
      ncol = 8
    )
    semComparisonmodel1 <- rbind(g_mr_sem_comparisonmodel_1, semComparisonmodel1)
    semComparisonmodel2 <- rbind(g_mr_sem_comparisonmodel_2, semComparisonmodel2)
    colnames(semComparisonmodel1) <- "Model"
    colnames(semComparisonmodel2) <- c("-2ML", "-2RL", "Scaled Deviance", "df", "cAIC","CTI","TLI","RMSEA")  
    g_mr_sem_comparisonmodel_1 <<- semComparisonmodel1
    g_mr_sem_comparisonmodel_2 <<- semComparisonmodel2
    
    SemSummary$mr_semComparisonmodel1=semComparisonmodel1
    SemSummary$mr_semComparisonmodel2=semComparisonmodel2
    
    if (!is.null(input$mr_sem_g_comparison) && input$mr_sem_g_comparison) {
      appendTab(inputId = "mr_sem_resulttabset",
                tabPanel(
                  "Model Comparison",
                  br(),
                  h4("History of Models"),
                  tableOutput("mr_semComparisonmodel1"),
                  tableOutput("mr_semComparisonmodel2")
                )
      )
    }
    # SemSummary$Summary1 <- SemSummary$FIT[1:2]
    # SemSummary$Summary1 <- SemSummary$FIT[1:2]
    # SemSummary$Summary1 <- SemSummary$FIT[1:2]
    # HSEM
    # Factor_list <- vector(mode = "list", length = input$mr_sem_g_factor_count)
    # Sem_list <- vector(mode = "list", length = input$mr_sem_g_factor_count)
    # 
    # for(k in 1:input$mr_sem_g_factor_count) {
    #     FactorModel <- hsemmodeling(
    #         Model = "factor",
    #         LinPred = as.formula(input[[paste0("mr_sem_f_model_", k)]])
    #     )
    #     
    #     SemModel <- hsemmodeling(
    #         Model = "sem",
    #         LinPred = as.formula(input[[paste0("mr_sem_s_model_", k)]])
    #     )
    #     
    #     Factor_list[[k]] <- FactorModel
    #     Sem_list[[k]] <- SemModel
    # }
    # 
    # fittedModel <- hsemfit(
    #     RespDist = "gaussian",
    #     DataMain = data2,
    #     factorModel = Factor_list,
    #     semModel = Sem_list
    # )
    
    # Set Progressbar
    setProgress(1, detail = "Finish")
    return(SemSummary)
    
    }) # End Progressbar
})



output$mr_semComparisonmodel1<-renderTable({
  input$file
  input$resetData
  input$di_option_run
  input$mr_sem_g_run
  
  if (is.null(mr_sem_results()$mr_semComparisonmodel1))
    return()
  
  mr_sem_results()$mr_semComparisonmodel1
}, rownames = TRUE, bordered = TRUE, caption = "", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_semComparisonmodel2<-renderTable({
  input$file
  input$resetData
  input$di_option_run
  input$mr_sem_g_run
  
  if (is.null(mr_sem_results()$mr_semComparisonmodel2))
    return()
  
  mr_sem_results()$mr_semComparisonmodel2
}, rownames = TRUE, bordered = TRUE, caption = "", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)


# SEM Factor Names ####
# reactiveValues object for storing current data set.
mr_sem_g_factor_name <- reactiveValues(
    sem1 = "factor1", 
    sem2 = "factor2", 
    sem3 = "factor3", 
    sem4 = "factor4", 
    sem5 = "factor5", 
    sem6 = "factor6"
)

output$mr_sem_g_names <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run

    factorCount = mr_sem_g_factor_value()
    
    return(
        do.call(div, c(lapply(1:mr_sem_g_factor_value(), function(i) {
            textInput(
                paste0("mr_sem_g_factor_name_",i), 
                paste0("Type Name Factor",i),
                placeholder = paste0("Factor",i)
            )
        })))
    )
})

# Return the UI for a modal dialog with data selection input. If 'failed' is
# TRUE, then display a message that the previous value was invalid.
mr_sem_g_modal <- function(failed = FALSE) {
    modalDialog(
        title = "Change Factor Names",
        easyClose = TRUE,
        uiOutput("mr_sem_g_names"),
        footer = tagList(
            modalButton("Cancel"),
            actionButton("ok", "OK")
        )
    )
}

# Show modal when button is clicked.
observeEvent(input$mr_sem_g_changefactorname, {
    showModal(mr_sem_g_modal())
})

# When OK button is pressed, attempt to load the data set. If successful,
# remove the modal. If not show another modal, but this time with a failure
# message.
observeEvent(input$ok, {
    # Check that data object exists and is data frame.
    for (i in 1:mr_sem_g_factor_value()) {
        if (input[[paste0("mr_sem_g_factor_name_", i)]] == "")
            mr_sem_g_factor_name[[paste0("sem", i)]] <- paste0("factor", i)
        else 
            mr_sem_g_factor_name[[paste0("sem", i)]] <- input[[paste0("mr_sem_g_factor_name_", i)]]
    }
    removeModal()
})

# SEM Model ####

output$mr_sem_g_model <- renderUI({
    Model1 = NULL
    Model2 = NULL
    Model3 = NULL
    for (k in 1:mr_sem_g_factor_value()) {
        Model1 <- c(Model1, input[[paste0('mr_sem_f_model_',k)]])
        Model2 <- c(Model2, input[[paste0('mr_sem_s_model1_',k)]])
        Model3 <- c(Model3, input[[paste0('mr_sem_s_model2_',k)]])
    }
    
    comment1 = "# Measurement Model "
    comment2 = "# Structural Model "
    
    Model1 = paste0(Model1, collapse="\n")
    Model2 = paste0(Model2, collapse="\n")
    Model3 = paste0(Model3, collapse="\n")

    Model = paste0( Model1)
    
    Height = paste0(mr_sem_g_factor_value() * 40, "px")
    textAreaInput("mr_sem_g_model", " Measurement Model", value = Model, height = Height)
})


output$mr_sem_g_model_2 <- renderUI({
  Model1 = NULL
  Model2 = NULL
  Model3 = NULL
  for (k in 1:mr_sem_g_factor_value()) {
    Model1 <- c(Model1, input[[paste0('mr_sem_f_model_',k)]])
    Model2 <- c(Model2, input[[paste0('mr_sem_s_model1_',k)]])
    Model3 <- c(Model3, input[[paste0('mr_sem_s_model2_',k)]])
  }
  
  comment1 = "# Measurement Model "
  comment2 = "# Structural Model "
  
  Model1 = paste0(Model1, collapse="\n")
  Model2 = paste0(Model2, collapse="\n")
  Model3 = paste0(Model3, collapse="\n")
  
  Model = paste0(Model2, "\n", Model3)
  
  Height = paste0(mr_sem_g_factor_value() * 40, "px")
  textAreaInput("mr_sem_g_model_2", "Structural Model", value = Model, height = Height)
})

output$mr_sem_g_dist<-renderUI({
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
    "mr_sem_g_dist",
    "Distribution", 
    choices = dist, 
    multiple = FALSE
  )
})

output$mr_sem_g_link<-renderUI({
  input$file
  input$resetData
  input$di_option_run
  
  link = NULL
  selection = NULL
  
  if(is.null(input$mr_sem_g_dist))
    return()
  
  if (input$mr_sem_g_dist == "gaussian"){
    selection ="identity"
    link = c("identity"="identity", "log"="log", "inverse"="inverse")
  } else if (input$mr_sem_g_dist == "binomial"){
    selection="logit"
    link = c("logit"="logit", "probit"="probit", "cauchit"="cauchit", "log"="log", "cloglog"="cloglog", "identity"="identity")
  } else if (input$mr_sem_g_dist == "poisson"){
    selection="log"
    link = c("log"="log", "identity"="identity", "sqrt"="sqrt")
  } else if (input$mr_sem_g_dist == "Gamma"){
    selection="log"
    link = c("inverse"="inverse", "identity"="identity", "log"="log")
  }
  selectInput("mr_sem_g_link", "Link Function", choices = link, selected = selection, multiple = FALSE)
})


output$mr_sem_g_comparison<-renderUI({
  input$file
  input$resetData
  input$di_option_run
  
  checkboxInput("mr_sem_g_comparison", "Model Comparison", value = FALSE)
})


# SEM Tab Panel ####
output$mr_sem_g_tabpanel <- renderUI({
    do.call(tabsetPanel, c(id = 'semtabp', lapply(1:mr_sem_g_factor_value(), function(i) {
        tabPanel(
            title = paste0('Factor', i),
            br(),
            uiOutput(paste0('mr_sem_r_accordion_', i))
        )
    })))
})

# SEM Components ####
lapply(1:6, function(k) { 
    
    # Accordion ####
    
    output[[paste0('mr_sem_r_accordion_',k)]]<-renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        semAccordion <- bs_accordion_sidebar(
            id = paste0("mr_sem_r_accordion", k),
            spec_side = c(width = 3, offset = 0),
            spec_main = c(width = 9, offset = 0)    
        )
        
        semAccordion <- semAccordion %>%
            bs_append(
                title_side = "Measurement Model",
                content_side = NULL,
                content_main = div(
                    h3(strong("Measurement Model")),
                    uiOutput(paste0("mr_sem_f_model_", k)),
                    uiOutput(paste0("mr_sem_f_resp_", k)),
                    uiOutput(paste0("mr_sem_f_variable_", k))
                )
            )
        
        semAccordion <- semAccordion %>%
        bs_append(
            title_side = "Structural Model",
            content_side = uiOutput(paste0("mr_sem_s_check_sem_", k)),
            content_main = div(
                h3(strong("Structural Model")),
                uiOutput(paste0("mr_sem_s_model1_", k)),
                uiOutput(paste0("mr_sem_s_model2_", k)),
                uiOutput(paste0("mr_sem_s_resp_", k)),
                uiOutput(paste0("mr_sem_s_variable_", k))
            )
        )
        
        div(
            semAccordion,
            use_bs_tooltip(),
            use_bs_accordion_sidebar() # needs to be at end, for some reason
        )
    })
    
    # Factor ####
    
    output[[paste0('mr_sem_f_model_',k)]] <- renderUI({
        input[[paste0('mr_sem_f_resp_',k)]]
        input[[paste0('mr_sem_f_variable_', k)]]$right
        
        responsePart = input[[paste0('mr_sem_f_resp_',k)]]
        variablePart = paste0(input[[paste0('mr_sem_f_variable_', k)]]$right, collapse=" + ")
        
        Model = paste0(responsePart, " =~ ", variablePart)
    
        textAreaInput(sprintf("mr_sem_f_model_%s", k), "Model",value = Model, height = "60px")
    })
    
    output[[paste0('mr_sem_f_resp_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        nameValue = mr_sem_g_factor_name[[paste0("sem", k)]]
        
        selectInput(
            sprintf("mr_sem_f_resp_%s", k),
            "Response Variable", 
            choices = as.list(c("", nameValue)),
            selected = nameValue,
            multiple = FALSE
        )
    })
    
    output[[paste0('mr_sem_f_variable_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        nameValue <- names(data2)
    
        chooserInput(
            sprintf("mr_sem_f_variable_%s", k), 
            "Variable", 
            "Selected", 
            nameValue, 
            c(), 
            size = 15, 
            multiple = TRUE, 
            width = 125
        )
    })
    
    # SEM ####
    
    output[[paste0('mr_sem_s_check_sem_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run

        checkboxInput(
            sprintf("mr_sem_s_check_sem_%s", k), 
            "Use", 
            value = FALSE
        )
    })
    
    output[[paste0('mr_sem_s_model1_',k)]] <- renderUI({
        input[[paste0('mr_sem_s_resp_',k)]]
        input[[paste0('mr_sem_s_variable_', k)]]$right
        
        if(is.null(input[[paste0('mr_sem_s_check_sem_',k)]]) || !input[[paste0('mr_sem_s_check_sem_',k)]])
            return()
        
        variableName = names(data2)
        
        responsePart = input[[paste0('mr_sem_s_resp_',k)]]
        variables = input[[paste0('mr_sem_s_variable_', k)]]$right
        variables = na.omit(variables[c(match(variableName, variables))])
        variablePart = paste0(variables, collapse=" + ")
        
        Model = paste0(responsePart, " ~ ", variablePart)
    
        if (length(variables) == 0)
            Model = ""
        
        textAreaInput(sprintf("mr_sem_s_model1_%s", k), "Model1", value = Model, height = "60px")
    })
    
    output[[paste0('mr_sem_s_model2_',k)]] <- renderUI({
        input[[paste0('mr_sem_s_resp_',k)]]
        input[[paste0('mr_sem_s_variable_', k)]]$right
        
        if(is.null(input[[paste0('mr_sem_s_check_sem_',k)]]) || !input[[paste0('mr_sem_s_check_sem_',k)]])
            return()
        
        factorName = NULL
        for(i in 1:mr_sem_g_factor_value())
            factorName <- c(factorName, mr_sem_g_factor_name[[paste0("sem", i)]])
        
        responsePart = input[[paste0('mr_sem_s_resp_',k)]]
        variables = input[[paste0('mr_sem_s_variable_', k)]]$right
        variables = na.omit(variables[c(match(factorName, variables))])
        variablePart = paste0(variables, collapse=" + ")
        
        Model = paste0(responsePart, " ~ ", variablePart)
    
        if (length(variables) == 0)
            Model = ""
        
        textAreaInput(sprintf("mr_sem_s_model2_%s", k), "Model2", value = Model, height = "60px")
    })
    
    output[[paste0('mr_sem_s_resp_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        if(is.null(input[[paste0('mr_sem_s_check_sem_',k)]]) || !input[[paste0('mr_sem_s_check_sem_',k)]])
            return()
        
        nameValue = mr_sem_g_factor_name[[paste0("sem", k)]]
        
        selectInput(
            sprintf("mr_sem_s_resp_%s", k),
            "Response Variable", 
            choices = as.list(c("", nameValue)),
            selected = nameValue,
            multiple = FALSE
        )
    })
    
    output[[paste0('mr_sem_s_variable_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        if(is.null(input[[paste0('mr_sem_s_check_sem_',k)]]) || !input[[paste0('mr_sem_s_check_sem_',k)]])
            return()
        
        nameValue <- NULL
        
        for(i in 1:mr_sem_g_factor_value()) {
            if (i == k)
                next;
            
            nameValue <- c(nameValue, mr_sem_g_factor_name[[paste0("sem", i)]])
        }
        nameValue <- c(nameValue, names(data2))
        
        chooserInput(
            sprintf("mr_sem_s_variable_%s", k), 
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

# SEM Results ####

output$mr_sem_r_summary <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$mr_sem_g_run
    
    if (g_resetresult == FALSE || is.null(mr_sem_results()$Summary))
        return()
  
    mr_sem_results()$Summary
}, rownames = FALSE, colnames = FALSE, bordered = TRUE, caption = "Summary", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)

output$mr_sem_r_pe1 <- renderTable({
  input$file
  input$resetData
  input$di_option_run
  input$mr_sem_g_run
  input$mr_sem_g_factor_count
  
  if (g_resetresult == FALSE || is.null(mr_sem_results()$PE))
    return()
  
  res<-summary(FittedModel,fit.measures=TRUE)
  expression<-subset(res$PE, op == "=~")
  n_resp<-nrow(expression)
  
  expression<-expression[,-c(2,4)]
  colnames(expression)<-c("","","Estimate","Std.Err","t-value","p-value")
  expression
  
}, rownames = TRUE, bordered = TRUE, caption = "factor loadings", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_sem_r_pe1_1 <- renderTable({
  input$file
  input$resetData
  input$di_option_run
  input$mr_sem_g_run
  input$mr_sem_g_factor_count
  if (g_resetresult == FALSE || is.null(mr_sem_results()$PE))
    return()
  true_k<-rep(0,input$mr_sem_g_factor_count)
  for (k in 1:input$mr_sem_g_factor_count) {
    if(is.null(input[[paste0('mr_sem_s_check_sem_',k)]]) || !input[[paste0('mr_sem_s_check_sem_',k)]]) {
      true_k[k]<-1
    }
  }
  if (sum(true_k)==input$mr_sem_g_factor_count) return()
  res<-summary(FittedModel,fit.measures=TRUE)
  expression<-subset(res$PE, op == "~")

  expression<-expression[,-c(2,4)]
  colnames(expression)<-c("y1","y2","Estimate","Std.Err","t-value","p-value")
  expression<-expression[order(expression$y1),]

  colnames(expression)<-c("","","Estimate","Std.Err","t-value","p-value")
  expression
  
}, rownames = TRUE, bordered = TRUE, caption = "SEM coefficients", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_sem_r_pe1_0_1 <- renderTable({
  input$file
  input$resetData
  input$di_option_run
  input$mr_sem_g_run
  input$mr_sem_g_factor_count
  input$mr_sem_g_factor_name
  
  if (g_resetresult == FALSE || is.null(mr_sem_results()$PE))
    return()
  true_k<-rep(0,input$mr_sem_g_factor_count)
  
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
  Sigma<-fitted(FittedModel)$cov[1:n_resp,1:n_resp]
  
  sd_theta=(n-1)*abs(log(det(Sigma))+tr(SS%*%solve(Sigma))-log(det(SS))-n_resp)
  caic=sd_theta+2*df
  cfi=res2$FIT['cfi']
  tli=res2$FIT['tli']
  rmsea=res2$FIT['rmsea']
  
#  rel1<-reliability(FittedModel)
#  length1<-ncol(rel1)
#  cron<-rel1[1,1:length1-1]
#  avevar<-rel1[4,1:length1-1]
#  CR<-CR(FittedModel)
  
#  measures<-rbind(cron,CR,avevar)
#  rownames(measures)<-c("Cronbach's alpha","Composite Reliablity","Average Varaince Extracted")
  
  sd1<-cbind(sd_theta,df,caic)
  rownames(sd1)<-""
  colnames(sd1)<-c("Scaled Deviance", "df", "cAIC")
  
  ctrm<-cbind(cfi,tli,rmsea)
  rownames(ctrm)<-""
  colnames(ctrm)<-c("CTI","TLI","RMSEA")
  
  selection<-cbind(sd_theta,df,caic,cfi,tli,rmsea)
  colnames(selection)<-c("Scaled Deviance", "df", "cAIC","CTI","TLI","RMSEA")
  rownames(selection)<-""

  measures
}, rownames = TRUE, bordered = TRUE, caption = "", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)


output$mr_sem_r_pe1_0_2 <- renderTable({
  input$file
  input$resetData
  input$di_option_run
  input$mr_sem_g_run
  input$mr_sem_g_factor_count
  input$mr_sem_g_factor_name
  
  if (g_resetresult == FALSE || is.null(mr_sem_results()$PE))
    return()
  true_k<-rep(0,input$mr_sem_g_factor_count)
  
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
  Sigma<-fitted(FittedModel)$cov[1:n_resp,1:n_resp]
  
  sd_theta=(n-1)*abs(log(det(Sigma))+tr(SS%*%solve(Sigma))-log(det(SS))-n_resp)
  caic=sd_theta+2*df
  cfi=res2$FIT['cfi']
  tli=res2$FIT['tli']
  rmsea=res2$FIT['rmsea']
  
#  rel1<-reliability(FittedModel)
#  length1<-ncol(rel1)
#  cron<-rel1[1,1:length1-1]
#  avevar<-rel1[4,1:length1-1]
#  CR<-CR(FittedModel)
  
#  measures<-rbind(cron,CR,avevar)
#  rownames(measures)<-c("Cronbach's alpha","Composite Reliablity","Average Varaince Extracted")
  
  sd1<-cbind(sd_theta,df,caic)
  rownames(sd1)<-""
  colnames(sd1)<-c("Scaled Deviance", "df", "cAIC")
  
  ctrm<-cbind(cfi,tli,rmsea)
  rownames(ctrm)<-""
  colnames(ctrm)<-c("CTI","TLI","RMSEA")
  
  selection<-cbind(sd_theta,df,caic,cfi,tli,rmsea)
  colnames(selection)<-c("Scaled Deviance", "df", "cAIC","CTI","TLI","RMSEA")
  rownames(selection)<-""
  
  # Comparison Model ####
  model11=paste0(input$mr_sem_g_model,"\n",input$mr_sem_g_model_2)
  semComparisonmodel1 = NULL
  semComparisonmodel2 = NULL
  semComparisonmodel1 <- matrix(
    c(
      model11
    ),
    ncol = 1
  )
  semComparisonmodel2 <- matrix(
    c(
      sd_theta,df,caic,cfi,tli,rmsea
    ),
    ncol = 6
  )
  semComparisonmodel1 <- rbind(g_mr_sem_comparisonmodel_1, semComparisonmodel1)
  semComparisonmodel2 <- rbind(g_mr_sem_comparisonmodel_2, semComparisonmodel2)
  colnames(semComparisonmodel1) <- ""
  colnames(semComparisonmodel2) <- c("Scaled Deviance", "df", "cAIC","CTI","TLI","RMSEA")  
  g_mr_sem_comparisonmodel_1 <<- semComparisonmodel1
  g_mr_sem_comparisonmodel_2 <<- semComparisonmodel2
  if (!is.null(input$mr_sem_g_comparison)) {
    return(semComparisonmodel2)
  }
}, rownames = TRUE, bordered = TRUE, caption = "", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_sem_r_pe1_1_1 <- renderTable({
  input$file
  input$resetData
  input$di_option_run
  input$mr_sem_g_run
  input$mr_sem_g_factor_count
  input$mr_sem_g_factor_name
  
  if (g_resetresult == FALSE || is.null(mr_sem_results()$PE))
    return()
  true_k<-rep(0,input$mr_sem_g_factor_count)
  
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
  Sigma<-fitted(FittedModel)$cov[1:n_resp,1:n_resp]
  
  sd_theta=(n-1)*abs(log(det(Sigma))+tr(SS%*%solve(Sigma))-log(det(SS))-n_resp)
  df=df
  caic=sd_theta+2*df
  cfi=res2$FIT['cfi']
  tli=res2$FIT['tli']
  rmsea=res2$FIT['rmsea']
  
#  rel1<-reliability(FittedModel)
#  length1<-ncol(rel1)
#  cron<-rel1[1,1:length1-1]
#  avevar<-rel1[4,1:length1-1]
#  CR<-CR(FittedModel)
  
#  measures<-rbind(cron,CR,avevar)
#  rownames(measures)<-c("Cronbach's alpha","Composite Reliablity","Average Varaince Extracted")
  
  sd1<-cbind(sd_theta,df,caic)
  rownames(sd1)<-""
  colnames(sd1)<-c("Scaled Deviance", "df", "cAIC")
  
  ctrm<-cbind(cfi,tli,rmsea)
  rownames(ctrm)<-""
  colnames(ctrm)<-c("CTI","TLI","RMSEA")
  
  measures
  
}, rownames = TRUE, bordered = TRUE, caption = "", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_sem_r_pe1_1_2 <- renderTable({
  input$file
  input$resetData
  input$di_option_run
  input$mr_sem_g_run
  input$mr_sem_g_factor_count
  input$mr_sem_g_factor_name
  
  if (g_resetresult == FALSE || is.null(mr_sem_results()$PE))
    return()
  true_k<-rep(0,input$mr_sem_g_factor_count)
  
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
  Sigma<-fitted(FittedModel)$cov[1:n_resp,1:n_resp]
  
  sd_theta=(n-1)*abs(log(det(Sigma))+tr(SS%*%solve(Sigma))-log(det(SS))-n_resp)
  caic=sd_theta+2*df
  cfi=res2$FIT['cfi']
  tli=res2$FIT['tli']
  rmsea=res2$FIT['rmsea']
  
#  rel1<-reliability(FittedModel)
#  length1<-ncol(rel1)
#  cron<-rel1[1,1:length1-1]
#  avevar<-rel1[4,1:length1-1]
#  CR<-CR(FittedModel)
  
#  measures<-rbind(cron,CR,avevar)
#  rownames(measures)<-c("Cronbach's alpha","Composite Reliablity","Average Varaince Extracted")
  
  sd1<-cbind(sd_theta,df,caic)
  rownames(sd1)<-""
  colnames(sd1)<-c("Scaled Deviance", "df", "cAIC")
  
  ctrm<-cbind(cfi,tli,rmsea)
  rownames(ctrm)<-""
  colnames(ctrm)<-c("CTI","TLI","RMSEA")
  
  selection<-cbind(sd_theta,df,caic,cfi,tli,rmsea)
  colnames(selection)<-c("Scaled Deviance", "df", "cAIC","CTI","TLI","RMSEA")
  rownames(selection)<-""
  
  selection
}, rownames = TRUE, bordered = TRUE, caption = "", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)


output$mr_sem_r_pe1_1_3 <- renderTable({
  input$file
  input$resetData
  input$di_option_run
  input$mr_sem_g_run
  input$mr_sem_g_factor_count
  input$mr_sem_g_factor_name
  
  if (g_resetresult == FALSE || is.null(mr_sem_results()$PE))
    return()
  true_k<-rep(0,input$mr_sem_g_factor_count)
  
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
  Sigma<-fitted(FittedModel)$cov[1:n_resp,1:n_resp]
  
  sd_theta=(n-1)*abs(log(det(Sigma))+tr(SS%*%solve(Sigma))-log(det(SS))-n_resp)
  caic=sd_theta+2*df
  cfi=res2$FIT['cfi']
  tli=res2$FIT['tli']
  rmsea=res2$FIT['rmsea']
  
#  rel1<-reliability(FittedModel)
#  length1<-ncol(rel1)
#  cron<-rel1[1,1:length1-1]
#  avevar<-rel1[4,1:length1-1]
#  CR<-CR(FittedModel)
  
#  measures<-rbind(cron,CR,avevar)
#  rownames(measures)<-c("Cronbach's alpha","Composite Reliablity","Average Varaince Extracted")
  
  sd1<-cbind(sd_theta,df,caic)
  rownames(sd1)<-""
  colnames(sd1)<-c("Scaled Deviance", "df", "cAIC")
  
  ctrm<-cbind(cfi,tli,rmsea)
  rownames(ctrm)<-""
  colnames(ctrm)<-c("CTI","TLI","RMSEA")
  
  ctrm
  
}, rownames = TRUE, bordered = TRUE, caption = "", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$mr_sem_r_pe2 <- renderTable({
  input$file
  input$resetData
  input$di_option_run
  input$mr_sem_g_run
  input$mr_sem_g_factor_count
  input$mr_sem_g_factor_name
  
  if (g_resetresult == FALSE || is.null(mr_sem_results()$PE))
    return()
  true_k<-rep(0,input$mr_sem_g_factor_count)
  
  if (sum(true_k)==input$mr_sem_g_factor_count) {
  
#  for (k in 1:input$mr_sem_g_factor_count) {
#    if(is.null(input[[paste0('mr_sem_s_check_sem_',k)]]) || !input[[paste0('mr_sem_s_check_sem_',k)]]) {
#      true_k[k]<-1
#    }
#  }
#  if (sum(true_k)<input$mr_sem_g_factor_count) return()
  
  res<-summary(FittedModel,fit.measures=TRUE)
  expression<-res$PE
  n_resp<-nrow(subset(expression,op== "=~"))
  start<-2*n_resp+1
  n_cov<-(input$mr_sem_g_factor_count+1)*input$mr_sem_g_factor_count/2
  end<-2*n_resp+n_cov
  expression<-expression[start:end,5]
  cov.matrix<-matrix(0,input$mr_sem_g_factor_count,input$mr_sem_g_factor_count)
  index<-1
  for (i in 1:input$mr_sem_g_factor_count) {
    for (j in 1:input$mr_sem_g_factor_count) {
      if (i==j) {
        cov.matrix[i,i]<-expression[index]
        index<-index+1
      }
    }
  }
  for (i in 1:input$mr_sem_g_factor_count) {
    for (j in 1:input$mr_sem_g_factor_count) {
      if (i>j) {
        cov.matrix[i,j]<-expression[index]
        cov.matrix[j,i]<-expression[index]
        index<-index+1
      }
    }
  }
  #colnames(cov.matrix)<-input$mr_sem_g_factor_name
  } else {
    res<-summary(FittedModel,fit.measures=TRUE)
    expression<-res$PE
    n_resp<-nrow(subset(expression,op== "=~"))
    expression<-subset(expression,op== "~~")
    print("expression")
    print(expression)
    cov.matrix<-matrix(0,input$mr_sem_g_factor_count,input$mr_sem_g_factor_count)
    for (i in 1:input$mr_sem_g_factor_count) cov.matrix[i,i]<-expression[n_resp+i,5]
  }
  cov.matrix
  
}, rownames = TRUE, bordered = TRUE, caption = "cov. matrix for factors", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 4)


output$mr_sem_r_pe3 <- renderTable({
  input$file
  input$resetData
  input$di_option_run
  input$mr_sem_g_run
  input$mr_sem_g_factor_count
  input$mr_sem_g_factor_name
  
  if (g_resetresult == FALSE || is.null(mr_sem_results()$PE))
    return()
  
  true_k<-rep(0,input$mr_sem_g_factor_count)
  
  for (k in 1:input$mr_sem_g_factor_count) {
    if(is.null(input[[paste0('mr_sem_s_check_sem_',k)]]) || !input[[paste0('mr_sem_s_check_sem_',k)]]) {
      true_k[k]<-1
      print(true_k)
    }
  }
#  if (sum(true_k)==input$mr_sem_g_factor_count) {
  if (sum(true_k)>=0) {
  
    n_resp<-length(mr_sem_results()$intercept)
    result<-matrix(0,n_resp,4)
    result[,1]<-mr_sem_results()$intercept
    result[,2]<-mr_sem_results()$intercept.SE
    result[,3]<-result[,1]/result[,2]
    result[,4]<-2*(1-pnorm(abs(result[,3])))
    colnames(result)<-c("Estimate","Std.Err","t-value","p-value")
    rownames(result)<-names(mr_sem_results()$intercept)
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


output$mr_sem_r_pe4 <- renderTable({
  input$file
  input$resetData
  input$di_option_run
  input$mr_sem_g_run
  input$mr_sem_g_factor_count
  input$mr_sem_g_factor_name
  input$mr_sem_g_dist
  
  if (g_resetresult == FALSE || is.null(mr_sem_results()$PE) || input$mr_sem_g_dist =="binomial" || input$mr_sem_g_dist == "poisson")
    return()
  
  true_k<-rep(0,input$mr_sem_g_factor_count)
  
  for (k in 1:input$mr_sem_g_factor_count) {
    if(is.null(input[[paste0('mr_sem_s_check_sem_',k)]]) || !input[[paste0('mr_sem_s_check_sem_',k)]]) {
      true_k[k]<-1
      print(true_k)
    }
  }
  if (sum(true_k)==input$mr_sem_g_factor_count) {
  
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


output$mr_sem_r_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$mr_sem_g_run
    
    if (g_resetresult == FALSE || is.null(mr_sem_results()$Plot1))
        return()

    local({
        output$mr_sem_r_downloadplot1 <- downloadPlot(
            mr_sem_results()$Plot1
        )
    })
    mr_sem_results()$Plot1
})

output$mr_sem_r_plot01 <- renderPlot({
  input$file
  input$resetData
  input$di_option_run
  input$mr_sem_g_run
  
  if (g_resetresult == FALSE || is.null(mr_sem_results()$Plot01))
    return()
  
  local({
    output$mr_sem_r_downloadplot01 <- downloadPlot(
      mr_sem_results()$Plot01
    )
  })
  mr_sem_results()$Plot01
})


output$mr_sem_r_plot02 <- renderPlot({
  input$file
  input$resetData
  input$di_option_run
  input$mr_sem_g_run
  
  if (g_resetresult == FALSE || is.null(mr_sem_results()$Plot02))
    return()
  
  local({
    output$mr_sem_r_downloadplot02 <- downloadPlot(
      mr_sem_results()$Plot02
    )
  })
  mr_sem_results()$Plot02
})


output$mr_sem_r_showplot01 <- renderUI({
  input$file
  input$resetData
  input$di_option_run
  input$mr_sem_g_run
  
  if (g_resetresult == FALSE || is.null(mr_sem_results()$Plot01))
    return()
  div(
    style = "position: relative; border: 1px solid #D3D3D3;",
    withSpinner(
      plotOutput("mr_sem_r_plot01", height="600px"),
      type = 1,
      color = "#2c3e50",
      size = 1.2
    ),
    div(
      style = "position: absolute; left:0.5em; bottom: 0.5em;",
      dropdown(
        downloadButton(outputId = "mr_sem_r_downloadplot01", label = "Download Plot"),
        size = "xs",
        icon = icon("download", class = "opt"),
        up = TRUE
      )
    )
  )    
  
})


output$mr_sem_r_showplot02 <- renderUI({
  input$file
  input$resetData
  input$di_option_run
  input$mr_sem_g_run
  
  if (g_resetresult == FALSE || is.null(mr_sem_results()$Plot02))
    return()
  div(
    style = "position: relative; border: 1px solid #D3D3D3;",
    withSpinner(
      plotOutput("mr_sem_r_plot02", height="600px"),
      type = 1,
      color = "#2c3e50",
      size = 1.2
    ),
    div(
      style = "position: absolute; left:0.5em; bottom: 0.5em;",
      dropdown(
        downloadButton(outputId = "mr_sem_r_downloadplot02", label = "Download Plot"),
        size = "xs",
        icon = icon("download", class = "opt"),
        up = TRUE
      )
    )
  )    
  
})

output$mr_sem_r_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$mr_sem_g_run

    if (g_resetresult == FALSE || is.null(mr_sem_results()$Plot1))
        return()
   
    div(
        style = "position: relative; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("mr_sem_r_plot1", height="600px"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left:0.5em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "mr_sem_r_downloadplot1", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})

output$mr_sem_r_showplot2 <- renderPlot({
  y<-c(StudentResidual<-mr_sem_results()$residual)
  x<-mu<-c(mr_sem_results()$muhat)[1:length(y)]
  par(mfrow=c(2,2))
  fit<- supsmu(x,y)
  plot(x, y, main="Residuals vs Fitted", xlab="scaled fitted values", ylab="Studentized Residual", cex=0.5) #plot data point
  lines(fit$x, fit$y) #plot smooth spline fit
  y<-abs(StudentResidual)
  fit<- supsmu(x,y)
  plot(x, y, main="|Residuals| vs Fitted",xlab="scaled fitted values", ylab="|Studentized Residual|", cex=0.5) #plot data point
  lines(fit$x, fit$y) #plot smooth spline fit
  qqnorm(StudentResidual,main="Normal Probability Plot")
  qqline(StudentResidual) # Normal probability plot
  hist(StudentResidual,main="Histogram of Student Redidual")  
})


output$mr_sem_r_showplot3 <- renderPlot({
  if (mr_sem_results()$n_sem>0) {
  y<-c(StudentResidual<-mr_sem_results()$residual2)
  x<-mu<-c(mr_sem_results()$muhat2)[1:length(y)]
  par(mfrow=c(1,2))
#  fit<- supsmu(x,y)
#  plot(x, y, main="Residuals vs Fitted", xlab="scaled fitted values", ylab="Studentized Residual", cex=0.5) #plot data point
#  lines(fit$x, fit$y) #plot smooth spline fit
#  y<-abs(StudentResidual)
#  fit<- supsmu(x,y)
#  plot(x, y, main="|Residuals| vs Fitted",xlab="scaled fitted values", ylab="|Studentized Residual|", cex=0.5) #plot data point
#  lines(fit$x, fit$y) #plot smooth spline fit
  qqnorm(StudentResidual,main="Normal Probability Plot")
  qqline(StudentResidual) # Normal probability plot
  hist(StudentResidual,main="Histogram of Student Redidual")  
  }
})

output$mr_sem_var1_select<-renderUI({
  nresp<-mr_sem_results()$nresp
  res<-mr_sem_results()$res
  expression<-subset(res$PE, op == "=~")
  data3<-data2
  #    muhat<-matrix(0,nrow(data3),n_resp)
  res3<-c("Overall",colnames(data3[,expression[,"rhs"][1:nresp]]))
  selectInput("mr_sem_var1_select","",choices =as.list(res3),selected = "Overall", multiple = FALSE)
})

output$mr_sem_r_showplot4 <- renderPlot({
  input$mr_sem_var1_select
  nresp<-mr_sem_results()$nresp
  samplesize<-mr_sem_results()$samplesize
  res<-mr_sem_results()$res
  expression<-subset(res$PE, op == "=~")
  data3<-data2
  #    muhat<-matrix(0,nrow(data3),n_resp)
  res3<-colnames(data3[,expression[,"rhs"][1:nresp]])
  
  res<-mr_sem_results()$res
  expression<-subset(res$PE, op == "=~")
  data3<-data2
  #    muhat<-matrix(0,nrow(data3),n_resp)
  res3<-colnames(data3[,expression[,"rhs"][1:nresp]])
  location<-0
  for (i in 1:nresp) {
    if (res3[i]==input$mr_sem_var1_select) location=i
  }
  if (location>=1) {
    start<-(location-1)*samplesize+1
    end<-location*samplesize
  } else {
    start<-1
    end<-length(mr_sem_results()$residual)
  }
  y<-c(StudentResidual<-mr_sem_results()$residual[start:end])
  x<-mu<-c(mr_sem_results()$muhat)[start:end]
  par(mfrow=c(2,2))
  fit<- supsmu(x,y)
  plot(x, y, main="Residuals vs Fitted", xlab="scaled fitted values", ylab="Studentized Residual", cex=0.5) #plot data point
  lines(fit$x, fit$y) #plot smooth spline fit
  y<-abs(StudentResidual)
  fit<- supsmu(x,y)
  plot(x, y, main="|Residuals| vs Fitted",xlab="scaled fitted values", ylab="|Studentized Residual|", cex=0.5) #plot data point
  lines(fit$x, fit$y) #plot smooth spline fit
  qqnorm(StudentResidual,main="Normal Probability Plot")
  qqline(StudentResidual) # Normal probability plot
  hist(StudentResidual,main="Histogram of Student Redidual")  
})


output$mr_sem_var2_select<-renderUI({
  nresp<-mr_sem_results()$nresp
  res<-mr_sem_results()$res
  nfactor<-mr_sem_results()$n_sem
  expression<-subset(res$PE, op == "=~")
  data3<-data2
  #    muhat<-matrix(0,nrow(data3),n_resp)
  res3<-c("Overall","factor1","factor2","factor3","factor4","factor5","factor6")
  end<-nfactor+1
  res3<-res3[1:end]
  selectInput("mr_sem_var2_select","", choices =as.list(res3),selected = "Overall",multiple = FALSE)
})

output$mr_sem_r_showplot5 <- renderPlot({
  n_sem=mr_sem_results()$n_sem
  if (n_sem>0) {
  input$mr_sem_var2_select
  nresp<-mr_sem_results()$nresp
  nfactor<-mr_sem_results()$nresp
  samplesize<-mr_sem_results()$samplesize
  res<-mr_sem_results()$res
  expression<-subset(res$PE, op == "=~")
  data3<-data2
  #    muhat<-matrix(0,nrow(data3),n_resp)
  res3<-colnames(data3[,expression[,"rhs"][1:nresp]])
  
  res<-mr_sem_results()$res
  expression<-subset(res$PE, op == "=~")
  data3<-data2
  #    muhat<-matrix(0,nrow(data3),n_resp)
  res3<-c("factor1","factor2","factor3","factor4","factor5","factor6")
  res3<-res3[1:n_sem]

  location<-0
  for (i in 1:n_sem) {
    if (res3[i]==input$mr_sem_var2_select) location=i
  }
  if (location>=1) {
    start<-(location-1)*samplesize+1
    end<-location*samplesize
  } else {
    start<-1
    end<-length(mr_sem_results()$residual2)
  }
  y<-c(StudentResidual<-mr_sem_results()$residual2[start:end])
  x<-mu<-c(mr_sem_results()$muhat2)[start:end]
  par(mfrow=c(1,2))
#  fit<- supsmu(x,y)
#  plot(x, y, main="Residuals vs Fitted", xlab="scaled fitted values", ylab="Studentized Residual", cex=0.5) #plot data point
#  lines(fit$x, fit$y) #plot smooth spline fit
#  y<-abs(StudentResidual)
#  fit<- supsmu(x,y)
#  plot(x, y, main="|Residuals| vs Fitted",xlab="scaled fitted values", ylab="|Studentized Residual|", cex=0.5) #plot data point
#  lines(fit$x, fit$y) #plot smooth spline fit
  qqnorm(StudentResidual,main="Normal Probability Plot")
  qqline(StudentResidual) # Normal probability plot
  hist(StudentResidual,main="Histogram of Student Redidual")  
  }
})


