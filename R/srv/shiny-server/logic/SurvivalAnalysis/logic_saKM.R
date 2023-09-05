# K-M Run ####

observeEvent(input$sa_km_run, {
    g_resetresult <<- TRUE
})

sa_km_results <- eventReactive(input$sa_km_run, {
    # ¦§ Warning & Notify ####
    if (input$sa_km_survivaltime == "" || input$sa_km_indicator == "") {
        showNotification("Please select survival time and indicator", type="warning")
        return()
    }
    
    if (input$sa_km_checkgroup == TRUE && input$sa_km_groupvariable == "") {
        showNotification("Please choose group variable", type="warning")
        return()
    }
    
    if (input$sa_km_check_initial == TRUE && input$sa_km_initialtime == "") {
        showNotification("Please choose initial", type="warning")
        return()
    }
    
    # Start Progressbar
    withProgress(message = 'K-M', style = "notification", value = 0, {
    
    removeTab(inputId = "sa_km_resulttabset", target = "Group K-M Estimates")
    removeTab(inputId = "sa_km_resulttabset", target = "Group K-M Curve")
    
    # ¦§ Title ####
    
    kmTitle <- input$sa_km_survivaltime
    if(input$sa_km_check_initial && !is.null(input$sa_km_initialtime) && input$sa_km_initialtime != "")
        kmTitle <- paste0(input$sa_km_survivaltime, ' - ', input$sa_km_initialtime)
    
    # ¦§ Create a survival object ####
    
    if(input$sa_km_check_initial && !is.null(input$sa_km_initialtime) && input$sa_km_initialtime != "")
        kmsurvival <<- with(
            data2, 
            Surv(
                time = data2[[input$sa_km_initialtime]],
                time2 = data2[[input$sa_km_survivaltime]], 
                event = data2[[input$sa_km_indicator]],
                type = "counting"
            )
        )
    else 
        kmsurvival <<- with(
            data2, 
            Surv(
                time = rep(0,length(data2[[input$sa_km_survivaltime]])),
                time2 = data2[[input$sa_km_survivaltime]],
                event = data2[[input$sa_km_indicator]],
                type = "counting"
            )
        )
    
    incProgress(0.5, detail = paste("Loading..."))
    
    kmresults <- survfit(kmsurvival ~ 1, data=data2, conf.type = "plain")
    
    if (input$sa_km_checkgroup == FALSE) {
        return(list(
            Title = kmTitle,
            Kmsurvfit = kmresults  
        ))
    }
    
    # ¦¦ Group Results ####
    # groupFormula <- as.formula(paste0("kmsurvival ~ ",input$sa_km_groupvariable))
    kmgroupresults <- survfit(
        formula = as.formula(paste("kmsurvival ~",input$sa_km_groupvariable)),
        data = data2, conf.type = "plain")
    
    kmgroupresults$call$formula <- as.formula(paste("kmsurvival ~",input$sa_km_groupvariable))
    
    logrankresults = NULL
    if(input$sa_km_checkgroup && input$sa_km_checklogrank) {
        if(input$sa_km_check_initial && !is.null(input$sa_km_initialtime) && input$sa_km_initialtime != "") {
            logrankresults = survdiff(
                formula = as.formula(
                    paste0(
                        "Surv(", 
                        input$sa_km_survivaltime, 
                        "-", 
                        input$sa_km_initialtime, 
                        ", ", 
                        input$sa_km_indicator, 
                        ") ~ ", 
                        input$sa_km_groupvariable
                    )
                ), 
                data = data2, 
                rho = 0
            ) 
            logrankresults$call$formula <- as.formula(
                paste0(
                    "Surv(", 
                    input$sa_km_survivaltime, 
                    "-", 
                    input$sa_km_initialtime, 
                    ", ", 
                    input$sa_km_indicator, 
                    ") ~ ", 
                    input$sa_km_groupvariable
                )
            )
        }
        else  {
            logrankresults = survdiff(
                formula = as.formula(
                    paste0(
                        "Surv(", 
                        input$sa_km_survivaltime, 
                        ", ", 
                        input$sa_km_indicator, 
                        ") ~ ", 
                        input$sa_km_groupvariable
                    )
                ), 
                data = data2, 
                rho = 0
            )  
            logrankresults$call$formula <- as.formula(
                paste0(
                    "Surv(", 
                    input$sa_km_survivaltime, 
                    ", ", 
                    input$sa_km_indicator, 
                    ") ~ ", 
                    input$sa_km_groupvariable
                )
            )
        }
    }
    
    groupvariable <- as.factor(data2[[input$sa_km_groupvariable]])
    grouplength <- length(levels(groupvariable))
    
    groupdata <- list()
    
    for(i in 1:grouplength) {
        groupvariable == levels(groupvariable)[i]
        tempdata <-data2[groupvariable == levels(groupvariable)[i],]

        kmgroupsurvival = NULL
        if(input$sa_km_check_initial && !is.null(input$sa_km_initialtime) && input$sa_km_initialtime != "")
            kmgroupsurvival = with(
                tempdata, 
                Surv(
                    time = tempdata[[input$sa_km_initialtime]],
                    time2 = tempdata[[input$sa_km_survivaltime]], 
                    event = tempdata[[input$sa_km_indicator]],
                    type = "counting"
                )
            )
        else 
            kmgroupsurvival = with(
                tempdata, 
                Surv(
                    time = rep(0,length(tempdata[[input$sa_km_survivaltime]])),
                    time2 = tempdata[[input$sa_km_survivaltime]],
                    event = tempdata[[input$sa_km_indicator]],
                    type = "counting"
                )
            )
        
        tempgroupresults <- survfit(kmgroupsurvival ~ 1, data=tempdata, conf.type = "plain")
        
        groupdata[[i]]<- tempgroupresults
        
    }
    names(groupdata) <- levels(groupvariable)
    
    appendTab(
        inputId = "sa_km_resulttabset",
        tabPanel(
            "Group K-M Estimates",
            br(),
            uiOutput("sa_km_select_group"),
            tableOutput("sa_km_groupkmestimates")
        )
    )
    
    appendTab(
        inputId = "sa_km_resulttabset",
        tabPanel(
            "Group K-M Curve",
            br(),
            uiOutput("sa_km_showgroupkmcurveplot1")
        )
    )
    
    incProgress(0.5, detail = paste("Loading..."))
    
    return(list(
        Title = kmTitle,
        Kmsurvfit = kmresults, 
        Kmgroupsurvfit = kmgroupresults,
        Kmlogrank = logrankresults,
        Groupdata = groupdata
    ))
    
    })
})

# K-M Components ####

output$sa_km_survivaltime <- renderUI({
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
            "sa_km_survivaltime",
            "Survival time",
            choices = as.list(c("", nameValue)), 
            selected = NULL, 
            multiple = FALSE
        ),
        bsTooltip(
            "sa_km_survivaltime", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$sa_km_check_initial<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    checkboxInput("sa_km_check_initial", "Initial Time", value = FALSE)
})

output$sa_km_initialtime <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$sa_km_check_initial) || !input$sa_km_check_initial)
        return()
        
    nameValue = names(select_if(data2, is.numeric))
    
    div(
        selectInput(
            "sa_km_initialtime",
            "Initial time",
            choices = as.list(c("", nameValue)), 
            selected = NULL, 
            multiple = FALSE
        ),
        bsTooltip(
            "sa_km_initialtime", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$sa_km_indicator <- renderUI({
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
            "sa_km_indicator",
            "Indicator",
            choices = as.list(c("", nameValue)),
            selected = NULL,
            multiple = FALSE
        ),
        bsTooltip(
            "sa_km_indicator", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$sa_km_checkgroup <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput(
        inputId = 'sa_km_checkgroup',
        label = 'Use Group',
        value = FALSE
    )
})

output$sa_km_groupvariable <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if (is.null(input$sa_km_checkgroup) || input$sa_km_checkgroup == FALSE)
        return()
    
    nameValue <- names(select_if(data2, is.factor))
    nameValue <- c(nameValue, names(select_if(data2, is.integer)))
    
    if (!is.null(input$sa_km_survivaltime) && is.integer(data2[[input$sa_km_survivaltime]])) {
        nameValue <- nameValue[-which(nameValue==input$sa_km_survivaltime)]
    }
    
    selectInput(
        "sa_km_groupvariable",
        "Groups",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = FALSE
    )
})

output$sa_km_checklogrank <- renderUI({
    input$file
    input$resetData
    input$di_option_run

    if(is.null(input$sa_km_checkgroup) || !input$sa_km_checkgroup)
        return()
    
    checkboxInput(
        inputId = 'sa_km_checklogrank',
        label = 'Log-Rank Test',
        value = FALSE
    )
})

# unfinished object
# output$sa_km_checkgehan <- renderUI({
#     input$file
#     input$resetData
#     
#     if(!is.null(input$sa_km_checkgroup) && input$sa_km_checkgroup)
#         checkboxInput(
#             inputId = 'sa_km_checkgehan',
#             label = 'Log Rank Test',
#             value = FALSE
#         )
# })
# 
# output$sa_km_checktaroneware <- renderUI({
#     input$file
#     input$resetData
#     
#     if(!is.null(input$sa_km_checkgroup) && input$sa_km_checkgroup)
#         checkboxInput(
#             inputId = 'sa_km_checktaroneware',
#             label = 'Log Rank Test',
#             value = FALSE
#         )
# })

# K-M Results 1 ####

output$sa_km_title <- renderUI({
    input$sa_km_run
    
    titleOutput = NULL
    if (g_resetresult == FALSE || is.null(sa_km_results()$Title))
        titleOutput = ""
    else
        titleOutput = sa_km_results()$Title #how to say?
    
    return(paste('Survival Time : ', titleOutput))
})

output$sa_km_kmresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_km_run
    
    if (g_resetresult == FALSE || is.null(sa_km_results()$Kmsurvfit))
        return()
    
    #if we use tempTable, it is transpose form
    tempTable<-summary(sa_km_results()$Kmsurvfit)$table
    resultTable<-Matrix(data=tempTable, nrow=1)
    colnames(resultTable) <- names(tempTable)
    
    return(as.matrix(resultTable))
}, bordered = TRUE, caption = "Data Summary", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 2)

output$sa_km_groupkmresults <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_km_run

    if (g_resetresult == FALSE || is.null(sa_km_results()$Kmgroupsurvfit))
        return()

    return(summary(sa_km_results()$Kmgroupsurvfit)$table)
}, bordered = TRUE, rownames = TRUE, caption = "Group Data Summary", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 2)

output$sa_km_logrankresults1 <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_km_run
    
    if (g_resetresult == FALSE || is.null(sa_km_results()$Kmgroupsurvfit) 
        || is.null(sa_km_results()$Kmlogrank))
        return()
    
    N <- sa_km_results()$Kmlogrank$n
    O <- sa_km_results()$Kmlogrank$obs
    E <- sa_km_results()$Kmlogrank$exp
    V <- sa_km_results()$Kmlogrank$var[1]
    logrankresultmatrix1 <- matrix(data = c(N, O, E, (O-E)^2/E,(O-E)^2/V ), ncol=5)
    colnames(logrankresultmatrix1) <- c("N", "Observed", "Expected", "(O-E)^2/E", "(O-E)^2/V")
    rownames(logrankresultmatrix1) <- names(sa_km_results()$Kmgroupsurvfit$strata)
    
    return(logrankresultmatrix1)
}, bordered = TRUE, caption = "Log-Rank Test Table", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 3)

output$sa_km_logrankresults2 <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_km_run
    
    if (g_resetresult == FALSE || is.null(sa_km_results()$Kmgroupsurvfit) 
        || is.null(sa_km_results()$Kmlogrank))
        return()
    
    chidf <- length(sa_km_results()$Kmgroupsurvfit$strata) - 1
    chisq <- sa_km_results()$Kmlogrank$chisq
    p_value <- 1 - pchisq(chisq, df = chidf)
    logrankresultmatrix2 <- matrix(data=c(chisq, chidf, p_value), nrow=1)
    colnames(logrankresultmatrix2) <- c("Chisq", "Degrees of freedom", "p-value")
    
    
    return(logrankresultmatrix2)
}, bordered = TRUE, caption = "Log-Rank Test Results", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 3)

# K-M Results 2 (Kaplan-Meier) ####

output$sa_km_kmestimates <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_km_run
    
    if (g_resetresult == FALSE || is.null(sa_km_results()$Kmsurvfit))
        return()
    
    resultEstimates<-matrix(
        data = c(
            sa_km_results()$Kmsurvfit$time, 
            sa_km_results()$Kmsurvfit$n.risk, 
            sa_km_results()$Kmsurvfit$n.event, 
            sa_km_results()$Kmsurvfit$n.censor, 
            sa_km_results()$Kmsurvfit$surv, 
            sa_km_results()$Kmsurvfit$std.err, 
            sa_km_results()$Kmsurvfit$lower, 
            sa_km_results()$Kmsurvfit$upper), ncol=8)
    colnames(resultEstimates)<-c("time", "n.risk", "n.event", "censored", "survival", "std.err", "lower 95% CI", "upper 95% CI")
    
    return(resultEstimates)
}, bordered = TRUE, caption = "Kaplan-Meier Estimates", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 3)

output$sa_km_kmcurveplot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$sa_km_run
    
    if (g_resetresult == FALSE || is.null(sa_km_results()$Kmsurvfit))
        return()
    
    local({
        output$sa_km_downloadkmcurveplot1 <- downloadSurvPlot({
            ggsurvplot(sa_km_results()$Kmsurvfit, conf.int = F, risk.table = TRUE)
        })  
    })
    
    ggsurvplot(sa_km_results()$Kmsurvfit, conf.int = F, risk.table = TRUE)
})

output$sa_km_showkmcurveplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$sa_km_run

    if(g_resetresult == FALSE || is.null(sa_km_results()$Kmsurvfit))
        return()
    
    div(
        style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("sa_km_kmcurveplot1", height="700px"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left:0.5em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "sa_km_downloadkmcurveplot1", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})

# K-M Results 3 (Kaplan-Meier Group) ####

output$sa_km_select_group <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$sa_km_run
    
    if (g_resetresult == FALSE || is.null(sa_km_results()$Kmgroupsurvfit))
        return()
    
    nameValue<-names(sa_km_results()$Groupdata)

    selectInput(
        "sa_km_select_group",
        "Choose Group",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = FALSE
    )
})

output$sa_km_groupkmestimates <-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_km_run
    input$sa_km_select_group
    
    if (g_resetresult == FALSE || is.null(sa_km_results()$Kmgroupsurvfit)
        || is.null(input$sa_km_select_group) || input$sa_km_select_group == ""
        || is.null(sa_km_results()$Groupdata) )
        return()
    
    groupdata <- sa_km_results()$Groupdata[[input$sa_km_select_group]]
    
    groupresultEstimates<-matrix(
        data = c(
            groupdata$time, 
            groupdata$n.risk, 
            groupdata$n.event, 
            groupdata$n.censor, 
            groupdata$surv, 
            groupdata$std.err, 
            groupdata$lower, 
            groupdata$upper), ncol=8)
    colnames(groupresultEstimates)<-c("time", "n.risk", "n.event", "censored", "survival", "std.err", "lower 95% CI", "upper 95% CI")
    
    return(groupresultEstimates)
}, bordered = TRUE, caption = "Group K-M Estimates", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 3)


output$sa_km_groupkmcurveplot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$sa_km_run
    
    if (g_resetresult == FALSE || is.null(sa_km_results()$Kmgroupsurvfit))
        return()
    
    local({
        output$sa_km_downloadgroupkmcurveplot1 <- downloadSurvPlot({
            ggsurvplot(sa_km_results()$Kmgroupsurvfit, conf.int = F, risk.table = TRUE)
        })   
    })
    
    ggsurvplot(sa_km_results()$Kmgroupsurvfit, conf.int = F, risk.table = TRUE)
})

output$sa_km_showgroupkmcurveplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$sa_km_run

    if(g_resetresult == FALSE || is.null(sa_km_results()$Kmgroupsurvfit))
        return()
    
    div(
        style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("sa_km_groupkmcurveplot1", height="700px"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left:0.5em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "sa_km_downloadgroupkmcurveplot1", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})

   
