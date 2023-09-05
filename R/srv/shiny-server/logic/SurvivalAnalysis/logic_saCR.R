# Competing Risk Run ####

observeEvent(input$sa_cr_g_run, {
    g_resetresult <<- TRUE
})

sa_cr_results <- eventReactive(input$sa_cr_g_run, {
        
    # Start Progressbar
    withProgress(message = 'Competing Risk Model', style = "notification", value = 0.0, {
        
        incProgress(0.5, detail = paste("Loading..."))
        
        fittedModel <<- cmpHL(
            as.formula(input$sa_cr_s1_model),
            as.formula(input$sa_cr_s2_model),
            data2
        )
        rownames(fittedModel$F.Est) <- c("Interesting Event", "Competing Event")
        colnames(fittedModel$F.Est) <- c("Estimate", "Std. Error", "t-value", "p-value")
        # Increase Progressbar
        incProgress(0.5, detail = paste("Loading..."))
        
        return(fittedModel)
    }) # End Progressbar
})

# Competing Risk UI (Accordion) ####

output$sa_cr_r_accordion<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    crAccordion <- bs_accordion_sidebar(
        id = "sa_cr_r_accordion",
        spec_side = c(width = 3, offset = 0),
        spec_main = c(width = 9, offset = 0)    
    )
    
    # ¦§ Interesting Event ####
    
    crAccordion <- crAccordion %>%
        bs_append(
            title_side = "Event 1",
            content_side = NULL,
            content_main = div(
                h3(strong("Interesting Event")),
                uiOutput("sa_cr_s1_survivaltime"),
                uiOutput("sa_cr_s1_check_initial"), 
                uiOutput("sa_cr_s1_initialtime"),
    
                uiOutput("sa_cr_s1_variable"),
                fluidRow(
                    column(
                        8,
                        uiOutput("sa_cr_s1_interaction")
                    ),
                    column(
                        4,
                        uiOutput("sa_cr_s1_interactionappend"),
                        style = "text-align:right; padding:15px"
                    )
                ),
                uiOutput("sa_cr_s1_rand"),
                uiOutput("sa_cr_s1_indicator"),
                uiOutput("sa_cr_s1_status"),
                uiOutput("sa_cr_s1_dist")
            )
        )
    
    # ¦¦ Competing Event ####
    
    crAccordion <- crAccordion %>%
    bs_append(
        title_side = "Event 2",
        content_side = NULL,
        content_main = div(
            h3(strong("Competing Event")),
            uiOutput("sa_cr_s2_survivaltime"),
            uiOutput("sa_cr_s2_check_initial"), 
            uiOutput("sa_cr_s2_initialtime"),

            uiOutput("sa_cr_s2_variable"),
            fluidRow(
                column(
                    8,
                    uiOutput("sa_cr_s2_interaction")
                ),
                column(
                    4,
                    uiOutput("sa_cr_s2_interactionappend"),
                    style = "text-align:right; padding:15px"
                )
            ),
            uiOutput("sa_cr_s2_rand"),
            uiOutput("sa_cr_s2_indicator"),
            uiOutput("sa_cr_s2_status"),
            uiOutput("sa_cr_s2_dist")
        )
    )
  
    div(
        crAccordion,
        use_bs_tooltip(),
        use_bs_accordion_sidebar() # needs to be at end, for some reason
    )
})

# Competing Risk Components ####

output$sa_cr_s1_model<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$sa_cr_s1_survivaltime
    input$sa_cr_s1_variable$right
    input$sa_cr_s1_indicator
    input$sa_cr_s1_initialtime
    
    if(any(is.null(input$sa_cr_s1_survivaltime), is.null(input$sa_cr_s1_variable), is.null(input$sa_cr_s1_indicator)))
        return()
    
    responsePart = NULL
    variablePart = paste0(input$sa_cr_s1_variable$right, collapse="+")
    randomPart = NULL
    modelEquation = NULL
    
    if (!is.null(input$sa_cr_s1_initialtime) && input$sa_cr_s1_initialtime != "") {
        responsePart = paste0("Surv(", input$sa_cr_s1_survivaltime," - ",input$sa_cr_s1_initialtime,", ",input$sa_cr_s1_indicator," == ", input$sa_cr_s1_status, ")")
    } else {
        responsePart = paste0("Surv(", input$sa_cr_s1_survivaltime,", ",input$sa_cr_s1_indicator," == ", input$sa_cr_s1_status, ")")
    }
    
    if(length(input$sa_cr_s1_variable$right) == 0) {
        variablePart = "1"
    }
    
    if(!is.null(input$sa_cr_s1_rand) && length(input$sa_cr_s1_rand) > 0) {
        randomPart = paste0(input$sa_cr_s1_rand, collapse = ") + (1|")
        randomPart = paste0(" + (1|", randomPart, ")")
    }
    
    modelEquation = paste0(responsePart, " ~ ", variablePart, randomPart)
    
    textAreaInput("sa_cr_s1_model","Model for Interesting Event",value = modelEquation, height = "60px")
})

output$sa_cr_s1_survivaltime<-renderUI({
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
            "sa_cr_s1_survivaltime",
            "Survival Time", 
            choices = as.list(c("",nameValue)),
            selected = NULL,
            multiple = FALSE
        ),
        bsTooltip(
            "sa_cr_s1_survivaltime", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$sa_cr_s1_check_initial<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("sa_cr_s1_check_initial", "Initial Time", value = FALSE)
})

output$sa_cr_s1_initialtime<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$sa_cr_s1_check_initial) || !input$sa_cr_s1_check_initial) 
        return()
    
    nameValue = names(select_if(data2, is.numeric))
    
    div(
        selectInput(
            "sa_cr_s1_initialtime",
            "Initial Time",
            choices = as.list(c("",nameValue)),
            multiple = FALSE
        ),
        bsTooltip(
            "sa_cr_s1_initialtime", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$sa_cr_s1_variable<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$sa_cr_s1_interactionappend
    
    nameValue <- c(names(data2), g_sa_cr_s1_interaction)
    
    chooserInput("sa_cr_s1_variable", "Variable", "Selected", nameValue, c(), size = 15, multiple = TRUE, width=150)
})

output$sa_cr_s1_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue <- names(data2)
    
    selectInput(
        "sa_cr_s1_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$sa_cr_s1_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    actionButton("sa_cr_s1_interactionappend", "Append")
})

observeEvent(input$sa_cr_s1_interactionappend, {
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    g_sa_cr_s1_interaction <<- c(g_sa_cr_s1_interaction, paste(input$sa_cr_s1_interaction, collapse=":"))
})

output$sa_cr_s1_rand <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$sa_cr_s1_randinteractionappend
    
    nameValue <- c(names(data2), g_sa_fm_randinteraction)
    
    selectInput(
        "sa_cr_s1_rand",
        "Random Effects",
        choices = as.list(nameValue), 
        multiple = TRUE
    )
})

output$sa_cr_s1_indicator<-renderUI({
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
            "sa_cr_s1_indicator",
            "Censoring Indicator",
            choices = as.list(c("",nameValue)),
            multiple = FALSE
        ),
        bsTooltip(
            "sa_cr_s1_indicator", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$sa_cr_s1_status <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    selectInput(
        "sa_cr_s1_status",
        "Censoring Value",
        choices = as.list(c("0" = "0", "1" = "1", "2" = "2")),
        selected = "1",
        multiple = FALSE
    )
})

output$sa_cr_s1_dist <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    selectInput(
        "sa_cr_s1_dist",
        "Distribution for Random Effect",
        choices = as.list(c("gaussian")), #, "gamma")),
        selected = "FM",
        multiple = FALSE
    )
})

output$sa_cr_s2_model<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$sa_cr_s2_survivaltime
    input$sa_cr_s2_variable$right
    input$sa_cr_s2_indicator
    input$sa_cr_s2_initialtime
    
    if(any(is.null(input$sa_cr_s2_survivaltime), is.null(input$sa_cr_s2_variable), is.null(input$sa_cr_s2_indicator)))
        return(textAreaInput("sa_cr_s2_model","Model for Competing Event",value = "", height = "60px"))
    
    responsePart = NULL
    variablePart = paste0(input$sa_cr_s2_variable$right, collapse="+")
    randomPart = NULL
    modelEquation = NULL
    
    if (!is.null(input$sa_cr_s2_initialtime) && input$sa_cr_s2_initialtime != "") {
        responsePart = paste0("Surv(", input$sa_cr_s2_survivaltime," - ",input$sa_cr_s2_initialtime,", ",input$sa_cr_s2_indicator," == ", input$sa_cr_s2_status, ")")
    } else {
        responsePart = paste0("Surv(", input$sa_cr_s2_survivaltime,", ",input$sa_cr_s2_indicator," == ", input$sa_cr_s2_status, ")")
    }
    
    if(length(input$sa_cr_s2_variable$right) == 0) {
        variablePart = "1"
    }
    
    if(!is.null(input$sa_cr_s2_rand) && length(input$sa_cr_s2_rand) > 0) {
        randomPart = paste0(input$sa_cr_s2_rand, collapse = ") + (1|")
        randomPart = paste0(" + (1|", randomPart, ")")
    }
    
    modelEquation = paste0(responsePart, " ~ ", variablePart, randomPart)
    
    textAreaInput("sa_cr_s2_model", "Model for Competing Event", value = modelEquation, height = "60px")
})

output$sa_cr_s2_survivaltime<-renderUI({
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
            "sa_cr_s2_survivaltime",
            "Survival Time", 
            choices = as.list(c("",nameValue)),
            selected = NULL,
            multiple = FALSE
        ),
        bsTooltip(
            "sa_cr_s2_survivaltime", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$sa_cr_s2_check_initial<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("sa_cr_s2_check_initial", "Initial Time", value = FALSE)
})

output$sa_cr_s2_initialtime<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$sa_cr_s2_check_initial) || !input$sa_cr_s2_check_initial) 
        return()
    
    nameValue = names(select_if(data2, is.numeric))
    
    div(
        selectInput(
            "sa_cr_s2_initialtime",
            "Initial Time",
            choices = as.list(c("",nameValue)),
            multiple = FALSE
        ),
        bsTooltip(
            "sa_cr_s2_initialtime", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$sa_cr_s2_variable<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$sa_cr_s2_interactionappend
    
    nameValue <- c(names(data2), g_sa_cr_s2_interaction)
    
    chooserInput("sa_cr_s2_variable", "Variable", "Selected", nameValue, c(), size = 15, multiple = TRUE, width=150)
})


output$sa_cr_s2_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue <- names(data2)
    
    selectInput(
        "sa_cr_s2_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$sa_cr_s2_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    actionButton("sa_cr_s2_interactionappend", "Append")
})

observeEvent(input$sa_cr_s2_interactionappend, {
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    g_sa_cr_s2_interaction <<- c(g_sa_cr_s2_interaction, paste(input$sa_cr_s2_interaction, collapse=":"))
})

output$sa_cr_s2_rand <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$sa_cr_s2_randinteractionappend
    
    nameValue <- c(names(data2), g_sa_fm_randinteraction)
    
    selectInput(
        "sa_cr_s2_rand",
        "Random Effects",
        choices = as.list(nameValue), 
        multiple = TRUE
    )
})

output$sa_cr_s2_indicator<-renderUI({
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
            "sa_cr_s2_indicator",
            "Censoring Indicator",
            choices = as.list(c("",nameValue)),
            multiple = FALSE
        ),
        bsTooltip(
            "sa_cr_s2_indicator", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$sa_cr_s2_status <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    selectInput(
        "sa_cr_s2_status",
        "Censoring Value",
        choices = as.list(c("0" = "0", "1" = "1", "2" = "2")),
        selected = "1",
        multiple = FALSE
    )
})

output$sa_cr_s2_dist <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    selectInput(
        "sa_cr_s2_dist",
        "Distribution for Random Effect",
        choices = as.list(c("gaussian")),#, "gamma")),
        selected = "FM",
        multiple = FALSE
    )
})

# Competing Risk Results ####

output$sa_cr_r_coefficients <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_cr_g_run
    
    if (g_resetresult == FALSE || is.null(sa_cr_results()$F.Est))
        return()

    sa_cr_results()$F.Est
}, rownames = TRUE, bordered = TRUE, caption = "Coefficients", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$sa_cr_r_estimates <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$sa_cr_g_run
    
    if (g_resetresult == FALSE || is.null(sa_cr_results()$D.Est))
        return()
    
    estimates <- matrix(
        c(
            log(sa_cr_results()$D.Est[,"alpha_h"]), 
            sa_cr_results()$D.Est[,"rho_h"]
        ), ncol = 2
    )
    colnames(estimates) <- c("Logarithm of Variance of Random Effect", "Shared Parameter")
    estimates
}, rownames = FALSE, bordered = TRUE, caption = "Estimates", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

# Competing Risk Plots ####

sa_cr_r_reactiveplot1 <- reactive({
    if(g_resetresult == FALSE || is.null(sa_cr_results()))
        return()
    
    res <- sa_cr_results()
    p <- res$p
    q <- res$q
    var <- diag(res$Hinv)[(p + 1):(p + q)]
    SE <- sqrt(var)
    v_h <- res$v_h
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
      ggtitle("Histogram of Frailty Effects") +
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5))
    
    lb <- sort(lb)
    ub <- sort(ub)
    randCount <- 1:q
        
    # Plot 3
    ciTable <- cbind(lb, ub, randCount)
    plot3 <- ggplot(data.frame(cbind(v_h, randCount)), aes(randCount, v_h)) + 
        geom_point() + 
        geom_line() +
        geom_segment(data = ciTable, aes(x = randCount, y = lb, xend = randCount, yend = ub)) +
        labs(x = "Frailty Number", y = "Estimated Frailty Effects") + 
        ggtitle("95% CI for Frailty Effect") +
        theme_bw() + 
        theme(plot.title = element_text(hjust = 0.5))
    
    return(ggarrange(ggarrange(plot1, plot2, ncol = 2), plot3, nrow = 2) )
})

output$sa_cr_r_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$sa_cr_g_run
    
    local({
        output$sa_cr_r_downloadplot1 <- downloadPlot({
            sa_cr_r_reactiveplot1()
        })     
    })
    
    sa_cr_r_reactiveplot1()
})

output$sa_cr_r_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$sa_cr_g_run

    if(g_resetresult == FALSE || is.null(sa_cr_results()))
        return()
    
    div(
        style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
        withSpinner(
            plotOutput("sa_cr_r_plot1", height="600px"),
            type = 1,
            color = "#2c3e50",
            size = 1.2
        ),
        div(
            style = "position: absolute; left:0.5em; bottom: 0.5em;",
            dropdown(
                downloadButton(outputId = "sa_cr_r_downloadplot1", label = "Download Plot"),
                size = "xs",
                icon = icon("download", class = "opt"),
                up = TRUE
            )
        )
    )
})

 
        