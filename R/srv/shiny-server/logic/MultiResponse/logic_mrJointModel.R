# Joint Model Components ####

output$mr_joint_g_respslider <- renderUI({
    sliderInput(
        inputId = sprintf("mr_joint_g_resp_count"),
        label = "Number of Response Variables",
        value = 1,
        min = 1,
        max = 6
    )
})

output$mr_joint_g_corr_structure<-renderUI({
    selectInput(
        "mr_joint_g_corr_structure",
        "Correlation Structure",
        choices = c("shared", "independent"),
        selected = "shared",
        multiple = FALSE
    )
})

mr_joint_g_resp_value <- reactive({
    if (is.null(input$mr_joint_g_resp_count))
        return(1)
    else
        return(input$mr_joint_g_resp_count)
})

output$mr_joint_g_file <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    fileInput("mr_joint_g_file","Data for Survival Analysis", multiple = TRUE)
})

source('logic/MultiResponse/JointModel/logic_mrJoint_s.R', local = T)
source('logic/MultiResponse/JointModel/logic_mrJoint_m.R', local = T)
source('logic/MultiResponse/JointModel/logic_mrJoint_p.R', local = T)
source('logic/MultiResponse/JointModel/logic_mrJoint_l.R', local = T)

# Joint Model Run ####

observeEvent(input$mr_joint_g_run, {
    g_resetresult <<- TRUE
})

mr_joint_results <- eventReactive(input$mr_joint_g_run, {
    if (input$mr_joint_g_resp_count != 1) {
        showNotification("Joint Model supports only one response.", type="warning")
        return()
    }
        
    # Start Progressbar
    withProgress(message = 'Joint Model', style = "notification", value = 0, {
        # ¦§ Variable Declaration ####
        fittedModel_jm <- NULL
        jm1 <- NULL
        jm2 <- NULL
        jm3 <- NULL
        
        # ¦§ Joint Modeling 1 ####
        jm1 <- jointmodeling(
            Model = "mean",
            RespDist = input$mr_joint_m_dist_1,
            Link = input$mr_joint_m_link_1,
            LinPred = as.formula(input$mr_joint_m_model_1),
            RandDist = "gaussian"
        )
        # ¦§ Joint Modeling 2 ####
        jm2 <- jointmodeling(
            Model = "mean",
            RespDist = input$mr_joint_s1_dist,
            Link = "log",
            LinPred = as.formula(input$mr_joint_s1_model),
            RandDist = "gaussian"
        )
        
        # Increase Progressbar
        incProgress(0.5, detail = paste("Loading..."))
        
        data_surv <- read.table(
            file = input$mr_joint_g_file$datapath,
            sep = input$sep,
            header = input$header,
            stringsAsFactors = TRUE
        )
        
        # ¦¦ Joint Modeling 3 ####
        if (!is.null(input$mr_joint_s2_check_competing) && input$mr_joint_s2_check_competing) {
            jm3 <- jointmodeling(
                Model = "mean",
                RespDist = input$mr_joint_s2_dist,
                Link = "log",
                LinPred = as.formula(input$mr_joint_s2_model),
                RandDist = "gaussian"
            )
            fittedModel_jm <- jmfit3(jm1, jm2, jm3, data2, data_surv = data_surv, Maxiter = 1)
        } else {
            fittedModel_jm <- jmfit(jm1, jm2, data2, data_surv = data_surv, Maxiter = 1, structure = input$mr_joint_g_corr_structure)
        }
        # Increase Progressbar
        incProgress(0.5, detail = paste("Loading..."))
        fittedModel_jm1 <<- fittedModel_jm
        fittedModel_jm$option$RespCount = input$mr_joint_g_resp_count
        return(fittedModel_jm)
    
    }) # End Progressbar
})

# Joint Model tabpanel ####
output$mr_joint_g_tabpanel <- renderUI({
    do.call(tabsetPanel, c(id = 'mr_joint_tabpanel', lapply(0:mr_joint_g_resp_value(), function(i) {
        if (i == 0) {
            tabPanel(
                title = paste0('Survival'),
                br(),
                uiOutput("mr_joint_s1_model"),
                uiOutput("mr_joint_s2_model"),
                uiOutput(paste0("mr_joint_r_accordion_", i))
            )
        } else {
            tabPanel(
                title = paste0('Response', i),
                br(),
                uiOutput(paste0("mr_joint_m_model_", i)),
                uiOutput(paste0("mr_joint_p_model_", i)),
                uiOutput(paste0("mr_joint_l_model_", i)),
                uiOutput(paste0("mr_joint_r_accordion_", i))
            )
        }
    })))
})

# Joint Model Accordion ####

lapply(0:6, function(k) {  
    output[[paste0('mr_joint_r_accordion_',k)]]<-renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$dm_mv_run
        input$dm_ct_run
        input$dm_sr_run
        input$dm_md_run
        
        jointAccordion <- bs_accordion_sidebar(
            id = paste0("mr_joint_r_accordion_", k),
            spec_side = c(width = 3, offset = 0),
            spec_main = c(width = 9, offset = 0)    
        )
        
        if (k == 0) {
            # ¦§ Survival ####
            jointAccordion <- jointAccordion %>%
            bs_append(
                title_side = "Event 1",
                content_side = NULL,
                content_main = div(
                    h3(strong("Interesting Event")),
                    uiOutput("mr_joint_s1_survivaltime"),
                    uiOutput("mr_joint_s1_check_initial"), 
                    uiOutput("mr_joint_s1_initialtime"),
                    uiOutput("mr_joint_s1_variable"),
                    fluidRow(
                        column(
                            8,
                            uiOutput("mr_joint_s1_interaction")
                        ),
                        column(
                            4,
                            uiOutput("mr_joint_s1_interactionappend"),
                            style = "text-align:right; padding:15px"
                        )
                    ),
                    uiOutput("mr_joint_s1_rand"),
                    uiOutput("mr_joint_s1_indicator"),
                    uiOutput("mr_joint_s1_status"),
                    uiOutput("mr_joint_s1_dist")
                )
            )
            
            jointAccordion <- jointAccordion %>%
            bs_append(
                title_side = "Event 2",
                content_side = uiOutput("mr_joint_s2_check_competing"),
                content_main = div(
                    h3(strong("Competing Event")),
                    
                    uiOutput("mr_joint_s2_survivaltime"),
                    uiOutput("mr_joint_s2_check_initial"), 
                    uiOutput("mr_joint_s2_initialtime"),
    
                    uiOutput("mr_joint_s2_variable"),
                    fluidRow(
                        column(
                            8,
                            uiOutput("mr_joint_s2_interaction")
                        ),
                        column(
                            4,
                            uiOutput("mr_joint_s2_interactionappend"),
                            style = "text-align:right; padding:15px"
                        )
                    ),
                    uiOutput("mr_joint_s2_rand"),
                    uiOutput("mr_joint_s2_indicator"),
                    uiOutput("mr_joint_s2_status"),
                    uiOutput("mr_joint_s2_dist")
                )
            )
        } else {
            # ¦§ Mean ####
            jointAccordion <- jointAccordion %>%
            bs_append(
                title_side = "Mean",
                content_side = NULL,
                content_main = div(
                    h3(strong("Model for Mean")),
                    uiOutput(paste0("mr_joint_m_resp_", k)),
                    uiOutput(paste0("mr_joint_m_variable_", k)),
                    fluidRow(
                        column(
                            8,
                            uiOutput(paste0("mr_joint_m_interaction_", k))
                        ),
                        column(
                            4,
                            uiOutput(paste0("mr_joint_m_interactionappend_", k)),
                            style = "text-align:right; padding:15px"
                        )
                    ),
                    uiOutput(paste0("mr_joint_m_rand_", k)),
                    uiOutput(paste0("mr_joint_m_check_randinteraction_", k)),
                    fluidRow(
                        column(
                            8,
                            uiOutput(paste0("mr_joint_m_randinteraction_", k))
                        ),
                        column(
                            4,
                            uiOutput(paste0("mr_joint_m_randinteractionappend_", k)),
                            style = "text-align:right; padding:15px"
                        )
                    ), #new components
                    uiOutput(paste0("mr_joint_m_randfamily_", k)),
                    hr(style = "border-color: #2C3E50;"),
                    uiOutput(paste0("mr_joint_m_dist_", k)),
                    # uiOutput(paste0("mr_joint_m_check_binomd_", k)),
                    # uiOutput(paste0("mr_joint_m_binomd_", k)),
                    uiOutput(paste0("mr_joint_m_link_", k)),
                    uiOutput(paste0("mr_joint_m_check_nointercept_", k)),
                    uiOutput(paste0("mr_joint_m_check_offset_", k)),
                    uiOutput(paste0("mr_joint_m_offset_", k)),
                    uiOutput(paste0("mr_joint_m_factor_", k)),
                    hr(style = "border-color: #2C3E50;"),
                    uiOutput(paste0("mr_joint_m_check_rcodes_", k))
                )
            )
        
            # ¦§ Phi ####
            jointAccordion <- jointAccordion %>%
            bs_append(
                title_side = "Phi",
                content_side = uiOutput(paste0("mr_joint_p_check_phi_", k)),
                content_main = div(
                    h3(strong("Model for phi")),
                    uiOutput(paste0("mr_joint_p_resp_", k)),
                    uiOutput(paste0("mr_joint_p_variable_", k)),
                    fluidRow(
                        column(
                            8,
                            uiOutput(paste0("mr_joint_p_interaction_", k))
                        ),
                        column(
                            4,
                            uiOutput(paste0("mr_joint_p_interactionappend_", k)),
                            style = "text-align:right; padding:15px"
                        )
                    ),
                    uiOutput(paste0("mr_joint_p_rand_", k)),
                    uiOutput(paste0("mr_joint_p_check_randinteraction_", k)),
                    fluidRow(
                        column(
                            8,
                            uiOutput(paste0("mr_joint_p_randinteraction_", k))
                        ),
                        column(
                            4,
                            uiOutput(paste0("mr_joint_p_randinteractionappend_", k)),
                            style = "text-align:right; padding:15px"
                        )
                    ), #new components
                    uiOutput(paste0("mr_joint_p_randfamily_", k)),
                    uiOutput(paste0("mr_joint_p_link_", k)),
                    uiOutput(paste0("mr_joint_p_check_offset_", k)),
                    uiOutput(paste0("mr_joint_p_offset_", k))
                )
            )
        
            # ¦¦ Lambda ####
            jointAccordion <- jointAccordion %>%
            bs_append(
                title_side = "Lambda",
                content_side = uiOutput(paste0("mr_joint_l_check_lambda_", k)),
                content_main = div(
                    h3(strong("Model for Lambda")),
                    uiOutput(paste0("mr_joint_l_resp_", k)),
                    uiOutput(paste0("mr_joint_l_variable_", k)),
                    fluidRow(
                        column(
                            8,
                            uiOutput(paste0("mr_joint_l_interaction_", k))
                        ),
                        column(
                            4,
                            uiOutput(paste0("mr_joint_l_interactionappend_", k)),
                            style = "text-align:right; padding:15px"
                        )
                    ),
                    uiOutput(paste0("mr_joint_l_rand_", k)),
                    uiOutput(paste0("mr_joint_l_check_randinteraction_", k)),
                    fluidRow(
                        column(
                            8,
                            uiOutput(paste0("mr_joint_l_randinteraction_", k))
                        ),
                        column(
                            4,
                            uiOutput(paste0("mr_joint_l_randinteractionappend_", k)),
                            style = "text-align:right; padding:15px"
                        )
                    ), #new components
                    uiOutput(paste0("mr_joint_l_randfamily_", k)),
                    uiOutput(paste0("mr_joint_l_link_", k)),
                    uiOutput(paste0("mr_joint_l_check_offset_", k)),
                    uiOutput(paste0("mr_joint_l_offset_", k))
                )
            )
        }
        
        div(
            jointAccordion,
            use_bs_tooltip(),
            use_bs_accordion_sidebar() # needs to be at end, for some reason
        )
    })
})

# Joint Model output ####
lapply(1:6, function(k) {
    output[[paste0('mr_joint_r_coefficients_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_joint_g_run
        
        if (g_resetresult == FALSE || is.null(mr_joint_results()$F.Est))
            return()
      
        mr_joint_results()$F.Est
    }, rownames = TRUE, bordered = TRUE, caption = "Coefficients", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
    
    output[[paste0('mr_joint_r_estimates_',k)]] <- renderTable({
        input$file
        input$resetData
        input$di_option_run
        input$mr_joint_g_run
        
        if (g_resetresult == FALSE || is.null(mr_joint_results()$D.Est))
            return()
      
        mr_joint_results()$D.Est
    }, rownames = TRUE, bordered = TRUE, caption = "Estimates", spacing = "m",
    caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)
})

# Joint Model Plots ####
mr_joint_m_reactiveplot1 <- reactive({
    if(g_resetresult == FALSE || is.null(mr_joint_results()))
        return()
    
    res <- mr_joint_results()
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

lapply(1:6, function(k) {
    
    output[[paste0('mr_joint_m_plot1_',k)]] <- renderPlot({
        input$file
        input$resetData
        input$di_option_run
        input$mr_joint_g_run
        
        local({
            output[[paste0('mr_joint_m_downloadplot1_',k)]] <- downloadPlot({
                mr_joint_m_reactiveplot1()
            })
        })
        
        mr_joint_m_reactiveplot1()
    })
    
    output[[paste0('mr_joint_m_showplot1_',k)]] <- renderUI({
        input$file
        input$resetData
        input$di_option_run
        input$mr_joint_g_run
        
        if (g_resetresult == FALSE || is.null(mr_joint_results()))
            return()
        
        div(
            h4("Model Checking Plots for Mean"),
            div(
                style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
                withSpinner(
                    plotOutput(paste0('mr_joint_m_plot1_', k), height = "600px"),
                    type = 1,
                    color = "#2c3e50",
                    size = 1.2
                ),
                div(
                    style = "position: absolute; left:0.5em; bottom: 0.5em;",
                    dropdown(
                        downloadButton(outputId = paste0('mr_joint_m_downloadplot1_', k), label = "Download Plot"),
                        size = "xs",
                        icon = icon("download", class = "opt"),
                        up = TRUE
                    )
                )
            ),
            br()
        )
    })
})