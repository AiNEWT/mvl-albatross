tabPanel(
    "Frailty Model", 
    sidebarLayout(
        sidebarPanel(
            fluidRow(
                column(
                    9, 
                    h3(strong("Frailty Model"))
                ),
                column(
                    3, 
                    actionButton("sa_fm_g_run", "Run", icon = icon("play")),
                    style = "text-align:right; padding:15px;"
                )
            ),
            hr(),
            uiOutput("sa_fm_s_model"), 
            uiOutput("sa_fm_s_survivaltime"),
            uiOutput("sa_fm_s_check_initial"), 
            uiOutput("sa_fm_s_initialtime"),
    
            uiOutput("sa_fm_s_variable"),
            fluidRow(
                column(
                    8,
                    uiOutput("sa_fm_s_interaction")
                ),
                column(
                    4,
                    uiOutput("sa_fm_s_interactionappend"),
                    style = "text-align:right; padding:15px"
                )
            ),
            
            uiOutput("sa_fm_s_rand"),
            uiOutput("sa_fm_s_check_randinteraction"),
            fluidRow(
                column(
                    8,
                    uiOutput("sa_fm_s_randinteraction")
                ),
                column(
                    4,
                    uiOutput("sa_fm_s_randinteractionappend"),
                    style = "text-align:right; padding:15px"
                )
            ),
            
            uiOutput("sa_fm_s_randfamily"),
            uiOutput("sa_fm_s_indicator"), 
            uiOutput("sa_fm_s_check_group"),
            
            # I don't know what it means in frailty model
            # hr(),
            # uiOutput("sa_fm_s_check_robustse"),
            # uiOutput("sa_fm_s_check_confint"),
            # uiOutput("sa_fm_s_check_exp"),
            
            hr(),
            uiOutput("sa_fm_s_check_comparison"),
            uiOutput("sa_fm_s_check_rcodes"),
            
            br(),
            h3(strong("Additional Settings")),
            h5(helpText("Order of Laplace Approximation for Likelihood(mean) and Restricted Likelihood(Dispersion)")),
            uiOutput("sa_fm_a_mord"),
            uiOutput("sa_fm_a_dord"),
            
            width = 3
        ), #sidebarPanel finished
        mainPanel(
            h3(uiOutput("sa_fm_r_title")),
            hr(),
            tabsetPanel(
                id = "sa_fm_resulttabset",
                tabPanel(
                    "Model Summary",
                    uiOutput("sa_fm_r_mainsummary"), 
                    uiOutput("sa_fm_fix_coef"),
                    uiOutput("sa_fm_rand_coef"),
                    uiOutput("sa_fm_r_likelihood"),
                    uiOutput("sa_fm_r_aic"),
                    uiOutput("sa_fm_r_comparisonmodel"),
                    uiOutput("sa_fm_r_rcodes")
                )
            ),
            width = 9
        )
    ) #sidebarLayout finished
) #tabPanel finished