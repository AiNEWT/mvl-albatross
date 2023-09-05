tabPanel(
    "Linear Model",
    sidebarLayout(
        sidebarPanel(
            fluidRow(
                column(
                    8,
                    h3(strong("Linear Model"))
                ),
                column(
                    4,
                    actionButton("rg_lm_g_run", "Run", icon = icon("play")),
                    style = "text-align:right; padding:15px"
                )
            ),
            uiOutput("rg_lm_m_model"),
            uiOutput("rg_lm_m_resp"),
            uiOutput("rg_lm_m_variable"),
            div(style = "padding:10px"),
            fluidRow(
                column(
                    8,
                    uiOutput("rg_lm_m_interaction")
                ),
                column(
                    4,
                    uiOutput("rg_lm_m_interactionappend"),
                    style = "text-align:right; padding:15px"
                )
            ),
            uiOutput("rg_lm_m_check_nointercept"),
            uiOutput("rg_lm_m_check_offset"),
            uiOutput("rg_lm_m_offset"),
            
            uiOutput("rg_lm_m_check_margins"),
            uiOutput("rg_lm_m_margins1"),
            uiOutput("rg_lm_m_margins2"),
            
            hr(style = "border-color: #2C3E50;"),
            uiOutput("rg_lm_m_check_vif"),
            uiOutput("rg_lm_m_check_robustse"),
            uiOutput("rg_lm_m_check_confint"),
            uiOutput("rg_lm_m_check_bptest"),
            
            hr(style = "border-color: #2C3E50;"),
            uiOutput("rg_lm_m_check_comparison"),
            uiOutput("rg_lm_m_check_rcodes"),
            width = 3
        ),
        mainPanel(
            h3(textOutput("rg_lm_r_model")),
            hr(),
            tabsetPanel(
                tabPanel(
                    "Model Summary",
                    br(),
                    fluidRow(
                        column(
                            8,
                            tableOutput("rg_lm_r_anova"),
                            tableOutput("rg_lm_r_coefficients")
                        ),
                        column(
                            4,
                            tableOutput("rg_lm_r_modelsummary")
                        )
                    ),
                    tableOutput("rg_lm_r_margins"),
                    hr(),
                    tableOutput("rg_lm_r_comparisonmodel"),
                    tableOutput("rg_lm_r_rcodes")
                ),
                tabPanel(
                    "Model Checking Plot",
                    br(),
                    uiOutput("rg_lm_r_showplot1")
                ),
                tabPanel(
                    "Prediction",
                    br(),
                    tableOutput("rg_lm_r_prediction")
                )
            ), width = 9
        )
    )
)