tabPanel(
    "HSEM",
    sidebarLayout(
        sidebarPanel(
            h3(strong("HSEM")),
            fluidRow(
                column(
                    6,
                    uiOutput("mr_dsem2_g_respslider")
#                   uiOutput("mr_dsem2_g_corr_structure_1")
                ),
                column(
                    6,
                    fluidRow(
                        column(
                            10
#                           uiOutput("mr_dsem2_g_corr_structure"),
#                           uiOutput("mr_dsem2_g_corr_structure_2")
                        ),
                        column(
                            4,
                            splitLayout(style = "padding: 9px;"),
                            actionButton("mr_dsem2_g_run", "Run", icon = icon("play")),
                            style = "text-align:center;"
                        )   
                    ),
                    uiOutput(paste0("mr_dsem2_m_check_comparison_", 1))   
                )
            ),
            # uiOutput("mr_sem_resp_select"),
            uiOutput("mr_dsem2_g_tabpanel"),
            width = 5
        ),
        mainPanel(
            h3("HSEM Summary"),
            hr(),
            tabsetPanel(
                id = "mr_dsem2_resulttabset",
                tabPanel(
                    "Model Summary",
                    br(),
                    mainPanel(
                        uiOutput("mr_dsem2_var1_select"),
                        uiOutput("mr_dsem2_g_summary"),
                        width = 12
                    )
                ),
                tabPanel(
                    "Model Checking Plot",
                    br(),
                    mainPanel(
                        uiOutput("mr_dsem2_var2_select"),
                        uiOutput("mr_dsem2_g_mcplots"),
                        width = 12
                    )
                ),
                tabPanel(
                    "Diagram",
                    br(),
                    h4("Correlation for Variables"),
                    uiOutput("mr_dsem2_r_showplot01"),
                    h4("Correlation for Factors"),
                    uiOutput("mr_dsem2_r_showplot02"),
                    h4("Path Diagram"),
                    h4("Model for Mean"),
                    uiOutput("mr_dsem2_r_showplot03"),
                    uiOutput("mr_dsem2_r_showplot04"),
#                    h4("Structural Model"),
#                    uiOutput("mr_dsem2_r_showplot05"),
                    width = 12
                )
            ),
            width = 7
        )
    )
)
