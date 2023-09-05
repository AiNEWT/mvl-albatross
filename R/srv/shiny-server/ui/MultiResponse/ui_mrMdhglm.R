tabPanel(
    "MDHGLM",
    sidebarLayout(
        sidebarPanel(
            h3(strong("MDHGLM")),
            fluidRow(
                column(
                    6,
                    uiOutput("mr_mdhglm_g_respslider")
                ),
                column(
                    6,
                    fluidRow(
                        column(
                        8,
                        uiOutput("mr_mdhglm_g_corr_structure")
                        ),
                        column(
                            4,
                            splitLayout(style = "padding: 9px;"),
                            actionButton("mr_mdhglm_g_run", "Run", icon = icon("play")),
                            style = "text-align:center;"
                        )   
                    ),
                    uiOutput(paste0("mr_mdhglm_m_check_comparison_", 1))   
                )
            ),
            uiOutput("mr_mdhglm_g_tabpanel"),
            width = 5
        ),
        mainPanel(
            h3("MDHGLM Summary"),
            hr(),
            tabsetPanel(
                id = "mr_mdhglm_resulttabset",
                tabPanel(
                    "Model Summary",
                    br(),
                    mainPanel(
                        uiOutput("mr_mdhglm_g_summary"),
                        width = 12
                    )
                )
            ),
            width = 7
        )
    )
)