tabPanel(
    "- Shared",
    sidebarLayout(
        sidebarPanel(
            h3(strong("MDHGLM Shared")),
            fluidRow(
                column(
                    6,
                    uiOutput("mr_shared_g_respslider")
                ),
                column(
                    4,
                    uiOutput("mr_shared_g_corr_structure")
                ),
                column(
                    2,
                    splitLayout(style = "padding: 9px;"),
                    actionButton("mr_shared_g_run", "Run", icon = icon("play")),
                    style = "text-align:center;"
                )
            ),
            uiOutput("mr_shared_g_tabpanel"),
            width = 5
        ),
        mainPanel(
            h3("MDHGLM Shared Summary"),
            hr(),
            uiOutput("mr_shared_r_betacoeff"),
            uiOutput("mr_shared_r_phicoeff"),
            uiOutput("mr_shared_r_lambdacoeff"),
            uiOutput("mr_shared_r_sharedparameter"),
            uiOutput("mr_shared_r_likelihood"),
            width = 7
        )
    )
)