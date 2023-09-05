tabPanel(
    "Factor Model",
    sidebarLayout(
        sidebarPanel(
            h3(strong("Factor Model")),
            fluidRow(
                column(
                    6,
                    uiOutput("mr_factor_g_factorslider")
                ),
                column(
                    4,
                    splitLayout(style = "padding: 9px;"),
                    actionButton("mr_factor_g_changefactorname", "Change Factor Names", icon = icon("sync")),
                    style = "text-align:center;"
                ),
                column(
                    2,
                    splitLayout(style = "padding: 9px;"),
                    actionButton("mr_factor_g_run", "Run", icon = icon("play")),
                    style = "text-align:center;"
                )
            ),
            uiOutput("mr_factor_g_model"),
            uiOutput("mr_factor_g_tabpanel"),
            uiOutput("mr_factor_g_dist"),
            uiOutput("mr_factor_g_link"),
            width = 5
        ),
        mainPanel(
            h3("Factor Model Summary"),
            hr(),
            uiOutput("mr_factor_g_summary"),
            width = 7
        )
    )
)