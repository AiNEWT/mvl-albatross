tabPanel(
    "- Factor Model",
    sidebarLayout(
        sidebarPanel(
            h3(strong("Factor Model")),
            fluidRow(
                column(
                    6,
                    uiOutput("mr_factor_g_respslider")
                ),
                column(
                    4,
                    uiOutput("mr_factor_g_corr_structure")
                ),
                column(
                    2,
                    splitLayout(style = "padding: 9px;"),
                    actionButton("mr_factor_g_run", "Run", icon = icon("play")),
                    style = "text-align:center;"
                )
            ),
            fluidRow(
                column(
                    6,
                    uiOutput("mr_factor_a_order")
                ),
                column(
                    6,
                    uiOutput("mr_factor_a_reml")
                )
            ),
            uiOutput("mr_factor_g_tabpanel"),
            width = 5
        ),
        mainPanel(
            h3("Factor Model Summary"),
            hr(),
            uiOutput("mr_factor_r_lambda"),
            uiOutput("mr_factor_r_gamma"),
            uiOutput("mr_factor_r_beta"),
            uiOutput("mr_factor_r_likelihood"),
            width = 7
        )
    )
)