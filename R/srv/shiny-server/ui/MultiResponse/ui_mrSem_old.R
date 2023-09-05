tabPanel(
    "SEM",
    sidebarLayout(
        sidebarPanel(
            h3(strong("SEM")),
            fluidRow(
                column(
                    6,
                    uiOutput("mr_sem_g_factorslider")
                ),
                column(
                    4,
                    splitLayout(style = "padding: 9px;"),
                    actionButton("mr_sem_g_changefactorname", "Change Factor Names", icon = icon("sync")),
                    style = "text-align:center;"
                ),
                column(
                    2,
                    splitLayout(style = "padding: 9px;"),
                    actionButton("mr_sem_g_run", "Run", icon = icon("play")),
                    style = "text-align:center;"
                )
            ),
            uiOutput("mr_sem_g_model"),
            uiOutput("mr_sem_g_tabpanel"),
            uiOutput("mr_sem_g_dist"),
            uiOutput("mr_sem_g_link"),
            width = 5
        ),
        mainPanel(
            h3("SEM Summary"),
            hr(),
            uiOutput("mr_sem_g_summary"),
            width = 7
        )
    )
)