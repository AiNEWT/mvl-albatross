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
            uiOutput("mr_sem_g_model_2"),
            uiOutput("mr_sem_g_tabpanel"),
            uiOutput("mr_sem_g_dist"),
            uiOutput("mr_sem_g_link"),
            uiOutput("mr_sem_g_comparison"),
            width = 5
        ),
        mainPanel(
            h3("SEM Summary"),
            hr(),
            tabsetPanel(
            #    tabPanel(
            #        "Summary",
            #        br(),
            #        uiOutput("mr_sem_r_summary")
            #    ),
                id = "mr_sem_resulttabset",
                tabPanel(
                    "Model Summary",
                    br(),
                    fluidRow(
                        column(
                            8,
            #                uiOutput("mr_sem_r_pe1_0_1"),
            #                uiOutput("mr_sem_r_pe1_0_2"),
            #               uiOutput("mr_sem_r_pe1_1_1"),
                            uiOutput("mr_sem_r_pe1_1_2"),
            #               uiOutput("mr_sem_r_pe1_1_3"),
                            uiOutput("mr_sem_r_pe1"),
                            uiOutput("mr_sem_r_pe1_1"),
                            uiOutput("mr_sem_r_pe2"),
                            uiOutput("mr_sem_r_pe3"),
                            uiOutput("mr_sem_r_pe4")
                        )
                    )
                ),
                tabPanel(
                    "Model Checking Plot",
                    br(),
                    h4("Model Checking Plots for Measurement Models"),
        #            plotOutput("mr_sem_r_showplot2"),
                    uiOutput("mr_sem_var1_select"),
                    plotOutput("mr_sem_r_showplot4"),
                    h4("Model Checking Plots for Structural Models"),
        #           plotOutput("mr_sem_r_showplot3"),
                    uiOutput("mr_sem_var2_select"),
                    plotOutput("mr_sem_r_showplot5")
                ),
                tabPanel(
                    "Diagram",
                    br(),
                    h4("Correlation for Variables"),
                    uiOutput("mr_sem_r_showplot01"),
                    h4("Correlation for Factors"),
                    uiOutput("mr_sem_r_showplot02"),
                    h4("Path Diagram"),
                    uiOutput("mr_sem_r_showplot1")
                )
            ),
            width = 7
        )
    )
)
