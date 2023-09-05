tabPanel(
    "Joint Model",
    sidebarLayout(
        sidebarPanel(
            h3(strong("Joint Model")),
            fluidRow(
                column(
                    6,
                    uiOutput("mr_joint_g_respslider")
                ),
                column(
                    4,
                    uiOutput("mr_joint_g_corr_structure")
                ),
                column(
                    2,
                    splitLayout(style = "padding: 9px;"),
                    actionButton("mr_joint_g_run", "Run", icon = icon("play")),
                    style = "text-align:center;"
                )
            ),
            uiOutput("mr_joint_g_file"),
            uiOutput("mr_joint_g_tabpanel"),
            width = 5
        ),
        mainPanel(
            h3("Joint Model Summary"),
            hr(),
            tabsetPanel(
                tabPanel(
                    "Model Summary",
                    br(),
                    mainPanel(
                        uiOutput(paste0('mr_joint_r_coefficients_', 1)),
                        uiOutput(paste0('mr_joint_r_estimates_', 1)),
                        width = 12
                    )
                ), 
                tabPanel(
                    "Random Effect Inferences",
                    br(),
                    mainPanel(
                        uiOutput(paste0('mr_joint_m_showplot1_', 1)),
                        width = 12
                    )
                )
            ),
            width = 7
        )
    )
)