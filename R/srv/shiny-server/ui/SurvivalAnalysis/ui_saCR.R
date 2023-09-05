tabPanel(
    "Competing Risk Model", 
    sidebarPanel(
        fluidRow(
            column(
                9,
                h3(strong("Competing Risk"))
            ),
            column(
                3,
                actionButton("sa_cr_g_run", "Run", icon = icon("play")),
                style = "text-align:right; padding:15px;"
            )
        ),
        uiOutput("sa_cr_s1_model"), 
        uiOutput("sa_cr_s2_model"),
        uiOutput("sa_cr_r_accordion"),
        width = 4
    ), #sidebarPanel finished
    mainPanel(
        h3("Competing Risk Model Summary"),
        hr(),
        tabsetPanel(
            tabPanel(
                "Model Summary",
                tableOutput("sa_cr_r_coefficients"),
                tableOutput("sa_cr_r_estimates")
            ),
            tabPanel(
                "Random Effect Inferences",
                hr(),
                uiOutput("sa_cr_r_showplot1")
            )
        ),
        width = 8
    )
) #tabPanel finished