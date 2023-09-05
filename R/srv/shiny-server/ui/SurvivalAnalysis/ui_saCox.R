tabPanel(
    "Cox Model",
    sidebarLayout(
        sidebarPanel(
            fluidRow(
                column(
                    9,
                    h3(strong("Cox Model"))
                ),
                column(
                    3,
                    actionButton("sa_cox_run", "Run", icon = icon("play")),
                    style = "text-align:right; padding:15px;"
                )
            ),
            hr(),
            uiOutput("sa_cox_model"),
            uiOutput("sa_cox_survivaltime"),
            uiOutput("sa_cox_check_initial"),
            uiOutput("sa_cox_initialtime"),
            
            uiOutput("sa_cox_variable"),
            fluidRow(
                column(
                    8,
                    uiOutput("sa_cox_interaction")
                ),
                column(
                    4,
                    uiOutput("sa_cox_interactionappend"),
                    style = "text-align:right; padding:15px"
                )
            ),
             
            uiOutput("sa_cox_indicator"),
            
            uiOutput("sa_cox_ties"),   
            uiOutput("sa_cox_check_robustse"),
            width = 3

        ),
        mainPanel(
            h3(uiOutput("sa_cox_title")),
            hr(),
            h3(uiOutput("sa_cox_nullmodelNotification")),
            fluidRow(
                column(
                    3,
                    tableOutput("sa_cox_numberresults"),
                    tableOutput("sa_cox_concordanceresults")
                ),
                column(
                    3,
                    tableOutput("sa_cox_robustresults")
                ),
                column(
                    6,
                    tableOutput("sa_cox_testresults")                   
                )
            ),
            tableOutput("sa_cox_coefresults"),
            tableOutput("sa_cox_confintresults"),
            tableOutput("sa_cox_zphresults")            
        )
    )
)