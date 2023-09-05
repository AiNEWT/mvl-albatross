tabPanel(
    "Kaplan-Meier Estimator",
    sidebarPanel(
        fluidRow(
            column(
                7,
                h3(strong("Kaplan-Meier Estimator"))
            ),
            column(
                5,
                actionButton("sa_km_run", "Run", icon = icon("play")),
                style = "text-align:right; padding:15px"
            )
        ),
        uiOutput("sa_km_survivaltime"),
        uiOutput("sa_km_check_initial"),
        uiOutput("sa_km_initialtime"),
        uiOutput("sa_km_indicator"),
        uiOutput("sa_km_checkgroup"),
        uiOutput("sa_km_groupvariable"),
        uiOutput("sa_km_checklogrank"),
        # uiOutput("sa_km_checkgehan"),   # unfinished
        # uiOutput("sa_km_checktaroneware")     # unfinished
        width = 3
    ),
    mainPanel(
        h3(uiOutput("sa_km_title")),
        hr(),
        tabsetPanel(
            id = "sa_km_resulttabset",
            tabPanel(
                "Data Summary",
                br(),
                tableOutput("sa_km_kmresults"), 
                tableOutput("sa_km_groupkmresults"),
                tableOutput("sa_km_logrankresults1"),
                tableOutput("sa_km_logrankresults2")
                # tableOutput("sa_km_gehanresults"),        # unfinished
                # tableOutput("sa_km_taronewareresults")        # unfinished
            ),
            tabPanel(
                "Kaplan-Meier Estimate",
                br(),
                tableOutput("sa_km_kmestimates")
            ),
            tabPanel(
                "Kaplan-Meier Curve",
                br(),
                uiOutput("sa_km_showkmcurveplot1")
            )
        ),
        width = 9
    )
)