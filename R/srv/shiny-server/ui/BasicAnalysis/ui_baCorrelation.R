tabPanel(
    "Correlation Analysis",
    sidebarPanel(
        fluidRow(
            column(
                7,
                h3(strong("Correlation Analysis"))
            ),
            column(
                5,
                actionButton("ba_co_run", "Run", icon = icon("play")),
                style = "text-align:right; padding:15px"
            )
        ),
        uiOutput("ba_co_selectvarname1"),
        div(style = "padding:10px"),
        uiOutput("ba_co_cofftype"),
        uiOutput("ba_co_plottype"),
        uiOutput("ba_co_bins")
        # uiOutput("ba_co_check_formula")
    ),
    mainPanel(
        h3(uiOutput("ba_co_title")),
        hr(),
        tabsetPanel(
            tabPanel(
                "Correlation Results",
                br(),
                tableOutput("ba_co_coresults"),
                tableOutput("ba_co_nresults"),
                tableOutput("ba_co_presults")
            ),
            tabPanel(
                "Correlation Plots",
                br(),
                uiOutput("ba_co_showplot1")
            )
        )
    )
)