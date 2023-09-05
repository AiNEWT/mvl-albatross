tabPanel(
    "ANOVA",
    sidebarPanel(
        fluidRow(
            column(
                7,
                h3(strong("ANOVA"))
            ),
            column(
                5,
                actionButton("ba_an_run", "Run", icon = icon("play")),
                style = "text-align:right; padding:15px"
            )
        ),
        uiOutput("ba_an_selectvarname1"),
        uiOutput("ba_an_selectvarname2"),
        fluidRow(
            column(
                8,
                uiOutput("ba_an_interaction")
            ),
            column(
                4,
                uiOutput("ba_an_interactionappend"),
                style = "text-align:right; padding:15px"
            )
        ),
#        uiOutput("ba_an_missing"),
        uiOutput("ba_an_shapirocheck"),
        uiOutput("ba_an_bpcheck"),
        uiOutput("ba_an_kwcheck")
    ),
    mainPanel(
        h3(uiOutput("ba_an_title")),
        hr(),
        tabsetPanel(
            tabPanel(
                "Model Summary",
                br(),
                uiOutput("ba_an_anovaresults"),
                uiOutput("ba_an_specificresults"),
                uiOutput("ba_an_shapiroresults"),
                uiOutput("ba_an_bpresults"),
                uiOutput("ba_an_kwresults")
            ),
            tabPanel(
                "Model Checking Plots",
                br(),
                uiOutput("ba_an_showplot1"),
                uiOutput("ba_an_showplot2")
            )
        )
    )
)