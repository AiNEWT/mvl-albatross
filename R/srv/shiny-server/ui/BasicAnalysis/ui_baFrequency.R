tabPanel(
    "Frequency Analysis",
    sidebarPanel(
        fluidRow(
            column(
                7,
                h3(strong("Frequency Analysis"))
            ),
            column(
                5,
                actionButton("ba_fa_run", "Run", icon = icon("play")),
                style = "text-align:right; padding:15px"
            )
        ),
        uiOutput("ba_fa_selectvarname1"),
        uiOutput("ba_fa_checkcolumn"),
        uiOutput("ba_fa_selectvarname2"),
        uiOutput("ba_fa_rpcheck"),
        uiOutput("ba_fa_cpcheck"),
        uiOutput("ba_fa_pcheck"),
        uiOutput("ba_fa_cscheck"),
        uiOutput("ba_fa_fecheck"), 
        uiOutput("ba_fa_fecheck_option"), 
        uiOutput("ba_fa_mncheck")
    ),
    mainPanel(
        h3(uiOutput("ba_fa_title")),
        hr(),
        tabsetPanel(
            tabPanel(
                "Data Summary",
                br(),
                uiOutput("ba_fa_faresults"),
                uiOutput("ba_fa_rpresults"),
                uiOutput("ba_fa_cpresults"),
                uiOutput("ba_fa_presults"),
                uiOutput("ba_fa_csresults"),
                uiOutput("ba_fa_feresults"),
                uiOutput("ba_fa_mnresults")
            ),
            tabPanel(
                "Chart",
                br(),
                uiOutput("ba_fa_showplot1")
            )
        )
    )
)