tabPanel(
    "Descriptive Statistics",
    sidebarPanel(
        fluidRow(
            column(
                7,
                h3(strong("Descriptive Statistics"))
            ),
            column(
                5,
                # actionButton("ba_ds_help", "Help", icon = icon("question")),
                actionButton("ba_ds_run", "Run", icon = icon("play")),
                style = "text-align:right; padding:15px"
            )
        ),
        uiOutput("ba_ds_selectvarname1"),
        uiOutput("ba_ds_check_percentile"),
        uiOutput("ba_ds_select_percentile"),
        uiOutput("ba_ds_checkgroup"),
        uiOutput("ba_ds_selectvarname2"),
        uiOutput("ba_ds_check_groupboxplot")
    ),
    mainPanel(
        h3(uiOutput("ba_ds_title")),
        hr(),
        tabsetPanel(
            id = "ba_ds_resulttabset",
            tabPanel(
                "Data Summary",
                br(),
                tableOutput("ba_ds_variableresults"),
                tableOutput("ba_ds_percentileresults"),
                tableOutput("ba_ds_groupresults"), 
                tableOutput("ba_ds_grouppercentileresults")
            ),
            tabPanel(
                "Variable Histogram",
                br(),
                uiOutput("ba_ds_bins"),
                uiOutput("ba_ds_showvariablehist")
            ), 
            tabPanel(
                "Box-Plot",
                br(),
                uiOutput("ba_ds_showboxplot")
            )
        )
    )
)