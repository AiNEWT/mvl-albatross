tabPanel(
    "Make Interval Variable",
    value = "dm_mvi_tab",
    fluidRow(
        column(
            7,
            h3(strong("Make Interval Variable"))
        ),
        column(
            5,
            actionButton("dm_mvi_run", "Run", icon = icon("play")),
            style = "text-align:right; padding:15px"
        )
    ),
    uiOutput("dm_mvi_newvariablename"),
#    uiOutput("dm_mvi_expression"),
    uiOutput("dm_mvi_selectvarname"),
 #   uiOutput("dm_mvi_check_tools"),
    uiOutput("dm_mvi_numberinterval"),
    uiOutput("dm_mvi_selectmethod")
)
