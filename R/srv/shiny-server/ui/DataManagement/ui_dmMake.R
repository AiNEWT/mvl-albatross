tabPanel(
    "Make Variable",
    value = "dm_mv_tab",
    fluidRow(
        column(
            7,
            h3(strong("Make Variable"))
        ),
        column(
            5,
            actionButton("dm_mv_run", "Run", icon = icon("play")),
            style = "text-align:right; padding:15px"
        )
    ),
    uiOutput("dm_mv_newvariablename"),
    uiOutput("dm_mv_expression"),
    uiOutput("dm_mv_check_tools"),
    uiOutput("dm_mv_selectvarname"),
    uiOutput("dm_mv_selectarith")
)
