tabPanel(
    "Rename variable",
    value = "dm_rn_tab",
    fluidRow(
        column(
            6,
            h3(strong("Rename variable"))
        ),
        column(
            6,
#            actionButton("dm_rn_preview", "Preview", icon = icon("eye")),
            actionButton("dm_rn_run", "Run", icon = icon("play")),
            style = "text-align:right; padding:15px"
        )
    ),
    uiOutput("dm_rn_selectvarname1"),
    uiOutput("dm_rn_expression")
)