tabPanel(
    "Delete variable",
    value = "dm_del_tab",
    fluidRow(
        column(
            6,
            h3(strong("Delete variable"))
        ),
        column(
            6,
            actionButton("dm_del_run", "Run", icon = icon("play")),
            style = "text-align:right; padding:15px"
        )
    ),
    uiOutput("dm_del_selectvarname1")
)