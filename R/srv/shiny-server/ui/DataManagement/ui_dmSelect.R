tabPanel(
    "Subset",
    value = "dm_sr_tab",
    fluidRow(
        column(
            6,
            h3(strong("Subset"))
        ),
        column(
            6,
            actionButton("dm_sr_preview", "Preview", icon = icon("eye")),
            actionButton("dm_sr_run", "Run", icon = icon("play")),
            style = "text-align:right; padding:15px"
        )
    ),
    uiOutput("dm_sr_expression"),
    uiOutput("dm_sr_check_tools"),
    uiOutput("dm_sr_selectvarname1"),
    uiOutput("dm_sr_selectvarname2"),
    uiOutput("dm_sr_selectarith")
)