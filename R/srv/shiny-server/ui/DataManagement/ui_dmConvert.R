tabPanel(
    "Convert Data Type",
    value = "dm_ct_tab",
    fluidRow(
        column(
            7,
            h3(strong("Convert Data Type"))
        ),
        column(
            5,
            actionButton("dm_ct_run", "Run", icon = icon("play")),
            style = "text-align:right; padding:15px"
        )
    ),
    uiOutput("dm_ct_selecttype"),
    uiOutput("dm_ct_variables"),
    width=6
)