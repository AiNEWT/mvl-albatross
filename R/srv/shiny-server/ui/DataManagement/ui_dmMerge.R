tabPanel(
    "Merge Datasets",
    fluidRow(
        column(
            6,
            h3(strong("Merge Datasets"))
        ),
        column(
            6,
            actionButton("dm_md_preview", "Preview", icon = icon("eye")),
            actionButton("dm_md_run", "Run", icon = icon("play")),
            style = "text-align:right; padding:15px"
        )
    ),
    fileInput("md_file","Upload Second Data", multiple = TRUE),
    uiOutput("dm_md_sheetindex"),
    checkboxInput(inputId = 'md_header', label = 'Header', value = TRUE),
    selectInput(
        "md_sep",
        "Separator",
        choices = c(
            Comma = ',',
            Semicolon = ';',
            Tab = '\t',
            Space = ''
        ),
        multiple = FALSE
    ),
    uiOutput("dm_md_selectvarname1"),
    uiOutput("dm_md_selectvarname2"),
    uiOutput("dm_md_mainall"),
    uiOutput("dm_md_secondall")
)
