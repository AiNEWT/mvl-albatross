tabPanel(
    "Data Import",
    introjsUI(),
    tags$head(
        tags$style(
            HTML(
            "caption {
                color: #2C3E50;
                font-size: 130%;
            }"
            )
        )
    ),
    fileInput("file", "Upload File", multiple = FALSE, width = '30%',
              accept = c(".csv",".txt",".xlsx",".xls",".RData", ".RDS", ".SAV"), 
              placeholder = "*.csv,*.txt,*.xls,*.xlsx,*.RData,*.RDS,*.SAV"),
    uiOutput("di_option"),
    downloadButton("downloadData", "Download Data"),
    downloadButton("download_rdata", "Download RData"),
    actionButton('resetData', 'Reset Data', icon = icon("redo")),
    uiOutput("editUI")
)
