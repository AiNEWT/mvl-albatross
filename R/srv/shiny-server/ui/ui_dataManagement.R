navbarMenu(
"Data Management",
    tabPanel(
        "Data Handling",
        sidebarPanel(
            tabsetPanel(
                id = "dm_inputtabset",
                source("ui/DataManagement/ui_dmMake.R", local = TRUE)[[1]],
                source("ui/DataManagement/ui_dmMakeInterval.R", local = TRUE)[[1]],
                source('ui/DataManagement/ui_dmConvert.R', local = TRUE)[[1]],
                source('ui/DataManagement/ui_dmSelect.R', local = TRUE)[[1]],
                source('ui/DataManagement/ui_dmRename.R', local = TRUE)[[1]],
                source('ui/DataManagement/ui_dmDelete.R', local = TRUE)[[1]]
            )
        ),
        mainPanel(
            h3(uiOutput("dm_datatitle")),
            hr(),
            tabsetPanel(
                id = "dm_resulttabset",
                tabPanel(
                    "Data Frame",
                    value = "dm_dftab",
                    br(),
                    tableOutput("dm_dataframe")
                ),
                tabPanel(
                    "Data Attributes",
                    value = "dm_datab",
                    br(),
                    tableOutput("dm_dataattribute")
                ),
                tabPanel(
                    "Select Rows Preview",
                    value = "dm_srptab",
                    br(),
                    tableOutput("dm_selectrowspreview")
                )
            )
        )
    ),
    tabPanel(
        "Merge Dataset",
        sidebarPanel(
            tabsetPanel(
                source('ui/DataManagement/ui_dmMerge.R', local = TRUE)[[1]]
            )
        ),
        mainPanel(
            h3(uiOutput("dm_mergetitle")),
            hr(),
            tabsetPanel(
                tabPanel(
                    "Data Frame",
                    br(),
                    tableOutput("dm_firstdataframe"),
                    hr(),
                    tableOutput("dm_seconddataframe")
                ),
                tabPanel(
                    "Merge Dataset Preview",
                    br(),
                    tableOutput("dm_mergepreview")
                )
            )
        )
    )
)
