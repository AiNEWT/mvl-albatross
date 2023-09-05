tabPanel(
    "t-test",
    sidebarPanel(
        fluidRow(
            column(
                7,
                h3(strong("t-test"))
            ),
            column(
                5,
                actionButton("ba_tt_run", "Run", icon = icon("play")),
                style = "text-align:right; padding:15px"
            )
        ),
        uiOutput("ba_tt_type"),
        uiOutput("ba_tt_selectvarname1"),
        uiOutput("ba_tt_selectvarname2"),
        uiOutput("ba_tt_groupvarname1"),
        uiOutput("ba_tt_groupvarname2"),
        uiOutput("ba_tt_mu"),
        uiOutput("ba_tt_equal"),
        uiOutput("ba_tt_alternative"),
        uiOutput("ba_tt_level"),
        uiOutput("ba_tt_shapirocheck"),
        uiOutput("ba_tt_levenecheck"),
        uiOutput("ba_tt_wilcoxcheck")
    ),
    mainPanel(
        h3(uiOutput("ba_tt_title")),
        hr(),
        tabsetPanel(
            id = "ba_tt_resulttabset",
            tabPanel(
                "Model Summary",
                br(),
                uiOutput("ba_tt_hypotheses"),
                uiOutput("ba_tt_descresults"),
                uiOutput("ba_tt_ttestresults"),
                uiOutput("ba_tt_shapiroresults"),
                uiOutput("ba_tt_leveneresults"),
                uiOutput("ba_tt_wilcoxresults")
            ),
            tabPanel(
                "Data View",
                br(),
                uiOutput("ba_tt_showplot3")
            )
        )
    )
)
