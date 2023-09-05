tabPanel(
    "Linear Mixed Model", 
    sidebarPanel(
        fluidRow(
            column(
                9, 
                h3(strong("LMM"))
            ),
            column(
                3, 
                actionButton("re_lmm_g_run", "Run", icon = icon("play")),
                style = "text-align:right; padding:15px"
            )
        ),
        uiOutput("re_lmm_m_model"),
        uiOutput("re_lmm_r_accordion"),
        width = 4
    ), #sidebarPanel closed
    mainPanel(
        h3("LMM Summary"),
        hr(),
        tabsetPanel(
            id = "re_lmm_resulttabset",
            tabPanel(
                "Model Summary",
                br(),
                uiOutput("re_lmm_r_mainsummary"),
                uiOutput("re_lmm_m_coeff"),
                uiOutput("re_lmm_l_coeff"),
                uiOutput("re_lmm_p_coeff"),
                uiOutput("re_lmm_r_likelihood"),
                hr(),
                tableOutput("re_lmm_r_comparisonmodel1"),
                tableOutput("re_lmm_r_comparisonmodel2"),
                tableOutput("re_lmm_r_rcodes")
            ),
            tabPanel(
                "Model Checking Plot",
                br(),
                uiOutput("re_lmm_m_showplot1")
            ),
            tabPanel(
                "Prediction", 
                br(),
                uiOutput("re_lmm_r_check_95mu"),
                uiOutput("re_lmm_r_box_95mu"),
                verbatimTextOutput("re_lmm_r_prediction"),
                
                uiOutput("re_lmm_r_check_calculator"),
                uiOutput("re_lmm_r_box_calculator"),
                verbatimTextOutput("re_lmm_r_calculator")
            )
        ),
        width = 8
    )
) #tabPanel closed