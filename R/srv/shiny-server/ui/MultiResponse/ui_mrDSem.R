tabPanel(
    "- Dynamic SEM", 
    sidebarPanel(
        fluidRow(
            column(
                9,
                h3(strong("Dynamic SEM"))
            ),
            column(
                3,
                actionButton("re_dsem_g_run", "Run", icon = icon("play")),
                style = "text-align:right; padding:15px;"
            )
        ),
        uiOutput("re_dsem_m_model"),
        uiOutput("re_dsem_p_model"),
        uiOutput("re_dsem_l_model"),
        uiOutput("re_dsem_r_accordion"),
        width = 4
    ), #sidebarPanel closed
    mainPanel(
        h3("Dynamic SEM Summary"),
        hr(),
        tabsetPanel(
            id = "re_dsem_resulttabset",
            tabPanel(
                "Model Summary",
                br(),
                uiOutput("re_dsem_r_mainsummary"),
                uiOutput("re_dsem_m_coeff"),
                uiOutput("re_dsem_l_coeff"),
                uiOutput("re_dsem_l_taucoeff"),
                uiOutput("re_dsem_p_coeff"),
                uiOutput("re_dsem_p_alphacoeff"),
                uiOutput("re_dsem_r_likelihood"),
                hr(),
                uiOutput("re_dsem_r_rcodes")
            ),
            tabPanel(
                "Model Checking Plot",
                br(),
                uiOutput("re_dsem_m_showplot1"),
                uiOutput("re_dsem_p_showplot1"),
                uiOutput("re_dsem_l_showplot1")
            ),
            tabPanel(
                "Prediction", 
                br(),
                uiOutput("re_dsem_r_check_95mu"),
                uiOutput("re_dsem_r_box_95mu"),
                uiOutput("re_dsem_r_check_95phi"),
                uiOutput("re_dsem_r_box_95phi"),
                verbatimTextOutput("re_dsem_r_prediction"),
                
                uiOutput("re_dsem_r_check_calculator"),
                uiOutput("re_dsem_r_box_calculator"),
                verbatimTextOutput("re_dsem_r_calculator")
            )
        ),
        width = 8
    )
) #tabPanel closed