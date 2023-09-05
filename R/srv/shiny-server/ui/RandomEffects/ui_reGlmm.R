tabPanel(
    "Generalized Linear Mixed Model", 
    sidebarPanel(
        fluidRow(
            column(
                9, 
                h3(strong("GLMM"))
            ),
            column(
                3, 
                actionButton("re_glmm_g_run", "Run", icon = icon("play")),
                style = "text-align:right; padding:15px"
            )
        ),
        uiOutput("re_glmm_m_model"),
        uiOutput("re_glmm_r_accordion"),
        width = 4
    ), #sidebarPanel closed
    mainPanel(
        h3("GLMM Summary"),
        hr(),
        tabsetPanel(
            id = "re_glmm_resulttabset",
            tabPanel(
                "Model Summary",
                br(),
                uiOutput("re_glmm_r_mainsummary"),
                uiOutput("re_glmm_m_coeff"),
                uiOutput("re_glmm_l_coeff"),
                uiOutput("re_glmm_p_coeff"),
                uiOutput("re_glmm_r_likelihood"),
                hr(),
                tableOutput("re_glmm_r_comparisonmodel1"),
                tableOutput("re_glmm_r_comparisonmodel2"),
                tableOutput("re_glmm_r_rcodes")
            ),
            tabPanel(
                "Model Checking Plot",
                br(),
                uiOutput("re_glmm_m_showplot1")
            ),
            tabPanel(
                "Prediction", 
                br(),
                uiOutput("re_glmm_r_check_95mu"),
                uiOutput("re_glmm_r_box_95mu"),
                verbatimTextOutput("re_glmm_r_prediction"),
                
                uiOutput("re_glmm_r_check_calculator"),
                uiOutput("re_glmm_r_box_calculator"),
                verbatimTextOutput("re_glmm_r_calculator")
            )
        ),
        width = 8
    )
) #tabPanel closed