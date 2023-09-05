tabPanel(
    "Hierarchical GLM", 
    sidebarPanel(
        fluidRow(
            column(
                9, 
                h3(strong("HGLM"))
            ),
            column(
                3, 
                actionButton("re_hglm_g_run", "Run", icon = icon("play")),
                style = "text-align:right; padding:15px"
            )
        ),
        uiOutput("re_hglm_m_model"),
        uiOutput("re_hglm_p_model"),
        uiOutput("re_hglm_l_model"),
        uiOutput("re_hglm_r_accordion"),
        width = 4
    ), #sidebarPanel closed
    mainPanel(
        h3("HGLM Summary"),
        hr(),
        tabsetPanel(
            id = "re_hglm_resulttabset",
            tabPanel(
                "Model Summary",
                br(),
                uiOutput("re_hglm_r_mainsummary"),
                uiOutput("re_hglm_m_coeff"),
                uiOutput("re_hglm_l_coeff"),
                uiOutput("re_hglm_p_coeff"),
                uiOutput("re_hglm_r_likelihood"),
                hr(),
                uiOutput("re_hglm_r_rcodes")
            ),
            tabPanel(
                "Model Checking Plot",
                br(),
                uiOutput("re_hglm_m_showplot1"),
                uiOutput("re_hglm_p_showplot1"),
                uiOutput("re_hglm_l_showplot1")
            ),
            tabPanel(
                "Prediction", 
                br(),
                uiOutput("re_hglm_r_check_95mu"),
                uiOutput("re_hglm_r_box_95mu"),
                uiOutput("re_hglm_r_check_95phi"),
                uiOutput("re_hglm_r_box_95phi"),
                verbatimTextOutput("re_hglm_r_prediction"),
                
                uiOutput("re_hglm_r_check_calculator"),
                uiOutput("re_hglm_r_box_calculator"),
                verbatimTextOutput("re_hglm_r_calculator")
            )
        ),
        width = 8
    )
) #tabPanel closed