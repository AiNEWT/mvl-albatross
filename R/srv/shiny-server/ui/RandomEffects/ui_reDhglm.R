tabPanel(
    "Double HGLM", 
    sidebarPanel(
        fluidRow(
            column(
                9,
                h3(strong("DHGLM"))
            ),
            column(
                3,
                actionButton("re_dhglm_g_run", "Run", icon = icon("play")),
                style = "text-align:right; padding:15px;"
            )
        ),
        uiOutput("re_dhglm_m_model"),
        uiOutput("re_dhglm_p_model"),
        uiOutput("re_dhglm_l_model"),
        uiOutput("re_dhglm_r_accordion"),
        width = 4
    ), #sidebarPanel closed
    mainPanel(
        h3("DHGLM Summary"),
        hr(),
        tabsetPanel(
            id = "re_dhglm_resulttabset",
            tabPanel(
                "Model Summary",
                br(),
                uiOutput("re_dhglm_r_mainsummary"),
                uiOutput("re_dhglm_m_coeff"),
                uiOutput("re_dhglm_l_coeff"),
                uiOutput("re_dhglm_l_taucoeff"),
                uiOutput("re_dhglm_p_coeff"),
                uiOutput("re_dhglm_p_alphacoeff"),
                uiOutput("re_dhglm_r_likelihood"),
                hr(),
                uiOutput("re_dhglm_r_rcodes")
            ),
            tabPanel(
                "Model Checking Plot",
                br(),
                uiOutput("re_dhglm_m_showplot1"),
                uiOutput("re_dhglm_p_showplot1"),
                uiOutput("re_dhglm_l_showplot1")
            ),
            tabPanel(
                "Prediction", 
                br(),
                uiOutput("re_dhglm_r_check_95mu"),
                uiOutput("re_dhglm_r_box_95mu"),
                uiOutput("re_dhglm_r_check_95phi"),
                uiOutput("re_dhglm_r_box_95phi"),
                verbatimTextOutput("re_dhglm_r_prediction"),
                
                uiOutput("re_dhglm_r_check_calculator"),
                uiOutput("re_dhglm_r_box_calculator"),
                verbatimTextOutput("re_dhglm_r_calculator")
            )
        ),
        width = 8
    )
) #tabPanel closed