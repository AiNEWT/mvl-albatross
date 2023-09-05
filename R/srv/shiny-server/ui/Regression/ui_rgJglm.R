tabPanel(
    "Joint GLM", 
    sidebarPanel(
        fluidRow(
            column(
                9,
                h3(strong("Joint GLM"))
            ),
            column(
                3,
                actionButton("rg_jglm_g_run", "Run", icon = icon("play")),
                style = "text-align:right; padding:15px"
            )
        ),
        uiOutput("rg_jglm_m_model"),
        uiOutput("rg_jglm_p_model"),
        uiOutput("rg_jglm_r_accordion"),
        width = 4
    ), #sidebarPanel closed
    mainPanel(
        h3("Joint GLM Summary"),
        hr(),
        tabsetPanel(
            id = "rg_jglm_resulttabset",
            tabPanel(
                "Model Summary",
                br(),
                uiOutput("rg_jglm_r_mainsummary"),
                uiOutput("rg_jglm_m_coeff"),
                uiOutput("rg_jglm_p_coeff"),
                uiOutput("rg_jglm_r_likelihood"),
                hr(),
                tableOutput("rg_jglm_r_comparisonmodel1"),
                tableOutput("rg_jglm_r_comparisonmodel2"),
                tableOutput("rg_jglm_r_rcodes")
            ),
            tabPanel(
                "Model Checking Plot",
                br(),
                uiOutput("rg_jglm_m_showplot1"),
                uiOutput("rg_jglm_p_showplot1")
            ),
            tabPanel(
                "Prediction", 
                br(),
                uiOutput("rg_jglm_r_check_95mu"),
                uiOutput("rg_jglm_r_box_95mu"),
                uiOutput("rg_jglm_r_check_95phi"),
                uiOutput("rg_jglm_r_box_95phi"),
                verbatimTextOutput("rg_jglm_r_prediction"),
                
                uiOutput("rg_jglm_r_check_calculator"),
                uiOutput("rg_jglm_r_box_calculator"),
                verbatimTextOutput("rg_jglm_r_calculator")
            )
        ),
        width = 8
    )
) #tabPanel closed