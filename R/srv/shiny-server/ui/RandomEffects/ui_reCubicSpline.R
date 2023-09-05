tabPanel(
    "- Cubic Spline", 
    sidebarPanel(
        fluidRow(
            column(
                9, 
                h3(strong("Cubic Spline"))
            ),
            column(
                3, 
                actionButton("re_cs_g_run", "Run"),
                style = "text-align:right; padding:15px"
            )
        ), 
        
        uiOutput("re_cs_m_resp"), 
        uiOutput("re_cs_m_variable"), 
        
        uiOutput("re_cs_m_select_method"), 
        uiOutput("re_cs_m_select_specific"), 
        
        uiOutput("re_cs_p_check_joint"),
        
        # uiOutput("re_cs_p_variable"),
        uiOutput("re_cs_p_select_method"),
        uiOutput("re_cs_p_select_specific"),
        
        # uiOutput("re_cs_m_dist"),
        # uiOutput("re_cs_m_link"),

        # uiOutput("re_cs_m_check_rcodes")


        width = 4
    ), #sidebarPanel closed
    mainPanel(
        h3("Cubic Spline Results"),
        hr(),
        uiOutput("re_cs_m_showplot"),
        uiOutput("re_cs_p_showplot"),
        # tabsetPanel(
        #     tabPanel(
        #         
        #         
        #     ), 
        #     tabPanel(
        #         plotOutput("re_cs_m_plot_residual", height="600px")
        #     ), 
        # ),
        width = 8
    )
) #tabPanel closed