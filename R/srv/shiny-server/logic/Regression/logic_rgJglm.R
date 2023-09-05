# JGLM Run ####

observeEvent(input$rg_jglm_g_run, {
    g_resetresult <<- TRUE
})

rg_jglm_results <- eventReactive(input$rg_jglm_g_run, {
    # ¦§ Warning & Notify ####
    # if user doesn't input some essential part, we wil return null value and show notification
    if (input$rg_jglm_m_resp == "")  {
      showNotification("Please choose response", type="warning")
      return()
    }
    
    if(input$rg_jglm_m_check_offset && input$rg_jglm_m_offset == "") {
      showNotification("Please choose offset variable", type="warning")
      return()
    }
    
    if(!is.null(input$rg_jglm_m_check_binomd) && input$rg_jglm_m_check_binomd && input$rg_jglm_m_binomd == "") {
      showNotification("Please choose binomial denominator", type="warning")
      return()
    }
    
    # Start Progressbar
    withProgress(message = 'Joint GLM', style = "notification", value = 0.1, {
        
    # ¦§ Variable Declaration 1 ####
    # variable declaration for which is used widely in rg_jglm_results
    dataLength = nrow(data2)
    
    # variable declaration for which is used in making MM(MeanModel), DM(DispersionModel)
    meanLink = input$rg_jglm_m_link
    phiLink = "log"
    
    if (!is.null(input$rg_jglm_p_check_phi) && input$rg_jglm_p_check_phi)
        phiLink = input$rg_jglm_p_link
    
    meanFormula = input$rg_jglm_m_model
    phiFormula = "constant"
    
    if (!is.null(input$rg_jglm_p_check_phi) && input$rg_jglm_p_check_phi)
        phiFormula = input$rg_jglm_p_model
    
    meanOffsetVariable = NULL
    phiOffsetVariable = NULL
    
    if (input$rg_jglm_m_check_offset)
        meanOffsetVariable = data2[[input$rg_jglm_m_offset]]
    if (!is.null(input$rg_jglm_p_check_offset) && input$rg_jglm_p_check_offset)
        phiOffsetVariable = data2[[input$rg_jglm_p_offset]]
    
    # ¦§ MeanModel (MM) ####
    MM = DHGLMMODELING(
        Model = "mean",
        Link = meanLink,
        LinPred = as.formula(meanFormula),
        Offset = meanOffsetVariable
    )
    
    # ¦§ DispersionModel (DM) ####
    if (phiFormula == "constant") 
        DM = DHGLMMODELING(
            Model = "dispersion",
            Link = "log"
        ) 
    else 
        DM = DHGLMMODELING(
            Model = "dispersion",
            Link = phiLink,
            LinPred = as.formula(phiFormula),
            Offset = phiOffsetVariable
        )
    
    incProgress(0.4, detail = paste("Loading..."))
    
    # ¦§ Variable Declaration 2 ####
    # Variable Declaration for dhglmfit
    if (is.null(input$rg_jglm_a_mord))
        mord = 1
    else
        mord = as.numeric(input$rg_jglm_a_mord)
    
    if (is.null(input$rg_jglm_a_dord))
        dord = 1
    else
        dord = as.numeric(input$rg_jglm_a_dord)
    
    if (is.null(input$rg_jglm_a_dispersion))
        Dmethod = "deviance"
    else
        Dmethod = input$rg_jglm_a_dispersion
    
    if (is.null(input$rg_jglm_a_reml))
        reml = TRUE
    else 
        reml = FALSE
    
    BinomialDen = NULL
    if (!is.null(input$rg_jglm_m_check_binomd) && input$rg_jglm_m_check_binomd) 
      BinomialDen = data2[[input$rg_jglm_m_binomd]]
    
    betafix = NULL
    
    # ¦§ Fitted Model ####
    
    fittedModel <<- dhglmfit(
        RespDist = isolate(input$rg_jglm_m_dist),
        BinomialDen = BinomialDen,
        DataMain = data2,
        MeanModel = MM,
        DispersionModel = DM,
        Dmethod = Dmethod,
        mord = mord,
        dord = dord,
        REML = reml
    )
    
    incProgress(0.3, detail = paste("Loading..."))

    fittedModel$option = list(
        RespDist = isolate(input$rg_jglm_m_dist),
        BinomialDen = BinomialDen,
        DataMain = data2,
        MeanModel = MM,
        DispersionModel = DM,
        Dmethod = Dmethod,
        mord = mord,
        dord = dord,
        REML = reml
    ) 

    # ¦§ Description ####
    
    modelDesc = matrix(c(
        input$rg_jglm_m_model,
        "constant",
        input$rg_jglm_m_link,
        phiLink,
        input$rg_jglm_m_dist,
        "gaussian"
    ),nrow=2, ncol=3)
    
    if (!is.null(input$rg_jglm_p_check_phi) && input$rg_jglm_p_check_phi) {
        modelDesc[2, 1] <- input$rg_jglm_p_model
        modelDesc[2, 2] <- input$rg_jglm_p_link
    }
    colnames(modelDesc) <- c("Model", "Link", "Dist")
    rownames(modelDesc) <- c("Mean", "Phi")
    fittedModel$modelDesc <- modelDesc
    
    # ¦§ Likelihood ####
    likeli_coeff <- matrix(fittedModel$likeli_coeff, ncol = 5)
    colnames(likeli_coeff) <- c("-2ML", "-2RL", "cAIC", "Scaled Deviance", "df")
    fittedModel$likeli_coeff = likeli_coeff
    
    # ¦§ Comparison Model ####
    jglmComparisonmodel1 = NULL
    jglmComparisonmodel2 = NULL
    
    if (!is.null(input$rg_jglm_m_check_comparison) && input$rg_jglm_m_check_comparison) {
        jglmComparisonmodel1 <- matrix(c(
            input$rg_jglm_m_model,
            input$rg_jglm_m_link,
            input$rg_jglm_m_dist,
            phiFormula
        ), ncol = 4)
      
        jglmComparisonmodel2 <- matrix(c(
            fittedModel$likeli_coeff[1,1],
            fittedModel$likeli_coeff[1,2],
            fittedModel$likeli_coeff[1,3],
            fittedModel$likeli_coeff[1,5],
            round(fittedModel$n.mean, digits=0),
            round(fittedModel$n.disp, digits=0)
        ),ncol = 6)
      
        jglmComparisonmodel1 <- rbind(g_rg_jglm_r_comparisonmodel_1, jglmComparisonmodel1)
        jglmComparisonmodel2 <- rbind(g_rg_jglm_r_comparisonmodel_2, jglmComparisonmodel2)
      
        colnames(jglmComparisonmodel1) <- c("Model", "Link", "Dist", "Phi Model")
        colnames(jglmComparisonmodel2) <- c("-2ML", "-2RL", "cAIC", "df", "n.mean", "n.disp")
        
        g_rg_jglm_r_comparisonmodel_1 <<- jglmComparisonmodel1
        g_rg_jglm_r_comparisonmodel_2 <<- jglmComparisonmodel2
        
        fittedModel$option$Comparisonmodel1 = jglmComparisonmodel1
        fittedModel$option$Comparisonmodel2 = jglmComparisonmodel2
    }
    
    # ¦¦ R Codes ####
    
    if (!is.null(input$rg_jglm_m_check_rcodes) && input$rg_jglm_m_check_rcodes) {
        Rcode1 <- NULL
        Rcode2 <- NULL
        Rcode3 <- NULL
        offsetRcode = NULL
        if (!is.null(meanOffsetVariable))
            offsetRcode <- paste0(", Offset = data[[\"", input$rg_jglm_m_offset, "\"]]")
      
        Rcode1 <- paste0(
            "MM <- DHGLMMODELING(Model = \"mean\", Link = \"",
            input$rg_jglm_m_link,
            "\", LinPred = ",
            meanFormula,
            offsetRcode,
            ")"
        )
        if(phiFormula == "constant") 
            Rcode2 <- "DM <- DHGLMMODELING(Model = \"dispersion\", Link = \"log\")"
        else
            Rcode2 <- paste0(
                "DM <- DHGLMMODELING(Model = \"dispersion\", Link = \"", phiLink, "\", LinPred = ",
                phiFormula,
                ")"
            )
      
        Rcode3 <- paste0(
            "dhglmfit(RespDist = \"", 
            input$rg_jglm_m_dist,
            "\", DataMain = data, MeanModel = MM, DispersionModel = DM)"
        )
        
        Rcodes <- matrix(c(Rcode1, Rcode2, Rcode3), ncol = 1)
        colnames(Rcodes) <- "Call"
        rownames(Rcodes) <- c("MM", "DM", "dhglmfit")
        fittedModel$option$Rcodes = Rcodes
    }
    
    setProgress(1, detail = "Finish")
    return(fittedModel)
    
    }) # End Progressbar
})

# JGLM UI (Accordion) ####

output$rg_jglm_r_accordion<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    jglmAccordion <- bs_accordion_sidebar(
        id = "rg_jglm_r_accordion",
        spec_side = c(width = 3, offset = 0),
        spec_main = c(width = 9, offset = 0)    
    )
    
    # ¦§ Mean ####
    
    jglmAccordion <- jglmAccordion %>%
    bs_append(
        title_side = "Mean",
        content_side = NULL,
        content_main = div(
            fluidRow(
                column(
                    9, 
                    h3(strong("Model for Mean"))
                ),
                column(
                    3,
                    style = "text-align:right; padding:15px"
                )
            ),
            uiOutput("rg_jglm_m_resp"), 
            uiOutput("rg_jglm_m_variable"), 
            fluidRow(
                column(
                    8,
                    uiOutput("rg_jglm_m_interaction") 
                ),
                column(
                    4,
                    uiOutput("rg_jglm_m_interactionappend"),
                    style = "text-align:right; padding:15px"
                )
            ),
            uiOutput("rg_jglm_m_dist"),
            uiOutput("rg_jglm_m_check_binomd"),
            uiOutput("rg_jglm_m_binomd"),
            uiOutput("rg_jglm_m_link"),
            
            hr(style = "border-color: #2C3E50;"),

            uiOutput("rg_jglm_m_check_nointercept"),
            uiOutput("rg_jglm_m_check_offset"),
            uiOutput("rg_jglm_m_offset"),
            
            # They don't work
            # uiOutput("rg_jglm_m_check_vif"),
            # uiOutput("rg_jglm_m_check_robustse"),
            # uiOutput("rg_jglm_m_check_confint"),
            # uiOutput("rg_jglm_m_check_exp"),
            
            hr(style = "border-color: #2C3E50;"),
            uiOutput("rg_jglm_m_check_comparison"),
            uiOutput("rg_jglm_m_check_rcodes")
        )
    )
    
    # ¦§ Phi ####
    
    jglmAccordion <- jglmAccordion %>%
    bs_append(
        title_side = "Phi",
        content_side = uiOutput("rg_jglm_p_check_phi"),
        content_main = div(
            h3(strong("Model for Phi")),
            
            uiOutput("rg_jglm_p_resp"),
            uiOutput("rg_jglm_p_variable"),
            fluidRow(
                column(
                    8,
                    uiOutput("rg_jglm_p_interaction")
                ),
                column(
                    4,
                    uiOutput("rg_jglm_p_interactionappend"),
                    style = "text-align:right; padding:15px"
                )
            ),

            uiOutput("rg_jglm_p_link"),
            
            uiOutput("rg_jglm_p_check_nointercept"),
            uiOutput("rg_jglm_p_check_offset"),
            uiOutput("rg_jglm_p_offset")
        )
    )
    
    # ¦¦ Additional Settings ####
    
    jglmAccordion <- jglmAccordion %>%
    bs_append(
        title_side = "Setting",
        content_side = NULL,
        content_main = div(
            h3(strong("Additional Settings")),
            uiOutput("rg_jglm_a_dispersion"),
            uiOutput("rg_jglm_a_reml"),
            uiOutput("rg_jglm_a_mord"),
            uiOutput("rg_jglm_a_dord"),
            h5(helpText("Order of Laplace Approximation for Likelihood(mean) and Restricted Likelihood(Dispersion)"))
        )
    )
    
    div(
        jglmAccordion,
        use_bs_tooltip(),
        use_bs_accordion_sidebar() # needs to be at end, for some reason
    )
})

# JGLM Components####
# ¦§ Mean ####

output$rg_jglm_m_model<-renderUI({ #algorithm not good
    input$rg_jglm_m_resp
    input$rg_jglm_m_variable$right
    input$rg_jglm_m_check_nointercept
    input$rg_jglm_m_binomd
    
    if(any(is.null(input$rg_jglm_m_resp), is.null(input$rg_jglm_m_variable)))
        return()
    
    responsePart = input$rg_jglm_m_resp
    variablePart = paste(input$rg_jglm_m_variable$right, collapse="+")
    meanModel = NULL
    
    if(!is.null(input$rg_jglm_m_check_nointercept) && !input$rg_jglm_m_check_nointercept) {
        if(variablePart == "")
            meanModel = paste(responsePart, "~", "1")
        else
            meanModel = paste(responsePart, "~", variablePart)
    }
    
    if(!is.null(input$rg_jglm_m_check_nointercept) && input$rg_jglm_m_check_nointercept) {
        if(variablePart == "")
            meanModel = paste(responsePart, "~", "-1")
        else
            meanModel = paste(responsePart, "~", "-1+", variablePart)
    }
   
    textAreaInput("rg_jglm_m_model", "Model for Mean", value = meanModel, height = "60px")
})

output$rg_jglm_m_resp <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue = names(select_if(data2, is.numeric))
    
    div(
        selectInput(
            "rg_jglm_m_resp",
            "Response Variable",
            choices = as.list(c("", nameValue)), 
            selected = NULL, 
            multiple = FALSE
        ),
        bsTooltip(
            "rg_jglm_m_resp", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$rg_jglm_m_variable <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$rg_jglm_m_interactionappend
    
    nameValue <- c(names(data2), g_rg_jglm_m_interaction)
    
    chooserInput("rg_jglm_m_variable", "Variable", "Selected", nameValue, c(), size = 15, multiple = TRUE, width = 150)
})

output$rg_jglm_m_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue <- names(data2)
    
    selectInput(
        "rg_jglm_m_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$rg_jglm_m_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    actionButton("rg_jglm_m_interactionappend", "Append")
})

observeEvent(input$rg_jglm_m_interactionappend, {
    if(length(input$rg_jglm_m_interaction) > 1)
        g_rg_jglm_m_interaction <<- c(g_rg_jglm_m_interaction, paste(input$rg_jglm_m_interaction, collapse=":"))
})

output$rg_jglm_m_dist<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    dist = c(
        "gaussian" = "gaussian",
        "binomial" = "binomial",
        "poisson" = "poisson",
        "gamma" = "gamma"
    )
    
    selectInput(
        "rg_jglm_m_dist",
        "Distribution", 
        choices = dist, 
        multiple = FALSE
    )
})

output$rg_jglm_m_check_binomd<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(any(is.null(input$rg_jglm_m_dist), input$rg_jglm_m_dist != "binomial"))
        return()
        
    checkboxInput("rg_jglm_m_check_binomd", "Binomial Denominator", value = FALSE)
})

output$rg_jglm_m_binomd<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$rg_jglm_m_check_binomd) || !(input$rg_jglm_m_check_binomd))
        return()
        
    nameValue = names(select_if(data2, is.numeric))
    
    selectInput(
        "rg_jglm_m_binomd",
        "Binomial Denominator",
        choices = as.list(c("", nameValue)),
        multiple = FALSE
    )
})

output$rg_jglm_m_link <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$rg_jglm_m_dist))
        return()
    
    selection = "identity"
    
    if (input$rg_jglm_m_dist == "binomial")
        selection = "logit"
    else if (input$rg_jglm_m_dist == "poisson")
        selection = "log"
    else if (input$rg_jglm_m_dist == "gamma")
        selection = "log"
    
    link = c(
        "identity" = "identity",
        "log" = "log",
        "logit" = "logit",
        "probit" = "probit",
        "cloglog" = "cloglog",
        "inverse" = "inverse"
    )
    
    selectInput(
        "rg_jglm_m_link",
        "Link Function",
        choices = link,
        selected = selection,
        multiple = FALSE
    )
})

output$rg_jglm_m_check_nointercept<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_jglm_m_check_nointercept", "No Intercept Model", value = FALSE)
})

output$rg_jglm_m_check_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_jglm_m_check_offset", "Offset Variable", value = FALSE)
})

output$rg_jglm_m_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    nameValue = names(select_if(data2, is.numeric))
    
    if(is.null(input$rg_jglm_m_check_offset) || !input$rg_jglm_m_check_offset)
        return()
       
    div(
        selectInput(
            "rg_jglm_m_offset",
            "Offset Variable",
            choices = as.list(c("", nameValue)),
            selected = NULL,
            multiple = FALSE
        ),
        bsTooltip(
            "rg_jglm_m_offset", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$rg_jglm_m_check_comparison<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_jglm_m_check_comparison", "model comparison", value = FALSE)
})

output$rg_jglm_m_check_rcodes<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_jglm_m_check_rcodes", "R Codes", value = FALSE)
})

# ¦§ Phi ####

output$rg_jglm_p_check_phi<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("rg_jglm_p_check_phi", strong("Use"), value = FALSE)
})

output$rg_jglm_p_model<-renderUI({ #algorithm not good
    input$rg_jglm_p_resp
    input$rg_jglm_p_variable$right
    input$rg_jglm_p_check_nointercept
    
    if(is.null(input$rg_jglm_p_check_phi) || !input$rg_jglm_p_check_phi)
        return()
    
    if(any(is.null(input$rg_jglm_p_resp), is.null(input$rg_jglm_p_variable)))
        return()
    
    responsePart = input$rg_jglm_p_resp
    variablePart = paste(input$rg_jglm_p_variable$right, collapse="+")
    phiModel = NULL
    
    if(!is.null(input$rg_jglm_p_check_nointercept) && !input$rg_jglm_p_check_nointercept) {
        if(variablePart == "")
            phiModel = paste(responsePart, "~", "1")
        else
            phiModel = paste(responsePart, "~", variablePart)
    }
    
    if(!is.null(input$rg_jglm_p_check_nointercept) && input$rg_jglm_p_check_nointercept) {
        if(variablePart == "")
            phiModel = paste(responsePart, "~", "-1")
        else
            phiModel = paste(responsePart, "~", "-1+", variablePart)
    }
    
    textAreaInput("rg_jglm_p_model", "Model for Phi", value = phiModel, height = "60px")
})

output$rg_jglm_p_resp <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$rg_jglm_p_check_phi) || !input$rg_jglm_p_check_phi)
        return()
    
    nameValue = c("phi" = "phi")
    
    selectInput(
        "rg_jglm_p_resp",
        "Residual Variance",
        choices = nameValue,
        selected = NULL, 
        multiple = FALSE
    )
})

output$rg_jglm_p_variable <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$rg_jglm_p_interactionappend
    
    if(is.null(input$rg_jglm_p_check_phi) || !input$rg_jglm_p_check_phi)
        return()
    
    nameValue <- c(names(data2), g_rg_jglm_p_interaction)
    
    chooserInput("rg_jglm_p_variable", "Variable", "Selected", nameValue, c(), size = 15, multiple = TRUE, width = 150)
})

output$rg_jglm_p_interaction <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    
    if(is.null(input$rg_jglm_p_check_phi) || !input$rg_jglm_p_check_phi)
        return()
    
    nameValue <- names(data2)
    
    selectInput(
        "rg_jglm_p_interaction",
        "Make Interaction Variable",
        choices = as.list(c("", nameValue)), 
        selected = NULL, 
        multiple = TRUE
    )
})

output$rg_jglm_p_interactionappend <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$rg_jglm_p_check_phi) || !input$rg_jglm_p_check_phi)
        return()
    
    actionButton("rg_jglm_p_interactionappend", "Append")
})

observeEvent(input$rg_jglm_p_interactionappend, {
    if(length(input$rg_jglm_p_interaction) > 1)
        g_rg_jglm_p_interaction <<- c(g_rg_jglm_p_interaction, paste(input$rg_jglm_p_interaction, collapse=":"))
})

output$rg_jglm_p_link <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$rg_jglm_p_check_phi) || !input$rg_jglm_p_check_phi)
        return()

    link = c(
        "log" = "log",
        "inverse" = "inverse",
        "identity" = "identity")
    
    selectInput(
        "rg_jglm_p_link",
        "Link Function",
        choices = link,
        multiple = FALSE
    )
})

output$rg_jglm_p_check_nointercept<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_p_check_phi
    
    if(is.null(input$rg_jglm_p_check_phi) || !input$rg_jglm_p_check_phi)
        return()
    
    checkboxInput("rg_jglm_p_check_nointercept", "No Intercept Model", value = FALSE)
})

output$rg_jglm_p_check_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_p_check_phi
    
    if(is.null(input$rg_jglm_p_check_phi) || !input$rg_jglm_p_check_phi)
        return()
    
    checkboxInput("rg_jglm_p_check_offset", "Offset Variable", value = FALSE)
})

output$rg_jglm_p_offset<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$dm_mv_run
    input$dm_ct_run
    input$dm_sr_run
    input$dm_md_run
    input$rg_jglm_p_check_phi
    
    if(is.null(input$rg_jglm_p_check_phi) || !input$rg_jglm_p_check_phi)
        return()
    
    if (is.null(input$rg_jglm_p_check_offset) || !input$rg_jglm_p_check_offset)
        return()
    
    nameValue = names(select_if(data2, is.numeric))
    
    div(
        selectInput(
            "rg_jglm_p_offset",
            "Offset Variable",
            choices = as.list(c("", nameValue)),
            selected = NULL,
            multiple = FALSE
        ),
        bsTooltip(
            "rg_jglm_p_offset", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

# ¦¦ Additional Settings ####

output$rg_jglm_a_dispersion<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    dispersionmethods = c(
        "deviance" = "deviance", 
        "Pearson" = "Pearson"
    )
    
    selectInput(
        "rg_jglm_a_dispersion",
        "Method of Estimating Phi", 
        choices = dispersionmethods, 
        multiple = FALSE
    )
})

output$rg_jglm_a_reml <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    selectInput(
        "rg_jglm_a_reml",
        "REML or ML",
        choice = c("REML" = TRUE, "ML" = FALSE),
        multiple = FALSE
    )
})

# output$rg_jglm_m_check_orthogonal <- renderUI({
#     input$file
#     input$resetData
#     
#     selectInput(
#         "rg_jglm_m_check_orthogonal",
#         "SE with orthogonality or not",
#         choice = c(TRUE, FALSE), 
#         selected = FALSE,
#         multiple = FALSE
#     )
# })

output$rg_jglm_a_mord <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    selectInput(
        "rg_jglm_a_mord",
        "Order for Mean Model",
        choice = c(1, 0),
        multiple = FALSE
    )
})

output$rg_jglm_a_dord <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    selectInput(
        "rg_jglm_a_dord",
        "Order for Dispersion Model",
        choice = c(1, 2),
        multiple = FALSE
    )
})

# JGLM Results ####

output$rg_jglm_r_mainsummary <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_jglm_results()$modelDesc))
        return()
  
    rg_jglm_results()$modelDesc
}, rownames = TRUE, bordered = TRUE, caption = "Model Description", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_jglm_m_coeff <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_jglm_results()$beta_coeff))
        return()
  
    rg_jglm_results()$beta_coeff
}, rownames = TRUE, bordered = TRUE, caption = "Estimate from Mean Model", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_jglm_p_coeff <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_jglm_results()$phi_coeff))
        return()
  
    rg_jglm_results()$phi_coeff
}, rownames = TRUE, bordered = TRUE, caption = "Estimate from Dispersion Model", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_jglm_r_likelihood <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_jglm_results()$likeli_coeff))
        return()
  
    rg_jglm_results()$likeli_coeff
}, rownames = FALSE, bordered = TRUE, caption = "Likelihood", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_jglm_r_comparisonmodel1<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_jglm_results()$option$Comparisonmodel1))
        return()
  
    rg_jglm_results()$option$Comparisonmodel1
}, rownames = TRUE, bordered = TRUE, caption = "Comparison of Mean & Phi Model", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_jglm_r_comparisonmodel2<-renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_jglm_results()$option$Comparisonmodel2))
        return()
  
    rg_jglm_results()$option$Comparisonmodel2
}, rownames = TRUE, bordered = TRUE, caption = "Comparison of Likelihood", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

output$rg_jglm_r_rcodes <- renderTable({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_g_run
    
    if (g_resetresult == FALSE || is.null(rg_jglm_results()$option$Rcodes))
        return()
  
    rg_jglm_results()$option$Rcodes
}, rownames = TRUE, bordered = TRUE, caption = "R Codes", spacing = "m",
caption.placement = getOption("xtable.caption.placement", "top"), digits = 5)

# JGLM Plot ####

output$rg_jglm_m_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_g_run
    
    if(g_resetresult == FALSE || is.null(rg_jglm_results()))
        return()
    
    local({
        output$rg_jglm_m_downloadplot1 <- downloadPlot({
            res<-rg_jglm_results()
            ggplotdhglm(res, type = "mean")
        })
    })
    
    res<-rg_jglm_results()
    ggplotdhglm(res, type = "mean")
})

output$rg_jglm_m_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_g_run

    if(g_resetresult == FALSE || is.null(rg_jglm_results()))
        return()
    
    div(
        h4("Model Checking Plots for Mean"),
        div(
            style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
            withSpinner(
                plotOutput("rg_jglm_m_plot1", height="600px"),
                type = 1,
                color = "#2c3e50",
                size = 1.2
            ),
            div(
                style = "position: absolute; left:0.5em; bottom: 0.5em;",
                dropdown(
                    downloadButton(outputId = "rg_jglm_m_downloadplot1", label = "Download Plot"),
                    size = "xs",
                    icon = icon("download", class = "opt"),
                    up = TRUE
                )
            )
        ),
        br()
    )
})

output$rg_jglm_p_plot1 <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_g_run
    
    if(g_resetresult == FALSE || is.null(rg_jglm_results()$resid_phi))
        return()
    
    local({
        output$rg_jglm_p_downloadplot1 <- downloadPlot({
            res<-rg_jglm_results()
            ggplotdhglm(res, type = "phi")
        })
    })
    
    res<-rg_jglm_results()
    ggplotdhglm(res, type = "phi")
})

output$rg_jglm_p_showplot1 <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_g_run

    if(g_resetresult == FALSE || is.null(rg_jglm_results()$resid_phi))
        return()
    
    div(
        h4("Model Checking Plots for Phi"),
        div(
            style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
            
            withSpinner(
                plotOutput("rg_jglm_p_plot1", height="600px"),
                type = 1,
                color = "#2c3e50",
                size = 1.2
            ),
            div(
                style = "position: absolute; left:0.5em; bottom: 0.5em;",
                dropdown(
                    downloadButton(outputId = "rg_jglm_p_downloadplot1", label = "Download Plot"),
                    size = "xs",
                    icon = icon("download", class = "opt"),
                    up = TRUE
                )
            )
        ),
        br()
    )
})

# JGLM Prediction ####

output$rg_jglm_r_check_95mu <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_g_run
    
    checkboxInput("rg_jglm_r_check_95mu", "95% Confidence Interval for mu", value = FALSE)
})

output$rg_jglm_r_box_95mu <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_g_run
    
    if(is.null(input$rg_jglm_r_check_95mu) || !input$rg_jglm_r_check_95mu)
        return()
    
    muLength = nrow(rg_jglm_results()$beta_coeff)
    if(rownames(rg_jglm_results()$beta_coeff)[1] == "(Intercept)")
        muLength = muLength - 1 
    
    textAreaInput("rg_jglm_r_box_95mu", "Type mean Covariates", value = paste(rep(1, muLength), collapse = ","))
})

output$rg_jglm_r_check_95phi <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_g_run
    
    if(is.null(rg_jglm_results()$modelDesc) || rg_jglm_results()$modelDesc[2,1] == "phi ~ 1"
       || rg_jglm_results()$modelDesc[2,1] == "phi~1")
        return()
    
    checkboxInput("rg_jglm_r_check_95phi", "95% Confidence Interval for phi", value = FALSE)
})

output$rg_jglm_r_box_95phi <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_g_run
    
    if(is.null(input$rg_jglm_r_check_95phi) || !input$rg_jglm_r_check_95phi)
        return()
    
    phiLength = nrow(rg_jglm_results()$phi_coeff)
    if(rownames(rg_jglm_results()$phi_coeff)[1] == "(Intercept)")
        phiLength = phiLength - 1 
    
    textAreaInput("rg_jglm_r_box_95phi", "Type phi Covariates", value = paste(rep(1, phiLength), collapse = ","))
})

output$rg_jglm_r_prediction <- renderPrint({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_g_run
    
    if(!is.null(input$rg_jglm_r_check_95mu) && input$rg_jglm_r_check_95mu && !is.null(rg_jglm_results()$beta_coeff)) {
        
        muLength = nrow(rg_jglm_results()$beta_coeff)
        if(rownames(rg_jglm_results()$beta_coeff)[1] == "(Intercept)")
            muLength = muLength - 1 
        
        muCovariates <- as.numeric(strsplit(input$rg_jglm_r_box_95mu, ",", fixed = TRUE)[[1]])
        if(muLength >= length(muCovariates))
            muCovariates <- c(muCovariates, rep(0,muLength - length(muCovariates)))
        else muCovariates <- muCovariates[1:muLength]
        
        if(rownames(rg_jglm_results()$beta_coeff)[1] == "(Intercept)")
            muCovariates<-c(1,muCovariates)
            
        muCovMatrix <- matrix(data = muCovariates, ncol=1)
        rownames(muCovMatrix) <- rownames(rg_jglm_results()$beta_coeff)
        print("Covariates for mean Predictor")
        print(t(muCovMatrix))
        muPrediction <- sum(as.vector(muCovMatrix) * rg_jglm_results()$beta_coeff[,1])
        muSE <- sqrt(t(muCovMatrix) %*% rg_jglm_results()$cov_mu %*% muCovMatrix)
        LL <- muPrediction - 1.96 * muSE
        UL <- muPrediction + 1.96 * muSE
        print(paste0("Lower Limit for Confidence Intervel of mean predictor : ", LL))
        print(paste0("Upper Limit for Confidence Intervel of mean predictor : ", UL))
    }
    
    if(!is.null(input$rg_jglm_r_check_95phi) && input$rg_jglm_r_check_95phi 
       && !is.null(rg_jglm_results()$phi_coeff)) {
        phiLength = nrow(rg_jglm_results()$phi_coeff)
        if(rownames(rg_jglm_results()$phi_coeff)[1] == "(Intercept)")
            phiLength = phiLength - 1
        
        phiCovariates <- as.numeric(strsplit(input$rg_jglm_r_box_95phi, ",", fixed = TRUE)[[1]])
        if(phiLength >= length(phiCovariates))
            phiCovariates <- c(phiCovariates, rep(0,phiLength - length(phiCovariates)))
        else phiCovariates <- phiCovariates[1:phiLength]
        
        if(rownames(rg_jglm_results()$phi_coeff)[1] == "(Intercept)")
            phiCovariates<-c(1,phiCovariates)
        
        phiCovMatrix <- matrix(data = phiCovariates, ncol=1)
        rownames(phiCovMatrix) <- rownames(rg_jglm_results()$phi_coeff)
        print("Covariates for phi Predictor")
        print(t(phiCovMatrix))
        phiPrediction <- sum(as.vector(phiCovMatrix) * rg_jglm_results()$phi_coeff[,1])
        phiSE <- sqrt(t(phiCovMatrix) %*% rg_jglm_results()$cov_phi %*% phiCovMatrix)
        LL <- phiPrediction - 1.96 * phiSE
        UL <- phiPrediction + 1.96 * phiSE
        print(paste0("Lower Limit for Confidence Intervel of phi predictor : ", LL))
        print(paste0("Upper Limit for Confidence Intervel of phi predictor : ", UL))
    }
})


output$rg_jglm_r_check_calculator <-renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_g_run
    
    checkboxInput("rg_jglm_r_check_calculator", "Calculator for data", value = FALSE)
})

output$rg_jglm_r_box_calculator <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_g_run
    
    if(is.null(input$rg_jglm_r_check_calculator) || !input$rg_jglm_r_check_calculator)
        return()
    
    textAreaInput("rg_jglm_r_box_calculator", "Type Arithmetic Form (e.g. sum(log(y))", value = paste0("sum(", input$rg_jglm_m_resp, ")"))
})

output$rg_jglm_r_calculator <- renderPrint({
    input$file
    input$resetData
    input$di_option_run
    input$rg_jglm_g_run
    
    if(is.null(input$rg_jglm_r_check_calculator) || !input$rg_jglm_r_check_calculator 
       || is.null(input$rg_jglm_r_box_calculator))
        return()
    
    calcEquation <- input$rg_jglm_r_box_calculator
    
    quo_var <- rlang::parse_expr(quo_name(enquo(calcEquation)))
    calcResults <- dplyr::mutate(data2, new1 = !!quo_var)
    
    print("Calculator Results")
    print(calcResults[1,ncol(calcResults)])
})
