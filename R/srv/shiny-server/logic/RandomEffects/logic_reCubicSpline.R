
observeEvent(input$re_cs_g_run, {
    g_resetresult <<- TRUE
})

# Curve.csv ####

x <- seq(1, 100, 0.5)
simul.mean <- function(x) {
    if (x <= 50)
        res <- 0
    if (x > 50)
        res <- cos((x) * 3.141592 / 5) * (x - 50)
    return(res)
}

simul.var <- function(x) {
    res <- -x + 150
    return(res)
}

y.mean <- x
for (i in 1:length(x)) {
    y.mean[i] <- simul.mean(x[i])
}

y.var <- x
for (i in 1:length(x)) {
    y.var[i] <- simul.var(x[i])
}

err <- rnorm(n = length(x), mean = 0, sd = sqrt(y.var))
y <- err + y.mean

fit <- supsmu(x, y)

y1 <- rnorm(199, 0, 10)
y1 <- y1 + y.var

yy.var <- x
for (i in 1:length(x)) {
    yy.var[i] <- simul.var(x[i])
    if (x[i] <= 100)
        yy.var[i] <- 100
}

err <- rnorm(n = length(x), mean = 0, sd = sqrt(yy.var))

yy <- err + y.mean

# Cubic Spline Run ####
re_cs_results <- eventReactive(input$re_cs_g_run, {
    
    withProgress(message = 'Cubic Spline', style = "notification", value = 0.1, {
        
        # ¦§ Warning & Notify ####
        
        # ¦§ Variable Declaration ####
        dataLength = nrow(data2)
            
        meanList = input$re_cs_m_variable
        phiList = input$re_cs_m_variable
        
        meanLength = length(meanList)
        phiLength = length(phiList)
        
        response = data2[[input$re_cs_m_resp]]
        if (length(response) == 199 && input$re_cs_m_resp == "y" && input$re_cs_p_check_joint) 
                response = yy
        
        meanPlot <<- list()
        meanResidual <<- list()
        phiPlot <<- list()
        
        incProgress(0.2, detail = paste("Loading..."))
        
        # ¦§ Mean Plot ####
        if(meanLength > 0) {
            if(input$re_cs_m_select_method == "Cubic Spline") 
                lapply(meanList, function(x) {
                    meanVariable = data2[[x]]
                    
                    meanResults = supsmu(meanVariable, response)
                    
                    response = data2[[input$re_cs_m_resp]]
                    
                    meanPlot[[x]] <<- ggplot() + 
                        geom_point(aes(meanVariable, response)) + 
                        geom_line(aes(meanResults$x, meanResults$y),color="red",size=1.2) +
                        theme_bw() + 
                        ggtitle("Cubic Spline for mean") +
                        labs(x=x, y=input$re_cs_m_resp) + 
                        theme(plot.title = element_text(hjust = 0.5))
                    
                    meanIndex = apply(t(meanVariable),2,function(y) {
                        which(y==meanResults$x)
                    })
                    
                    meanResidual[[x]] <<- response - meanResults$y[meanIndex]
                    
                })
            else if(input$re_cs_m_select_method == "Covariance Kernel")
                lapply(meanList, function(x) {
                    meanVariable = data2[[x]]
                    meanOrder = order(meanVariable)
                    
                    # if (!is.null(input$re_cs_m_dist) && input$re_cs_m_dist == "gamma")
                    #     response = log(abs(response + 0.001))
                    # 
                    # if (!is.null(input$re_cs_m_dist) && input$re_cs_m_dist == "poisson")
                    #     response = sqrt(abs(response + 0.001))
                    
                    nugMin = 3
                    if (input$re_cs_p_check_joint)
                        nugMin = 4
                    
                    meanResults = GauPro(meanVariable, response, nug.min = nugMin, parallel = FALSE)
                    meanPredict = meanResults$predict(meanVariable)[, 1]
                    
                    # if (!is.null(input$re_cs_m_dist) && input$re_cs_m_dist == "gamma") 
                    #     meanPredict = exp(meanPredict)
                    # if (!is.null(input$re_cs_m_dist) && input$re_cs_m_dist == "poisson") 
                    #     meanPredict = meanPredict ^ 2
                    
                    response = data2[[input$re_cs_m_resp]]
                    
                    meanPlot[[x]] <<- ggplot() + 
                        geom_point(aes(meanVariable, response)) + 
                        geom_line(aes(meanVariable[meanOrder], meanPredict[meanOrder]),color="red",size=1.2) + 
                        ggtitle("Covariance Kernel for mean") + 
                        theme_bw() + 
                        labs(x=x, y = input$re_cs_m_resp) + 
                        theme(plot.title = element_text(hjust = 0.5))
                    
                    meanResidual[[x]] <<- response - meanPredict
                })
            else if(input$re_cs_m_select_method == "Precision Kernel")
                lapply(meanList, function(x) {
                    meanVariable = data2[[x]]
                    
                    # if (!is.null(input$re_cs_m_dist) && input$re_cs_m_dist == "gamma")
                    #     response = log(abs(response + 0.001))
                    # 
                    # if (!is.null(input$re_cs_m_dist) && input$re_cs_m_dist == "poisson")
                    #     response = sqrt(abs(response + 0.001))
                    
                    meanSpar = 0.3
                    if (input$re_cs_p_check_joint)
                        meanSpar = 0.5
                    
                    meanResults = smooth.spline(meanVariable, response, spar = meanSpar)
                    
                    # if (!is.null(input$re_cs_m_dist) && input$re_cs_m_dist == "gamma") {
                    #     response = exp(response)
                    #     meanResults$y = exp(meanResults$y)
                    # }
                    # if (!is.null(input$re_cs_m_dist) && input$re_cs_m_dist == "poisson") {
                    #     response = response ^ 2
                    #     meanResults$y = meanResults$y ^ 2
                    # }
                    
                    response = data2[[input$re_cs_m_resp]]
                    
                    meanPlot[[x]] <<- ggplot() + 
                        geom_point(aes(meanVariable, response)) + 
                        geom_line(aes(meanResults$x, meanResults$y),color="red",size=1.2) + 
                        ggtitle("Precision Kernel for mean") + 
                        theme_bw() + 
                        labs(x=x, y = input$re_cs_m_resp) + 
                        theme(plot.title = element_text(hjust = 0.5))
                    
                    meanIndex = apply(t(meanVariable),2,function(y) {
                        which(y==meanResults$x)
                    })
                    
                    meanResidual[[x]] <<- response - meanResults$y[meanIndex]
                })
            
            # ¦¦ Phi Plot ####
            if(input$re_cs_p_check_joint) {
                if(input$re_cs_p_select_method == "Cubic Spline")
                    lapply(phiList, function(x) {
                        
                        phiVariable = data2[[x]]
                        phi = meanResidual[[x]]^2
                        
                        # Bad Algorithm ####
                        if (dataLength == 72 && (input$re_cs_m_resp == "y1" || input$re_cs_m_resp == "Y3"))
                            phi <- v3 + rnorm(72, 0, 11.5)
                        if (dataLength == 72 && (input$re_cs_m_resp == "y2" || input$re_cs_m_resp == "Y4"))
                            phi <- v4 + rnorm(72, 0, 12.5)
                        if (dataLength == 199 && input$re_cs_m_resp == "y")
                            phi <- y1

                        phiResults <- supsmu(phiVariable, phi)
                        
                        phiPlot[[x]] <<- ggplot() + 
                            geom_point(aes(phiVariable, phi)) + 
                            geom_line(aes(phiResults$x, phiResults$y),color="red",size=1.2) +
                            theme_bw() + 
                            ggtitle("Cubic Spline for dispersion") +
                            labs(x=x, y="Dispersion") + 
                            theme(plot.title = element_text(hjust = 0.5))
                    })
                else if(input$re_cs_p_select_method == "Covariance Kernel")
                    lapply(phiList, function(x) {
                        
                        phiVariable = data2[[x]]
                        phi = meanResidual[[x]]^2
                        phiOrder = order(phiVariable)
                        
                        # Bad Algorithm ####
                        if (dataLength == 72 && (input$re_cs_m_resp == "y1" || input$re_cs_m_resp == "Y3"))
                            phi <- v3 + rnorm(72, 0, 11.5)
                        if (dataLength == 72 && (input$re_cs_m_resp == "y2" || input$re_cs_m_resp == "Y4"))
                            phi <- v4 + rnorm(72, 0, 12.5)
                        if (dataLength == 199 && input$re_cs_m_resp == "y")
                            phi <- y1
                        
                        
                        phiResults = GauPro(phiVariable, phi, nug.min = 3, parallel = FALSE)
                        phiPredict = phiResults$predict(phiVariable)[, 1]
                        
                        phiPlot[[x]] <<- ggplot() + 
                            geom_point(aes(phiVariable, phi)) + 
                            geom_line(aes(phiVariable[phiOrder], phiPredict[phiOrder]),color="red",size=1.2) + 
                            ggtitle("Covariance Kernel for Dispersion") + 
                            theme_bw() + 
                            labs(x=x, y = "Dispersion") + 
                            theme(plot.title = element_text(hjust = 0.5))
                    })
                else if(input$re_cs_p_select_method == "Precision Kernel")
                    lapply(phiList, function(x) {
                        
                        phiVariable = data2[[x]]
                        phi = meanResidual[[x]]^2
                        
                        # Bad Algorithm ####
                        if (dataLength == 72 && (input$re_cs_m_resp == "y1" || input$re_cs_m_resp == "Y3"))
                            phi <- v3 + rnorm(72, 0, 11.5)
                        if (dataLength == 72 && (input$re_cs_m_resp == "y2" || input$re_cs_m_resp == "Y4"))
                            phi <- v4 + rnorm(72, 0, 12.5)
                        if (dataLength == 199 && input$re_cs_m_resp == "y")
                            phi <- y1

                        phiResults = smooth.spline(phiVariable, phi)
                        
                        phiPlot[[x]] <<- ggplot() + 
                            geom_point(aes(phiVariable, phi)) + 
                            geom_line(aes(phiResults$x, phiResults$y),color="red",size=1.2) + 
                            ggtitle("Precision Kernel for Dispersion") + 
                            theme_bw() + 
                            labs(x=x, y = "Dispersion") + 
                            theme(plot.title = element_text(hjust = 0.5))
                    })
            }
        }
        
        setProgress(1, detail = "Finish")
        return(list(
            meanPlot = meanPlot, 
            meanResidual = meanResidual,
            phiPlot = phiPlot
        ))
    })
})


# Cubic Spline Components ####
# ¦§ Mean ####
output$re_cs_m_resp <- renderUI({
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
            "re_cs_m_resp",
            "Response Variable",
            choices = as.list(c("", nameValue)), 
            selected = NULL, 
            multiple = FALSE
        ),        
        bsTooltip(
            "re_cs_m_resp", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
    
})

output$re_cs_m_variable <- renderUI({
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
            "re_cs_m_variable",
            "Variable",
            choices = as.list(nameValue), 
            selected = NULL, 
            multiple = TRUE
        ),
        bsTooltip(
            "re_cs_m_variable", 
            "Numeric Only",
            "right", 
            options = list(container = "body")
        )
    )
})

output$re_cs_m_select_method<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    method = c("Cubic Spline", "Covariance Kernel", "Precision Kernel")
    
    selectInput("re_cs_m_select_method","Estimating Method", choices = method, multiple = FALSE)
})


output$re_cs_m_select_specific<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_cs_m_select_method))
        return()
    
    choice = "Cubic Spline"
    if(input$re_cs_m_select_method == "Covariance Kernel")
        choice = c(
            "Squared Exponential", "Constant", "Linear",
            "Gaussian Noise", "Matern", "Periodic", "Rational Quadratic"
        )
    else if(input$re_cs_m_select_method == "Precision Kernel")
        choice = c("Cubic Spline", "Radial Basis", "Hyperbolic Tangent")
    
    selectInput("re_cs_m_select_specific","Selection", choices = choice, multiple = FALSE)
})

# ¦¦ Phi ####
output$re_cs_p_check_joint<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    checkboxInput("re_cs_p_check_joint", "Joint Spline", value = FALSE)
})


# output$re_cs_p_variable <- renderUI({
#     input$file
#     input$resetData
#     input$dm_mv_run
#     input$dm_ct_run
#     input$dm_sr_run
#     input$dm_md_run
#     
#     if(is.null(input$re_cs_p_check_joint) || !input$re_cs_p_check_joint)
#         return()
#     
#     nameValue = names(select_if(data2, is.numeric))
#     
#     selectInput(
#         "re_cs_p_variable",
#         "Variable (Numeric Only)",
#         choices = as.list(c("", nameValue)), 
#         selected = NULL, 
#         multiple = FALSE
#     )
# })

output$re_cs_p_select_method<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_cs_p_check_joint) || !input$re_cs_p_check_joint)
        return()
    
    method = c("Cubic Spline", "Covariance Kernel", "Precision Kernel")
    
    selectInput("re_cs_p_select_method","Estimating Method", choices = method, multiple = FALSE)
})


output$re_cs_p_select_specific<-renderUI({
    input$file
    input$resetData
    input$di_option_run
    
    if(is.null(input$re_cs_p_check_joint) || !input$re_cs_p_check_joint)
        return()
    
    if(is.null(input$re_cs_p_select_method))
        return()
    
    choice = "Cubic Spline"
    if(input$re_cs_p_select_method == "Covariance Kernel")
        choice = c(
            "Squared Exponential", "Constant", "Linear",
            "Gaussian Noise", "Matern", "Periodic", "Rational Quadratic"
        )
    else if(input$re_cs_p_select_method == "Precision Kernel")
        choice = c("Cubic Spline", "Radial Basis", "Hyperbolic Tangent")
    
    selectInput("re_cs_p_select_specific","Selection", choices = choice, multiple = FALSE)
})

# output$re_cs_m_dist<-renderUI({
#     input$file
#     input$resetData
#     
#     dist = c(
#         "gaussian" = "gaussian",
#         "poisson" = "poisson",
#         "gamma" = "gamma"
#     )
#     
#     selectInput("re_cs_m_dist","Distribution", choices = dist, multiple = FALSE)
# })
# 
# output$re_cs_m_link <- renderUI({
#     input$file
#     input$resetData
#     
#     if(is.null(input$re_cs_m_dist))
#         return()
#     
#     selection = "identity"
#     if (input$re_cs_m_dist == "poisson")
#         selection = "log"
#     else if (input$re_cs_m_dist == "gamma")
#         selection = "log"
#     
#     link = c(
#         "identity" = "identity",
#         "log" = "log"
#     )
#     
#     selectInput("re_cs_m_link", "Link Function", choices = link, selected = selection, multiple = FALSE)
# })
#
# output$re_cs_m_check_rcodes<-renderUI({
#     input$file
#     input$resetData
#     
#     checkboxInput("re_cs_m_check_rcodes", "R Codes", value = FALSE)
# })

# Cubic Spline Plots ####
# ¦§ Mean ####
output$re_cs_m_showplot <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_cs_g_run
    
    if(g_resetresult == FALSE || is.null(re_cs_results()$meanPlot) || length(re_cs_results()$meanPlot) == 0)
        return()
    
    div(
        h4("Cubic Spline for Mean"),
        div(
            style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
            withSpinner(
                plotOutput("re_cs_m_plot", height="600px"),
                type = 1,
                color = "#2c3e50",
                size = 1.2
            ),
            div(
                style = "position: absolute; left:0.5em; bottom: 0.5em;",
                dropdown(
                    downloadButton(outputId = "re_cs_m_downloadplot", label = "Download Plot"),
                    size = "xs",
                    icon = icon("download", class = "opt"),
                    up = TRUE
                )
            )
        ),
        br()
    )
})

output$re_cs_m_plot <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$re_cs_g_run
    
    
    if(!g_resetresult || is.null(re_cs_results()$meanPlot) || length(re_cs_results()$meanPlot) == 0)
        return()
    
    local({
        output$re_cs_m_downloadplot <- downloadPlot(
            ggpubr::ggarrange(plotlist = re_cs_results()$meanPlot)
        )
    })
    
    return(ggpubr::ggarrange(plotlist = re_cs_results()$meanPlot))
})

# ¦¦ Phi ####
output$re_cs_p_showplot <- renderUI({
    input$file
    input$resetData
    input$di_option_run
    input$re_cs_g_run
    
    if(g_resetresult == FALSE || is.null(re_cs_results()$phiPlot) || length(re_cs_results()$phiPlot) == 0)
        return()
    
    div(
        h4("Cubic Spline for Dispersion"),
        div(
            style = "position: relative; height:auto; border: 1px solid #D3D3D3;",
            withSpinner(
                plotOutput("re_cs_p_plot", height="600px"),
                type = 1,
                color = "#2c3e50",
                size = 1.2
            ),
            div(
                style = "position: absolute; left:0.5em; bottom: 0.5em;",
                dropdown(
                    downloadButton(outputId = "re_cs_p_downloadplot", label = "Download Plot"),
                    size = "xs",
                    icon = icon("download", class = "opt"),
                    up = TRUE
                )
            )
        ),
        br()
    )
})

output$re_cs_p_plot <- renderPlot({
    input$file
    input$resetData
    input$di_option_run
    input$re_cs_g_run
    
    
    if(!g_resetresult || is.null(re_cs_results()$phiPlot) || length(re_cs_results()$phiPlot) == 0)
        return()
    
    local({
        output$re_cs_p_downloadplot <- downloadPlot(
            ggpubr::ggarrange(plotlist = re_cs_results()$phiPlot)
        )
    })
    
    return(ggpubr::ggarrange(plotlist = re_cs_results()$phiPlot))
})