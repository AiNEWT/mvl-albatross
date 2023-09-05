ggplotglm <- function (OUTPUT) {
    
    regression <- OUTPUT
    
    x <- mu <- regression$fitted.values
    y <- StudentResidual <- rstandard(regression)
    fit <- data.frame(supsmu(x,y))
    
    plot1 <- ggplot(data.frame(cbind(x, y)), aes(x, y)) + 
      geom_point() + 
      geom_line(data = fit, aes(x, y)) + 
      labs(x = "Scaled Fitted Values", y = "Studentized Residual") + 
      ggtitle("Residuals vs Fitted") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    
    # Plot 2
    y <- abs(StudentResidual)
    fit <- data.frame(supsmu(x,y))
    plot2 <- ggplot(data.frame(cbind(x, y)), aes(x, y)) + 
      geom_point() + 
      geom_line(data = fit, aes(x, y)) + 
      labs(x = "Scaled Fitted Values", y = "|Studentized Residual|") + 
      ggtitle("|Residuals| vs Fitted") + 
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5))
    
    # Plot 3
    y <- StudentResidual <- rstandard(regression)
    plot3 <- ggplot(data = data.frame(y), aes(sample=y)) + 
      stat_qq() + 
      stat_qq_line(linetype = 2) + 
      labs(x = "Theoretical Quantiles", y = "Standardized Residual") + 
      ggtitle("Normal Probability Plot") + 
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5))
      
    # Plot 4
    plot4 <- ggplot(data.frame(y), aes(y)) + geom_histogram(binwidth=0.5) +
      labs(x = "Studentized Residual", y = "Frequency") +
      ggtitle("Histogram of Studentized Residual") +
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5))
    
    Plots <- ggpubr::ggarrange(plot1, plot2, plot3, plot4)
    
    return(Plots)
}

ggplotdhglm <- function (OUTPUT, type = "mean", spatial = "") {
    if (type == "ci" || type == "ci1" || type == "ci2" || type == "ci3" || type == "ci4" || type == "phici" || type == "lambdaci") {
        plotNumber = 0
        if (type == "ci" || type == "phici" || type == "lambdaci")
            plotNumber = 0
        else if (type == "ci1")
            plotNumber = 1
        else if (type == "ci2")
            plotNumber = 2
        else if (type == "ci3")
            plotNumber = 3
        else if (type == "ci4")
            plotNumber = 4
        
        res <- OUTPUT
        if (spatial == "Matern")
            q <- length(res$matern$v_h)
        else if (spatial == "IAR" || spatial == "MRF")
            q <- length(fittedModel$v_h)
        else
            q <- res$q
        
        qcum <- c(0,cumsum(q))
        if (plotNumber > 0) {
            v_h <- res$v_h[(qcum[plotNumber]+1):(qcum[plotNumber+1])]
            SE <- sqrt(v_h ^ 2 / 2)
        } else {
            if (spatial == "Matern")
                v_h <- res$matern$v_h
            else if (type == "phici")
                v_h <- res$phi_v_h
            else if (type == "lambdaci")
                v_h <- res$lambda_v_h
            else
                v_h <- res$v_h
            SE <- sqrt(v_h ^ 2 / 4) ### I don't know why SE is sqrt(v_h ^ 2 / 4), and why divided by 4 or 2?
        }
        lb <- v_h - 1.96 * SE # lower bound
        ub <- v_h + 1.96 * SE # upper bound
        v_h <- sort(v_h)
        lb <- sort(lb)
        ub <- sort(ub)
        if (plotNumber > 0) randCount <- 1:q[plotNumber]
        else randCount <- 1:q

        # Plot
        ciTable <- cbind(lb, ub, randCount)
        ciPlot <- ggplot(data.frame(cbind(v_h, randCount)), aes(randCount, v_h)) + 
            geom_point() + 
            geom_line() +
            geom_segment(data = ciTable, aes(x = randCount, y = lb, xend = randCount, yend = ub)) +
            labs(x = "Number", y = "Estimated Random Effects") + 
            ggtitle("95% CI for Random Effect") +
            theme_bw() + 
            theme(plot.title = element_text(hjust = 0.5))
        
        return(ciPlot)
    }
    
    mu <- NULL
    StudentResidual <- NULL
    q <- OUTPUT$q
    
    if (spatial == "Matern") {
        v_h <- OUTPUT$matern$v_h / sd(OUTPUT$matern$v_h)
        if (type == "v") {
            StudentResidual <- v_h
        } else if (type == "mean") {
            StudentResidual <- as.vector(OUTPUT$matern$std_dev_res)
            
            for (i in 1:length(StudentResidual))
                if (is.nan(StudentResidual[i]))
                    StudentResidual[i] <- 0
            
            mu <- as.vector(OUTPUT$matern$fv)
            StudentResidual <- StudentResidual
        }
    } else if (spatial == "IAR" || spatial == "MRF") {
        v_h <- OUTPUT$v_h / sd(OUTPUT$v_h)
        if (type == "v") {
            StudentResidual <- v_h
        } else if (type == "mean") {
            StudentResidual <- as.vector(OUTPUT$mean_residual)
            
            for (i in 1:length(StudentResidual))
                if (is.nan(StudentResidual[i]))
                    StudentResidual[i] <- 0
                
            mu <- as.vector(OUTPUT$mu)
            StudentResidual <- StudentResidual
        }
    } else if (type=="mean") {
	    mu<-OUTPUT[[7]]
	    StudentResidual<-OUTPUT[[1]]
    } else if (type=="phi") {
	    mu<-OUTPUT[[4]]
	    StudentResidual<-OUTPUT[[3]]
    } else if (type=="lambda") {
	    mu<-OUTPUT[[6]]
	    StudentResidual<-OUTPUT[[5]]
    } else if (type=="lambda1") {
	    mu<-OUTPUT[[6]][1:60]
	    StudentResidual<-OUTPUT[[5]][1:60]
    } else if (type=="lambda2") {
	    mu<-OUTPUT[[6]][61:120]
	    StudentResidual<-OUTPUT[[5]][61:120]
    } else if (type=="v") {
	    mu<-OUTPUT[[24]][1:60]
	    StudentResidual<-OUTPUT$sv_h
    } else if (type=="v1") {
        temp1<-q[1]
	    mu<-OUTPUT[[24]][1:temp1]
	    StudentResidual<-OUTPUT$sv_h[1:temp1]
    } else if (type=="v2") {
        temp2<-q[1]+1
        temp3<-q[1]+q[2]
	    mu<-OUTPUT[[24]][temp2:temp3]
	    StudentResidual<-OUTPUT$sv_h[temp2:temp3]
    } else if (type=="v3") {
        temp2<-q[1]+q[2]+1
        temp3<-q[1]+q[2]+q[3]
	    mu<-OUTPUT[[24]][temp2:temp3]
	    StudentResidual<-OUTPUT$sv_h[temp2:temp3]
    } else if (type=="v4") {
        temp2<-q[1]+q[2]+q[3]+1
        temp3<-q[1]+q[2]+q[3]+q[4]
	    mu<-OUTPUT[[24]][temp2:temp3]
	    StudentResidual<-OUTPUT$sv_h[temp2:temp3]
    } else if (type=="phiv") {
	    mu<-OUTPUT[[4]]
	    StudentResidual<-OUTPUT$phi_sv_h
    } else if (type=="lambdav") {
	    mu<-OUTPUT[[6]]
	    StudentResidual<-OUTPUT$lambda_sv_h
    } else if (type=="alpha") {
	    mu<-OUTPUT[[32]]
	    StudentResidual<-OUTPUT[[31]]
    }
    
    x <- mu
    y <- StudentResidual
    
    # Plot 3
    plot3 <- ggplot(data = data.frame(y), aes(sample=y)) + 
      stat_qq() + 
      stat_qq_line(linetype = 2) + 
      labs(x = "Theoretical Quantiles", y = "Standardized Residual") + 
      ggtitle("Normal Probability Plot") + 
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5))
              
    # Plot 4
    plot4 <- ggplot(data.frame(y), aes(y)) + geom_histogram(binwidth=0.5) +
      labs(x = "Studentized Residual", y = "Frequency") +
      ggtitle("Histogram of Studentized Residual") +
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5))
         
#    if ((type == "mean" || type == "phi" || type == "lambda") && nlevels(as.factor(x)) > 1) {
    if ((type == "mean" || type == "phi") && nlevels(as.factor(x)) > 1) {
        fit<- data.frame(supsmu(x, y))
        
        plot1 <- ggplot(data.frame(cbind(x, y)), aes(x, y)) + 
          geom_point() + 
          geom_line(data = fit, aes(x, y)) + 
          labs(x = "Scaled Fitted Values", y = "Studentized Residual") + 
          ggtitle("Residuals vs Fitted") +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5))
        
        # Plot 2
        absy <- abs(StudentResidual)
        absfit <- data.frame(supsmu(x,absy))
        plot2 <- ggplot(data.frame(cbind(x, absy)), aes(x, absy)) + 
          geom_point() + 
          geom_line(data = absfit, aes(x, y)) + 
          labs(x = "Scaled Fitted Values", y = "|Studentized Residual|") + 
          ggtitle("|Residuals| vs Fitted") + 
          theme_bw() + 
          theme(plot.title = element_text(hjust = 0.5))
        return(ggpubr::ggarrange(plot1, plot2, plot3, plot4))
    } else {
        return(ggpubr::ggarrange(plot3, plot4))
    }
}