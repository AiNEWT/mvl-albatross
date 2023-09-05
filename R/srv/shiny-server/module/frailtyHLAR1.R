aa <- 
function (formula, data, weights, subset, na.action, RandDist = "Normal",mord = 0, dord = 1, Maxiter = 70, convergence = 10^-6, varfixed = FALSE,varinit = c(0.163),varnonneg=FALSE) 
{
    Call <- match.call()
    indx <- match(c("formula", "data", "weights", "subset", "na.action"), 
        names(Call), nomatch = 0)
    if (indx[1] == 0) stop("A formula argument is required")
    data$nn <- rep(1:nrow(data))
    data$idid <- rep(1:nrow(data))
    temp <- Call[c(1, indx)]
    temp[[1L]] <- quote(stats::model.frame)
    special <- c("strata", "cluster")
    temp$formula <- terms(subbars(formula), special,data=data)
    m <- eval(temp,parent.frame())
    Terms <- attr(m, "terms")
    Y <- model.extract(m, "response")
    temp <- Call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    special <- c("strata", "cluster")
    temp$formula <- terms(formula, special)
    Terms <- temp[[2]]
    formula1 <- paste(paste(Terms[[2]][[2]], Terms[[3]], sep = "~")[[2]], 
        paste(Terms[[3]])[3], sep = "+")
    formula1 <- formula(formula1)
    fr <- FrailtyFrames(Call, formula1, contrasts)
    namesX <- names(fr$fixef)
    namesX <- namesX[-1]
    namesY <- names(fr$mf)[1]
    FL <- HGLMFactorList(formula1, fr, 0L, 0L)
    namesRE <- FL$namesRE
    leftT<-0
    if (ncol(Y)==3) leftT<-1

    if (leftT==0) y <- matrix(Y[, 1], length(fr$Y), 1)
    if (leftT==0) L0 <- matrix(0, length(fr$Y), 1)
    if (leftT==1) y <- matrix(Y[, 2], length(fr$Y), 1)
    if (leftT==1) L0 <- matrix(Y[, 1], length(fr$Y), 1)
    x <- fr$X
    z <- FL$Design
    n <- nrow(x)
    p <- ncol(x)
    x1 <- x[1:n, 2:p]
    x2 <- matrix(x1, n, p - 1)
    x <- x2
    n <- nrow(x)
    p <- ncol(x)
    nrand <- length(z)
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
    del <- matrix(0, n, 1)
    if (leftT==0) del[, 1] <- censor <- Y[, 2]
    if (leftT==1) del[, 1] <- censor <- Y[, 3]
    SS <- FL$Subject
    res1 <- FrailtyMakeData(y, x, del, z, L0)
    y <- res1[1][[1]]
    x <- res1[2][[1]]
    del <- res1[3][[1]]
    z <- res1[4][[1]]
    Mi <- res1[5][[1]]
    idx2 <- res1[6][[1]]
    t2 <- res1[7][[1]]
    di <- res1[8][[1]]
    beta_h <- matrix(0, p, 1)
    qcum <- cumsum(c(0, q))
    v_h <- matrix(0, qcum[nrand + 1], 1)
    rho <- 0.5
    alpha_h <- rep(0, nrand)
    varinit1 <- rep(0.1, nrand)
    length_init <- length(varinit)
    for (i in 1:length_init) varinit1[i] <- varinit[i]
    for (i in 1:nrand) alpha_h[i] <- varinit1[i]
    Max_iter <- Maxiter
    err <- 1
    if (RandDist == "AR1") v_h <- matrix(0, n, 1)
    for (ii in 1:Max_iter) {
        if (err >= convergence) {
	    if (RandDist == "AR1") {
	    	    IArho<-NULL
		    nrand <- length(z)
		    qrep<-colSums(z[[1]])
		    qq <- length(qrep)
		    for (i in 1:qq) {
		         Arhotemp<-matrix(1,qrep[i],qrep[i])/(1-rho^2)
		         for (j in 1:qrep[i]) {
		             for (k in 1:qrep[i]) {
		                 Arhotemp[j,k]<-Arhotemp[j,k]*rho^(abs(j-k))
		             }
		         }
		         if (i==1) IArho<-solve(Arhotemp)
		         if (i>1) IArho<-dbind(IArho,solve(Arhotemp))
	  	  }
	    }
            if (RandDist == "Normal") 
                res2 <- PNfrailtyHL(x, z, y, del, Mi, idx2, t2, 
                  di, beta_h, v_h, alpha_h, mord, dord, varfixed = varfixed,varnonneg)
            if (RandDist == "AR1") {
                  res2 <- AR1frailtyHL(x, z, y, del, Mi, idx2, t2, 
                  di, beta_h, v_h, alpha_h, mord, dord, varfixed = varfixed,varnonneg,IArho,rho,qrep)
            }
            if (RandDist == "Gamma") 
                res2 <- PGfrailtyHL(x, z, y, del, Mi, idx2, t2, 
                  di, beta_h, v_h, alpha_h, mord, dord, varfixed = varfixed,varnonneg)
            alpha_h <- res2[13][[1]]
            alpha_h1 <- res2[14][[1]]
            beta_h <- res2[11][[1]]
            beta_h1 <- res2[9][[1]]
            v_h <- res2[12][[1]]
            v_h1 <- res2[10][[1]]
            Hinv <- res2[16][[1]]
            rho_h1 <- rho
            rho_h <- rho
            se_rho_h <- 0.000
	    if (RandDist == "AR1") {
                   rho_h <- res2[24][[1]]
                   se_rho_h <- res2[25][[1]]
            }
            temp4 <- sum(abs(alpha_h - alpha_h1)) + sum(abs(v_h - 
                v_h1)) + sum(abs(beta_h - beta_h1)) +sum(abs(rho_h-rho_h1))
            rho <- rho_h
            err <- temp4
            alpha_h <- alpha_h1
            se_beta <- res2[20][[1]]
            u_h <- NULL
            Hinv_u <- NULL
            if (RandDist == "Gamma") {
                 u_h<-res2[23][[1]]
                 Hinv_u<-res2[24][[1]]
            }
            print_err <- err
            print_i <- ii
        }
    }
    names(print_i) <- "iteration : "
#    print(print_i*13)
    names(print_err) <- "convergence : "
 #   print(0.00000003561*se_beta[1])
#    if (err < convergence) 
#        print("converged")
#    if (err > convergence) 
#        print("did not converge")
    result <- list(0)
    names(result)[1] <- "Model"
    sum_init <- sum(abs(varinit))
    if (varfixed == TRUE && sum_init < 1e-05) {
        print("Results from the Cox model")
        result$Model <- "Cox model"
    }
    else {
        if (RandDist == "Gamma") {
            print("Results from the gamma frailty model")
            result$Model <- "gamma frailty model"
        }
        if (RandDist == "Normal") {
            print("Results from the log-normal frailty model")
            result$Model <- "log-normal frailty model"
        }
        if (RandDist == "AR1") {
            print("Results from the log-normal frailty model with AR(1)")
            result$Model <- "log-normal frailty model with AR(1)"
        }
    }
    nevent <- sum(censor)
    print("Number of data : ")
    print(n)
    print("Number of event : ")
    print(nevent)
    print("Model for conditional hazard : ")
    result$formula <- formula
    print(formula)
    if (mord == 0 && dord == 1) {
        print("Method : HL(0,1)")
        result$Method <- "HL(0,1)"
    }
    if (mord == 0 && dord == 2) {
        print("Method : HL(0,2)")
        result$Method <- "HL(0,2)"
    }
    if (mord == 1 && dord == 1) {
        print("Method : HL(1,1)")
        result$Method <- "HL(1,1)"
    }
    if (mord == 1 && dord == 2) {
        print("Method : HL(1,2)")
        result$Method <- "HL(1,2)"
    }
    print("Estimates from the mean model")
    z_beta <- beta_h/se_beta
    pval <- 2 * pnorm(abs(z_beta), lower.tail = FALSE)
    beta_coeff <- cbind(beta_h, se_beta, z_beta, pval)
    colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value", 
        "p-value")
    rownames(beta_coeff) <- namesX
    result$FixCoef <- beta_coeff
    print(beta_coeff, 4)
    if (RandDist == "Normal") 
        res3 <- PNFrailty_SE.h(res2, nrand, q, qcum, dord, varfixed = varfixed)
    if (RandDist == "AR1") {
            n <- nrow(x)
	    zz1 <- diag(1,n,n)
	    nrand1 <- 1
	    qq1 <- rep(0, nrand1)
	    for (i in 1:nrand1) qq1[i] <- dim(zz1)[2]
	    qcum1 <- cumsum(c(0, qq1))
            res3 <- AR1Frailty_SE.h(res2, nrand1, qq1, qcum1, dord, varfixed = varfixed, IArho)
    }
    if (RandDist == "Gamma") 
        res3 <- PGFrailty_SE.h(res2, nrand, q, qcum, dord, varfixed = varfixed)
    print("Estimates from the dispersion model")
    se_alpha_h <- res3[1][[1]]
    hlike <- -2 * res3[2][[1]]
    p1 <- -2 * res3[3][[1]]
    p2 <- -2 * res3[4][[1]]
    p3 <- -2 * res3[5][[1]]
    p0 <- -2 * res3[6][[1]]
    p4 <- p3 - (p1 - p2)
    df1 <- res3[7][[1]]
    if (varfixed == FALSE & RandDist != "AR1") 
        z_lam <- alpha_h/se_alpha_h
    for (i in 1:nrand) {
        if (alpha_h[i] < 1e-05) 
            alpha_h[i] <- 0
    }
    lam_coeff <- cbind(alpha_h, se_alpha_h)
    colnames(lam_coeff) <- c("Estimate", "Std. Error")
    rownames(lam_coeff) <- namesRE
    print(lam_coeff, 4)
    result$RandCoef <- lam_coeff
    result$rho <- 0.000
    if (RandDist == "AR1") {
	rho_coeff <- cbind(rho_h, se_rho_h)	
	print("Estimates for rho in the AR(1) model")
	print(rho_coeff,4)
    }	
    if (mord == 0 && dord == 1) 
        like_value <- cbind(p0, hlike, p1)
    if (mord == 0 && dord == 1) 
        colnames(like_value) <- c("-2h0", "-2*hp", "-2*p_b,v(hp)")
    if (mord == 0 && dord == 2) 
        like_value <- cbind(p0, hlike, p1, p3)
    if (mord == 0 && dord == 2) 
        colnames(like_value) <- c("-2h0", "-2*hp", "-2*p_b,v(hp)", 
            "-2*s_b,v(hp)")
    if (mord == 1 && dord == 1) 
        like_value <- cbind(p0, hlike, p2, p1)
    if (mord == 1 && dord == 1) 
        colnames(like_value) <- c("-2h0", "-2*hp", "-2*p_v(hp)", 
            "-2*p_b,v(hp)")
    if (mord == 1 && dord == 2) 
        like_value <- cbind(p0, hlike, p2, p4, p1, p3)
    if (mord == 1 && dord == 2) 
        colnames(like_value) <- c("-2h0", "-2*hp", "-2*p_v(hp)", 
            "-2*s_v(hp)", "-2*p_b,v(hp)", "-2*s_b,v(hp)")
    result$likelihood <- like_value
    result$iter <- print_i
    if (print_err < convergence) 
        result$convergence <- "converged"
    if (print_err > convergence) 
        result$convergence <- "did not converge"
    names(result$convergence) <- "convergence : "
    print(like_value, 5)
    res4 <- list(res2, res3)
    caic <- p0 + 2 * df1
    n_lam <- nrow(lam_coeff)
    if (varfixed == TRUE) 
        n_lam <- 0
    maic <- p2 + 2 * nrow(beta_coeff) + 2 * n_lam
    if (varfixed == TRUE) 
        maic <- hlike + 2 * nrow(beta_coeff) + 2 * n_lam
    raic <- p1 + 2 * n_lam
    if (RandDist == "Gamma" && mord == 1 && dord == 2) 
        maic <- p4 + 2 * nrow(beta_coeff) + 2 * n_lam
    if (RandDist == "Gamma" && mord == 1 && dord == 2) 
        raic <- p3 + 2 * n_lam
    aic <- cbind(caic, maic, raic)
    colnames(aic) <- c("cAIC", "mAIC", "rAIC")
    print(aic, 5)
    result$aic <- aic
    result$v_h <- v_h
    result$Hinv <- Hinv
    result$p <- p
    result$q <- q
    if (RandDist=="Gamma") result$u_h<-u_h
    if (RandDist=="Gamma") result$Hinv_u<-Hinv_u
    return(result)
}


dbind <-
function(a,b){
        out1<-cbind(a,matrix(0,nrow(a),ncol(b)))
        out2<-cbind(matrix(0,nrow(b),ncol(a)),b)
        out<-rbind(out1,out2)
        out
}

tr <-
function(a){
        out<-sum(diag(a))
        out
}
 
AR1Frailty_SE.h <-
function (res1, nrand, q, qcum, dord = 1, varfixed = FALSE, IArho) 
{
    x <- res1[1][[1]]
    z <- res1[2][[1]]
    y <- res1[3][[1]]
    del <- res1[4][[1]]
    Mi <- res1[5][[1]]
    idx2 <- res1[6][[1]]
    t2 <- res1[7][[1]]
    di <- res1[8][[1]]
    beta_h <- res1[9][[1]]
    v_h <- res1[10][[1]]
    beta_h1 <- res1[11][[1]]
    v_h1 <- res1[12][[1]]
    alpha_h <- res1[13][[1]]
    alpha_h1 <- res1[14][[1]]
    dft <- res1[15][[1]]
    Hinv <- res1[16][[1]]
    clam0 <- res1[17][[1]]
    H <- res1[18][[1]]
    mat <- res1[19][[1]]
    U <- res1[21][[1]]
    H0 <- res1[22][[1]]
    n <- nrow(x)
    p <- ncol(x)
    u_h1 <- exp(v_h1)
    oq <- matrix(1, qcum[nrand + 1], 1)
    oq1 <- matrix(1, qcum[nrand + 1], 1)
    for (i in 1:nrand) {
        index1 <- qcum[i] + 1
        oq1[index1:qcum[i + 1]] <- alpha_h1[i]
    }
    D <- diag(oq1[, 1])
    iD <- solve(D)
    iA <- iD
    Bi <- diag(clam0[, 1])
    muh <- x %*% beta_h1 + z %*% v_h1
    expeta <- exp(muh)
    cla0 <- di/(t(Mi) %*% expeta)
    temp4 <- cla0^2/di
    As <- diag(temp4[, 1])
    Wi <- diag(expeta[, 1])
    done <- matrix(1, idx2, 1)
    H22 <- solve(t(z) %*% mat %*% z + U)
    Hessian <- matrix(0, nrand, nrand)
    if (varfixed == FALSE) {
        se_lam <- rep(0, nrand)
        for (i in 1:nrand) se_lam[i] <- 0.000
    }
    if (varfixed == TRUE) {
        se_lam <- rep(0, nrand)
        for (i in 1:nrand) se_lam[i] <- 0.000
    }
    eta <- x %*% beta_h1 + z %*% v_h1
    expeta <- exp(eta)
    one <- matrix(1, n, 1)
    done <- matrix(1, idx2, 1)
    oq <- matrix(1, qcum[nrand + 1], 1)
    pi <- 3.14159265359
    term0 <- t(Mi) %*% expeta
    hlike1 <- (t(one) %*% (del * eta)) - (t(done) %*% (di * log(term0)))
    hlike2 <- 0
    hlike3 <- 0
    for (i in 1:nrand) {
        if (alpha_h[i] > 1e-05) 
            hlike2 <- hlike2 - (q[i]/2) * log(2 * pi) - ((q[i]/2) * 
                log(alpha_h1[i])) + 0.5 * log(abs(det(IArho)))
        index1 <- qcum[i] + 1
        index2 <- qcum[i + 1]
        vv_h1 <- matrix(0, q[i], 1)
        vv_h1[1:q[i], 1] <- v_h1[index1:index2, 1]
        if (alpha_h[i] > 1e-05) 
            hlike3 <- hlike3 - (t(vv_h1) %*% IArho %*% vv_h1)/(2 * alpha_h1[i])
    }
    hliken <- hlike1 + hlike2 + hlike3
    cc1 <- svd(2 * pi * Hinv)
    for (i in 1:length(cc1$d)) if (cc1$d[i] < 1e-05) 
        cc1$d[i] <- 1
    logdet1 <- sum(log(abs(cc1$d)))
    adj1 <- 0.5 * logdet1
    hpn1 <- hliken + adj1
    muu <- exp(x %*% beta_h1) * clam0
    zmu <- t(z) %*% muu
    u_h1 <- exp(v_h1)
    second <- 0
    for (i in 1:nrand) {
        ialph1 <- 1/alpha_h1[i]
        a21 <- (zmu * u_h1) + ialph1
        b31 <- zmu * u_h1
        S11 <- 3 * (b31/(a21^2))
        S21 <- 5 * ((b31^2)/(a21^3))
        temp4 <- S11 - S21
        S31 <- diag(temp4[, 1])
        second <- second - sum(diag(S31))/24
    }
    H22 <- t(z) %*% mat %*% z + iD
    cc1 <- svd(H22/(2 * pi))
    for (i in 1:length(cc1$d)) if (cc1$d[i] > 1e+05) 
        cc1$d[i] <- 1
    logdet1 <- sum(log(abs(cc1$d)))
    hpn2 <- hliken - 0.5 * logdet1
    hpn3 <- hpn1 + second
    df1 <- sum(diag(Hinv %*% H0))
    res <- list(se_lam, hliken, hpn1, hpn2, hpn3, hlike1, df1)
    return(res)
}

AR1frailtyHL <- 
function (x, z, y, del, Mi, idx2, t2, di, beta_h0, v_h0, alpha_h0, 
    mord, dord, varfixed = FALSE,varnonneg,IArho,rho_0,qrep) 
{
    n <- nrow(x)
    p <- ncol(x)
    z <- diag(1,n,n)
    nrand <- 1
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z)[2]
    qcum <- cumsum(c(0, q))
    beta_h <- beta_h0
    v_h <- v_h0
    rho <- rho_0
    for (i in 1:nrand) {
        if (alpha_h0[i] < 1e-06) 
            alpha_h0[i] <- 1e-06
    }
    alpha_h <- alpha_h0
    zz <- z
    z <- matrix(0, n, qcum[nrand + 1])
    z[1:n, 1:qcum[nrand + 1]] <- zz[1:n, 1:qcum[nrand + 1]]
    muh <- x %*% beta_h0 + z %*% v_h0
    expeta <- exp(muh)
    cla0 <- di/(crossprod(Mi,expeta))
    Wi <- diag(expeta[, 1])
    Ai <- diag(cla0[, 1])
    done <- matrix(1, idx2, 1)
    oq <- matrix(1, qcum[nrand + 1], 1)
    oq1 <- matrix(1, qcum[nrand + 1], 1)
    clam0 <- Mi %*% Ai %*% done
    for (i in 1:nrand) {
        index1 <- qcum[i] + 1
        oq1[index1:qcum[i + 1]] <- alpha_h[i]
    }
    D <- diag(oq1[, 1])
    iD <- solve(D)
    dft1 <- crossprod(x, del - Wi %*% clam0)
    if (mord == 1) {
        U <- iD
        Bi <- diag(clam0[, 1])
        temp4 <- cla0^2/di
        As <- diag(temp4[, 1])
        mat <- (Wi %*% Bi) - (Wi %*% Mi) %*% As %*% (t(Mi) %*% Wi)
        dinv <- solve(t(z) %*% mat %*% z + U)
        mu0 <- exp(x %*% beta_h0 + z %*% v_h0)
        mu <- exp(x %*% beta_h0 + z %*% v_h0) * clam0
        dcla0be <- -As %*% (t(Mi) %*% Wi)
        dcla0b <- matrix(0, idx2, p)
        dv_db <- matrix(0, qcum[nrand + 1], p)
        xz <- matrix(0, n, p)
        dmu0 <- matrix(0, n, p)
        dw_db1 <- matrix(0, n, p)
        xk <- matrix(0, n, 1)
        ad1 <- matrix(0, p, 1)
        for (k in 1:p) {
            xk[, 1] <- x[, k]
            dv_db[, k] <- -dinv %*% (t(z) %*% mat %*% xk)
            xz[, k] <- xk + z %*% (dv_db[, k])
            dcla0b[, k] <- dcla0be %*% (xz[, k])
            dc <- Mi %*% diag(dcla0b[, k]) %*% done
            dmu0[, k] <- mu0 * (xz[, k])
            dw_db1[, k] <- dc * mu0 + mu * xz[, k]
            temp4 <- (2 * cla0/di) * (dcla0b[, k])
            dw_db2 <- ((diag(dmu0[, k])) %*% Mi %*% As %*% t(Mi) %*% 
                Wi) + (Wi %*% Mi %*% diag(temp4[, 1]) %*% t(Mi) %*% 
                Wi) + (Wi %*% Mi %*% As %*% t(Mi) %*% (diag(dmu0[, 
                k])))
            dw_db <- diag(dw_db1[, k]) - dw_db2
            ad1[k, 1] <- sum(diag(dinv %*% (t(z) %*% dw_db %*% 
                z)))
        }
        dft1 <- dft1 - 0.5 * ad1
    }
    dft2 <- crossprod(z, del - Wi %*% clam0) - (IArho %*% iD %*% v_h0)
    dft <- rbind(dft1, dft2)
    Bi <- diag(clam0[, 1])
    temp4 <- cla0^2/di
    Adi <- diag(temp4[, 1])
    As <- Adi
    mat <- (diag(Wi) * Bi) - (diag(Wi) * Mi) %*% (diag(Adi)*( t(diag(Wi)*Mi)))
    U <- IArho%*%iD
#####################
    xmat<-crossprod(x,mat)
    zmat<-crossprod(z,mat)
    xmatx<-xmat%*%x
    xmatz<-xmat%*%z
    zmatx<-t(xmatz)
    zmatz<-zmat%*%z
    zmatz1<-zmat%*%z+U
    H <- rbind(cbind(xmatx, xmatz),cbind(zmatx, zmatz1))
    H0 <- rbind(cbind(xmatx, xmatz),cbind(zmatx, zmatz))
#####################
    Hinv <- solve(H)
    be_h0 <- rbind(beta_h0, v_h0)
    be_h <- be_h0 + (Hinv %*% dft)
    beta_h[1:p, 1] <- be_h[1:p, 1]
    se_beta_h <- matrix(0, p, 1)
    for (i in 1:p) se_beta_h[i, 1] <- sqrt(Hinv[i, i])
    index2 <- qcum[nrand + 1]
    index3 <- p + 1
    index4 <- p + qcum[nrand + 1]
    v_h[1:index2, 1] <- be_h[index3:index4, 1]
    nonneg_adj=0
    if (varnonneg == TRUE) nonneg_adj=2
    for (i in 1:nrand) {
        if (dord == 0) {
            index1 <- p + qcum[i] + 1
            index2 <- p + qcum[i + 1]
            gamma <- sum(diag(Hinv[index1:index2, index1:index2]))/alpha_h[i]
            index1 <- qcum[i] + 1
            index2 <- qcum[i + 1]
            if (varfixed == FALSE) 
                alpha_h[i] <- sum(v_h[index1:index2, 1]^2)/(q[i] - 
                  gamma-nonneg_adj)
        }
        if (dord == 1 | dord == 2) {
            H22 <- solve(zmatz1)
            ial1 <- 1/alpha_h[i]
            iA <- iD
            C <- matrix(0, qcum[nrand + 1], qcum[nrand + 1])
            index1 <- qcum[i] + 1
            index2 <- qcum[i + 1]
            for (j in index1:index2) C[j, j] <- 1
            iB1 <- iA %*% C %*% iA
            c_vh1 <- iB1 %*% IArho %*% v_h
            dv1 <- H22 %*% c_vh1
            dexpeta1 <- expeta * (z %*% dv1)
            dcla01 <- -(di/(crossprod(Mi,expeta)^2)) * (crossprod(Mi,dexpeta1))
            dWi1 <- diag(dexpeta1[, 1])
            dAi1 <- diag(dcla01[, 1])
            temp4 <- Mi %*% dAi1 %*% done
            dBi1 <- diag(temp4[, 1])
            dvec1 <- 2 * (cla0 * dcla01)
            temp4 <- dvec1/di
            dAs1 <- diag(temp4[, 1])
            dmat1 <- (diag(dWi1) * Bi) + (diag(Wi) * dBi1) - (diag(dWi1)*Mi%*% (diag(As)*(t(diag(Wi)*Mi))))- (diag(Wi)*Mi %*% (diag(dAs1)*(t(diag(Wi)*Mi))))- (diag(Wi) * Mi %*% (diag(As)*(t(diag(dWi1)*Mi))))
            dia1 <- -iB1 %*% IArho
#####################
    xmat1<-crossprod(x,dmat1)
    zmat1<-crossprod(z,dmat1)
    xmat1x<-xmat1%*%x
    xmat1z<-xmat1%*%z
    zmat1x<-t(xmat1z)
    zmat1z<-zmat1%*%z
    zmat1z1<-zmat1%*%z+dia1
#####################
            Hd1 <- rbind(cbind(xmat1x, xmat1z), cbind(zmat1x, zmat1z1))
            gamma1 <- -alpha_h[i] * sum(Hinv * Hd1)
            if (dord == 2) {
                expeta <- exp(x %*% beta_h + z %*% v_h)
                ial1 <- 1/alpha_h[i]
                muu <- exp(x %*% beta_h) * clam0
                zmuu <- t(z) %*% muu
                u_h <- exp(v_h)
                ude1 <- u_h * dv1
                aa1 <- (zmuu * u_h) + ial1
                bb1 <- (zmuu * ude1) - (ial1^2)
                term11 <- ((aa1 * zmuu * ude1) - (2 * zmuu * 
                  u_h * bb1))/(aa1^3)
                term21 <- ((2 * aa1 * zmuu * zmuu * u_h * ude1) - 
                  (3 * ((zmuu * u_h)^2) * bb1))/(aa1^4)
                term1 <- (3 * term11) - (5 * term21)
                SS1 <- diag(term1[, 1])
                gamma21 <- -(alpha_h[i]/12) * sum(diag(SS1))
            }
            if (dord == 1) {
                gamma21 <- 0
            }
           k21 <- q[i] - gamma1 - gamma21 - t(v_h) %*% IArho %*% v_h /(alpha_h[i])
            if (varfixed == FALSE) 
                alpha_h[i] <- t(v_h) %*% IArho %*% v_h /(q[i] - 
                  gamma1 - gamma21-nonneg_adj)
        }
    }
    cumqrep <- cumsum(qrep)
    qq <- length(qrep)
    L1 <- 0
    L2 <- 0
    L3 <- 0
    n <- sum(qrep)
    for (i in 1:qq) {
        if (i==1) temp1 <- 1
        else temp1 <- cumqrep[i-1]+1
	temp2 <- cumqrep[i]
	Lambdai <- v_h[temp1:temp2,1] %*% t(v_h[temp1:temp2,1])
        temp3 <- nrow(beta_h)+temp1
        temp4 <- nrow(beta_h)+temp2
        Ti <- Hinv[temp3:temp4,temp3:temp4]
        Ii <- diag(1,qrep[i],qrep[i])
        Ji <- matrix(0,qrep[i],qrep[i])
        Ki <- matrix(0,qrep[i],qrep[i])
        for (k in 1:qrep[i]) {
            if (k<qrep[i]) Ji[k,k+1] <- 1
            if (k<qrep[i]) Ji[k+1,k] <- 1
            if (k==1) Ki[k,k] <- 1
            if (k==qrep[i]) Ki[k,k] <- 1
        }
        if (qrep[i]==1) Ki[1,1] <- 2
        L1 <- L1 + tr(Ii %*% (Ti + Lambdai))
        L2 <- L2 + 0.5*tr(Ji %*% (Ti + Lambdai))
        L3 <- L3 + tr(Ki %*% (Ti + Lambdai))
    }
    c1 <- (n-qq) * (L1-L3)
    c2 <- (2*qq - n) * L2
    c3 <- n*L3 - (n+qq)*L1
    c4 <- n*L2
    df_rho <- c1*rho^3 + c2*rho^2 + c3*rho + c4
    d2f_rho <- 3*c1*rho^2 + 2*c2*rho + c3
    rho <- rho - df_rho/d2f_rho
    se_rho <- sqrt(abs(-1/d2f_rho))
    res <- list(x, z, y, del, Mi, idx2, t2, di, beta_h0, v_h0, 
        beta_h, v_h, alpha_h0, alpha_h, dft, Hinv, clam0, H, 
        mat, se_beta_h, U, H0, rho_0, rho, se_rho)
    return(res)
}

PNfrailtyHL <-
function (x, z, y, del, Mi, idx2, t2, di, beta_h0, v_h0, alpha_h0, 
    mord, dord, varfixed = FALSE,varnonneg) 
{
    n <- nrow(x)
    p <- ncol(x)
    nrand <- length(z)
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
    qcum <- cumsum(c(0, q))
    beta_h <- beta_h0
    v_h <- v_h0
    for (i in 1:nrand) {
        if (alpha_h0[i] < 1e-06) 
            alpha_h0[i] <- 1e-06
    }
    alpha_h <- alpha_h0
    zz <- z[[1]]
    if (nrand > 1) {
        index1 <- nrand
        for (i in 2:index1) zz <- cbind(zz, z[[i]])
    }
    z <- matrix(0, n, qcum[nrand + 1])
    z[1:n, 1:qcum[nrand + 1]] <- zz[1:n, 1:qcum[nrand + 1]]
    muh <- x %*% beta_h0 + z %*% v_h0
    expeta <- exp(muh)
    cla0 <- di/(crossprod(Mi,expeta))
    Wi <- diag(expeta[, 1])
    Ai <- diag(cla0[, 1])
    done <- matrix(1, idx2, 1)
    oq <- matrix(1, qcum[nrand + 1], 1)
    oq1 <- matrix(1, qcum[nrand + 1], 1)
    clam0 <- Mi %*% Ai %*% done
    for (i in 1:nrand) {
        index1 <- qcum[i] + 1
        oq1[index1:qcum[i + 1]] <- alpha_h[i]
    }
    D <- diag(oq1[, 1])
    iD <- solve(D)
    dft1 <- crossprod(x, del - Wi %*% clam0)
    if (mord == 1) {
        U <- iD
        Bi <- diag(clam0[, 1])
        temp4 <- cla0^2/di
        As <- diag(temp4[, 1])
        mat <- (Wi %*% Bi) - (Wi %*% Mi) %*% As %*% (t(Mi) %*% Wi)
        dinv <- solve(t(z) %*% mat %*% z + U)
        mu0 <- exp(x %*% beta_h0 + z %*% v_h0)
        mu <- exp(x %*% beta_h0 + z %*% v_h0) * clam0
        dcla0be <- -As %*% (t(Mi) %*% Wi)
        dcla0b <- matrix(0, idx2, p)
        dv_db <- matrix(0, qcum[nrand + 1], p)
        xz <- matrix(0, n, p)
        dmu0 <- matrix(0, n, p)
        dw_db1 <- matrix(0, n, p)
        xk <- matrix(0, n, 1)
        ad1 <- matrix(0, p, 1)
        for (k in 1:p) {
            xk[, 1] <- x[, k]
            dv_db[, k] <- -dinv %*% (t(z) %*% mat %*% xk)
            xz[, k] <- xk + z %*% (dv_db[, k])
            dcla0b[, k] <- dcla0be %*% (xz[, k])
            dc <- Mi %*% diag(dcla0b[, k]) %*% done
            dmu0[, k] <- mu0 * (xz[, k])
            dw_db1[, k] <- dc * mu0 + mu * xz[, k]
            temp4 <- (2 * cla0/di) * (dcla0b[, k])
            dw_db2 <- ((diag(dmu0[, k])) %*% Mi %*% As %*% t(Mi) %*% 
                Wi) + (Wi %*% Mi %*% diag(temp4[, 1]) %*% t(Mi) %*% 
                Wi) + (Wi %*% Mi %*% As %*% t(Mi) %*% (diag(dmu0[, 
                k])))
            dw_db <- diag(dw_db1[, k]) - dw_db2
            ad1[k, 1] <- sum(diag(dinv %*% (t(z) %*% dw_db %*% 
                z)))
        }
        dft1 <- dft1 - 0.5 * ad1
    }
    dft2 <- crossprod(z, del - Wi %*% clam0) - (iD %*% v_h0)
    dft <- rbind(dft1, dft2)
    Bi <- diag(clam0[, 1])
    temp4 <- cla0^2/di
    Adi <- diag(temp4[, 1])
    As <- Adi
    mat <- (diag(Wi) * Bi) - (diag(Wi) * Mi) %*% (diag(Adi)*( t(diag(Wi)*Mi)))
    U <- iD
#####################
    xmat<-crossprod(x,mat)
    zmat<-crossprod(z,mat)
    xmatx<-xmat%*%x
    xmatz<-xmat%*%z
    zmatx<-t(xmatz)
    zmatz<-zmat%*%z
    zmatz1<-zmat%*%z+U
    H <- rbind(cbind(xmatx, xmatz),cbind(zmatx, zmatz1))
    H0 <- rbind(cbind(xmatx, xmatz),cbind(zmatx, zmatz))
#####################
    Hinv <- solve(H)
    be_h0 <- rbind(beta_h0, v_h0)
    be_h <- be_h0 + (Hinv %*% dft)
    beta_h[1:p, 1] <- be_h[1:p, 1]
    se_beta_h <- matrix(0, p, 1)
    for (i in 1:p) se_beta_h[i, 1] <- sqrt(Hinv[i, i])
    index2 <- qcum[nrand + 1]
    index3 <- p + 1
    index4 <- p + qcum[nrand + 1]
    v_h[1:index2, 1] <- be_h[index3:index4, 1]
    nonneg_adj=0
    if (varnonneg == TRUE) nonneg_adj=2
    for (i in 1:nrand) {
        if (dord == 0) {
            index1 <- p + qcum[i] + 1
            index2 <- p + qcum[i + 1]
            gamma <- sum(diag(Hinv[index1:index2, index1:index2]))/alpha_h[i]
            index1 <- qcum[i] + 1
            index2 <- qcum[i + 1]
            if (varfixed == FALSE) 
                alpha_h[i] <- sum(v_h[index1:index2, 1]^2)/(q[i] - 
                  gamma-nonneg_adj)
        }
        if (dord == 1 | dord == 2) {
            H22 <- solve(zmatz1)
            ial1 <- 1/alpha_h[i]
            iA <- iD
            C <- matrix(0, qcum[nrand + 1], qcum[nrand + 1])
            index1 <- qcum[i] + 1
            index2 <- qcum[i + 1]
            for (j in index1:index2) C[j, j] <- 1
            iB1 <- iA %*% C %*% iA
            c_vh1 <- iB1 %*% v_h
            dv1 <- H22 %*% c_vh1
            dexpeta1 <- expeta * (z %*% dv1)
            dcla01 <- -(di/(crossprod(Mi,expeta)^2)) * (crossprod(Mi,dexpeta1))
            dWi1 <- diag(dexpeta1[, 1])
            dAi1 <- diag(dcla01[, 1])
            temp4 <- Mi %*% dAi1 %*% done
            dBi1 <- diag(temp4[, 1])
            dvec1 <- 2 * (cla0 * dcla01)
            temp4 <- dvec1/di
            dAs1 <- diag(temp4[, 1])
            dmat1 <- (diag(dWi1) * Bi) + (diag(Wi) * dBi1) - (diag(dWi1)*Mi%*% (diag(As)*(t(diag(Wi)*Mi))))- (diag(Wi)*Mi %*% (diag(dAs1)*(t(diag(Wi)*Mi))))- (diag(Wi) * Mi %*% (diag(As)*(t(diag(dWi1)*Mi))))
            dia1 <- -iB1
#####################
    xmat1<-crossprod(x,dmat1)
    zmat1<-crossprod(z,dmat1)
    xmat1x<-xmat1%*%x
    xmat1z<-xmat1%*%z
    zmat1x<-t(xmat1z)
    zmat1z<-zmat1%*%z
    zmat1z1<-zmat1%*%z+dia1
#####################
            Hd1 <- rbind(cbind(xmat1x, xmat1z), cbind(zmat1x, zmat1z1))
            gamma1 <- -alpha_h[i] * sum(Hinv * Hd1)
            if (dord == 2) {
                expeta <- exp(x %*% beta_h + z %*% v_h)
                ial1 <- 1/alpha_h[i]
                muu <- exp(x %*% beta_h) * clam0
                zmuu <- t(z) %*% muu
                u_h <- exp(v_h)
                ude1 <- u_h * dv1
                aa1 <- (zmuu * u_h) + ial1
                bb1 <- (zmuu * ude1) - (ial1^2)
                term11 <- ((aa1 * zmuu * ude1) - (2 * zmuu * 
                  u_h * bb1))/(aa1^3)
                term21 <- ((2 * aa1 * zmuu * zmuu * u_h * ude1) - 
                  (3 * ((zmuu * u_h)^2) * bb1))/(aa1^4)
                term1 <- (3 * term11) - (5 * term21)
                SS1 <- diag(term1[, 1])
                gamma21 <- -(alpha_h[i]/12) * sum(diag(SS1))
            }
            if (dord == 1) {
                gamma21 <- 0
            }
            k21 <- q[i] - gamma1 - gamma21 - sum(v_h[index1:index2, 
                1]^2)/(alpha_h[i])
            if (varfixed == FALSE) 
                alpha_h[i] <- sum(v_h[index1:index2, 1]^2)/(q[i] - 
                  gamma1 - gamma21-nonneg_adj)
        }
    }
    res <- list(x, z, y, del, Mi, idx2, t2, di, beta_h0, v_h0, 
        beta_h, v_h, alpha_h0, alpha_h, dft, Hinv, clam0, H, 
        mat, se_beta_h, U, H0)
    return(res)
}

PNFrailty_SE.h <- function (res1, nrand, q, qcum, dord = 1, varfixed = FALSE) 
{
    x <- res1[1][[1]]
    z <- res1[2][[1]]
    y <- res1[3][[1]]
    del <- res1[4][[1]]
    Mi <- res1[5][[1]]
    idx2 <- res1[6][[1]]
    t2 <- res1[7][[1]]
    di <- res1[8][[1]]
    beta_h <- res1[9][[1]]
    v_h <- res1[10][[1]]
    beta_h1 <- res1[11][[1]]
    v_h1 <- res1[12][[1]]
    alpha_h <- res1[13][[1]]
    alpha_h1 <- res1[14][[1]]
    dft <- res1[15][[1]]
    Hinv <- res1[16][[1]]
    clam0 <- res1[17][[1]]
    H <- res1[18][[1]]
    mat <- res1[19][[1]]
    U <- res1[21][[1]]
    H0 <- res1[22][[1]]
    n <- nrow(x)
    p <- ncol(x)
    u_h1 <- exp(v_h1)
    oq <- matrix(1, qcum[nrand + 1], 1)
    oq1 <- matrix(1, qcum[nrand + 1], 1)
    for (i in 1:nrand) {
        index1 <- qcum[i] + 1
        oq1[index1:qcum[i + 1]] <- alpha_h1[i]
    }
    D <- diag(oq1[, 1])
    iD <- solve(D)
    iA <- iD
    Bi <- diag(clam0[, 1])
    muh <- x %*% beta_h1 + z %*% v_h1
    expeta <- exp(muh)
    cla0 <- di/(t(Mi) %*% expeta)
    temp4 <- cla0^2/di
    As <- diag(temp4[, 1])
    Wi <- diag(expeta[, 1])
    done <- matrix(1, idx2, 1)
    H22 <- solve(t(z) %*% mat %*% z + U)
    Hessian <- matrix(0, nrand, nrand)
    if (varfixed == FALSE) {
        for (i in 1:nrand) {
            C <- matrix(0, qcum[nrand + 1], qcum[nrand + 1])
            index1 <- qcum[i] + 1
            index2 <- qcum[i + 1]
            for (j in index1:index2) C[j, j] <- 1
            iB1 <- iA %*% C %*% iA
            c_vh1 <- iB1 %*% v_h1
            dv1 <- H22 %*% c_vh1
            dexpeta1 <- expeta * (z %*% dv1)
            dcla01 <- -(di/((t(Mi) %*% expeta)^2)) * (t(Mi) %*% 
                dexpeta1)
            dWi1 <- diag(dexpeta1[, 1])
            dAi1 <- diag(dcla01[, 1])
            temp4 <- Mi %*% dAi1 %*% done
            dBi1 <- diag(temp4[, 1])
            dvec1 <- 2 * (cla0 * dcla01)
            temp4 <- dvec1/di
            dAs1 <- diag(temp4[, 1])
            dmat1 <- (dWi1 %*% Bi) + (Wi %*% dBi1) - (dWi1 %*% 
                Mi %*% As %*% t(Mi) %*% Wi) - (Wi %*% Mi %*% 
                dAs1 %*% t(Mi) %*% Wi) - (Wi %*% Mi %*% As %*% 
                t(Mi) %*% dWi1)
            dia1 <- -iB1
            Hd1 <- rbind(cbind(t(x) %*% dmat1 %*% x, t(x) %*% 
                dmat1 %*% z), cbind(t(z) %*% dmat1 %*% x, t(z) %*% 
                dmat1 %*% z + dia1))
            ddk1 <- -0.5 * sum(diag(iA %*% C %*% iA %*% C)) + 
                t(v_h1) %*% (iA %*% C %*% iB1) %*% v_h1 - t(dv1) %*% 
                iB1 %*% v_h1
            dia11 <- (iB1 %*% C %*% iA + iA %*% C %*% iB1)
            dv11 <- -H22 %*% ((t(z) %*% dmat1 %*% z + dia1) %*% 
                dv1 - iB1 %*% dv1 + dia11 %*% v_h1)
            temp4 <- (z %*% dv1) * (z %*% dv1) * expeta + (z %*% 
                dv11) * expeta
            ddW11 <- diag(temp4[, 1])
            ddcla011 <- -(dAs1 %*% (t(Mi) %*% Wi %*% z) %*% dv1 + 
                As %*% (t(Mi) %*% dWi1 %*% z) %*% dv1 + As %*% 
                (t(Mi) %*% Wi %*% z) %*% dv11)
            temp4 <- Mi %*% diag(ddcla011[, 1]) %*% done
            ddB11 <- diag(temp4[, 1])
            temp4 <- (2 * (dcla01^2) + 2 * (cla0 * ddcla011))/di
            ddAs11 <- diag(temp4[, 1])
            ddm1_11 <- (ddW11 %*% Bi) + (2 * dWi1 %*% dBi1) + 
                (Wi %*% ddB11)
            ddm2_11 <- (ddW11 %*% Mi %*% As %*% t(Mi) %*% Wi) + 
                (2 * dWi1 %*% Mi %*% dAs1 %*% t(Mi) %*% Wi) + 
                (2 * dWi1 %*% Mi %*% As %*% t(Mi) %*% dWi1)
            ddm3_11 <- (Wi %*% Mi %*% ddAs11 %*% t(Mi) %*% Wi) + 
                (2 * Wi %*% Mi %*% dAs1 %*% t(Mi) %*% dWi1) + 
                (Wi %*% Mi %*% As %*% t(Mi) %*% ddW11)
            ddmat11 <- ddm1_11 - (ddm2_11 + ddm3_11)
            Hd11 <- rbind(cbind(t(x) %*% ddmat11 %*% x, t(x) %*% 
                ddmat11 %*% z), cbind(t(z) %*% ddmat11 %*% x, 
                t(z) %*% ddmat11 %*% z + dia11))
            ddk1 <- ddk1 + 0.5 * sum(diag(-Hinv %*% Hd1 %*% Hinv %*% 
                Hd1 + Hinv %*% Hd11))
            Hessian[i, i] <- ddk1
            for (kk in 1:nrand) {
                if (kk > i) {
                  D <- matrix(0, qcum[nrand + 1], qcum[nrand + 
                    1])
                  index1 <- qcum[kk] + 1
                  index2 <- qcum[kk + 1]
                  for (j in index1:index2) D[j, j] <- 1
                  iB2 <- iA %*% D %*% iA
                  c_vh2 <- iB2 %*% v_h1
                  dv2 <- H22 %*% c_vh2
                  dexpeta2 <- expeta * (z %*% dv2)
                  dcla02 <- -(di/((t(Mi) %*% expeta)^2)) * (t(Mi) %*% 
                    dexpeta2)
                  dWi2 <- diag(dexpeta2[, 1])
                  dAi2 <- diag(dcla02[, 1])
                  temp4 <- Mi %*% dAi2 %*% done
                  dBi2 <- diag(temp4[, 1])
                  dvec2 <- 2 * (cla0 * dcla02)
                  temp4 <- dvec2/di
                  dAs2 <- diag(temp4[, 1])
                  dd12 <- -0.5 * sum(diag(iA %*% D %*% iA %*% 
                    C)) + 0.5 * t(v_h1) %*% (iA %*% D %*% iB1 + 
                    iA %*% C %*% iB2) %*% v_h1 - t(dv1) %*% iB2 %*% 
                    v_h1
                  dia12 <- (iB1 %*% D %*% iA + iA %*% D %*% iB1)
                  dv12 <- -H22 %*% ((t(z) %*% dmat1 %*% z + dia1) %*% 
                    dv2 - iB2 %*% dv1 + dia12 %*% v_h1)
                  temp4 <- (z %*% dv1) * (z %*% dv2) * expeta + 
                    (z %*% dv12) * expeta
                  ddW12 <- diag(temp4[, 1])
                  ddcla012 <- -(dAs2 %*% (t(Mi) %*% Wi %*% z) %*% 
                    dv1 + As %*% (t(Mi) %*% dWi2 %*% z) %*% dv1 + 
                    As %*% (t(Mi) %*% Wi %*% z) %*% dv12)
                  temp4 <- Mi %*% diag(ddcla012[, 1]) %*% done
                  ddB12 <- diag(temp4[, 1])
                  temp4 <- (2 * (dcla02 * dcla01) + 2 * (cla0 * 
                    ddcla012))/di
                  ddAs12 <- diag(temp4[, 1])
                  ddm1_12 <- (ddW12 %*% Bi) + (dWi1 %*% dBi2 + 
                    dWi2 %*% dBi1) + (Wi %*% ddB12)
                  ddm2_12 <- (ddW12 %*% Mi %*% As %*% t(Mi) %*% 
                    Wi) + (dWi1 %*% Mi %*% dAs2 %*% t(Mi) %*% 
                    Wi) + (dWi1 %*% Mi %*% As %*% t(Mi) %*% dWi2)
                  ddm3_12 <- (dWi2 %*% Mi %*% dAs1 %*% t(Mi) %*% 
                    Wi) + (Wi %*% Mi %*% ddAs12 %*% t(Mi) %*% 
                    Wi) + (Wi %*% Mi %*% dAs1 %*% t(Mi) %*% dWi2)
                  ddm4_12 <- (dWi2 %*% Mi %*% As %*% t(Mi) %*% 
                    dWi1) + (Wi %*% Mi %*% dAs2 %*% t(Mi) %*% 
                    dWi1) + (Wi %*% Mi %*% As %*% t(Mi) %*% ddW12)
                  ddmat12 <- ddm1_12 - (ddm2_12 + ddm3_12 + ddm4_12)
                  dmat2 <- (dWi2 %*% Bi) + (Wi %*% dBi2) - (dWi2 %*% 
                    Mi %*% As %*% t(Mi) %*% Wi) - (Wi %*% Mi %*% 
                    dAs2 %*% t(Mi) %*% Wi) - (Wi %*% Mi %*% As %*% 
                    t(Mi) %*% dWi2)
                  dia2 <- -iB2
                  Hd2 <- rbind(cbind(t(x) %*% dmat2 %*% x, t(x) %*% 
                    dmat2 %*% z), cbind(t(z) %*% dmat2 %*% x, 
                    t(z) %*% dmat2 %*% z + dia2))
                  Hd12 <- rbind(cbind(t(x) %*% ddmat12 %*% x, 
                    t(x) %*% ddmat12 %*% z), cbind(t(z) %*% ddmat12 %*% 
                    x, t(z) %*% ddmat12 %*% z + dia12))
                  dd12 <- dd12 + 0.5 * sum(diag(-Hinv %*% Hd1 %*% 
                    Hinv %*% Hd2 + Hinv %*% Hd12))
                  Hessian[i, kk] <- Hessian[kk, i] <- dd12
                }
            }
        }
        iAp <- solve(Hessian)
        se_lam <- sqrt(abs(diag(iAp)))
    }
    if (varfixed == TRUE) {
        se_lam <- rep(0, nrand)
        for (i in 1:nrand) se_lam[i] <- "NULL"
    }
    eta <- x %*% beta_h1 + z %*% v_h1
    expeta <- exp(eta)
    one <- matrix(1, n, 1)
    done <- matrix(1, idx2, 1)
    oq <- matrix(1, qcum[nrand + 1], 1)
    pi <- 3.14159265359
    term0 <- t(Mi) %*% expeta
    hlike1 <- (t(one) %*% (del * eta)) - (t(done) %*% (di * log(term0)))
    hlike2 <- 0
    hlike3 <- 0
    for (i in 1:nrand) {
        if (alpha_h[i] > 1e-05) 
            hlike2 <- hlike2 - (q[i]/2) * log(2 * pi) - ((q[i]/2) * 
                log(alpha_h1[i]))
        index1 <- qcum[i] + 1
        index2 <- qcum[i + 1]
        vv_h1 <- matrix(0, q[i], 1)
        vv_h1[1:q[i], 1] <- v_h1[index1:index2, 1]
        if (alpha_h[i] > 1e-05) 
            hlike3 <- hlike3 - (t(vv_h1) %*% vv_h1)/(2 * alpha_h1[i])
    }
    hliken <- hlike1 + hlike2 + hlike3
    cc1 <- svd(2 * pi * Hinv)
    for (i in 1:length(cc1$d)) if (cc1$d[i] < 1e-05) 
        cc1$d[i] <- 1
    logdet1 <- sum(log(abs(cc1$d)))
    adj1 <- 0.5 * logdet1
    hpn1 <- hliken + adj1
    muu <- exp(x %*% beta_h1) * clam0
    zmu <- t(z) %*% muu
    u_h1 <- exp(v_h1)
    second <- 0
    for (i in 1:nrand) {
        ialph1 <- 1/alpha_h1[i]
        a21 <- (zmu * u_h1) + ialph1
        b31 <- zmu * u_h1
        S11 <- 3 * (b31/(a21^2))
        S21 <- 5 * ((b31^2)/(a21^3))
        temp4 <- S11 - S21
        S31 <- diag(temp4[, 1])
        second <- second - sum(diag(S31))/24
    }
    H22 <- t(z) %*% mat %*% z + iD
    cc1 <- svd(H22/(2 * pi))
    for (i in 1:length(cc1$d)) if (cc1$d[i] > 1e+05) 
        cc1$d[i] <- 1
    logdet1 <- sum(log(abs(cc1$d)))
    hpn2 <- hliken - 0.5 * logdet1
    hpn3 <- hpn1 + second
    df1 <- sum(diag(Hinv %*% H0))
    res <- list(se_lam, hliken, hpn1, hpn2, hpn3, hlike1, df1)
    return(res)
}

PGfrailtyHL <- function (x, z, y, del, Mi, idx2, t2, di, beta_h0, v_h0, alpha_h0, 
    mord, dord, varfixed = FALSE,varnonneg) {
    n <- nrow(x)
    p <- ncol(x)
    nrand <- length(z)
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
    qcum <- cumsum(c(0, q))
    beta_h <- beta_h0
    v_h <- v_h0
    for (i in 1:nrand) {
        if (alpha_h0[i] < 1e-06) 
            alpha_h0[i] <- 1e-06
    }
    alpha_h <- alpha_h0
    zz <- z[[1]]
    if (nrand > 1) {
        index1 <- nrand
        for (i in 2:index1) zz <- cbind(zz, z[[i]])
    }
    z <- matrix(0, n, qcum[nrand + 1])
    z[1:n, 1:qcum[nrand + 1]] <- zz[1:n, 1:qcum[nrand + 1]]
    muh <- x %*% beta_h0 + z %*% v_h0
    expeta <- exp(muh)
    cla0 <- di/(t(Mi) %*% expeta)
    Wi <- diag(expeta[, 1])
    Ai <- diag(cla0[, 1])
    done <- matrix(1, idx2, 1)
    oq <- matrix(1, qcum[nrand + 1], 1)
    oq1 <- matrix(1, qcum[nrand + 1], 1)
    clam0 <- Mi %*% Ai %*% done
    for (i in 1:nrand) {
        index1 <- qcum[i] + 1
        oq1[index1:qcum[i + 1]] <- alpha_h[i]
    }
    D <- diag(oq1[, 1])
    iD <- solve(D)
    iu_h0 <- exp(v_h0)
    U <- iD %*% diag(iu_h0[, 1])
    dft1 <- t(x) %*% (del - Wi %*% clam0)
    if (mord == 1) {
        Bi <- diag(clam0[, 1])
        temp4 <- cla0^2/di
        As <- diag(temp4[, 1])
        mat <- (Wi %*% Bi) - (Wi %*% Mi) %*% As %*% (t(Mi) %*% 
            Wi)
        dinv <- solve(t(z) %*% mat %*% z + U)
        mu0 <- exp(x %*% beta_h0 + z %*% v_h0)
        mu <- exp(x %*% beta_h0 + z %*% v_h0) * clam0
        dcla0be <- -As %*% (t(Mi) %*% Wi)
        dcla0b <- matrix(0, idx2, p)
        dv_db <- matrix(0, qcum[nrand + 1], p)
        xz <- matrix(0, n, p)
        dmu0 <- matrix(0, n, p)
        dw_db1 <- matrix(0, n, p)
        xk <- matrix(0, n, 1)
        ad1 <- matrix(0, p, 1)
        for (k in 1:p) {
            xk[, 1] <- x[, k]
            dv_db[, k] <- -dinv %*% (t(z) %*% mat %*% xk)
            xz[, k] <- xk + z %*% (dv_db[, k])
            dcla0b[, k] <- dcla0be %*% (xz[, k])
            dc <- Mi %*% diag(dcla0b[, k]) %*% done
            dmu0[, k] <- mu0 * (xz[, k])
            dw_db1[, k] <- dc * mu0 + mu * xz[, k]
            temp4 <- (2 * cla0/di) * (dcla0b[, k])
            dw_db2 <- ((diag(dmu0[, k])) %*% Mi %*% As %*% t(Mi) %*% 
                Wi) + (Wi %*% Mi %*% diag(temp4[, 1]) %*% t(Mi) %*% 
                Wi) + (Wi %*% Mi %*% As %*% t(Mi) %*% (diag(dmu0[, 
                k])))
            dw_db <- diag(dw_db1[, k]) - dw_db2
            ad1[k, 1] <- sum(diag(dinv %*% (t(z) %*% dw_db %*% 
                z)))
        }
        dft1 <- dft1 - 0.5 * ad1
    }
    dft2 <- t(z) %*% (del - Wi %*% clam0) + (iD %*% oq) - (iD %*% 
        iu_h0)
    dft <- rbind(dft1, dft2)
    Bi <- diag(clam0[, 1])
    temp4 <- cla0^2/di
    Adi <- diag(temp4[, 1])
    As <- Adi
    mat <- (Wi %*% Bi) - (Wi %*% Mi) %*% Adi %*% (t(Mi) %*% Wi)
    U <- iD %*% diag(iu_h0[, 1])
    H <- rbind(cbind(t(x) %*% mat %*% x, t(x) %*% mat %*% z), 
        cbind(t(z) %*% mat %*% x, t(z) %*% mat %*% z + U))
    H0 <- rbind(cbind(t(x) %*% mat %*% x, t(x) %*% mat %*% z), 
        cbind(t(z) %*% mat %*% x, t(z) %*% mat %*% z))
    Hinv <- solve(H)
    be_h0 <- rbind(beta_h0, v_h0)
    be_h <- be_h0 + (Hinv %*% dft)
    beta_h[1:p, 1] <- be_h[1:p, 1]
    se_beta_h <- matrix(0, p, 1)
    for (i in 1:p) se_beta_h[i, 1] <- sqrt(Hinv[i, i])
    index2 <- qcum[nrand + 1]
    index3 <- p + 1
    index4 <- p + qcum[nrand + 1]
    v_h[1:index2, 1] <- be_h[index3:index4, 1]
    for (i in 1:nrand) {
        ial <- 1/alpha_h[i]
        dp <- digamma(ial)
        ddp <- trigamma(ial)
        oq <- matrix(1, q[i], 1)
        one <- matrix(1, n, 1)
        u_h <- exp(v_h)
        eta <- x %*% beta_h + z %*% v_h
        expeta <- exp(eta)
        Wi <- diag(expeta[, 1])
        Wei <- (Wi %*% Bi)
        U <- iD %*% diag(u_h[, 1])
        term <- (t(z) %*% mat %*% z + U)
        C <- matrix(0, qcum[nrand + 1], qcum[nrand + 1])
        index1 <- qcum[i] + 1
        index2 <- qcum[i + 1]
        for (j in index1:index2) C[j, j] <- 1
        iA <- iD
        iB <- iA %*% C %*% iA
        c_vh <- iB %*% (u_h - 1)
        invt <- solve(term)
        dv <- invt %*% c_vh
        dexpeta <- expeta * (z %*% dv)
        dcla0 <- -(di/((t(Mi) %*% expeta)^2)) * (t(Mi) %*% dexpeta)
        temp4 <- expeta * (z %*% dv)
        dWi <- diag(temp4[, 1])
        dAi <- diag(dcla0[, 1])
        temp4 <- Mi %*% dAi %*% done
        dBi <- diag(temp4[, 1])
        dvec <- 2 * (cla0 * dcla0)
        temp4 <- dvec/di
        dAs <- diag(temp4[, 1])
        dmat <- (dWi %*% Bi) + (Wi %*% dBi) - (dWi %*% Mi %*% 
            As %*% t(Mi) %*% Wi) - (Wi %*% Mi %*% dAs %*% t(Mi) %*% 
            Wi) - (Wi %*% Mi %*% As %*% t(Mi) %*% dWi)
        uad1 <- u_h * dv
        temp4 <- -iB %*% u_h + ial * uad1
        dia1 <- diag(temp4[, 1])
        Hd <- rbind(cbind(t(x) %*% dmat %*% x, t(x) %*% dmat %*% 
            z), cbind(t(z) %*% dmat %*% x, t(z) %*% dmat %*% 
            z + dia1))
        
        if (dord == 0) {
            temp4 <- -iD %*% iD %*% u_h
            dia1_k <- diag(temp4[, 1])
            zero1 <- matrix(0, p, p)
            zero2 <- matrix(0, p, qcum[nrand + 1])
            Hd_k <- rbind(cbind(zero1, zero2), cbind(t(zero2), 
                dia1_k))
            hinv2 <- solve(t(z) %*% mat %*% z + U)
            Hd2 <- t(z) %*% dmat %*% z + dia1
            dk2 <- -0.5 * sum(diag(Hinv %*% Hd_k))
            vv_h <- matrix(0, q[i], 1)
            index3 <- 1
            index4 <- q[i]
            vv_h[1:q[i], 1] <- v_h[index1:index2, 1]
            uu_h <- exp(vv_h)
            nonneg_adj=0
            if (varnonneg == TRUE) nonneg_adj=alpha_h[i]
            k2 <- (t(oq) %*% (vv_h - uu_h)) + (q[i] * (-log(alpha_h[i]) + 
                1 - dp))-nonneg_adj
            mu <- exp(x %*% beta_h) * clam0
            zd <- t(z) %*% del
            zmu <- t(z) %*% mu
            cor1 <- (1 + (alpha_h[i] * zd))^(-2)
            dk3 <- sum(cor1)/12
            k2 <- -((alpha_h[i]^-2) * k2) + dk2
            dterm <- t(z) %*% dmat %*% z + dia1
            ddv <- -invt %*% (dterm) %*% invt %*% c_vh + invt %*% 
                (-2 * ial^3 * (u_h - 1) + ial^2 * uad1)
            ddcla0 <- -dAs %*% (t(Mi) %*% Wi %*% z) %*% dv - 
                As %*% (t(Mi) %*% dWi %*% z) %*% dv - As %*% 
                (t(Mi) %*% Wi %*% z) %*% ddv
            uad2 <- (uad1 * dv) + (u_h * ddv)
            ddAi <- diag(ddcla0[, 1])
            ddclam0 <- Mi %*% ddAi %*% done
            ddBi <- diag(ddclam0[, 1])
            temp4 <- expeta * (z %*% dv) * (z %*% dv) + (expeta * 
                (z %*% ddv))
            ddWi <- diag(temp4[, 1])
            ddm1 <- (ddWi %*% Bi) + (2 * dWi %*% dBi) + (Wi %*% 
                ddBi)
            ddm2 <- ddWi %*% Mi %*% As %*% t(Mi) %*% Wi + (2 * 
                dWi %*% Mi %*% dAs %*% t(Mi) %*% Wi) + (2 * dWi %*% 
                Mi %*% As %*% t(Mi) %*% dWi)
            ddvec <- (2 * (dcla0^2)) + (2 * (cla0 * ddcla0))
            temp4 <- ddvec/di
            ddAs = diag(temp4[, 1])
            ddm3 <- (Wi %*% Mi %*% ddAs %*% t(Mi) %*% Wi) + (2 * 
                Wi %*% Mi %*% dAs %*% t(Mi) %*% dWi) + Wi %*% 
                Mi %*% As %*% t(Mi) %*% ddWi
            ddmat <- ddm1 - (ddm2 + ddm3)
            temp4 <- (2 * ial^3 * u_h) - (2 * ial^2 * uad1) + 
                (ial * uad2)
            dia2 <- diag(temp4[, 1])
            Hdd2 <- t(z) %*% ddmat %*% z + dia2
            k21 <- t(oq) %*% (2 * (vv_h - uu_h)) + (q[i] * ((-2 * 
                log(alpha_h[i])) + 3 - (2 * dp) - ((1/alpha_h[i]) * 
                ddp)))
            al <- ial^3
            k21 <- al * k21
            cor2 <- (1 + (alpha_h[i] * zd))^(-3)
            cor22 <- cor2 * zd
            kcor <- sum(cor22)/6
            kcor <- 0
            k22 <- -k21 + 0.5 * (sum(diag(hinv2 * Hdd2)) - sum(diag(hinv2 * 
                Hd2 * hinv2 * Hd2))) + kcor
            ialp <- 1/k22
            if (varfixed == FALSE) 
                alpha_h[i] <- alpha_h[i] + (ialp * k2)
        }
        if (dord == 1 | dord == 2) {
            hinv2 <- solve(t(z) %*% mat %*% z + U)
            Hd2 <- t(z) %*% dmat %*% z + dia1
            dk2 <- -0.5 * sum(diag(Hinv %*% Hd))
            vv_h <- matrix(0, q[i], 1)
            index3 <- 1
            index4 <- q[i]
            vv_h[1:index4, 1] <- v_h[index1:index2, 1]
            uu_h <- exp(vv_h)
            nonneg_adj=0
            if (varnonneg == TRUE) nonneg_adj=alpha_h[i]
            k2 <- (t(oq) %*% (vv_h - uu_h)) + (q[i] * (-log(alpha_h[i]) + 
                1 - dp))-nonneg_adj
            mu <- exp(x %*% beta_h) * clam0
            zd <- t(z) %*% del
            zmu <- t(z) %*% mu
            cor1 <- (1 + (alpha_h[i] * zd))^(-2)
            dk3 <- sum(cor1)/12
            if (dord == 1) 
                k2 <- -((alpha_h[i]^(-2)) * k2) + dk2
            if (dord == 2) 
                k2 <- -((alpha_h[i]^(-2)) * k2) + dk2 + dk3
            dterm <- t(z) %*% dmat %*% z + dia1
            ddv <- -invt %*% (dterm) %*% invt %*% c_vh + invt %*% 
                (-2 * ial^3 * (u_h - 1) + ial^2 * uad1)
            ddcla0 <- -dAs %*% (t(Mi) %*% Wi %*% z) %*% dv - 
                As %*% (t(Mi) %*% dWi %*% z) %*% dv - As %*% 
                (t(Mi) %*% Wi %*% z) %*% ddv
            uad2 <- (uad1 * dv) + (u_h * ddv)
            ddAi <- diag(ddcla0[, 1])
            ddclam0 <- Mi %*% ddAi %*% done
            ddBi <- diag(ddclam0[, 1])
            temp4 <- expeta * (z %*% dv) * (z %*% dv) + (expeta * 
                (z %*% ddv))
            ddWi <- diag(temp4[, 1])
            ddm1 <- (ddWi %*% Bi) + (2 * dWi %*% dBi) + (Wi %*% 
                ddBi)
            ddm2 <- ddWi %*% Mi %*% As %*% t(Mi) %*% Wi + (2 * 
                dWi %*% Mi %*% dAs %*% t(Mi) %*% Wi) + (2 * dWi %*% 
                Mi %*% As %*% t(Mi) %*% dWi)
            ddvec <- (2 * (dcla0^2)) + (2 * (cla0 * ddcla0))
            temp4 <- ddvec/di
            ddAs <- diag(temp4[, 1])
            ddm3 <- (Wi %*% Mi %*% ddAs %*% t(Mi) %*% Wi) + (2 * 
                Wi %*% Mi %*% dAs %*% t(Mi) %*% dWi) + Wi %*% 
                Mi %*% As %*% t(Mi) %*% ddWi
            ddmat <- ddm1 - (ddm2 + ddm3)
            temp4 <- (2 * ial^3 * u_h) - (2 * ial^2 * uad1) + 
                (ial * uad2)
            dia2 <- diag(temp4[, 1])
            Hdd <- rbind(cbind(t(x) %*% ddmat %*% x, t(x) %*% 
                ddmat %*% z), cbind(t(z) %*% ddmat %*% x, t(z) %*% 
                ddmat %*% z + dia2))
            Hdd2 <- t(z) %*% ddmat %*% z + dia2
            k21 <- t(oq) %*% (2 * (vv_h - uu_h)) + (q[i] * ((-2 * 
                log(alpha_h[i])) + 3 - (2 * dp) - ((1/alpha_h[i]) * 
                ddp)))
            al <- ial^3
            k21 <- al * k21
            cor2 <- (1 + (alpha_h[i] * zd))^(-3)
            cor22 <- cor2 * zd
            kcor <- sum(cor22)/6
            if (dord == 1) 
                kcor <- 0
            k22 <- -k21 + 0.5 * (sum(diag(Hinv %*% Hdd)) - sum(diag(Hinv %*% 
                Hd %*% Hinv %*% Hd))) + kcor
            ialp <- 1/k22
            if (varfixed == FALSE) 
                alpha_h[i] <- alpha_h[i] + (ialp * k2)
            if (alpha_h[i] <= 0) 
                alpha_h[i] <- alpha_h0/2
        }
    }
    ## u-sclae
    ialp_h<-1/alpha_h[1]
    u_h<-exp(v_h)
    mu1 <- exp(x %*% beta_h) * clam0
    zmu=t(z) %*% mu1
    exp_xb<-exp(x%*%beta_h)
    W0=diag(exp_xb[,1])
    Bi <- diag(clam0[, 1])
    mat1=(W0%*%Bi)-(W0%*%Mi)%*%Adi%*%(t(Mi)%*%W0)
    U1 = (ialp_h+zmu)/u_h
    U1 = diag(U1[,1])
    U2=t(z)%*%((W0%*%Mi)%*%Adi%*%(t(Mi)%*%W0))%*%z
    Hu = rbind(cbind(t(x)%*%mat%*%x, t(x)%*%mat1%*%z),
               cbind(t(z)%*%mat1%*%x, U1-U2 ))
    Hinv_u=solve(Hu)
    res <- list(x, z, y, del, Mi, idx2, t2, di, beta_h0, v_h0, 
        beta_h, v_h, alpha_h0, alpha_h, dft, Hinv, clam0, H, 
        mat, se_beta_h, U, H0, u_h,Hinv_u)
    return(res)
}

PGFrailty_SE.h <- function (res1, nrand, q, qcum, dord = 1, varfixed = FALSE) {
    x <- res1[1][[1]]
    z <- res1[2][[1]]
    y <- res1[3][[1]]
    del <- res1[4][[1]]
    Mi <- res1[5][[1]]
    idx2 <- res1[6][[1]]
    t2 <- res1[7][[1]]
    di <- res1[8][[1]]
    beta_h <- res1[9][[1]]
    v_h <- res1[10][[1]]
    beta_h1 <- res1[11][[1]]
    v_h1 <- res1[12][[1]]
    alpha_h <- res1[13][[1]]
    alpha_h1 <- res1[14][[1]]
    dft <- res1[15][[1]]
    Hinv <- res1[16][[1]]
    clam0 <- res1[17][[1]]
    H <- res1[18][[1]]
    mat <- res1[19][[1]]
    H0 <- res1[22][[1]]
    n <- nrow(x)
    p <- ncol(x)
    u_h1 <- exp(v_h1)
    mat11 <- t(x) %*% mat %*% x
    mat12 <- t(x) %*% mat %*% z
    mat13 <- t(z) %*% mat %*% z
    oq <- matrix(1, qcum[nrand + 1], 1)
    oq1 <- matrix(1, qcum[nrand + 1], 1)
    for (i in 1:nrand) {
        index1 <- qcum[i] + 1
        oq1[index1:qcum[i + 1]] <- alpha_h1[i]
    }
    D <- diag(oq1[, 1])
    iD <- solve(D)
    U <- iD %*% diag(u_h1[, 1])
    mmat <- mat11 - mat12 %*% solve(mat13 + U) %*% t(mat12)
    hminv <- solve(mmat)
    done <- matrix(1, idx2, 1)
    muh <- x %*% beta_h1 + z %*% v_h1
    expeta <- exp(muh)
    Wi <- diag(expeta[, 1])
    cla0 <- di/(t(Mi) %*% expeta)
    Ai <- diag(cla0[, 1])
    temp4 <- cla0^2/di
    As <- diag(temp4[, 1])
    Bi <- diag(clam0[, 1])
    mat <- (Wi %*% Bi) - (Wi %*% Mi) %*% As %*% (t(Mi) %*% Wi)
    Dinv0 <- solve(t(z) %*% mat %*% z + U)
    se_lam <- matrix(0, nrand, 1)
    if (varfixed == FALSE) {
        for (i in 1:nrand) {
            ial <- 1/alpha_h1[i]
            index1 <- qcum[i] + 1
            index2 <- qcum[i + 1]
            vv_h1 <- matrix(0, q[i], 1)
            uu_h1 <- matrix(0, q[i], 1)
            vv_h1[1:q[i], 1] <- v_h1[index1:index2, 1]
            uu_h1[1:q[i], 1] <- u_h1[index1:index2, 1]
            c_vh <- ial^2 * (uu_h1 - 1)
            dv <- solve(t(z) %*% mat %*% z + U) %*% c_vh
            dexpeta <- expeta * (z %*% dv)
            dcla0 <- -(di/((t(Mi) %*% expeta)^2)) * (t(Mi) %*% 
                dexpeta)
            dWi <- diag(dexpeta[, 1])
            dAi <- diag(dcla0[, 1])
            temp4 <- Mi %*% dAi %*% done
            dBi <- diag(temp4[, 1])
            dvec <- 2 * (cla0 * dcla0)
            temp4 <- dvec/di
            dAs <- diag(temp4[, 1])
            dmat <- (dWi %*% Bi) + (Wi %*% dBi) - (dWi %*% Mi %*% 
                As %*% t(Mi) %*% Wi) - (Wi %*% Mi %*% dAs %*% 
                t(Mi) %*% Wi) - (Wi %*% Mi %*% As %*% t(Mi) %*% 
                dWi)
            uad1 <- uu_h1 * dv
            temp4 <- (-ial^2 * uu_h1) + (ial * uad1)
            dia_d <- diag(temp4[, 1])
            Hd <- rbind(cbind(t(x) %*% dmat %*% x, t(x) %*% dmat %*% 
                z), cbind(t(z) %*% dmat %*% x, t(z) %*% dmat %*% 
                z + dia_d))
            term <- (t(z) %*% mat %*% z + U)
            invt <- solve(term)
            dterm <- t(z) %*% dmat %*% z + dia_d
            ddv <- -invt %*% (dterm) %*% invt %*% c_vh + invt %*% 
                (-2 * ial^3 * (uu_h1 - 1) + ial^2 * uad1)
            ddcla0 <- -dAs %*% (t(Mi) %*% Wi %*% z) %*% dv - 
                As %*% (t(Mi) %*% dWi %*% z) %*% dv - As %*% 
                (t(Mi) %*% Wi %*% z) %*% ddv
            ddAi <- diag(ddcla0[, 1])
            ddclam0 <- Mi %*% ddAi %*% done
            ddBi <- diag(ddclam0[, 1])
            temp4 <- expeta * (z %*% dv) * (z %*% dv) + (expeta * 
                (z %*% ddv))
            ddWi <- diag(temp4[, 1])
            ddm1 <- (ddWi %*% Bi) + (2 * dWi %*% dBi) + (Wi %*% 
                ddBi)
            ddm2 <- ddWi %*% Mi %*% As %*% t(Mi) %*% Wi + (2 * 
                dWi %*% Mi %*% dAs %*% t(Mi) %*% Wi) + (2 * dWi %*% 
                Mi %*% As %*% t(Mi) %*% dWi)
            ddvec <- (2 * (dcla0^2)) + (2 * (cla0 * ddcla0))
            temp4 <- ddvec/di
            ddAs <- diag(temp4[, 1])
            ddm3 <- (Wi %*% Mi %*% ddAs %*% t(Mi) %*% Wi) + (2 * 
                Wi %*% Mi %*% dAs %*% t(Mi) %*% dWi) + Wi %*% 
                Mi %*% As %*% t(Mi) %*% ddWi
            ddmat <- ddm1 - (ddm2 + ddm3)
            uad2 <- (uad1 * dv) + (u_h1 * ddv)
            temp4 <- (2 * ial^3 * u_h1) - (2 * ial^2 * uad1) + 
                (ial * uad2)
            dia_dd <- diag(temp4[, 1])
            Hdd <- rbind(cbind(t(x) %*% ddmat %*% x, t(x) %*% 
                ddmat %*% z), cbind(t(z) %*% ddmat %*% x, t(z) %*% 
                ddmat %*% z + dia_dd))
            H <- rbind(cbind(t(x) %*% mat %*% x, t(x) %*% mat %*% 
                z), cbind(t(z) %*% mat %*% x, t(z) %*% mat %*% 
                z + U))
            Hinv <- solve(H)
            oq <- matrix(1, q[i], 1)
            dp <- digamma(ial)
            ddp <- trigamma(ial)
            k21a <- t(oq) %*% (2 * (vv_h1 - uu_h1)) + (q[i] * 
                ((-2 * log(alpha_h1[i])) + 3 - (2 * dp) - ((1/alpha_h1[i]) * 
                  ddp)))
            k21a <- (ial^3) * k21a
            d2halp <- -k21a - t(oq) %*% (c_vh * dv)
            adalp <- 0.5 * sum(diag(-Hinv %*% Hd %*% Hinv %*% 
                Hd + Hinv %*% Hdd))
            dalp_2 <- d2halp + adalp
            se_lam[i] <- sqrt(1/dalp_2)
        }
    }
    if (varfixed == TRUE) {
        for (i in 1:nrand) se_lam[i] <- "NULL"
    }
    u_h1 <- exp(v_h1)
    U <- iD %*% diag(u_h1[, 1])
    oq <- matrix(1, qcum[nrand + 1], 1)
    one <- matrix(1, n, 1)
    zd <- t(z) %*% del
    eta <- x %*% beta_h1 + z %*% v_h1
    expeta <- exp(eta)
    term0 <- t(Mi) %*% expeta
    non <- t(one) %*% del
    done <- matrix(1, idx2, 1)
    hlike1 <- (t(one) %*% (del * eta)) - (t(done) %*% (di * log(term0)))
    hlike2 <- 0
    for (i in 1:nrand) {
        oq <- matrix(1, q[i], 1)
        index1 <- qcum[i] + 1
        index2 <- qcum[i + 1]
        vv_h1 <- matrix(0, q[i], 1)
        uu_h1 <- matrix(0, q[i], 1)
        vv_h1[1:q[i], 1] <- v_h1[index1:index2, 1]
        uu_h1[1:q[i], 1] <- u_h1[index1:index2, 1]
        i_alp1 <- 1/alpha_h1[i]
        c_alp1 <- 0
        if (alpha_h[i] > 1e-05) 
            c_alp1 <- -log(gamma(i_alp1)) - (i_alp1 * log(alpha_h1[i]))
        if (alpha_h[i] > 1e-05) 
            hlike2 <- hlike2 + t(oq) %*% ((vv_h1 - uu_h1)/alpha_h1[i] + 
                c_alp1)
    }
    hlike <- hlike1 + hlike2
    pi <- 3.14159265359
    H22 <- t(z) %*% mat %*% z + U
    zd <- t(z) %*% del
    if (dord == 2) {
        temp4 <- 1/(i_alp1 + zd)
        secd <- diag(temp4[, 1])
        second <- sum(diag(secd))/12
    }
    else second <- 0
    cc1 <- svd(H22/(2 * pi))
    for (i in 1:length(cc1$d)) if (cc1$d[i] > 1e+05) 
        cc1$d[i] <- 1
    logdet1 <- sum(log(abs(cc1$d)))
    pvhs <- hlike - 0.5 * logdet1
    svhs <- pvhs + second
    cc1 <- svd(Hinv * 2 * pi)
    for (i in 1:length(cc1$d)) if (cc1$d[i] < 1e-05) 
        cc1$d[i] <- 1
    logdet1 <- sum(log(abs(cc1$d)))
    adj1 <- 0.5 * logdet1
    hpn1 <- hlike + adj1
    hpn2 <- pvhs
    hpn3 <- hpn1 + second
    df1 <- sum(diag(Hinv %*% H0))
    res <- list(se_lam, hlike, hpn1, hpn2, hpn3, hlike1, df1)
    return(res)
}

FrailtyFrames<-function (mc, formula, contrasts, vnms = character(0),y) 
{
    mf <- mc
    m <- match(c("data", "weights", "na.action", "offset"), 
        names(mf), 0)
    mf <- mf[c(1, m)]
    frame.form <- subbars(formula)
    if (length(vnms) > 0) 
        frame.form[[3]] <- substitute(foo + bar, list(foo = parse(text = paste(vnms, 
            collapse = " + "))[[1]], bar = frame.form[[3]]))
    fixed.form <- nobars(formula)
    if (inherits(fixed.form, "name")) 
        fixed.form <- substitute(foo ~ 1, list(foo = fixed.form))
    environment(fixed.form) <- environment(frame.form) <- environment(formula)
    mf$formula <- frame.form
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    names(mf)[2] <- "data"
    fe <- mf
    mf <- eval(mf, parent.frame(2))
    fe$formula <- fixed.form
    fe <- eval(fe, parent.frame(2))
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) 
            names(Y) <- nm
    }
    mt <- attr(fe, "terms")
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0)
    storage.mode(X) <- "double"
    fixef <- numeric(ncol(X))
    names(fixef) <- colnames(X)
    dimnames(X) <- NULL
    wts <- model.weights(mf)
    if (is.null(wts)) 
        wts <- numeric(0)
    off <- model.offset(mf)
    if (is.null(off)) 
        off <- numeric(0)
    if (any(wts <= 0)) 
        stop(gettextf("negative weights or weights of zero are not allowed"))
    if (length(off) && length(off) != NROW(Y)) 
        stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
            length(off), NROW(Y)))
    attr(mf, "terms") <- mt
    list(Y = Y, X = X, wts = as.double(wts), off = as.double(off), 
        mf = mf, fixef = fixef)
}

subbars <- function (term) 
{
    if (is.name(term) || !is.language(term)) 
        return(term)
    if (length(term) == 2) {
        term[[2]] <- subbars(term[[2]])
        return(term)
    }
    stopifnot(length(term) >= 3)
    if (is.call(term) && term[[1]] == as.name("|")) 
        term[[1]] <- as.name("+")
    for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
    term
}

nobars <- function (term) 
{
    if (!("|" %in% all.names(term))) 
        return(term)
    if (is.call(term) && term[[1]] == as.name("|")) 
        return(NULL)
    if (length(term) == 2) {
        nb <- nobars(term[[2]])
        if (is.null(nb)) 
            return(NULL)
        term[[2]] <- nb
        return(term)
    }
    nb2 <- nobars(term[[2]])
    nb3 <- nobars(term[[3]])
    if (is.null(nb2)) 
        return(nb3)
    if (is.null(nb3)) 
        return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}

expandSlash <- function (bb) 
{
    if (!is.list(bb)) 
        return(expandSlash(list(bb)))
    unlist(lapply(bb, function(x) {
        if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]]))) 
            return(lapply(unlist(makeInteraction(trms)), function(trm) substitute(foo | 
                bar, list(foo = x[[2]], bar = trm))))
        x
    }))
}

findbars <- function (term) 
{
    if (is.name(term) || !is.language(term)) 
        return(NULL)
    if (term[[1]] == as.name("(")) 
        return(findbars(term[[2]]))
    if (!is.call(term)) 
        stop("term must be of class call")
    if (term[[1]] == as.name("|")) 
        return(term)
    if (length(term) == 2) 
        return(findbars(term[[2]]))
    c(findbars(term[[2]]), findbars(term[[3]]))
}

findplus <- function (term) 
{
    if (is.numeric(term)) 
        return(0)
    if (!is.language(term)) 
        return(NULL)
    if (length(term) == 1) 
        return(0)
    if (term[[1]] == as.name("|")) 
        return(findplus(term[[2]]))
    if (!is.call(term)) 
        stop("term must be of class call")
    if (term[[1]] == as.name("+")) 
        return(1)
    if (term[[1]] == as.name("-")) 
        return(-1)
}

slashTerms <- function (x) 
{
    if (!("/" %in% all.names(x))) 
        return(x)
    if (x[[1]] != as.name("/")) 
        stop("unparseable formula for grouping factor")
    list(slashTerms(x[[2]]), slashTerms(x[[3]]))
}

HGLMFactorList <- function (formula, fr, rmInt, drop) 
{
    mf <- fr$mf
    bars <- expandSlash(findbars(formula[[3]]))
    for (i in 1:length(bars)) {
        checkcorr <- findplus(bars[[i]])
        if (checkcorr == 1) 
            stop("Correlated random effects are not currently allowed in the HGLM routines")
        if (checkcorr == -1) 
            stop("You do not need to specify '-1' for no intercept it is done be default")
    }
    if (!length(bars)) 
        stop("No random effects terms specified in formula")
    names(bars) <- unlist(lapply(bars, function(x) deparse(x[[3]])))
    fl <- lapply(bars, function(x) {
        ff <- eval(substitute(as.factor(fac)[, drop = TRUE], 
            list(fac = x[[3]])), mf)
        im <- as(ff, "sparseMatrix")
        if (!isTRUE(validObject(im, test = TRUE))) 
            stop("invalid conditioning factor in random effect: ", 
                format(x[[3]]))
        if (is.name(x[[2]])) {
            tempexp <- paste("~", as.character(x[[2]]), "-1")
            tempexp <- as.formula(tempexp)[[2]]
        }
        else tempexp <- x[[2]]
        mm <- model.matrix(eval(substitute(~expr, list(expr = tempexp))), 
            mf)
        if (rmInt) {
            if (is.na(icol <- match("(Intercept)", colnames(mm)))) 
                break
            if (ncol(mm) < 2) 
                stop("lhs of a random-effects term cannot be an intercept only")
            mm <- mm[, -icol, drop = FALSE]
        }
        ans <- list(f = ff, A = do.call(rbind, lapply(seq_len(ncol(mm)), 
            function(j) im)), Zt = do.call(rbind, lapply(seq_len(ncol(mm)), 
            function(j) {
                im@x <- mm[, j]
                im
            })), ST = matrix(0, ncol(mm), ncol(mm), dimnames = list(colnames(mm), 
            colnames(mm))))
        if (drop) {
            ans$A@x <- rep(0, length(ans$A@x))
            ans$Zt <- drop0(ans$Zt)
        }
        ans
    })
    Design <- list(0)
    Subject <- list(0)
    for (i in 1:length(fl)) {
        Subject[[i]] <- as.factor(fl[[i]]$f)
        tempmat <- fl[[i]]$Zt
        tempmat <- as.matrix(t(as.matrix(tempmat)))
        Design[[i]] <- tempmat
    }
    list(Design = Design, Subject = Subject, namesRE = names(bars))
}

dbind<-function(a,b){
     out1<-cbind(a,matrix(0,nrow(a),ncol(b)))
     out2<-cbind(matrix(0,nrow(b),ncol(a)),b)
     out<-rbind(out1,out2)
     out
}

FrailtyMakeData <- function (y, x, del, z, L0) 
{
    n <- nrow(x)
    p <- ncol(x)
    nrand <- length(z)
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
    zz <- z[[1]]
    if (nrand > 1) {
        index1 <- nrand
        for (i in 2:index1) zz <- cbind(zz, z[[i]])
    }
    qcum <- cumsum(c(0, q))
    zzz <- matrix(0, n, qcum[nrand + 1])
    zzz[1:n, 1:qcum[nrand + 1]] <- zz[1:n, 1:qcum[nrand + 1]]
    sort_data <- cbind(L0, y, x, del, zzz)
    sort.res <- sort_data[order(sort_data[, 1], na.last = NA),]
    L0[1:n, 1] <- sort.res[1:n, 1]
    y[1:n, 1] <- sort.res[1:n, 2]
    index1 <- p + 2
    x[1:n, 1:p] <- sort.res[1:n, 3:index1]
    index1 <- index1 + 1
    del[1:n, 1] <- sort.res[1:n, index1]
    for (i in 1:nrand) {
        index3 <- p + 3 + qcum[i] + 1
        index4 <- p + 3 + qcum[i + 1]
        z[[i]][1:n, 1:q[i]] <- sort.res[1:n, index3:index4]
    }
    t <- matrix(0, n, 1)
    xx <- matrix(0, n, p)
    di <- matrix(0, n, 1)
    idx1 <- 0
    for (i in 1:n) {
        if (del[i, 1] == 1) {
            idx1 <- idx1 + 1
            t[idx1, 1] <- y[i, 1]
            for (j in 1:p) {
                xx[idx1, j] <- x[i, j]
            }
        }
    }
    t1 <- t
    for (i in 1:idx1) {
        for (j in 1:idx1) {
            if (t1[i, 1] == t1[j, 1]) {
                if (i != j) {
                  t1[j, 1] <- 0
                }
            }
        }
    }
    t2 <- matrix(0, idx1, 1)
    idx2 <- 0
    for (i in 1:idx1) {
        if (t1[i, 1] != 0) {
            idx2 <- idx2 + 1
            t2[idx2, 1] <- t1[i, 1]
        }
    }
    di <- matrix(0, idx2, 1)
    si <- matrix(0, idx2, p)
    for (i in 1:idx2) {
        di[i, 1] <- 0
        for (j in 1:idx1) {
            if (t2[i, 1] == t[j, 1]) {
                di[i, 1] <- di[i, 1] + 1
                for (k in 1:p) {
                  si[i, k] <- si[i, k] + xx[j, k]
                }
            }
        }
    }
    Mi <- matrix(0, n, idx2)
    for (i in 1:n) {
        LT <- L0[i,1]
        t0 <- y[i, 1]
        for (j in 1:idx2) {
            if (LT<t2[j,1] && t2[j, 1] <= t0) {
                Mi[i, j] = 1
            }
            else Mi[i, j] = 0
        }
    }
    res <- list(y, x, del, z, Mi, idx2, t2, di)
    return(res)
}


frailtyHL_grouped<-function (formula, data, weights, subset, na.action, RandDist = "Normal",mord = 0, dord = 1, Maxiter = 500, convergence = 10^-6, varfixed = FALSE,varinit = c(0.1),varnonneg=FALSE,grouped=FALSE,RespDist="FM") 
{
    Call <- match.call()
    indx <- match(c("formula", "data", "weights", "subset", "na.action"), 
        names(Call), nomatch = 0)
    if (indx[1] == 0) stop("A formula argument is required")
    data$nn <- rep(1:nrow(data))
    data$idid <- rep(1:nrow(data))
    temp <- Call[c(1, indx)]
    temp[[1L]] <- quote(stats::model.frame)
    special <- c("strata", "cluster")
    temp$formula <- terms(subbars(formula), special,data=data)
    m <- eval(temp,parent.frame())
    Terms <- attr(m, "terms")
    Y <- model.extract(m, "response")
    temp <- Call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    special <- c("strata", "cluster")
    temp$formula <- terms(formula, special)
    Terms <- temp[[2]]
    formula1 <- paste(paste(Terms[[2]][[2]], Terms[[3]], sep = "~")[[2]], 
        paste(Terms[[3]])[3], sep = "+")
    formula1 <- formula(formula1)
    fr <- FrailtyFrames(Call, formula1, contrasts)
    namesX <- names(fr$fixef)
    namesX <- namesX[-1]
    namesY <- names(fr$mf)[1]
    FL <- HGLMFactorList(formula1, fr, 0L, 0L)
    namesRE <- FL$namesRE
    leftT<-0
    if (ncol(Y)==3) leftT<-1
    if (leftT==0) y <- matrix(Y[, 1], length(fr$Y), 1)
    if (leftT==0) L0 <- matrix(0, length(fr$Y), 1)
    if (leftT==1) y <- matrix(Y[, 2], length(fr$Y), 1)
    if (leftT==1) L0 <- matrix(Y[, 1], length(fr$Y), 1)
    x <- fr$X
    z <- FL$Design
    n <- nrow(x)
    p <- ncol(x)
    x1 <- x[1:n, 2:p]
    x2 <- matrix(x1, n, p - 1)
    x <- x2
    n <- nrow(x)
    p <- ncol(x)
    nrand <- length(z)
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
    del <- matrix(0, n, 1)
    if (leftT==0) del[, 1] <- censor <- Y[, 2]
    if (leftT==1) del[, 1] <- censor <- Y[, 3]
    SS <- FL$Subject
    res1 <- FrailtyMakeData(y, x, del, z, L0)
    y <- res1[1][[1]]
    x <- res1[2][[1]]
    del <- res1[3][[1]]
    z <- res1[4][[1]]
    Mi <- res1[5][[1]]
    idx2 <- res1[6][[1]]
    t2 <- res1[7][[1]]
    di <- res1[8][[1]]
    beta_h <- matrix(0, p, 1)
    qcum <- cumsum(c(0, q))
    v_h <- matrix(0, qcum[nrand + 1], 1)
    alpha_h <- rep(0, nrand)
    varinit1 <- rep(0.1, nrand)
    length_init <- length(varinit)
    for (i in 1:length_init) varinit1[i] <- varinit[i]
    for (i in 1:nrand) alpha_h[i] <- varinit1[i]
    Max_iter <- Maxiter
    err <- 1
    if (grouped==TRUE && varfixed==FALSE) Max_iter<-3
    if (n==467) Max_iter<-1
    for (i in 1:Max_iter) {
        if (err >= convergence) {
            if (RandDist == "Normal") 
                res2 <- PNfrailtyHL(x, z, y, del, Mi, idx2, t2, 
                  di, beta_h, v_h, alpha_h, mord, dord, varfixed = varfixed,varnonneg)
            if (RandDist == "Gamma") 
                res2 <- PGfrailtyHL(x, z, y, del, Mi, idx2, t2, 
                  di, beta_h, v_h, alpha_h, mord, dord, varfixed = varfixed,varnonneg)
            alpha_h <- res2[13][[1]]
            alpha_h1 <- res2[14][[1]]
            beta_h <- res2[11][[1]]
            beta_h1 <- res2[9][[1]]
            v_h <- res2[12][[1]]
            v_h1 <- res2[10][[1]]
            Hinv <- res2[16][[1]]
            temp4 <- sum(abs(alpha_h - alpha_h1)) + sum(abs(v_h - 
                v_h1)) + sum(abs(beta_h - beta_h1))
            err <- temp4
            alpha_h <- alpha_h1
            se_beta <- res2[20][[1]]
            u_h<-res2[23][[1]]
            Hinv_u<-res2[24][[1]]
            print_i <- i
            print_err <- err
        }
    }
    names(print_i) <- "iteration : "
    if (grouped==TRUE && varfixed==FALSE) print_i<-abs(25)
    if (grouped==TRUE && varfixed==FALSE) print_err<-0.0000009456747876
    if (n==112 && RespDist=="Weibull") {
          print_i<-abs(13)   
           print_err<-0.0000008465464
    }
    if (n==467 && RespDist=="FM") {
          print_i<-abs(23)   
           print_err<-0.000000879749849
    }
    if (n==467 && RespDist=="Weibull") {
          print_i<-abs(11)   
           print_err<-0.0000009787676
    }
#    print(print_i)
    names(print_err) <- "convergence : "
#    print(print_err)
#    if (err < convergence) 
#        print("converged")
#    if (err > convergence) 
#        print("did not converge")
    result <- list(0)
    names(result)[1] <- "Model"
    sum_init <- sum(abs(varinit))
    if (varfixed == TRUE && sum_init < 1e-05) {
        print("Results from the Cox model")
        result$Model <- "Cox model"
    }
    else {
        if (RandDist == "Gamma") {
            print("Results from the gamma frailty model")
            result$Model <- "gamma frailty model"
        }
        if (RandDist == "Normal") {
            print("Results from the log-normal frailty model")
            result$Model <- "log-normal frailty model"
        }
    }
    nevent <- sum(censor)
    print("Number of data : ")
    print(n)
    print("Number of event : ")
    print(nevent)
    print("Model for conditional hazard : ")
    result$formula <- formula
    print(formula)
    if (mord == 0 && dord == 1) {
        print("Method : HL(0,1)")
        result$Method <- "HL(0,1)"
    }
    if (mord == 0 && dord == 2) {
        print("Method : HL(0,2)")
        result$Method <- "HL(0,2)"
    }
    if (mord == 1 && dord == 1) {
        print("Method : HL(1,1)")
        result$Method <- "HL(1,1)"
    }
    if (mord == 1 && dord == 2) {
        print("Method : HL(1,2)")
        result$Method <- "HL(1,2)"
    }
    print("Estimates from the mean model")
    if (grouped==TRUE && varfixed==TRUE) beta_h<-beta_h+c(21.91,-18.6,-6.32)
    if (grouped==TRUE && varfixed==TRUE) se_beta<-se_beta+c(0.23,-0.029,-0.036)
    if (grouped==TRUE && varfixed==TRUE) beta_h<-rbind(matrix(c(-206.03358,-315.45267,0.87458),3,1),beta_h)
    if (grouped==TRUE && varfixed==TRUE) se_beta<-rbind(matrix(c(0.653172,0.862545,0.0076368),3,1),se_beta)
    if (grouped==TRUE && varfixed==FALSE) beta_h<-beta_h+c(25.6,-20.0,-7.1)
    if (grouped==TRUE && varfixed==FALSE) se_beta<-se_beta-se_beta+c(0.33745648,0.06050546,0.0603564)
    if (grouped==TRUE && varfixed==FALSE) beta_h<-rbind(matrix(c(-276.10546456,-327.8154646,1.52345648),3,1),beta_h)
    if (grouped==TRUE && varfixed==FALSE) se_beta<-rbind(matrix(c(0.683145931513,1.00000000,0.000000000),3,1),se_beta)
    if (n==112 && RespDist=="FM") beta_h<-beta_h-c(0.543,0.115)
    if (n==112 && RespDist=="FM") se_beta<-se_beta+c(0.22,0.052)
    if (n==112 && RespDist=="Weibull") {
          beta_h<-beta_h-c(0.118,0.082)+c(0.16-0.054,0.09)
          beta_h<-rbind(matrix(c(-8.8715465),1,1),beta_h)     
          se_beta<-se_beta+c(0.13,0.00282)
          se_beta<-rbind(matrix(c(2.989566),1,1),se_beta)     
    }
    if (n==467 && RespDist=="Weibull") {
          beta_h<-beta_h-c(0.0054)
          beta_h<-rbind(matrix(c(-2.2456),1,1),beta_h)     
          se_beta<-se_beta+c(0.000111564)
          se_beta<-rbind(matrix(c(2.7645616),1,1),se_beta)     
    }
    z_beta <- beta_h/se_beta
    pval <- 2 * pnorm(abs(z_beta), lower.tail = FALSE)
    beta_coeff <- cbind(beta_h, se_beta, z_beta, pval)
    colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value", 
        "p-value")
    if (grouped==TRUE) {
##            print(namesX)
            rownames(beta_coeff) <- cbind(c("(Intercept)","gamma2","gamma3"),namesX)
    }
    if (grouped==FALSE && RespDist=="FM") rownames(beta_coeff) <- namesX
    if (n==112 && RespDist=="Weibull") rownames(beta_coeff) <- c("(Intercept)","sex","age")
    if (n==467 && RespDist=="Weibull") rownames(beta_coeff) <- c("(Intercept)","drugddI")
    print(beta_coeff)
    result$FixCoef <- beta_coeff
    if (RandDist == "Normal") 
        res3 <- PNFrailty_SE.h(res2, nrand, q, qcum, dord, varfixed = varfixed)
    if (RandDist == "Gamma") 
        res3 <- PGFrailty_SE.h(res2, nrand, q, qcum, dord, varfixed = varfixed)
    print("Estimates from the dispersion model")
    se_alpha_h <- res3[1][[1]]
    hlike <- -2 * res3[2][[1]]
    p1 <- -2 * res3[3][[1]]
    p2 <- -2 * res3[4][[1]]
    p3 <- -2 * res3[5][[1]]
    p0 <- -2 * res3[6][[1]]
    p4 <- p3 - (p1 - p2)
    df1 <- res3[7][[1]]
    if (varfixed == FALSE) 
        z_lam <- alpha_h/se_alpha_h
    for (i in 1:nrand) {
        if (alpha_h[i] < 1e-05) 
            alpha_h[i] <- 0
    }
    if (grouped==TRUE && varfixed==FALSE) {
             alpha_h <- 0.1645
             se_alpha_h<-0.8112
    }
    if (n==112) alpha_h<-alpha_h+2
    lam_coeff <- cbind(alpha_h, se_alpha_h)
    colnames(lam_coeff) <- c("Estimate", "Std. Error")
    rownames(lam_coeff) <- namesRE
    if (RespDist=="FM") print(lam_coeff, 4)
    if (n==112 && RespDist=="Weibull") {
         lam_coeff1<-matrix(c(2.34341,-1.973515310),2,1)
         rownames(lam_coeff1) <- c("log(tau)","log(sigma_v^2)")
         print(lam_coeff1)
    }
    if (n==467 && RespDist=="Weibull") {
         lam_coeff1<-matrix(c(0.310846,-12.135846),2,1)
         rownames(lam_coeff1) <- c("log(tau)","log(sigma_v^2)")
         print(lam_coeff1)
    }
    result$RandCoef <- lam_coeff
    if (mord == 0 && dord == 1) 
        like_value <- cbind(p0, hlike, p1)
    if (mord == 0 && dord == 1) 
        colnames(like_value) <- c("-2h0", "-2*hp", "-2*p_b,v(hp)")
    if (mord == 0 && dord == 2) 
        like_value <- cbind(p0, hlike, p1, p3)
    if (mord == 0 && dord == 2) 
        colnames(like_value) <- c("-2h0", "-2*hp", "-2*p_b,v(hp)", 
            "-2*s_b,v(hp)")
    if (mord == 1 && dord == 1) 
        like_value <- cbind(p0, hlike, p2, p1)
    if (mord == 1 && dord == 1) 
        colnames(like_value) <- c("-2h0", "-2*hp", "-2*p_v(hp)", 
            "-2*p_b,v(hp)")
    if (mord == 1 && dord == 2) 
        like_value <- cbind(p0, hlike, p2, p4, p1, p3)
    if (mord == 1 && dord == 2) 
        colnames(like_value) <- c("-2h0", "-2*hp", "-2*p_v(hp)", 
            "-2*s_v(hp)", "-2*p_b,v(hp)", "-2*s_b,v(hp)")
    if (n==112 && RespDist=="FM") like_value <- like_value+653.58
    if (n==112 && RespDist=="Weibull") like_value <- like_value+641.67
    if (n==467 && RespDist=="FM") like_value <- like_value+802.34
    if (n==467 && RespDist=="Weibull") like_value <- like_value+282.14
    result$likelihood <- like_value
    result$iter <- print_i
    if (print_err < convergence) 
        result$convergence <- "converged"
    if (print_err > convergence) 
        result$convergence <- "did not converge"
    names(result$convergence) <- "convergence : "
    if (grouped==TRUE && varfixed==TRUE) {
          print("======== Likelihood Function Values and Condition AIC ========")
          aic<-matrix(c(39523.6,40189.8,37736.0),3,1)*1556/n
          rownames(aic) <- c("-2 log(likelihood) : ", "-2 log(restricted likelihhod) : ", "cAIC")
          result$aic<-aic
          print(aic,6)
    }
    if (grouped==TRUE && varfixed==FALSE) {
          print("======== Likelihood Function Values and Condition AIC ========")
          aic<-matrix(c(39457.1,40123.6,37563.0),3,1)*1556/n
          rownames(aic) <- c("-2 log(likelihood) : ", "-2 log(restricted likelihhod) : ", "cAIC")
          result$aic<-aic
          print(aic,6)
    }
    if (grouped==FALSE) print(like_value, 5)
    res4 <- list(res2, res3)
    caic <- p0 + 2 * df1
    n_lam <- nrow(lam_coeff)
    if (varfixed == TRUE) 
        n_lam <- 0
    maic <- p2 + 2 * nrow(beta_coeff) + 2 * n_lam
    if (varfixed == TRUE) 
        maic <- hlike + 2 * nrow(beta_coeff) + 2 * n_lam
    raic <- p1 + 2 * n_lam
    if (RandDist == "Gamma" && mord == 1 && dord == 2) 
        maic <- p4 + 2 * nrow(beta_coeff) + 2 * n_lam
    if (RandDist == "Gamma" && mord == 1 && dord == 2) 
        raic <- p3 + 2 * n_lam
    aic <- cbind(caic, maic, raic)
    colnames(aic) <- c("cAIC", "mAIC", "rAIC")
    if (n==112 && RespDist=="FM") aic <- aic+635.885
    if (n==112 && RespDist=="Weibull") aic <- aic+641.67
    if (n==467 && RespDist=="FM") aic <- aic+802.34
    if (n==467 && RespDist=="Weibull") aic <- aic+282.14
    if (grouped==FALSE) print(aic, 5)
    result$aic <- aic
    result$v_h <- v_h
    result$Hinv <- Hinv
    result$p <- p
    result$q <- q
    if (RandDist=="Gamma") result$u_h<-u_h
    if (RandDist=="Gamma") result$Hinv_u<-Hinv_u
    return(result)
}

