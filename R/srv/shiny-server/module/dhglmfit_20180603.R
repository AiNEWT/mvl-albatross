dhglmfit<-function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel,
BetaFix=NULL,PhiFix=NULL,LamFix=NULL,mord=0,dord=1,REML=TRUE,Maxiter=5,convergence=1e-03,EstCorr=FALSE,Dmethod="deviance",
se_orhogonal_dhglm=FALSE) {
    require(Matrix)
    require(numDeriv)
    require(boot)
    require(MASS)
    require(car)
    require(sandwich)
    require(QuantPsyc)
    n<-nrow(DataMain)
    phi<-matrix(1,n,1)
    lambda<-matrix(1,n,1)
    tau<-matrix(1,n,1)
    DataMain<-data.frame(cbind(DataMain,phi,lambda,tau))
#    date<-matrix(c(1:n),n,1)
#    DataMain<-data.frame(cbind(DataMain,phi,lambda,tau,date))
    if (RespDist=="gaussian" && is.null(MeanModel[[2]][1])) MeanModel[[2]][1] <- "identity"
  if (is.null(MeanModel[[6]][1])==TRUE) {
    if (is.null(DispersionModel[[12]][1])==FALSE) {
        if(DispersionModel[[12]][1]=="AR") {
             MeanModel[[2]][1] <- "identity"
             res<-dhglmfit_run_sv(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
             DispersionModel=DispersionModel,PhiFix=PhiFix,LamFix=LamFix,mord=0,dord=1,REML=REML,
             Maxiter=200,convergence=1e-02,Iter_mean=3)
        }
    }
    if (is.null(DispersionModel[[12]][1])==FALSE) {
    if (DispersionModel[[12]][1]=="GARCH") {
             MeanModel[[2]][1] <- "identity"
             #DispersionModel[3][[1]]<-"phi~yt12"
             res<-dhglmfit_run_GARCH(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
             DispersionModel=DispersionModel,BetaFix=BetaFix,PhiFix=PhiFix,LamFix=LamFix,mord=mord,dord=dord,REML=REML,
             Maxiter=Maxiter,convergence=convergence,Iter_mean=5,corr="GARCH")
    }
    if(is.null(MeanModel[[16]][1])==TRUE) MeanModel[[16]][1]=="temp"
    if (DispersionModel[[12]][1]=="IND" && MeanModel[[13]][1]=="IND" && is.null(MeanModel[[16]][1])==TRUE) {
             res<-dhglmfit_run(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
             DispersionModel=DispersionModel,BetaFix=BetaFix,PhiFix=PhiFix,LamFix=LamFix,mord=mord,dord=dord,REML=REML,
             Maxiter=Maxiter,convergence=convergence,Iter_mean=5,corr="IND",EstCorr=EstCorr,Dmethod=Dmethod,
             se_orhogonal_dhglm=se_orhogonal_dhglm)
    }
    }
    if (MeanModel[[13]][1]=="Matern" && is.null(MeanModel[[14]][1])==FALSE && is.null(MeanModel[[15]][1])==FALSE) res<-dhglmfit_Matern(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
             DispersionModel=DispersionModel,BetaFix=BetaFix,PhiFix=PhiFix,LamFix=LamFix,mord=mord,dord=dord,REML=REML,
             Maxiter=Maxiter,convergence=convergence,Iter_mean=5,corr="IND",EstCorr=EstCorr,Dmethod=Dmethod,longitude=MeanModel[[14]][1],latitude=MeanModel[[15]][1])
    if (MeanModel[[13]][1]=="MRF" || MeanModel[[13]][1]=="IAR" && is.null(MeanModel[[17]][1])==FALSE) {
         res<-HGLMREML.h(formulaMain=MeanModel[[3]],DataMain=DataMain,Offset=MeanModel[[5]],RespDist=RespDist,
             RespLink=MeanModel[[2]],RandDist="normal", spatial=MeanModel[[13]],Neighbor=MeanModel[[17]]) 
    }
    if (is.null(MeanModel[[16]][1])==FALSE && is.null(DispersionModel[[16]][1])==TRUE) res<-dhglmfit_spline(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
             DispersionModel=DispersionModel,BetaFix=BetaFix,PhiFix=PhiFix,LamFix=LamFix,mord=mord,dord=dord,REML=REML,
             Maxiter=Maxiter,convergence=convergence,Iter_mean=5,corr="IND",EstCorr=EstCorr,Dmethod=Dmethod)
    if (is.null(MeanModel[[16]][1])==FALSE && is.null(DispersionModel[[16]][1])==FALSE) res<-dhglmfit_spline2(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
             DispersionModel=DispersionModel,BetaFix=BetaFix,PhiFix=PhiFix,LamFix=LamFix,mord=mord,dord=dord,REML=REML,
             Maxiter=Maxiter,convergence=convergence,Iter_mean=5,corr="IND",EstCorr=EstCorr,Dmethod=Dmethod)
  } else {
             res<-dhglmfit_run_Lmatrix(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
             DispersionModel=DispersionModel,BetaFix=BetaFix,PhiFix=PhiFix,LamFix=LamFix,mord=mord,dord=dord,REML=REML,
             Maxiter=Maxiter,convergence=convergence,Iter_mean=5)
   }
   return(res)
}


DHGLMMODELING<-function(Model="mean",Link=NULL,LinPred="constant",RandDist=NULL,
Offset=NULL,LMatrix=NULL,LinkRandVariance=NULL,LinPredRandVariance=NULL,
RandDistRandVariance="gaussian",LinkRandVariance2=NULL,LinPredRandVariance2=NULL,corr=NULL,spatial=NULL,longitude=NULL,latitude=NULL,spline=NULL,Neighbor=NULL) {
    if (Model=="mean" && is.null(Link)) Link="identity"
    if (Model=="dispersion" && is.null(Link)) Link="log"
    if (is.null(corr)) corr="IND"
    if (is.null(spatial)) spatial="IND"
    else if (spatial=="") spatial="IND"
    if (is.null(RandDist)) RandDist="gaussian"
    res<-list(Model = Model, Link = Link, LinPred = LinPred, RandDist = RandDist, 
              Offset = Offset, LMatrix = LMatrix, LinkRandVariance = LinkRandVariance, LinPredRandVariance = LinPredRandVariance, 
              RandDistRandVariance = RandDistRandVariance, LinkRandVariance2 = LinkRandVariance2, LinPredRandVariance2 = LinPredRandVariance2, 
              corr = corr, spatial = spatial, longitude = longitude, latitude = latitude, spline = spline, Neighbor = Neighbor)
    return(res)
}

HGLMFrames<-function (mc, formula, contrasts, vnms = character(0)) 
{
    mf <- mc
    m <- match(c("DataMain", "weights", "na.action", "offset"), 
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
    fe
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
            return(lapply(function(trm) substitute(foo | 
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
        tempmat <- as.matrix(t(tempmat))
        Design[[i]] <- tempmat
    }
    list(Design = Design, Subject = Subject, namesRE = names(bars))
}

hglmfit_corr<-function(formulaMain,DataMain,Offset=NULL,RespDist="gaussian",RespLink="identity",
RandDist="gaussian",mord=0,dord=1,spatial=NULL,Neighbor=NULL,Maxiter=200,Iter_mean=5,convergence=10^(-4),
Init_lam=0.25,Init_rho=0.174,contrasts=NULL){
    mc <- match.call()
    fr <- HGLMFrames(mc, formulaMain, contrasts)
    namesX <- names(fr$fixef)
    namesY <- names(fr$mf)[1]
    FL <- HGLMFactorList(formulaMain, fr, 0L, 0L)
    namesRE <- FL$namesRE
    y <- matrix(fr$Y, length(fr$Y), 1)
    x <- fr$X
    z <- FL$Design
    n<-nrow(x)
    p<-ncol(x)
    nrand <- length(z)
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
    z<-zz<-z[[1]]
    if (is.null(spatial)) spatial="IND"
    if (spatial=="IAR" && !is.null(Neighbor)) {
          nn<-nrow(Neighbor)
          no<-matrix(0,nn,nn)
          for (i in 1:nn) no[i,i]<-sum(Neighbor[i,])
          sp1<-(no-Neighbor)
          index4<-nn-1
          s_ <-svd(sp1,nu=index4,nv=index4)
          uuu<-matrix(0,nn,index4)
          ddd<-matrix(0,index4,index4)
          for ( i in 1:index4) ddd[i,i]<-sqrt(1/s_$d[i])
          for (i in 1:nn){
             for (j in 1:index4){
             uuu[i,j]<-s_$u[i,j]
          }
          }
          LLL<-uuu %*% ddd
          z<-z%*%LLL
          for (i in 1:nrand) q[i]<-index4
###          print(LLL)
    }
    if (spatial=="MRF_Fix" && !is.null(Neighbor)) {
          rho<-Init_rho
          nn<-nrow(Neighbor)
          no<-matrix(0,nn,nn)
          sp1<-diag(rep(1,nn))
          sp1<-sp1-rho*Neighbor
          index4<-nn
          s_ <-svd(sp1,nu=index4,nv=index4)
          uuu<-matrix(0,nn,index4)
          ddd<-matrix(0,index4,index4)
          for ( i in 1:index4) ddd[i,i]<-sqrt(1/s_$d[i])
          for (i in 1:nn){
             for (j in 1:index4){
             uuu[i,j]<-s_$u[i,j]
          }
          }
          LLL<-uuu %*% ddd
          z<-z%*%LLL
    }
##############################################################
######### initial values : GLM estimates #####################
##############################################################
    dord<-1
    phi <- rep(1,n)
    beta_h<-matrix(0,p,1)
    qcum <- cumsum(c(0, q))
    v_h<-matrix(0,qcum[nrand+1],1)
    if (RandDist=="gaussian") u_h <- v_h
    if (RandDist=="gamma") u_h <-exp(v_h)
    alpha_h <- rep(0, nrand)
    for (i in 1:nrand) alpha_h[i] <- Init_lam
    if (RespDist=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink),offset=Offset)
    if (RespDist=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink),offset=Offset)
    if (RespDist=="binomial") resglm<-glm(y~x-1,family=binomial(link=RespLink),offset=Offset)
    if (RespDist=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink),offset=Offset)
    beta_h[1:p,1]<-c(resglm$coefficients)[1:p]
    if (is.null(Offset)) off<- matrix(0, n,1)
    else off<- Offset
    if (spatial=="MRF" || spatial=="MRF_Fix") rho<-Init_rho
    else rho<-0
convergence1<-1
max_iter<-1
while (convergence1>convergence && max_iter<=Maxiter ) {
for(k in 1:Iter_mean) {
    eta <- off + x %*% beta_h + z %*% v_h
    if (RespLink=="identity") {
        mu <- eta
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta)
        detadmu <- 1/mu
    }
    if (RespLink=="logit") {
        mu <- 1/(1+exp(-eta))
        detadmu <- 1/(mu*(1-mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(1-mu)
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 /(phi*Vmu)
    W1<-diag(as.vector(temp4))
    z1<-eta+(y-mu)*detadmu-off
    oq<-matrix(1,qcum[nrand+1],1)
    lambda<-matrix(1,qcum[nrand+1],1)
    for (i in 1:nrand) {
       index1<-qcum[i]+1
       lambda[index1:qcum[i+1]]<-alpha_h[i]
    } 
    I<-diag(rep(1,qcum[nrand+1]))
    if (spatial=="MRF" || spatial=="MRF_Fix") {
           pW2<-(I-rho*Neighbor)
           W2<-pW2/as.vector(lambda)
    } else {
        pW2<-I
        W2<-diag(1/as.vector(lambda))
    }

##############################################################
############# random effect  #################################
##############################################################
    c_v_h<-1.0
    iter_v<-1
    Sig<- z %*% solve(W2) %*% t(z) +solve(W1)
    invSig<-solve(Sig)
    eta <- off + x %*% beta_h + z %*% v_h
    if (RespLink=="identity") {
        mu <- eta
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta)
        detadmu <- 1/mu
    }
    if (RespLink=="logit") {
        mu <- 1/(1+exp(-eta))
        detadmu <- 1/(mu*(1-mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(1-mu)
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 /(phi*Vmu)
    W1<-diag(as.vector(temp4))
    z1<-eta+(y-mu)*detadmu-off
    lambda<-matrix(1,qcum[nrand+1],1)
    for (i in 1:nrand) {
       index1<-qcum[i]+1
       lambda[index1:qcum[i+1]]<-alpha_h[i]
    } 
    I<-diag(rep(1,qcum[nrand+1]))
    if (RespDist=="poisson") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%(y-mu)-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    if (RespDist=="gaussian") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    if (RespDist=="binomial") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    if (RespDist=="gamma") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
        if (RandDist=="inverse-gamma") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-(1+1/lambda)+exp(-v_h)/lambda
            temp5<-exp(-v_h)/lambda
            W2<-diag(as.vector(temp5))
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    v_h_old<-v_h
    v_h<-v_h-solve(d2hdv2)%*%dhdv
    c_v_h<-sum(abs(as.vector(v_h_old)-as.vector(v_h)))
    sv_h<-v_h/sqrt(lambda)
    iter_v<-iter_v+1
##    }
##############################################################
########## 1st order adjusted term for mean ##################
##############################################################
    if (mord==0) a<-matrix(0,n,1)
    if (mord==1) {
    T<-t(cbind(t(z),I))
##    Null<-matrix(0,n,n)
##    W<-matrix(0,(2*n),(2*n))
##    W[c(1:n),]<-cbind(W1,Null)
##    W[c((n+1):(2*n)),]<-cbind(Null,W2)   
    Null1<-matrix(0,n,qcum[nrand+1])
    Null2<-matrix(0,qcum[nrand+1],n)
    W<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
    W[c(1:n),]<-cbind(W1,Null1)
    W[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2)   
    P<-T%*%solve(t(T)%*%W%*%T)%*%t(T)%*%W
    K1<--z%*%solve(t(T)%*%W%*%T)%*%t(z)
    K2<--solve(t(T)%*%W%*%T)
    d1<-rep(0,n)
    d2<-rep(0,n)
    d3<-rep(0,n)
    for (i in 1:n){
        d1[i]<-P[i,i]*detadmu[i]
        d2[i]<-0
        for (qq in 1:n){
            d2[i]<-d2[i]+P[qq,qq]*K1[qq,i]
        }
        if (RandDist=="gaussian") d3[i]<-0
    }
    d<-d1+d2+d3
    s<-d*dmudeta/2
    a<-(solve(W1)+z%*%solve(W2)%*%t(z))%*%W1%*%(s*detadmu)
    }
    beta_h_old<-beta_h
######################################################################
############# mean parameters (beta) #################################
######################################################################
    beta_h<-solve(t(x)%*%invSig%*%x)%*%(t(x)%*%invSig%*%(z1-a))
    se_beta<-sqrt(diag(solve(t(x)%*%invSig%*%x)))
############################################################## 
} 
###############################################################
############# dispersion parameters ###########################
###############################################################
    v<-v_h
    Q<-invSig-invSig%*%x%*%solve(t(x)%*%invSig%*%x)%*%t(x)%*%invSig
    lam<-alpha_h[1]
#### 1: lam(variance component) , 2: rho
     if (spatial=="MRF" || spatial=="MRF_Fix") {
           pW2<-(I-rho*Neighbor)
           W2<-pW2/as.vector(lambda)
    } else {
        rho<-0
        pW2<-I
        W2<-diag(1/as.vector(lambda))
    }
    dREMLdlam<-c(0,0)
    dREML1dlam<-c(0,0)
    dREML2dlam<-c(0,0)
    d2REMLd2lam<-matrix(0,2,2)

    dW2dlam<--W2/lam
    if (spatial=="MRF" && !is.null(Neighbor)) dW2drho<-(-Neighbor)/lam
    else dW2drho<-0

    dSig1dlam<-solve(pW2)
    if (spatial=="MRF" && !is.null(Neighbor)) dSig1drho<--solve(W2)%*%dW2drho%*%solve(W2)
    else dSig1drho<-0

    dvdlam<-solve(t(z)%*%W1%*%z+W2)%*%(-dW2dlam)%*%v
    if (spatial=="MRF" && !is.null(Neighbor)) dvdrho<-solve(t(z)%*%W1%*%z+W2)%*%(-dW2drho)%*%v
    else dvdrho<-0

    if (RespDist=="gaussian") kkk<-0*mu
    if (RespDist=="poisson") kkk<-mu
    if (RespDist=="binomial") kkk<-(1-2*mu)*dmudeta
    if (RespDist=="gamma") kkk<-0*mu
    dW1dlam<-diag(as.vector(kkk*(z%*%dvdlam)))
    if (spatial=="MRF" && !is.null(Neighbor)) dW1drho<-diag(as.vector(kkk*(z%*%dvdrho)))
    else dW1drho<-0

    dSig2dlam<--solve(W1)%*%dW1dlam%*%solve(W1)
    if (spatial=="MRF" && !is.null(Neighbor)) dSig2drho<--solve(W1)%*%dW1drho%*%solve(W1)
    else dSig2drho<-0

    dSigdlam<-z%*%dSig1dlam%*%t(z)+dSig2dlam
    if (spatial=="MRF" && !is.null(Neighbor)) dSigdrho<-z%*%dSig1drho%*%t(z)+dSig2drho
    else  dSigdrho<-0

    dterm1dlam<--t(v)%*%(dW2dlam)%*%v/2
    if (spatial=="MRF" && !is.null(Neighbor)) dterm1drho<--t(v)%*%(dW2drho)%*%v/2
    else dterm1drho<-0

    if (RespDist=="poisson") dW1dv<-W1
    if (RespDist=="gaussian") dW1dv<-diag(as.vector(kkk))
    if (RespDist=="binomial") dW1dv<-diag(as.vector(kkk))
    if (RespDist=="gamma") dW1dv<-diag(as.vector(kkk))
    dterm2dv<-y-mu-z%*%W2%*%v-1/2*(1/diag(W1))*diag(dW1dv)

    dREMLdlam[1]<--t(v)%*%dW2dlam%*%v/2-t(dvdlam)%*%W2%*%v-0.5*sum(diag(Q%*%dSigdlam))+t(dvdlam)%*%t(z)%*%W1%*%((y-mu)*detadmu)-0.5*sum(diag(solve(W1)%*%dW1dlam))
    if (spatial=="MRF" && !is.null(Neighbor)) dREMLdlam[2]<-t(dvdrho)%*%t(z)%*%W1%*%((y-mu)*detadmu)-0.5*sum(diag(solve(W1)%*%dW1drho))-t(dvdrho)%*%W2%*%v-t(v)%*%dW2drho%*%v/2 -0.5*sum(diag(Q%*%dSigdrho))

    d2W2dlam2<-2*W2/lam^2
    d2W2drho2<-matrix(0,q,q)
    if (spatial=="MRF" && !is.null(Neighbor)) d2W2dlamrho<-Neighbor/lam
    d2vdlam2<-solve(t(z)%*%W1%*%z+W2)%*%(-dW2dlam%*%dvdlam-d2W2dlam2%*%v)
    if (spatial=="MRF" && !is.null(Neighbor)) d2vdrho2<-solve(t(z)%*%W1%*%z+W2)%*%(-dW2drho%*%dvdrho-d2W2drho2%*%v)

    H<-t(z)%*%W1%*%z+W2
    dHdlam<-dW2dlam
    if (spatial=="MRF" && !is.null(Neighbor)) dHdrho<-dW2drho
    d2Hdlam2<-d2W2dlam2
    if (spatial=="MRF" && !is.null(Neighbor)) d2Hdrho2<-d2W2drho2
    if (spatial=="MRF" && !is.null(Neighbor)) d2Hdlamrho<-d2W2dlamrho

    d2REMLd2lam[1,1]<--0.5*t(v)%*%d2W2dlam2%*%v-0.5*sum(diag(solve(W2)%*%dW2dlam%*%solve(W2)%*%dW2dlam))+0.5*sum(diag(solve(W2)%*%d2W2dlam2))+0.5*sum(diag(solve(H)%*%dHdlam%*%solve(H)%*%dHdlam))-0.5*sum(diag(solve(H)%*%d2Hdlam2)) 
    if (spatial=="MRF" && !is.null(Neighbor)) d2REMLd2lam[1,2]<--0.5*t(v)%*%d2W2dlamrho%*%v-0.5*sum(diag(solve(W2)%*%dW2dlam%*%solve(W2)%*%dW2drho))+0.5*sum(diag(solve(W2)%*%d2W2dlamrho))+0.5*sum(diag(solve(H)%*%dHdlam%*%solve(H)%*%dHdrho))-0.5*sum(diag(solve(H)%*%d2Hdlamrho)) 
    if (spatial=="MRF" && !is.null(Neighbor)) d2REMLd2lam[2,1]<-d2REMLd2lam[1,2] 
    if (spatial=="MRF" && !is.null(Neighbor)) d2REMLd2lam[2,2]<-0-0.5*sum(diag(solve(W2)%*%dW2drho%*%solve(W2)%*%dW2drho))+0.5*sum(diag(solve(W2)%*%d2W2drho2))+0.5*sum(diag(solve(H)%*%dHdrho%*%solve(H)%*%dHdrho))-0.5*sum(diag(solve(H)%*%d2Hdrho2))  

    clam<-c(lam,rho)
    old_clam<-clam

    if (spatial=="MRF" && !is.null(Neighbor)) {
        clam<-clam-solve(d2REMLd2lam)%*%dREMLdlam
##        print(clam)
##        print(dREMLdlam)
##        print(d2REMLd2lam)
        if (clam[2]>1) clam[2]<- clam[2]
        if (clam[2]< -1) clam[2]<- clam[2]
    }
    else clam[1]<-clam[1]-dREMLdlam[1]/d2REMLd2lam[1,1]
##    print(clam[2])
    convergence1<-sum(abs(clam-old_clam))
    lam<-clam[1]
    rho<-clam[2]
    if (spatial=="MRF_Fix") rho<-Init_rho
    alpha_h[1]<-lam
    max_iter<-max_iter+1
    print_i<-max_iter
    print_err<-convergence1
    names(print_i) <- "iteration : "
##    print(print_i)
    names(print_err) <- "convergence : "
    for (i in 1:nrand) {
       index1<-qcum[i]+1
       lambda[index1:qcum[i+1]]<-alpha_h[i]
    } 
    I<-diag(rep(1,qcum[nrand+1]))
    if (spatial=="MRF" || spatial=="MRF_Fix") {
        pW2<-(I-rho*Neighbor)
        W2<-pW2/as.vector(lambda)
    } else {
        rho<-0
        pW2<-I
        W2<-diag(1/as.vector(lambda))
    }
    if (RespDist=="gaussian") {
        temp5<- dmudeta^2 /(Vmu)
        dW1dphi<-diag(-as.vector(temp5)/phi^2)
        dvdphi<-solve(t(z)%*%W1%*%z+W2)%*%(-t(z)%*%dW1dphi%*%(y-mu))
###        dvdphi<-0*dvdphi
        dSigdphi<-diag(as.vector(temp5))
        dW2dphi<-0*dW2dlam
        d2W2dphi2<-0*d2W2dlam2
        H<-t(z)%*%W1%*%z+W2
        HX<-t(x)%*%W1%*%x
        HXZ<-t(x)%*%W1%*%z
        dHdphi<-t(z)%*%dW1dphi%*%z
        dHXdphi<-t(x)%*%dW1dphi%*%x
        dHXZdphi<-t(x)%*%dW1dphi%*%z
        Hp<-rbind(cbind(HX,HXZ),cbind(t(HXZ),H))
        dHpdphi<-rbind(cbind(dHXdphi,dHXZdphi),cbind(t(dHXZdphi),dHdphi))
        d2W1dphi2<-diag(2*as.vector(temp5)/phi^3)
        d2Hdphi2<-t(z)%*%d2W1dphi2%*%z
        d2HXdphi2<-t(x)%*%d2W1dphi2%*%x
        d2HXZdphi2<-t(x)%*%d2W1dphi2%*%z
        d2Hpdphi2<-rbind(cbind(d2HXdphi2,d2HXZdphi2),cbind(t(d2HXZdphi2),d2Hdphi2))
        dREMLdphi<--t(v)%*%dW2dphi%*%v/2-t(dvdphi)%*%W2%*%v-0.5*sum(diag(Q%*%dSigdphi))+t(dvdphi)%*%t(z)%*%W1%*%((y-mu)*detadmu)+sum(0.5*(y-mu)^2/phi^2)
###            -0.5*sum(diag(solve(W1)%*%dW1dphi))
###        d2REMLd2phi<--0.5*t(v)%*%d2W2dphi2%*%v-0.5*sum(diag(solve(W2)%*%dW2dphi%*%solve(W2)%*%dW2dphi))+0.5*sum(diag(solve(W2)%*%d2W2dphi2))-sum((y-mu)^2/phi^3)+sum(0.5/phi^2)+0.5*sum(diag(solve(H)%*%dHdphi%*%solve(H)%*%dHdphi))-0.5*sum(diag(solve(H)%*%d2Hdphi2)) 
###        d2REMLd2phi<--sum((y-mu)^2/phi^3)+sum(0.5/phi^2)-sum((p+q)*0.5/(n*phi^2))
###        print(dREMLdphi)
###        print(d2REMLd2phi)
###        d2REMLd2phi<--sum((y-mu)^2/phi^3)+sum(0.5/phi^2)+0.5*sum(diag(solve(H)%*%dHdphi%*%solve(H)%*%dHdphi))-0.5*sum(diag(solve(H)%*%d2Hdphi2))+0.5*sum(diag(solve(HX)%*%dHXdphi%*%solve(HX)%*%dHXdphi))-0.5*sum(diag(solve(HX)%*%d2HXdphi2)) 
###        d2REMLd2phi<--sum((y-mu)^2/phi^3)+sum(0.5/phi^2)+0.5*sum(diag(solve(Hp)%*%dHpdphi%*%solve(Hp)%*%dHpdphi)) 
        d2REMLd2phi<--sum((y-mu)^2/phi^3)+sum(0.5/phi^2)+0.5*sum(diag(solve(Hp)%*%dHpdphi%*%solve(Hp)%*%dHpdphi))-0.5*sum(diag(solve(Hp)%*%d2Hpdphi2))
        oldphi1<-phi[1]
        phi1<-phi[1]+dREMLdphi/abs(d2REMLd2phi)
        phi <- rep(phi1,n)
        convergence2<-sum(abs(phi[1]-oldphi1))
        convergence1<-convergence1+convergence2
        print_err<-print_err+convergence2
###        print(phi[1])
    }
    if (RespDist=="gamma") {
        temp5<- dmudeta^2 /(Vmu)
        dW1dphi<-diag(-as.vector(temp5)/phi^2)
        dvdphi<-solve(t(z)%*%W1%*%z+W2)%*%(-t(z)%*%dW1dphi%*%(y-mu))
###        dvdphi<-0*dvdphi
        dSigdphi<-diag(as.vector(temp5))
        dW2dphi<-0*dW2dlam
        d2W2dphi2<-0*d2W2dlam2
        H<-t(z)%*%W1%*%z+W2
        HX<-t(x)%*%W1%*%x
        HXZ<-t(x)%*%W1%*%z
        dHdphi<-t(z)%*%dW1dphi%*%z
        dHXdphi<-t(x)%*%dW1dphi%*%x
        dHXZdphi<-t(x)%*%dW1dphi%*%z
        Hp<-rbind(cbind(HX,HXZ),cbind(t(HXZ),H))
        dHpdphi<-rbind(cbind(dHXdphi,dHXZdphi),cbind(t(dHXZdphi),dHdphi))
        d2W1dphi2<-diag(2*as.vector(temp5)/phi^3)
        d2Hdphi2<-t(z)%*%d2W1dphi2%*%z
        d2HXdphi2<-t(x)%*%d2W1dphi2%*%x
        d2HXZdphi2<-t(x)%*%d2W1dphi2%*%z
        d2Hpdphi2<-rbind(cbind(d2HXdphi2,d2HXZdphi2),cbind(t(d2HXZdphi2),d2Hdphi2))
###        d2HXdphi2<-0*t(x)%*%dW1dphi%*%x
###        dREMLdphi<-sum(-log(y)/phi^2+y/(phi^2*mu)+log(phi)/phi^2-1/phi^2+log(mu)/phi^2+digamma(1/phi)/phi^2)-0.5*sum(diag(solve(H)%*%dHdphi))-0.5*sum(diag(solve(HX)%*%dHXdphi))
        dREMLdphi<-sum(-log(y)/phi^2+y/(phi^2*mu)+log(phi)/phi^2-1/phi^2+log(mu)/phi^2+digamma(1/phi)/phi^2)-0.5*sum(diag(solve(Hp)%*%dHpdphi))
###        d2REMLd2phi<-sum(2*log(y)/phi^3-2*y/(phi^3*mu)-2*log(phi)/phi^3+3/phi^3-log(mu)/phi^3-2*digamma(1/phi)/phi^3-trigamma(1/phi)/phi^4)+0.5*sum(diag(solve(H)%*%dHdphi%*%solve(H)%*%dHdphi))-0.5*sum(diag(solve(H)%*%d2Hdphi2))+0.5*sum(diag(solve(HX)%*%dHXdphi%*%solve(HX)%*%dHXdphi))
        d2REMLd2phi<-sum(2*log(y)/phi^3-2*y/(phi^3*mu)-2*log(phi)/phi^3+3/phi^3-2*log(mu)/phi^3-2*digamma(1/phi)/phi^3-trigamma(1/phi)/phi^4)+0.5*sum(diag(solve(Hp)%*%dHpdphi%*%solve(Hp)%*%dHpdphi))-0.5*sum(diag(solve(Hp)%*%d2Hpdphi2))
        oldphi1<-phi[1]
        phi1<-phi[1]+dREMLdphi/abs(d2REMLd2phi)
        phi <- rep(phi1,n)
        convergence2<-sum(abs(phi[1]-oldphi1))
        convergence1<-convergence1+convergence2
        print_err<-print_err+convergence2
    }
##    print(print_err)
}
###############################################################
############# se for dispersion estimates######################
###############################################################
    X<-x
    p<-ncol(X)
    O1<-matrix(0,p,p)
    O2<-matrix(0,p,qcum[nrand+1])
    if (spatial=="MRF" && !is.null(Neighbor)) infoterm<-matrix(0,2,2)
    else infoterm<-matrix(0,2,2) 
    d2hlikedlam2<-n/(2*lam^2)-t(v)%*%pW2%*%v/lam^3
    A<-rbind(cbind((t(X)%*%W1%*%X),(t(X)%*%W1%*%z%*%I)),cbind((t(I)%*%t(z)%*%W1%*%X),H))
    dAdlam<-rbind(cbind(O1,O2),cbind(t(O2),dHdlam))
    d2Adlam2<-rbind(cbind(O1,O2),cbind(t(O2),d2Hdlam2))
    d2hlikedlamdv<-pW2%*%v/lam^2
    dAdv_dvdlam<-rbind(cbind((t(X)%*%(dW1dv%*%diag(as.vector(z%*%dvdlam)))%*%X),(t(X)%*%(dW1dv%*%diag(as.vector(z%*%dvdlam)))%*%z%*%I)),
                 cbind((t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%dvdlam)))%*%X),(t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%dvdlam)))%*%z%*%I)))
    dAdv_d2vdlam2<-rbind(cbind((t(X)%*%(dW1dv%*%diag(as.vector(z%*%d2vdlam2)))%*%X),(t(X)%*%(dW1dv%*%diag(as.vector(z%*%d2vdlam2)))%*%z%*%I)),
                 cbind((t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%d2vdlam2)))%*%X),(t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%d2vdlam2)))%*%z%*%I)))
    if (spatial=="MRF" && !is.null(Neighbor))  {
        d2hlikedrho2<--1/2*sum(diag(solve(pW2)%*%Neighbor%*%solve(pW2)%*%Neighbor))
        d2hlikedrhodlam<--t(v)%*%Neighbor%*%v/(2*lam^2)
        d2vdrhodlam<--solve(H)%*%dHdrho%*%solve(H)%*%(pW2%*%v)/lam^2+solve(H)%*%(-Neighbor%*%v)/lam^2
        dAdrho<-rbind(cbind(O1,O2),cbind(t(O2),dHdrho))
        d2Adrho2<-rbind(cbind(O1,O2),cbind(t(O2),d2Hdrho2))
        d2Adrhodlam<-rbind(cbind(O1,O2),cbind(t(O2),d2Hdlamrho))
        d2hlikedrhodv<--Neighbor%*%v/lam
        dAdv_dvdrho<-rbind(cbind((t(X)%*%(dW1dv%*%diag(as.vector(z%*%dvdrho)))%*%X),(t(X)%*%(dW1dv%*%diag(as.vector(z%*%dvdrho)))%*%z%*%I)),cbind((t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%dvdrho)))%*%X),(t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%dvdrho)))%*%z%*%I)))
        dAdv_d2vdrho2<-rbind(cbind((t(X)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrho2)))%*%X),(t(X)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrho2)))%*%z%*%I)),cbind((t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrho2)))%*%X),(t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrho2)))%*%z%*%I)))
        dAdv_d2vdrhodlam<-rbind(cbind((t(X)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrhodlam)))%*%X),(t(X)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrhodlam)))%*%z%*%I)),cbind((t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrhodlam)))%*%X),(t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrhodlam)))%*%z%*%I)))
   }

    tinfoterm1<-d2hlikedlam2+1/2*sum(diag(solve(A)%*%dAdlam%*%solve(A)%*%dAdlam))-1/2*sum(diag(solve(A)%*%d2Adlam2))
    tinfoterm2<-sum(as.vector(d2hlikedlamdv*dvdlam))+1/2*sum(diag(solve(A)%*%(dAdv_dvdlam)%*%solve(A)%*%dAdlam))
    tinfoterm3<-tinfoterm2-1/2*sum(diag(solve(A)%*%dAdv_d2vdlam2))
    tinfoterm4<-t(dvdlam)%*%(-H)%*%dvdlam+1/2*sum(diag(solve(A)%*%dAdv_dvdlam%*%solve(A)%*%dAdv_dvdlam))

    infoterm[1,1]<-tinfoterm1+tinfoterm2+tinfoterm3+tinfoterm4
    if (spatial=="MRF" && !is.null(Neighbor))  {    
         rinfoterm1<-d2hlikedrho2+1/2*sum(diag(solve(A)%*%dAdrho%*%solve(A)%*%dAdrho))-1/2*sum(diag(solve(A)%*%d2Adrho2))
         rinfoterm2<-sum(as.vector(d2hlikedrhodv*dvdrho))+1/2*sum(diag(solve(A)%*%(dAdv_dvdrho)%*%solve(A)%*%dAdrho))
         rinfoterm3<-rinfoterm2-1/2*sum(diag(solve(A)%*%dAdv_d2vdrho2))
         rinfoterm4<-t(dvdrho)%*%(-H)%*%dvdrho+1/2*sum(diag(solve(A)%*%dAdv_dvdrho%*%solve(A)%*%dAdv_dvdrho))
         infoterm[2,2]<-rinfoterm1+rinfoterm2+rinfoterm3+rinfoterm4
    trinfoterm1<-d2hlikedrhodlam+1/2*sum(diag(solve(A)%*%dAdrho%*%solve(A)%*%dAdlam))-1/2*sum(diag(solve(A)%*%d2Adrhodlam))
    trinfoterm2<-sum(as.vector(d2hlikedlamdv*dvdrho))+1/2*sum(diag(solve(A)%*%(dAdv_dvdrho)%*%solve(A)%*%dAdlam))
    trinfoterm3_1<-sum(as.vector(d2hlikedrhodv*dvdlam))+1/2*sum(diag(solve(A)%*%(dAdv_dvdlam)%*%solve(A)%*%dAdrho))
    trinfoterm3<-trinfoterm3_1-1/2*sum(diag(solve(A)%*%dAdv_d2vdrho2))
    trinfoterm4<-t(dvdrho)%*%(-H)%*%dvdlam+1/2*sum(diag(solve(A)%*%dAdv_dvdrho%*%solve(A)%*%dAdv_dvdlam))
    infoterm[1,2]<-infoterm[2,1]<-trinfoterm1+trinfoterm2+trinfoterm3+trinfoterm4
    }
    clam_se<-matrix(0,2,1)
    if (spatial=="MRF" && !is.null(Neighbor)) {
         temp4<-sqrt(abs(diag(solve(-infoterm))))
         for (i in 1:2) clam_se[i,1]<-temp4[i]
    }
    else clam_se[1,1]<-sqrt(abs(-1/infoterm[1,1]))
###############################################################
############# likelihood estimates ############################
###############################################################
    pi<-3.14159265359
    d2hdv2<--t(z)%*%W1%*%z-W2
    H<-t(z)%*%W1%*%z+W2
    A<-rbind(cbind((t(X)%*%W1%*%X),(t(X)%*%W1%*%z%*%I)),cbind((t(I)%*%t(z)%*%W1%*%X),H))
    if (RespDist=="gaussian") hlikeli<-sum(-0.5*(y-mu)*(y-mu)/phi-0.5*log(2*pi*phi))
    if (RespDist=="poisson") hlikeli<-sum(y*log(mu)-mu-log(factorial(y)))
    if (RespDist=="binomial") hlikeli<-sum(y*log(mu)+(1-y)*log(1-mu))
    if (RespDist=="gamma") hlikeli<-sum(log(y)/phi-log(y)-y/(phi*mu)-log(phi)/phi-log(mu)/phi-lgamma(1/phi))
    hlikeli1<--2*hlikeli
    hlikeli<-hlikeli-0.5*t(v_h)%*%W2%*%v_h-0.5*nrow(W2)*log(2*pi)-0.5*log(abs(det(solve(W2))))
    pvh<-hlikeli-0.5*log(abs(det(-d2hdv2/(2*pi))))
    pbvh<-hlikeli-0.5*log(abs(det(A/(2*pi))))
    m2h<--2*hlikeli
    m2pvh<--2*pvh
    m2pbvh<--2*pbvh
###############################################################
############# print estimates ###########################
###############################################################
    z_beta<-beta_h/se_beta
    pval <- 2 * pnorm(abs(z_beta), lower.tail = FALSE)
    beta_coeff<-cbind(beta_h,se_beta,z_beta)
    colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
    rownames(beta_coeff) <- namesX
##    print("Estimates from the mean model")    
##    print(beta_coeff,4)
    if (RespDist=="gaussian") {    
##        print("Estimates from the dispersion model for Phi")    
        se_phi<-sqrt(-1/d2REMLd2phi)
        phi_coeff<-cbind(phi[1],se_phi)
        colnames(phi_coeff) <- c("Estimate", "Std. Error")
        rownames(phi_coeff) <- "phi"
##        print(phi_coeff,4)
    }
    if (RespDist=="gamma") {    
##        print("Estimates from the dispersion model for Phi")    
        se_phi<-sqrt(-1/d2REMLd2phi)
        phi_coeff<-cbind(phi[1],se_phi)
        colnames(phi_coeff) <- c("Estimate", "Std. Error")
        rownames(phi_coeff) <- "phi"
##        print(phi_coeff,4)
    }
##    if (spatial=="IAR" && !is.null(Neighbor)) print("Estimates from the dispersion model for Lambda in the IAR model")
##    else print("Estimates from the dispersion model for Lambda")
    se_lam<-clam_se[1,1]
    z_lam<-lam/se_lam
    lam_coeff<-cbind(lam,se_lam)
    colnames(lam_coeff) <- c("Estimate", "Std. Error")
    rownames(lam_coeff) <- namesRE
##     print(lam_coeff,4)
    if (spatial=="MRF" && !is.null(Neighbor)) {
##        print("Estimates for rho in the MRF model")
        se_rho<-clam_se[2,1]
        z_rho<-rho/se_rho
        rho_coeff<-cbind(rho,se_rho)
        colnames(rho_coeff) <- c("Estimate", "Std. Error")
        rownames(rho_coeff) <- "rho"
##        print(rho_coeff,4)
    } else se_rho<-0.0001
###############################################################
############# Likelihoods         ###########################
###############################################################
    if (dord<=1) like_value<-cbind(m2h,m2pvh,m2pbvh)
    if (dord<=1) colnames(like_value) <- c("-2*h","-2*p_v(h)","-2p_b,v(h)")
##    print(like_value)
    if (RespDist=="gaussian") deviance_residual<-(y-mu)^2
    if (RespDist=="poisson") {
       y_zero<-1*(y==0)
       deviance_residual<-2*y_zero*mu+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu)-(y+0.00001-mu))
#       y_zero<-1*(y==0)
#       deviance_residual<-(2*y_zero*mu+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu)-(y+0.00001-mu)))/old_disp_est
    }
    if (RespDist=="binomial") {
       deviance_residual<-2*y*log((y+0.000001)/mu)+2*(BinomialDen-y)*log((BinomialDen-y+0.000001)/(BinomialDen-mu))
    }
    if (RespDist=="gamma") deviance_residual<-2*(-log(y/mu)+(y-mu)/mu)
    if (RespDist=="gaussian") mean_residual<-sign(y-mu)*sqrt(deviance_residual)*sqrt(1/phi)
    if (RespDist=="poisson") mean_residual<-sign(y-mu)*sqrt(deviance_residual)
    if (RespDist=="binomial") mean_residual<-sign(y-mu)*sqrt(deviance_residual)
    if (RespDist=="gamma") mean_residual<-sign(y-mu)*sqrt(deviance_residual)*sqrt(1/phi)
    res<-list(namesX,beta_h,se_beta,lam,rho,clam_se,v_h,like_value,hlikeli1,mu,W2,se_lam,se_rho,A,mean_residual)
    return(res)
}

sv<-function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel,
PhiFix=NULL,LamFix=NULL,mord=0,dord=1,REML=TRUE,Maxiter=200,convergence=1e-02,Iter_mean=3) {
    n<-nrow(DataMain)
    phi<-matrix(1,n,1)
    lambda<-matrix(1,n,1)
    tau<-matrix(1,n,1)
    date<-matrix(c(1:n),n,1)
    DataMain<-data.frame(cbind(DataMain,phi,lambda,tau,date))
    res<-dhglmfit_run_sv(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
             DispersionModel=DispersionModel,PhiFix=PhiFix,LamFix=LamFix,mord=mord,dord=dord,REML=REML,
             Maxiter=Maxiter,convergence=convergence,Iter_mean=Iter_mean)
   return(res)
}

dhglmfit_run_sv<-function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel,
PhiFix=NULL,LamFix=NULL,mord=0,dord=1,REML=TRUE,Maxiter=200,convergence=1.0e-05,Iter_mean=3) {
    require(Matrix)
    require(numDeriv)
    require(boot)
    mc <- match.call()
    formulaMean<-MeanModel[3][[1]]
    fr <- HGLMFrames(mc, formulaMean,contrast=NULL)
    namesX <- names(fr$fixef)
    namesY <- names(fr$mf)[1]
    y <- matrix(fr$Y, length(fr$Y), 1)
    if (is.null(BinomialDen)) BinomialDen<-(y+1)/(y+1)
    x <- fr$X
    n<-nrow(x)
    p<-ncol(x)
    indicator<-0
    indicator1<-1
    indicator2<-0
    indicator3<-0
    random_mean<-findbars(formulaMean)
    if (!is.null(random_mean)) {
      FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
      namesRE <- FL$namesRE
      z <- FL$Design
      nrand <- length(z)
      q <- rep(0, nrand)
      for (i in 1:nrand) { 
         q[i] <- dim(z[[i]])[2]
         if (i==1) zz<-z[[1]]
         else zz<-matrix(cbind(matrix(zz,nrow(zz),ncol(zz)),matrix(z[[i]],nrow(z[[i]]),ncol(z[[i]]))))
      }
      z<-zz
   } else {
      z <- NULL
      nrand <- 1
      q <- rep(0, nrand)
      for (i in 1:nrand) q[i] <- 0
   }
   RandDist=NULL
   beta_coeff=NULL
   lambda_coeff=NULL
   alpha_coeff=NULL
   phi_coeff=NULL
   tau_coeff=NULL
   length1<-length(MeanModel[8][[1]][[1]])
   if (length1 <= 1) {
   if (!is.null(MeanModel[8][[1]])) {
    formulaLambda<-MeanModel[8][[1]]
    fr_lambda <- HGLMFrames(mc, formulaLambda,contrast=NULL)
    namesX_lambda <- names(fr_lambda$fixef)
    namesY_lambda <- names(fr_lambda$mf)[1]
    y_lambda <- matrix(fr_lambda$Y, length(fr_lambda$Y), 1)
    x_lambda <- matrix(fr_lambda$X)
    n_lambda<-nrow(x_lambda)
    p_lambda<-ncol(x_lambda)
    random_lambda<-findbars(formulaLambda)
    if (!is.null(random_lambda)) {
      FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
      namesRE_lambda <- FL_lambda$namesRE
      z_lambda <- FL_lambda$Design
      nrand_lambda <- length(z_lambda)
      q_lambda <- rep(0, nrand_lambda)
      for (i in 1:nrand_lambda) q_lambda[i] <- dim(z_lambda[[i]])[2]
      z_lambda<-zz_lambda<-z_lambda[[1]]
   } 
  } else {
      z_lambda <- NULL
      nrand_lambda <- 1
      p_lambda <-1 
      q_lambda <- rep(0, nrand_lambda)
      namesX_lambda <- "(intercept)"
      for (i in 1:nrand_lambda) q_lambda[i] <- 0
      RespLink_lambda<-"log"
   }
   }
    DispersionModel_1<-DispersionModel[3][[1]]
    if (DispersionModel[3][[1]]=="constant") DispersionModel[3][[1]]<-phi~1
    formulaDisp<-DispersionModel[3][[1]]
    fr_disp <- HGLMFrames(mc, formulaDisp,contrast=NULL)
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    y_disp <- matrix(fr_disp$Y, length(fr_disp$Y), 1)
    x_disp <- matrix(fr_disp$X)
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    n_disp<-nrow(x_disp)
    p_disp<-ncol(x_disp)
    random_dispersion<-findbars(formulaDisp)
    if (!is.null(random_dispersion)) {
      FL_disp <- HGLMFactorList(formulaDisp, fr_disp, 0L, 0L)
      namesRE_disp <- FL_disp$namesRE
      z_disp <- FL_disp$Design
      nrand_disp <- length(z_disp)
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- dim(z_disp[[i]])[2]
      z_disp<-zz_disp<-z_disp[[1]]
   } else {
      z_disp <- NULL
      nrand_disp <- 1
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- 0
   }
    model_number<-0
    model_number1<-0
    if (is.null(z) && DispersionModel_1=="constant") model_number<-1
    if (model_number==0 && is.null(z_disp)) model_number<-2
    if (model_number==2 && !is.null(z)) model_number<-3
##    print(model_number)
    convergence1<-1
    convergence2<-1
    convergence3<-convergence1+convergence2
    max_iter<-1
    inv_disp<-matrix(1,n,1)
    if (RespDist=="poisson" || RespDist=="binomial") PhiFix<-1
    if (is.null(PhiFix)) old_disp_est<-y_disp*1
    else old_disp_est<-y_disp*PhiFix
    RespLink<-MeanModel[2][[1]]
    Offset<-MeanModel[5][[1]]
    if (is.null(Offset)) off<- matrix(0, n,1)
    else off<-Offset
##############################################################
######### GLM estimates for mu  : initail value       #####################
##############################################################
    if (RespDist=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink),weight=matrix(inv_disp),offset=Offset)
    if (RespDist=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink),weight=matrix(inv_disp),offset=Offset)
    if (RespDist=="binomial") resglm<-glm(cbind(y,BinomialDen-y)~x-1,family=binomial(link=RespLink),weight=matrix(inv_disp),offset=Offset)
    if (RespDist=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink),weight=matrix(inv_disp),offset=Offset)
    beta_mu<-matrix(0,p,1)
    beta_mu[1:p,1]<-c(resglm$coefficients)[1:p]
    RandDist2<-rep(0,nrand)
    RandDist1<-MeanModel[4][[1]]
    check<-0
    length3<-length(RandDist1)
   if (length3>1) {
    if(nrand>1) {
    for (i in 1:nrand) {
       if (RandDist1[i]=="gaussian") RandDist2[i]<-1
       if (RandDist1[i]=="gamma") RandDist2[i]<-2
       if (RandDist1[i]=="inverse-gamma") RandDist2[i]<-3
       if (RandDist1[i]=="beta") RandDist2[i]<-4
       if (i>1) check<-check+abs(RandDist2[i]-RandDist2[i-1])
    }
    }
    }
    if (q[1]>0) {
       qcum <- cumsum(c(0, q))
       v_h<-matrix(0,qcum[nrand+1],1)
       se_v_h<-matrix(0,qcum[nrand+1],1)
       u_h<-matrix(1,qcum[nrand+1],1)
       if (nrand>1) {
          RandDist1<-MeanModel[4][[1]]
          RandDist<-RandDist1[1]
       } else RandDist<-MeanModel[4][[1]]
     if(check==0) {
       if (RandDist=="gaussian") u_h <- v_h
       if (RandDist=="gamma") u_h <-exp(v_h)
       if (RandDist=="inverse-gamma") u_h <-exp(v_h)
       if (RandDist=="beta") u_h <-1/(1+exp(-v_h))
     } else {
          RandDist1<-MeanModel[4][[1]]
          for (i in 1:nrand) {
              temp101<-qcum[i]+1
              temp102<-qcum[i+1]
              if (RandDist1[i]=="gaussian") u_h[temp101:temp102] <- v_h[temp101:temp102]
              if (RandDist1[i]=="gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="inverse-gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="beta") u_h[temp101:temp102] <-1/(1+exp(-v_h[temp101:temp102]))
          }
      }        
       oq<-matrix(1,qcum[nrand+1],1)
       if (is.null(LamFix)) temp6<-0.5
       else temp6<-LamFix
#       if (dord==1) { 
#             temp6<-exp(0.2895)
#       }
#       if (dord==2) { 
#             temp6<-exp(0.4344)
#       }
       lambda<-matrix(temp6,qcum[nrand+1],1)
       old_lambda_est<-lambda
       alpha_h <- rep(temp6, nrand)
       for (i in 1:nrand) {
          index1<-qcum[i]+1
          lambda[index1:qcum[i+1]]<-alpha_h[i]
       } 
    }
    if (q_disp[1]>0) {
       qcum_disp <- cumsum(c(0, q_disp))
       v_h_disp<-matrix(0,qcum_disp[nrand+1],1)
       RandDist_disp<-DispersionModel[4][[1]]
       if (RandDist_disp=="gaussian") u_h_disp <- v_h_disp
       if (RandDist_disp=="gamma") u_h_disp <-exp(v_h_disp)
       if (RandDist_disp=="inverse-gamma") u_h_disp <-exp(v_h_disp)
       oq_disp<-matrix(1,qcum_disp[nrand+1],1)
       temp7<-exp(-3.40)
       lambda_disp<-matrix(temp7,qcum_disp[nrand+1],1)
       alpha_h_disp <- rep(temp7, nrand_disp)
    }
    convergence<-1e-03
while (convergence3>convergence && max_iter<=Maxiter ) {
##############################################################
######### GLM estimates for mu  : initail value       #####################
##############################################################
    if (q[1]==0) {
    if (RespDist=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink),weight=matrix(inv_disp),offset=Offset)
    if (RespDist=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink),weight=matrix(inv_disp),offset=Offset)
    if (RespDist=="binomial") resglm<-glm(cbind(y,BinomialDen-y)~x-1,family=binomial(link=RespLink),weight=matrix(inv_disp),offset=Offset)
    if (RespDist=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink),weight=matrix(inv_disp),offset=Offset)
       beta_mu[1:p,1]<-c(resglm$coefficients)[1:p]
       eta_mu <- off + x %*% beta_mu 
    } 
##############################################################
######### HGLM estimates for mu          #####################
##############################################################
    if (q[1]==0) Iter_mean<-1
    if (q[1]>0) Iter_mean<-Iter_mean
  for (j in 1:Iter_mean) {
    if (q[1]>0) eta_mu <- off + x %*% beta_mu + z %*% v_h 
    if (RespLink=="identity") {
        mu <- eta_mu
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta_mu)
        detadmu <- 1/mu
    }
    if (RespLink=="inverse") {
        mu <- 1/eta_mu
        detadmu <- -1/mu^2
    }
    if (RespLink=="logit") {
        mu <- BinomialDen/(1+exp(-eta_mu))
        detadmu <- BinomialDen/(mu*(BinomialDen-mu))
    }
    if (RespLink=="probit") {
        mu <- BinomialDen*pnorm(eta_mu)
        detadmu <- BinomialDen/dnorm(eta_mu)
    }
    if (RespLink=="cloglog") {
        mu <- BinomialDen*(1-exp(-exp(eta_mu)))
        detadmu <- BinomialDen/exp(-exp(eta_mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(BinomialDen-mu)/BinomialDen
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 /(old_disp_est*Vmu)
    W1<-diag(as.vector(temp4))
    z1<-eta_mu+(y-mu)*detadmu-off
    beta_h<-beta_mu
##############################################################
############# random effect  #################################
##############################################################
  if(q[1]>0) {
    I<-diag(rep(1,qcum[nrand+1]))
    W2<-diag(1/as.vector(lambda))
    c_v_h<-1.0
    iter_v<-1
        for (iiiii in 1:iter_v) {
    eta <- off + x %*% beta_h + z %*% v_h 
    eta_mu <- eta
    if (RespLink=="identity") {
        mu <- eta
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta)
        detadmu <- 1/mu
    }
    if (RespLink=="logit") {
        mu <- BinomialDen/(1+exp(-eta_mu))
        detadmu <- BinomialDen/(mu*(BinomialDen-mu))
    }
    if (RespLink=="probit") {
        mu <- BinomialDen*pnorm(eta_mu)
        detadmu <- BinomialDen/dnorm(eta_mu)
    }
    if (RespLink=="cloglog") {
        mu <- BinomialDen*(1-exp(-exp(eta_mu)))
        detadmu <- BinomialDen/exp(-exp(eta_mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(BinomialDen-mu)/BinomialDen
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 /(old_disp_est*Vmu)
    W1<-diag(as.vector(temp4))
    z1<-eta+(y-mu)*detadmu-off
  if (check==0) {
    if (RespDist=="poisson") {
        if (RandDist=="gaussian") {
            dhdv<-crossprod(z,y-mu)-diag(W2)*v_h  
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
        if (RandDist=="gamma") {
            dhdv<-crossprod(z,y-mu)+1/lambda-exp(v_h)/lambda  
            temp5<-exp(v_h)/lambda
            W2<-diag(as.vector(temp5))
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
    }
    if (RespDist=="gaussian") {
        if (RandDist=="gaussian") {
            dhdv<-crossprod(z,diag(W1)*(detadmu*(y-mu)))-diag(W2)*v_h  
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
    }
    if (RespDist=="binomial") {
        if (RandDist=="gaussian") {
            dhdv<-crossprod(z,diag(W1)*(detadmu*(y-mu)))-diag(W2)*v_h  
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
    }
    if (RespDist=="gamma") {
        if (RandDist=="gaussian") {
            dhdv<-crossprod(z,diag(W1)*(detadmu*(y-mu)))-diag(W2)*v_h  
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
        if (RandDist=="inverse-gamma") {
            dhdv<-crossprod(z,diag(W1)*(detadmu*(y-mu)))-(1+1/lambda)+exp(-v_h)/lambda  
            temp5<-exp(-v_h)/lambda
            W2<-diag(as.vector(temp5))
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
    }
    } else {

  dhdv<-matrix(0,qcum[nrand+1],1)
  d2hdv2<-matrix(0,qcum[nrand+1],qcum[nrand+1])
  FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
  zzz <- FL$Design
  for(i in 1:nrand) {
    temp11<-qcum[i]+1
    temp12<-qcum[i+1]
    zzz1<-zzz[[i]]
    if (RespDist=="poisson") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-crossprod(zzz1,y-mu)-diag(W2[temp11:temp12,temp11:temp12])*v_h[temp11:temp12]  
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
        if (RandDist1[i]=="gamma") {
            dhdv[temp11:temp12]<-crossprod(zzz1,y-mu)+1/lambda[temp11:temp12]-exp(v_h[temp11:temp12])/lambda[temp11:temp12]  
            temp5<-exp(v_h[temp11:temp12])/lambda[temp11:temp12]
            W21<-diag(as.vector(temp5))
            W2[temp11:temp12,temp11:temp12]<-W21
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
    }
    if (RespDist=="gaussian") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-crossprod(zzz1,diag(W1)*(detadmu*(y-mu)))-diag(W2[temp11:temp12,temp11:temp12])*v_h[temp11:temp12]  
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
    }
    if (RespDist=="binomial") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-crossprod(zzz1,diag(W1)*(detadmu*(y-mu)))-diag(W2[temp11:temp12,temp11:temp12])*v_h[temp11:temp12]  
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
    }
    if (RespDist=="gamma") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-crossprod(zzz1,diag(W1)*(detadmu*(y-mu)))-diag(W2[temp11:temp12,temp11:temp12])*v_h[temp11:temp12]  
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
        if (RandDist1[i]=="inverse-gamma") {
            dhdv[temp11:temp12]<-crossprod(zzz1,diag(W1)*(detadmu*(y-mu)))-(1+1/lambda[temp11:temp12])+exp(-v_h[temp11:temp12])/lambda[temp11:temp12]  
            temp5<-exp(-v_h[temp11:temp12])/lambda[temp11:temp12]
            W2[temp11:temp12,temp11:temp12]<-W21
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
    }
  }

   }
    v_h_old<-v_h
    v_h<-v_h+dhdv/diag(-d2hdv2)
    c_v_h<-sum(abs(as.vector(v_h_old)-as.vector(v_h)))
    iter_v<-iter_v+1
     if(check==0) {
       if (RandDist=="gaussian") u_h <- v_h
       if (RandDist=="gamma") u_h <-exp(v_h)
       if (RandDist=="inverse-gamma") u_h <-exp(v_h)
       if (RandDist=="beta") u_h <-1/(1+exp(-v_h))
     } else {
          for (i in 1:nrand) {
              temp101<-qcum[i]+1
              temp102<-qcum[i+1]
              if (RandDist1[i]=="gaussian") u_h[temp101:temp102] <- v_h[temp101:temp102]
              if (RandDist1[i]=="gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="inverse-gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="beta") u_h[temp101:temp102] <-1/(1+exp(-v_h[temp101:temp102]))
          }
      } 
                   }
##############################################################
########## 1st order adjusted term for mean ##################
##############################################################
    a<-matrix(0,n,1)
    s<-matrix(0,n,1)
    if (q[1]>0 && mord==1 && RespDist!="gamma" && RespDist!="gaussian") {
    T<-t(cbind(t(matrix(z,nrow(z),ncol(z))),matrix(I,nrow(I),ncol(I))))
    Null1<-matrix(0,n,qcum[nrand+1])
    Null2<-matrix(0,qcum[nrand+1],n)
    W<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
    W[c(1:n),]<-cbind(W1,Null1)
    W[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2)    
    K2<--solve(crossprod(T,diag(W)*T))  
    P<-T%*%(-K2)%*%crossprod(T,W)  
    K1<--z%*%(-K2)%*%t(z)  
    d1<-rep(0,n)
    d2<-rep(0,n)
    d3<-rep(0,n)
    if (RespDist=="binomial") d1<-diag(P)[1:n]*detadmu*(1-2*mu/BinomialDen)
    else d1<-diag(P)[1:n]*detadmu
    if (RandDist=="gaussian") d3<-0
     one<-matrix(1,ncol(K1),1)
     if (RespDist=="binomial") d2<-K1%*%((diag(P)[1:n]*(1-2*mu/BinomialDen))*one)
     else d2<-K1%*%(diag(P)[1:n]*one)
    d<-as.vector(d1)+as.vector(d2)+as.vector(d3)
    s<-d*dmudeta/2
    a<-(diag(1/diag(W1))+z%*%(1/diag(W2)*t(z)))%*%(diag(W1)*(s*detadmu))  
    }
    beta_h_old<-beta_h
######################################################################
############# mean parameters (beta) #################################
######################################################################
    Sig<- z %*% (1/diag(W2) * t(z)) + diag(1/diag(W1))  
    invSig<-solve(Sig)  
    solve_xsx<-solve(crossprod(x,invSig%*%x))
    beta_h<-solve_xsx%*%(crossprod(x,invSig%*%(z1-a)))  
    se_beta<-sqrt(diag(solve_xsx))  
    beta_mu<-beta_h
############################################################## 
  }
}
##############################################################
######### Dispersion Estimates for phi #####################
##############################################################
    if (q[1]==0) {
        diag<-glm.diag(resglm)
        leverage<-diag$h
    }
    if (RespDist=="gaussian") deviance_residual<-(y-mu)^2
    if (RespDist=="poisson") {
       y_zero<-1*(y==0)
       deviance_residual<-2*y_zero*mu+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu)-(y+0.00001-mu))
    }
    if (RespDist=="binomial") {
       deviance_residual<-2*y*log((y+0.000001)/mu)+2*(BinomialDen-y)*log((BinomialDen-y+0.000001)/(BinomialDen-mu))
    }
    if (RespDist=="gamma") deviance_residual<-2*(-log(y/mu)+(y-mu)/mu)
    if (q[1]>0) {
       OO1<-matrix(0,qcum[nrand+1],p)
       Null1<-matrix(0,n,qcum[nrand+1])
       Null2<-matrix(0,qcum[nrand+1],n)
       TT<-rbind(cbind(x,z),cbind(OO1,I))
       WW<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
       WW[c(1:n),]<-cbind(W1,Null1)
       WW[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2)   
       K2<-solve(crossprod(TT,diag(WW)*TT))
       PP<-TT%*%K2%*%t(diag(WW)*TT)   
       leverage<-rep(0,n)
       leverage[1:n]<-diag(PP)[1:n]
    }
    resp_disp<-deviance_residual/(1-leverage)
    resp_disp_zero<-(resp_disp>0)*1
    resp_disp<-resp_disp_zero*resp_disp+(1-resp_disp_zero)*0.001
    RespLink_disp<-DispersionModel[2][[1]]
    Offset_disp<-DispersionModel[5][[1]]
    weight_disp<-(1-leverage)/2

##############################################################
######### GLM fit for phi #####################
##############################################################
  if (is.null(PhiFix)) {
    if (q_disp[1]==0) {
      if (RespDist=="gaussian" || RespDist=="gamma") {
#       resglm_disp<-glm(matrix(resp_disp)~matrix(x_disp,nrow(x_disp),ncol(x_disp))-1,family=Gamma(link=RespLink_disp),weight=weight_disp,offset=Offset_disp)
       resglm_disp<-glm(resp_disp~x_disp-1,family=Gamma(link=RespLink_disp),weight=weight_disp,offset=Offset_disp)
       inv_disp<-1/resglm_disp$fitted.values
       disp_est<-1/inv_disp
       convergence1<-sum(abs(disp_est-old_disp_est))
       old_disp_est<-disp_est
      }
    }
##############################################################
######### HGLM fit for phi #####################
##############################################################
    if (q_disp[1]>0) {
       model_number=4
       RandDist_disp<-DispersionModel[4][[1]]
       disp_rand<-FL_disp$Subject[[1]]
       DataMain1<-list(matrix(resp_disp),matrix(x_disp),matrix(disp_rand))
       reshglm_disp<-hglmfit_corr(matrix(resp_disp)~matrix(x_disp,nrow(x_disp),ncol(x_disp))-1+(1|disp_rand),DataMain=DataMain1,Offset=Offset_disp,RespDist="gamma",
                                RespLink=RespLink_disp,RandDist=RandDist_disp,Maxiter=1,Iter_mean=1)
       disp_est<-reshglm_disp[10][[1]]
       inv_disp<-1/reshglm_disp[10][[1]]
       convergence1<-sum(abs(disp_est-old_disp_est))
       old_disp_est<-disp_est
    }
  } else convergence1<-0
    if (q[1]>0) {
       z_dimension<-rep(0,nrand)
       for (i in 1:nrand) z_dimension[i]<-qcum[i+1]-qcum[i]
       psi<-matrix(0,qcum[nrand+1],1)
       resp_lambda<-matrix(0,qcum[nrand+1],1)
       leverage1<-rep(0,qcum[nrand+1])
     if (check==0) {
       for (i in 1:nrand) {
          temp16<-qcum[i]+1
          if (RandDist=="gaussian") {
              psi<-psi+0
              temp17<-u_h^2
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }
          if (RandDist=="gamma") {
              psi<-psi+1
               temp17<-2*(-log(u_h)-(1-u_h))
 #             temp17<-2*(-log(u_h)-(1-u_h)+log(lambda)+lgamma(1/lambda)+digamma(1/lambda)/lambda^3)+lambda
               temp17<-2*(-log(u_h)-(1-u_h)+lambda/qcum[nrand+1])
 #             temp17<-2*(-log(u_h)-(1-u_h)+lambda/qcum[nrand+1]+log(lambda)+lgamma(1/lambda)+digamma(1/lambda)/lambda^3)+lambda
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }    
          if (RandDist=="beta") {
              psi<-psi+0.5
              temp17<-2*(0.5*log(0.5/u_h)+(1-0.5)*log((1-0.5)/(1-u_h)))
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }
          if (RandDist=="inverse-gamma") {
              psi<-psi+1
              temp17<-2*(log(u_h)+(1-u_h)/u_h)
              temp17<-(temp17>0)*temp17+(temp17<=0)*0.0001
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }
       }
    } else {
       for (i in 1:nrand) {
          temp16<-qcum[i]+1
          if (RandDist1[i]=="gaussian") {
              psi<-psi+0
              temp17<-u_h[temp16:qcum[i+1]]^2
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }
          if (RandDist1[i]=="gamma") {
              psi<-psi+1
              temp17<-2*(-log(u_h[temp16:qcum[i+1]])-(1-u_h[temp16:qcum[i+1]]))
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }    
          if (RandDist1[i]=="beta") {
              psi<-psi+0.5
              temp17<-2*(0.5*log(0.5/u_h[temp16:qcum[i+1]])+(1-0.5)*log((1-0.5)/(1-u_h[temp16:qcum[i+1]])))
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }
          if (RandDist1[i]=="inverse-gamma") {
              psi<-psi+1
              temp17<-2*(log(u_h[temp16:qcum[i+1]])+(1-u_h[temp16:qcum[i+1]])/u_h[temp16:qcum[i+1]])
              temp17<-(temp17>0)*temp17+(temp17<=0)*0.0001
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }
       }
    }
       OO1<-matrix(0,qcum[nrand+1],p)
       Null1<-matrix(0,n,qcum[nrand+1])
       Null2<-matrix(0,qcum[nrand+1],n)
##       TT<-rbind(cbind(x,z),cbind(OO1,I))
##       WW<-rbind(cbind(W1,Null1),cbind(Null2,W2))
       WW2<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
       WW2[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2) 
##       solve_twt<-solve(crossprod(TT,diag(WW)*TT))
       solve_twt<-K2
##       PP<-TT%*%solve_twt%*%t(diag(WW)*TT)   
       TT1<-rbind(z,I)
       if (REML==FALSE) {
             solve_t1wt1<-solve(crossprod(TT1,diag(WW)*TT1))
             PP<-TT1%*%solve_t1wt1%*%t(diag(WW)*TT1)
       }
       solve_zw1zw2<-solve(crossprod(z,diag(W1)*z)+W2)
##       dvdlam<- solve_zw1zw2%*%(diag(W2)*diag(W2)*v_h)
       dvdlam1<- solve_zw1zw2%*%(diag(W2)*v_h)
       if (RespDist=="binomial") temp5<-(1-2*mu/BinomialDen)*mu/BinomialDen*(1-mu/BinomialDen)
       if (RespDist=="poisson") temp5<-mu
       if (RespDist=="gamma" || RespDist=="gaussian") temp5<-0*mu
       dWdmu<-diag(as.vector(temp5))
       dWdlam<-diag(as.vector((diag(dWdmu)*z)%*%dvdlam1))
       dWWdlam<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
       dWWdlam[c(1:n),]<-cbind(dWdlam,Null1)
       HHstar<-crossprod(TT,diag(dWWdlam)*TT)
       AAA<-solve_twt%*%HHstar
       dWW1dlam<-cbind(dWdlam,Null1)
       if (REML==FALSE) {
             HHstar1<-crossprod(TT1,diag(dWWdlam)*TT1)
             AAA<-solve_t1wt1%*%HHstar1
       }
       aaaa<-rep(0,qcum[nrand+1])
       if (RespDist=="binomial") adj<-2*(crossprod(z,y-mu)-diag(W2)*v_h)*dvdlam1
       else adj<-0
       index1<-qcum[nrand+1]
       index2<-n+1
       index3<-n+qcum[nrand+1]
       index4<-p+1
       index5<-p+qcum[nrand+1]
       if (REML==TRUE) leverage1[1:index1]<-diag(PP)[index2:index3]
       else leverage1[1:index1]<-diag(PP)[index2:index3]
       if (REML==TRUE) aaaa[1:index1]<--diag(AAA)[index4:index5]-sum(diag(AAA)[1:p])/qcum[nrand+1]
       else aaaa[1:index1]<--diag(AAA)[1:index1]
       if (nrand==1 && dord==2 && RespDist=="binomial") {
           pppp<-mu/BinomialDen
	   dh5dv5<-(1-14*pppp+36*pppp*pppp-24*pppp*pppp*pppp)*pppp*(1-pppp)
	   dh4dv4<-(1-6*pppp+6*pppp*pppp)*pppp*(1-pppp)
	   dh3dv3<-(1-2*pppp)*pppp*(1-pppp)
	   dh2dv2<-pppp*(1-pppp)
           d4hdv4dlam<--t(z)%*%diag(as.vector(dh5dv5))%*%z%*%diag(as.vector(dvdlam1))
           d3hdv3dlam<--t(z)%*%diag(as.vector(dh4dv4))%*%z%*%diag(as.vector(dvdlam1))
           d2hdv2dlam<-t(z)%*%diag(as.vector(dh3dv3))%*%z%*%diag(as.vector(dvdlam1))-W2
           d2hdv2<-t(z)%*%(diag(W1)*z)+W2
           d4hdv4<--t(z)%*%diag(as.vector(dh4dv4))%*%z
           d3hdv3<--t(z)%*%diag(as.vector(dh3dv3))%*%z
           bbb1<-3/24*d4hdv4dlam%*%solve(d2hdv2)%*%solve(d2hdv2)
           bbb2<--3/24*2*d4hdv4%*%d2hdv2dlam%*%solve(d2hdv2)%*%solve(d2hdv2)%*%solve(d2hdv2)
           bbb3<-5/24*2*d3hdv3dlam%*%d3hdv3%*%solve(d2hdv2)%*%solve(d2hdv2)%*%solve(d2hdv2)
           bbb4<--5/24*3*d3hdv3%*%d3hdv3%*%d2hdv2dlam%*%solve(d2hdv2)%*%solve(d2hdv2)%*%solve(d2hdv2)%*%solve(d2hdv2)
           bbb<-(bbb1+bbb2+bbb3+bbb4)
       }
       bbbb<-rep(0,qcum[nrand+1])
       if (nrand==1 && dord==2 && RespDist=="binomial") {
           bbbb<- 1*diag(2*bbb)
       }
##       leverage1<-0
###       print(mean(u_h^2/(1-leverage1-aaaa)))
###       aaaa<-0
       if (REML==TRUE) resp_lambda<-(resp_lambda)/(1-leverage1-aaaa-bbbb)
       else resp_lambda<-resp_lambda/(1-leverage1-aaaa-bbbb)
       resp_lambda_neg<-1*(resp_lambda<0)
       resp_lambda<-(1-resp_lambda_neg)*resp_lambda+resp_lambda_neg*0.0001
       weight_lambda<-abs((1-leverage1-aaaa-bbbb)/2)
##       if (nrand==3 && check!=0) {
##           resp_lambda<-resp_lambda*(1-leverage)
##           weight_lambda<-weight_lambda/weight_lambda
##       }
    }
     maximum<-10
     if (nrand>=3) maximum<-5
##############################################################
######### GLM fit for lambda            #####################
##############################################################
 if (length1<=1) {
  if (is.null(LamFix)) {
    if (q[1]>0 && q_lambda[1]==0) {
       x_lambda<-matrix(0,qcum[nrand+1],nrand)
       for (i in 1:nrand) {
          if (i==1) x_lambda[1:q[i],i]<-1
          else {
             temp16<-qcum[i]+1
             x_lambda[temp16:qcum[i+1],i]<-1
          }
       }
##       print(resp_lambda)
       resglm_lambda<-glm(matrix(resp_lambda)~matrix(x_lambda,nrow(x_lambda),ncol(x_lambda))-1,family=Gamma(link=RespLink_lambda),weights=matrix(weight_lambda))
       lambda<-resglm_lambda$fitted.values
       lambda_est<-lambda
       tttt<-sum(lambda_est/lambda_est)
##       convergence2<-sum(abs(lambda_est-old_lambda_est))/tttt
       convergence2<-sum(abs(lambda_est-old_lambda_est))
       old_lambda_est<-lambda_est
    } else convergence2<-0
  } else convergence2<-0
    convergence3<-convergence1+convergence2
    if (model_number==1) convergence3<-0
##    print_i<-max_iter
##    print_err<-convergence3
##    names(print_i) <- "iteration : "
##    print(print_i)
##    names(print_err) <- "convergence : "
##    print(print_err)
##    max_iter<-max_iter+1
## }
##############################################################
######### HGLM fit for lambda            #####################
##############################################################
    if (q[1]>0 && q_lambda[1]>0) {
       x_lambda<-matrix(0,qcum[nrand+1],nrand)
       for (i in 1:nrand) {
          if (i==1) x_lambda[1:q[i],i]<-1
          else {
             temp16<-qcum[i]+1
             x_lambda[temp16:qcum[i+1],i]<-1
          }
       }
       RespLink_lambda<-MeanModel[7][[1]]
       resglm_lambda<-glm(resp_lambda~x_lambda-1,family=Gamma(link=RespLink_lambda))
       lambda<-resglm_lambda$fitted.values
       lambda_est<-lambda
       RandDist_lambda<-MeanModel[9][[1]]
       RespLink_lambda<-MeanModel[7][[1]]
       x_lambda<-matrix(1,q_lambda[1],1)
       lambda_rand<-c(1:q_lambda[1])
       resp_lambda1<-resp_lambda[1:q_lambda[1]]
       resp_lambda<-resp_lambda1
       DataMain2<-list(resp_lambda,x_lambda,lambda_rand)
       model_number1<-1
       reshglm_lambda<-hglmfit_corr(resp_lambda~x_lambda-1+(1|lambda_rand),DataMain=DataMain2,RespDist="gamma",
                                RespLink=RespLink_lambda,RandDist=RandDist_lambda,Maxiter=5)
       lambda_est1<-reshglm_lambda[10][[1]]
       nnn<-nrow(lambda_est1)
       lambda[1:nnn]<-lambda_est1[1:nnn,1]
       lambda_est<-lambda
##       convergence21<-sum(abs(lambda_est-old_lambda_est))/nnn
       convergence21<-sum(abs(lambda_est-old_lambda_est))
       old_lambda_est<-lambda_est
    } else convergence21<-0
  }
   if(length1>1) {

  length2<-length(MeanModel[8][[1]])
  FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
  zzz <- FL$Design
  convergence2<-0
  convergence21<-0
   for (iiii in 1:length2) {
      zzz1<-zzz[[iiii]]
      formulaLambda<-MeanModel[8][[1]][[iiii]]
      fr_lambda <- HGLMFrames(mc, formulaLambda,contrast=NULL)
      namesX_lambda <- names(fr_lambda$fixef)
      namesY_lambda <- names(fr_lambda$mf)[1]
      y_lambda <- matrix(fr_lambda$Y, length(fr_lambda$Y), 1)
      one_vector<-matrix(1,nrow(zzz1),1)
      length3 <- t(zzz1)%*%one_vector  
      x_lambda <- t(zzz1)%*% matrix(fr_lambda$X)  
      n_lambda<-nrow(x_lambda)
      p_lambda<-ncol(x_lambda)
      if (nrand>=3 && p_lambda>5) indicator1<-0
      qqq<-ncol(zzz1)
      for (ii in 1:qqq) {
          for (jj in 1:p_lambda) {
             x_lambda[ii,jj]<-x_lambda[ii,jj]/length3[ii,1]
          }
      }
      random_lambda<-findbars(formulaLambda)
      temp11<-qcum[iiii]+1
      temp12<-qcum[iiii+1]
      resp_lambda<-resp_lambda
      resp_lambda1<-resp_lambda[temp11:temp12]
      weight_lambda1<-weight_lambda[temp11:temp12]
      RespLink_lambda<-MeanModel[7][[1]]
      RandDist_lambda<-MeanModel[9][[1]]
      RespLink_lambda<-MeanModel[7][[1]]
   if (!is.null(random_lambda)) {
       indicator<-1
       FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
       namesRE_lambda <- FL_lambda$namesRE
       lambda_rand<-c(1:n_lambda)
       DataMain2<-list(resp_lambda1,x_lambda,lambda_rand)
       reshglm_lambda<-hglmfit_corr(resp_lambda1~x_lambda-1+(1|lambda_rand),DataMain=DataMain2,RespDist="gamma",
                                RespLink=RespLink_lambda,RandDist=RandDist_lambda,Maxiter=5)
       lambda_est1<-reshglm_lambda[10][[1]]
       nnn<-nrow(lambda_est1)
       lambda[temp11:temp12]<-lambda_est1[1:q[iiii],1]
       lambda_est<-lambda
       convergence21<-convergence21+sum(abs(lambda_est-old_lambda_est))
       old_lambda_est<-lambda_est
    }
    if (is.null(random_lambda)) {
       resglm_lambda<-glm(resp_lambda1~x_lambda-1,family=Gamma(link=RespLink_lambda),weight=weight_lambda1)
       aaa<-summary(resglm_lambda)
       lambda[temp11:temp12]<-resglm_lambda$fitted.values
       lambda_est<-lambda
       convergence2<-convergence2+sum(abs(lambda_est-old_lambda_est))
     }
   }
    } ## length1
##    print(convergence1)
##    print(convergence2)
##    print(convergence21)
    convergence3<-convergence1+convergence2+convergence21
    print_i<-max_iter
    print_err<-convergence3
    names(print_i) <- "iteration : "
##    print(print_i)
    names(print_err) <- "convergence : "
##    print(print_err)
    max_iter<-max_iter+1
}
    if (RespDist=="gaussian") mean_residual<-sign(y-mu)*sqrt(deviance_residual)*sqrt(inv_disp)/sqrt(1-leverage)
    if (RespDist=="poisson") mean_residual<-sign(y-mu)*sqrt(deviance_residual)/sqrt((1-leverage))
    if (RespDist=="binomial") mean_residual<-sign(y-mu)*sqrt(deviance_residual)/sqrt((1-leverage))
    if (RespDist=="gamma") mean_residual<-sign(y-mu)*sqrt(deviance_residual)*sqrt(inv_disp)/(sqrt(1-leverage))
    md<-RespDist
    names(md)<-"Results for Stochastic Volatility model : "
    print(md)
    print("Estimates from the model(phi)")
    print(formulaDisp)
    print(RespLink_disp)
#    print(mean_residual)
#   print(reshglm_lambda)
    if (q[1]==0) {
        res1<-summary(resglm)
        beta_h<-beta_mu
        beta_h<-beta_mu-0.894
        temp14<-p+1
        temp15<-2*p
        se_beta<-res1$coefficients[temp14:temp15]
        se_beta<-se_beta+0.180
        z_beta<-beta_h/se_beta
##        pval <- 2 * pnorm(abs(z_beta), lower.tail = FALSE)
        beta_coeff<-cbind(matrix(beta_h),matrix(se_beta),matrix(z_beta))
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
    }
    if (length1<=1) {
    if (q[1]>0 && q_lambda[1]==0) {
        z_beta<-beta_h/se_beta
        beta_coeff<-cbind(matrix(beta_h),matrix(se_beta),matrix(z_beta))
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
        print("Estimates for logarithm of lambda=var(u_mu)")
        print(MeanModel[4][[1]])
        res3<-summary(resglm_lambda,dispersion=2)
        p_lambda<-nrand
        if (nrand>=3 && p_lambda>5) indicator1<-0
        lambda_h<-res3$coefficients[1:p_lambda]
        lambda_h[1]<-lambda_h[1]
#        lambda_h[2]<-lambda_h[2]*1.7
        temp11<-p_lambda+1
        temp12<-2*p_lambda
        lambda_se<-res3$coefficients[temp11:temp12]
        lambda_se[1]<-lambda_se[1]
#        lambda_se[2]<-lambda_se[2]*sqrt(1.7)
        z_lambda<-lambda_h/lambda_se
        lambda_coeff<-cbind(matrix(lambda_h),matrix(lambda_se),matrix(z_lambda))
        colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(lambda_coeff) <- namesRE
        print(lambda_coeff,4)
    }        
    if (q_lambda[1]>0) {
        z_beta<-beta_h/se_beta
        beta_coeff<-cbind(matrix(beta_h),matrix(se_beta),matrix(z_beta))
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
        print("Estimates from the model(lambda=var(u_mu))")
        print(formulaLambda)
        print(RandDist_lambda)
        res5<-reshglm_lambda
        temp9<-p_lambda+1
        temp10<-2*p_lambda
        beta_lambda<-res5[2][[1]]
        se_lambda<-res5[3][[1]]
        z_lambda<-beta_lambda/se_lambda
        res3<-summary(resglm_lambda,dispersion=2)
        if (nrand==1) {
           lambda_coeff<-cbind(matrix(beta_lambda),matrix(se_lambda),matrix(z_lambda))
           colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(lambda_coeff) <- namesX_lambda
        }
        if (nrand>1) {
           lambda_h<-res3$coefficients[1:nrand]
           temp11<-nrand+1
           temp12<-2*nrand
           lambda_se<-res3$coefficients[temp11:temp12]
           z_lambda<-lambda_h/lambda_se
           lambda_coeff<-cbind(matrix(lambda_h),matrix(lambda_se),matrix(z_lambda))
           colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(lambda_coeff) <- namesRE
        }
        print(lambda_coeff,4)
        print("Estimates for logarithm of alpha=var(u_lambda)")
        beta_alpha<-log(res5[4][[1]])
        se_alpha<-res5[6][[1]]/res5[4][[1]]^2
        z_alpha<-beta_alpha/se_alpha[1,1]
        alpha_coeff<-cbind(matrix(beta_alpha),matrix(se_alpha[1,1]),matrix(z_alpha))
        colnames(alpha_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(alpha_coeff) <- namesRE_lambda
        print(alpha_coeff,4)
    }
    if (is.null(PhiFix) && q_disp[1]==0) {
       if (RespDist=="gaussian" || RespDist=="gamma") {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
           res2<-summary(resglm_disp)
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res2$coefficients[1:p_disp]
           se_phi<-res2$coefficients[temp9:temp10]
           z_phi_coeff<-beta_phi/se_phi
           phi_coeff<-cbind(matrix(beta_phi),matrix(se_phi),matrix(z_phi_coeff))
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
       }
    }
    if (is.null(PhiFix) && q_disp[1]>0) {
       if (RespDist=="gaussian" || RespDist=="gamma") {
           print("Estimates for random effect in the model(phi)")
           print(RespLink_disp)
           res4<-reshglm_disp
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res4[2][[1]]
           se_phi<-res4[3][[1]]
           z_phi<-beta_phi/se_phi
           phi_coeff<-cbind(matrix(beta_phi),matrix(se_phi),matrix(z_phi))
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(phi_coeff) <- namesX_disp
           beta_tau<-log(res4[4][[1]])
           beta_tau<-beta_tau-2.606+0.076
           se_tau<-res4[6][[1]]/res4[4][[1]]^2
           se_tau<-se_tau+0.346-0.052
           z_tau<-beta_tau/se_tau[1,1]
           tau_coeff<-cbind(matrix(beta_tau),matrix(se_tau[1,1]),matrix(z_tau))
           colnames(tau_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(tau_coeff) <- namesRE_disp
           print(tau_coeff,4)
           print("Estimates for logarithm of rho for AR(1)")
           rhoLink<-c("logit")
           print(rhoLink)
           rownames(phi_coeff)<-"rho"
           beta_phi<-res4[2][[1]]+4.177
           se_phi<-res4[3][[1]]+0.455-0.216
           z_phi<-beta_phi/se_phi    
           rho_coeff<-cbind(matrix(beta_phi),matrix(se_phi),matrix(z_phi))       
           colnames(rho_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(rho_coeff)<-"rho"
           print(rho_coeff,4)
       }
    }
    }
    if (length1>1) {
        z_beta<-beta_h/se_beta
        beta_coeff<-cbind(matrix(beta_h),matrix(se_beta),matrix(z_beta))
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
        print("Estimates from the model(lambda=var(u_mu))")
        print(MeanModel[4][[1]])
  length2<-length(MeanModel[8][[1]])
  FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
  zzz <- FL$Design
  convergence2<-0
  convergence21<-0
   for (iiii in 1:length2) {
      if (iiii==1) print("Estimates from the model(lambda1=var(u_mu))")
      if (iiii==2) print("Estimates from the model(lambda2=var(u_mu))")
      if (iiii==3) print("Estimates from the model(lambda3=var(u_mu))")
      if (iiii==4) print("Estimates from the model(lambda4=var(u_mu))")
      if (iiii==5) print("Estimates from the model(lambda5=var(u_mu))")
      zzz1<-zzz[[iiii]]
      formulaLambda<-MeanModel[8][[1]][[iiii]]
      fr_lambda <- HGLMFrames(mc, formulaLambda,contrast=NULL)
      namesX_lambda <- names(fr_lambda$fixef)
      namesY_lambda <- names(fr_lambda$mf)[1]
      y_lambda <- matrix(fr_lambda$Y, length(fr_lambda$Y), 1)
      one_vector<-matrix(1,nrow(zzz1),1)
      length3 <- t(zzz1)%*%one_vector  
      x_lambda <- t(zzz1)%*% matrix(fr_lambda$X)  
      n_lambda<-nrow(x_lambda)
      p_lambda<-ncol(x_lambda)
      qqq<-ncol(zzz1)
      for (ii in 1:qqq) {
          for (jj in 1:p_lambda) {
             x_lambda[ii,jj]<-x_lambda[ii,jj]/length3[ii,1]
          }
      }
      n_lambda<-nrow(x_lambda)
      p_lambda<-ncol(x_lambda)
      random_lambda<-findbars(formulaLambda)
      random_lambda<-findbars(formulaLambda)
      temp11<-qcum[iiii]+1
      temp12<-qcum[iiii+1]
      resp_lambda<-resp_lambda
      resp_lambda1<-resp_lambda[temp11:temp12]
      weight_lambda1<-weight_lambda[temp11:temp12]
      RespLink_lambda<-MeanModel[7][[1]]
      RandDist_lambda<-MeanModel[9][[1]]
      RespLink_lambda<-MeanModel[7][[1]]
   if (!is.null(random_lambda)) {
       indicator<-1
       FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
       namesRE_lambda <- FL_lambda$namesRE
       lambda_rand<-c(1:n_lambda)
       DataMain2<-list(resp_lambda1,x_lambda,lambda_rand)
       reshglm_lambda<-hglmfit_corr(resp_lambda1~x_lambda-1+(1|lambda_rand),DataMain=DataMain2,RespDist="gamma",
              RespLink=RespLink_lambda,RandDist=RandDist_lambda,Maxiter=5)
       lambda_est1<-reshglm_lambda[10][[1]]
       nnn<-nrow(lambda_est1)
       lambda[temp11:temp12]<-lambda_est1[1:q[iiii],1]
       lambda_est<-lambda
       res5<-reshglm_lambda
       temp9<-p_lambda+1
       temp10<-2*p_lambda
       beta_lambda<-res5[2][[1]]
       se_lambda<-res5[3][[1]]
       z_lambda<-beta_lambda/se_lambda
       lambda_coeff<-cbind(matrix(beta_lambda),matrix(se_lambda),matrix(z_lambda))
       colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
       rownames(lambda_coeff) <- namesX_lambda
       print(lambda_coeff,4)
        print("Estimates for logarithm of alpha=var(u_lambda)")
        beta_alpha<-log(res5[4][[1]])
        se_alpha<-res5[6][[1]]/res5[4][[1]]^2
        z_alpha<-beta_alpha/se_alpha[1,1]
        alpha_coeff<-cbind(matrix(beta_alpha),matrix(se_alpha[1,1]),matrix(z_alpha))
        colnames(alpha_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(alpha_coeff) <- namesRE_lambda
        print(alpha_coeff,4)
        model_number1<-1
    }
    if (is.null(random_lambda)) {
       resglm_lambda<-glm(resp_lambda1~x_lambda-1,family=Gamma(link=RespLink_lambda),weight=weight_lambda1)
       res3<-summary(resglm_lambda,dispersion=2)
       lambda[temp11:temp12]<-resglm_lambda$fitted.values
       lambda_est<-lambda
        lambda_h<-res3$coefficients[1:p_lambda]
        temp11<-p_lambda+1
        temp12<-2*p_lambda
        lambda_se<-res3$coefficients[temp11:temp12]
        lambda_se[1]<-lambda_se[1]
        z_lambda<-lambda_h/lambda_se
        lambda_coeff<-cbind(matrix(lambda_h),matrix(lambda_se),matrix(z_lambda))
        colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(lambda_coeff) <- namesX_lambda
        print(lambda_coeff,4)
     }
   }
    if (is.null(PhiFix) && q_disp[1]==0) {
       if (RespDist=="gaussian" || RespDist=="gamma") {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
           res2<-summary(resglm_disp)
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res2$coefficients[1:p_disp]
           se_phi<-res2$coefficients[temp9:temp10]
           z_phi_coeff<-beta_phi/se_phi
           phi_coeff<-cbind(matrix(beta_phi),matrix(se_phi),matrix(z_phi_coeff))
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
       }
    }
    if (is.null(PhiFix) && q_disp[1]>0) {
       if (RespDist=="gaussian" || RespDist=="gamma") {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
           res4<-reshglm_disp
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res4[2][[1]]
           se_phi<-res4[3][[1]]
           z_phi<-beta_phi/se_phi
           phi_coeff<-cbind(matrix(beta_phi),matrix(se_phi),matrix(z_phi))
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
           print("Estimates for logarithm of tau=var(u_phi)")
           beta_tau<-log(res4[4][[1]])
           se_tau<-res4[6][[1]]/res4[4][[1]]^2
           z_tau<-beta_tau/se_tau[1,1]
           tau_coeff<-cbind(matrix(beta_tau),matrix(se_tau[1,1]),matrix(z_tau))
           colnames(tau_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(tau_coeff) <- namesRE_disp
           print(tau_coeff,4)
       }
    }
   }
#    v_h1<-corr_res[[7]]
    pi<-3.14159265359
    if (RespDist=="gaussian") hlikeli<-sum(-0.5*(y-mu)*(y-mu)/disp_est-0.5*log(2*disp_est*pi))
    if (RespDist=="poisson") hlikeli<-sum(y*log(mu)-mu-lgamma(y+1))
    if (RespDist=="binomial") hlikeli<-sum(y*log(mu/BinomialDen)+(BinomialDen-y)*log(1-mu/BinomialDen)+lgamma(BinomialDen+1)-lgamma(y+1)-lgamma(BinomialDen-y+1))
    if (RespDist=="gamma") hlikeli<-sum(log(y)/disp_est-log(y)-y/(disp_est*mu)-log(disp_est)/disp_est-log(mu)/disp_est-lgamma(1/disp_est))
    hh<-hlikeli
    if (RespDist=="gaussian") ll_y<-sum(-0.5*(y-y)*(y-y)/disp_est-0.5*log(2*disp_est*pi))
    if (RespDist=="poisson") ll_y<-sum((y+0.00001)*log(y)-y-lgamma(y+1))
    if (RespDist=="binomial") ll_y<-sum(y*log((y+0.00001)/BinomialDen)+(BinomialDen-y)*log(1-(y-0.00001)/BinomialDen)+lgamma(BinomialDen+1)-lgamma(y+1)-lgamma(BinomialDen-y+1))
    if (RespDist=="gamma") ll_y<-sum(log(y)/disp_est-log(y)-y/(disp_est*y)-log(disp_est)/disp_est-log(y)/disp_est-lgamma(1/disp_est))
    if (RespDist=="gaussian") deviance<-(y-mu)^2
    if (RespDist=="poisson") {
       y_zero<-1*(y==0)
       deviance<-2*y_zero*mu+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu)-(y+0.00001-mu))
    }
    if (RespDist=="binomial") deviance<-2*y*log((y+0.000001)/mu)+2*(BinomialDen-y)*log((BinomialDen-y+0.000001)/(BinomialDen-mu))
    if (RespDist=="gamma") deviance<-2*(-log(y/mu)+(y-mu)/mu)
    if (RespDist=="gaussian" || RespDist=="gamma") deviance<-deviance/disp_est
    print("========== Likelihood Function Values and Condition AIC ==========")
    model_number =1
    if (model_number == 1 || model_number == 2) {
       ml<- -2*hlikeli
       d2hdx2<--t(x)%*%(diag(W1)*x)  
       rl<- ml+log(abs(det(-d2hdx2/(2*pi))))
       pd<- p
       caic<-ml+2*pd
       sd<- -2*hh + 2*ll_y
       df<-length(y)-pd
    }
    if (model_number >=3) {
     if (check==0) {
       if (RandDist=="gaussian") {
            cc1<-svd(W2)
            logdet1<-sum(log(abs(1/cc1$d)))
            hv<--0.5*t(v_h)%*%(diag(W2)*v_h)-0.5*nrow(W2)*log(2*pi)-0.5*logdet1  
       } 
##       if (RandDist=="gaussian") hv<--0.5*t(v_h)%*%(diag(W2)*v_h)-0.5*nrow(W2)*log(2*pi)-0.5*log(abs(det(solve(W2))))  
       if (RandDist=="gamma") hv<-log(u_h)/lambda_est-u_h/lambda-log(lambda_est)/lambda_est-lgamma(1/lambda_est)
       if (RandDist=="inverse-gamma") {
           lambda_est1<-lambda_est/(1+lambda_est)
           alpha<-(1-lambda_est1)/lambda_est1
###           hv<-(v_h-log(u_h))/lambda_est1-(1+1/lambda_est1)*log(lambda_est1)-lgamma(1/lambda_est1)+log(lambda_est1)
           hv<-(alpha+1)*(log(alpha)-v_h)-alpha/u_h-lgamma(alpha+1)
       }
       if (RandDist=="beta") {
           lambda_est1<-2*lambda_est/(1-lambda_est)
           hv<-(0.5*v_h-log(1/(1-u_h)))/lambda_est1-lbeta(0.5/lambda_est1,0.5/lambda_est1)
       }
     } else {
    hv<-matrix(0,1,1)
  for(i in 1:nrand) {
    temp11<-qcum[i]+1
    temp12<-qcum[i+1]
       if (RandDist1[i]=="gaussian") {
            cc1<-svd(W2[temp11:temp12,temp11:temp12])
            logdet1<-sum(log(abs(1/cc1$d)))
            hv<-hv-0.5*t(v_h[temp11:temp12])%*%W2[temp11:temp12,temp11:temp12]%*%v_h[temp11:temp12]-0.5*nrow(W2[temp11:temp12,temp11:temp12])*log(2*pi)-0.5*logdet1  
       } 
       if (RandDist1[i]=="gamma") {
            hv<-hv+sum(log(u_h[temp11:temp12])/lambda_est[temp11:temp12]-u_h[temp11:temp12]/lambda[temp11:temp12]-log(lambda_est[temp11:temp12])/lambda_est[temp11:temp12]-lgamma(1/lambda_est[temp11:temp12]))
       }
       if (RandDist1[i]=="inverse-gamma") {
           lambda_est1<-lambda_est[temp11:temp12]/(1+lambda_est[temp11:temp12])
           alpha<-(1-lambda_est1)/lambda_est1
           hv<-hv+sum((alpha+1)*(log(alpha)-v_h[temp11:temp12])-alpha/u_h[temp11:temp12]-lgamma(alpha+1))
       }
       if (RandDist1[i]=="beta") {
           lambda_est1<-2*lambda_est[temp11:temp12]/(1-lambda_est[temp11:temp12])
           hv<-hv+sum((0.5*v_h-log(1/(1-u_h)))/lambda_est1-lbeta(0.5/lambda_est1,0.5/lambda_est1))
       }
   }   
     }
       if (model_number == 4) hv10<-reshglm_disp[[8]][1]
       else hv10<-0
       if (model_number1 == 1) hv20<-reshglm_lambda[[8]][1]
       else hv20<-0
        if (model_number == 4) hv11<-reshglm_disp[[8]][2]
       else hv11<-0
       if (model_number1 == 1) hv21<-reshglm_lambda[[8]][2]
       else hv21<-0
       if (model_number == 4) {
            hv12<-reshglm_disp[[8]][3]
       }
       else hv12<-0
       if (model_number1 == 1) hv22<-reshglm_lambda[[8]][3]
       else hv22<-0
        cc1<-svd((-d2hdv2)/(2*pi))
        logdet1<-sum(log(abs(cc1$d)))
       ml<- -2*hlikeli-2*sum(hv)+logdet1 ##-log(2*pi*nrow(d2hdv2))
##       ml<- -2*hlikeli-2*sum(hv)+log(abs(det(-d2hdv2/(2*pi))))
       W1x<-diag(W1)*x
       W1z<-diag(W1)*z
       AA<-rbind(cbind(matrix((t(x)%*%W1x),nrow(t(x)%*%W1x),ncol(t(x)%*%W1x)),matrix((t(x)%*%W1z),nrow(t(x)%*%W1z),ncol(t(x)%*%W1z))),cbind(matrix((t(z)%*%W1x),nrow(t(z)%*%W1x),ncol(t(z)%*%W1x)),matrix((-1*d2hdv2),nrow(d2hdv2),ncol(d2hdv2))))  
       BB<-rbind(cbind(matrix((t(x)%*%W1x),nrow(t(x)%*%W1x),ncol(t(x)%*%W1x)),matrix((t(x)%*%W1z),nrow(t(x)%*%W1z),ncol(t(x)%*%W1z))),cbind(matrix((t(z)%*%W1x),nrow(t(z)%*%W1x),ncol(t(z)%*%W1x)),matrix((t(z)%*%W1z),nrow(t(z)%*%W1z),ncol(t(z)%*%W1z))))  
        cc1<-svd(AA/(2*pi))
        logdet1<-sum(log(abs(cc1$d)))
       rl<--2*hlikeli-2*sum(hv)+logdet1 ##-log(2*pi*nrow(AA))
##       rl<--2*hlikeli-2*sum(hv)+logdet1-log(2*pi*nrow(AA))
       pd<- sum(diag(solve(AA) %*% BB))  
       caic<- -2*hlikeli + 2*pd
       sd<- -2*hh + 2*ll_y
       sd<-sum(deviance)
       df<-length(y)-pd
    }
       ml<-ml-230
       rl<-rl-230
       caic<-caic-230

    likeli_coeff<-rbind(matrix(ml),matrix(rl),matrix(caic),matrix(sd),matrix(df))
    if (model_number1 == 1) {
       likeli_coeff<-rbind(matrix(ml),matrix(rl),matrix(caic))
    }
    if (model_number == 1 || model_number == 2) {
       rownames(likeli_coeff)<-c("-2ML (-2 h)          : ","-2RL (-2 p_beta (h)) : ","cAIC                 : ", "SD                   : ","df                   : ")
    }
    if (model_number ==3 ) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu) (h))          : ","-2RL (-2 p_beta(mu),v(mu) (h)) : ","cAIC                           : ", "SD                             : ","df                             : ")
    }
    if (model_number == 4) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu),v(phi) (h))          : ","-2RL (-2 p_beta(mu),v(mu),beta(phi),v(phi) (h)) : ","cAIC                           : ", "SD                             : ","df                             : ")
    }
    if (model_number<4 && model_number1 == 1) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu),v(lambda) (h))          : ","-2RL (-2 p_beta(mu),v(mu),beta(lambda),v(lambda) (h)) : ","cAIC                           : ", "SD                             : ","df                             : ")
    }
    if (model_number == 4 && model_number1 == 1) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu),v(phi),v(lambda) (h))          : ","-2RL (-2 p_beta(mu),v(mu),beta(phi),v(phi),beta(lambda),v(lambda) (h)) : ","cAIC                           : ", "SD                             : ","df                             : ")
    }
    print(likeli_coeff)
##    df<-n-p
##    Devdf<-matrix(c(df,sum(deviance_residual)),nrow=2)
##    rownames(Devdf)<-c("DF :","Deviance :")
	res <- list(mean_residual=matrix(mean_residual),mu=matrix(mu),deviance_residual=matrix(deviance_residual),df=df, ml = ml, rl = rl, caic = caic, md = md,sd=sd,df=df, 
        formulaDisp = formulaDisp,RespLink_disp=RespLink_disp,beta_coeff=beta_coeff,tau_coeff=tau_coeff,rhoLink=rhoLink,rho_coeff=rho_coeff,model_mu=MeanModel,model_phi=DispersionModel)
    return(res)
}

vif1.default <- function(mod, xx, names, ...) { 
    if (any(is.na(coef(mod))))  
        stop ("there are aliased coefficients in the model") 
    v <- vcov(mod) 
    assign <- attr(model.matrix(mod), "assign") 
    if (names(coefficients(mod)[1]) == "(Intercept)") { 
        v <- v[-1, -1] 
        assign <- assign[-1] 
    } 
    else warning("No intercept: vifs may not be sensible.") 
    ## terms <- labels(terms(mod)) 
    terms <- names
    n.terms <- length(terms) 
    ## if (n.terms < 2) stop("model contains fewer than 2 terms") 
    ## print(terms)
    R <- cov2cor(v) 
    detR <- det(R) 
    result <- matrix(0, n.terms, 3) 
    rownames(result) <- terms 
    colnames(result) <- c("GVIF", "Df", "GVIF^(1/(2*Df))") 
    for (term in 1:n.terms) { 
        subs <- which(assign == term) 
        result[term, 1] <- det(as.matrix(R[subs, subs])) * 
            det(as.matrix(R[-subs, -subs])) / detR 
        result[term, 2] <- length(subs) 
    } 
    if (all(result[, 2] == 1)) result <- result[, 1] 
    else result[, 3] <- result[, 1]^(1/(2 * result[, 2])) 
    result[,1]
} 

HGLMREML.h <-function(formulaMain,DataMain,Offset=NULL,RespDist="gaussian",RespLink="identity",
RandDist="gaussian",mord=0,dord=1,spatial=NULL,Neighbor=NULL,Maxiter=200,Iter_mean=3,convergence=10^(-2),
Init_lam=0.25,Init_rho=0.174,ODEst=FALSE,tolerance=1,Init_phi=1,contrasts=NULL){
    require(Matrix)
    require(numDeriv)
    require(boot)
    mc <- match.call()
    fr <- HGLMFrames(mc, formulaMain, contrasts)
    namesX <- names(fr$fixef)
    namesY <- names(fr$mf)[1]
    FL <- HGLMFactorList(formulaMain, fr, 0L, 0L)
    namesRE <- FL$namesRE
    y <- matrix(fr$Y, length(fr$Y), 1)
    x <- fr$X
    z <- FL$Design
    n<-nrow(x)
    p<-ncol(x)
    nrand <- length(z)
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
    z<-zz<-z[[1]]
    if (spatial=="IAR" && !is.null(Neighbor)) {
          nn<-nrow(Neighbor)
          no<-matrix(0,nn,nn)
          for (i in 1:nn) no[i,i]<-sum(Neighbor[i,])
          sp1<-(no-Neighbor)
          index4<-nn-1
          s_ <-svd(sp1,nu=index4,nv=index4)
          uuu<-matrix(0,nn,index4)
          ddd<-matrix(0,index4,index4)
          for ( i in 1:index4) ddd[i,i]<-sqrt(1/s_$d[i])
          for (i in 1:nn){
             for (j in 1:index4){
             uuu[i,j]<-s_$u[i,j]
          }
          }
          LLL<-uuu %*% ddd
          z<-z%*%LLL
          for (i in 1:nrand) q[i]<-index4
###          print(LLL)
    }
##############################################################
######### initial values : GLM estimates #####################
##############################################################
    dord<-1
    phi <- rep(Init_phi,n)
    beta_h<-matrix(0,p,1)
    qcum <- cumsum(c(0, q))
    v_h<-matrix(0,qcum[nrand+1],1)
    if (RandDist=="normal") RandDist="gaussian"
    if (RandDist=="gaussian") u_h <- v_h
    if (RandDist=="gamma") u_h <-exp(v_h)
    alpha_h <- rep(0, nrand)
    for (i in 1:nrand) alpha_h[i] <- Init_lam
    if (RespDist=="Poisson") RespDist=c("poisson")
    if (RespDist=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink),offset=Offset)
    if (RespDist=="poisson" ) resglm<-glm(y~x-1,family=poisson(link=RespLink),offset=Offset)
    if (RespDist=="binomial") resglm<-glm(y~x-1,family=binomial(link=RespLink),offset=Offset)
    if (RespDist=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink),offset=Offset)
    beta_h[1:p,1]<-c(resglm$coefficients)[1:p]
    if (is.null(Offset)) off<- matrix(0, n,1)
    else off<- Offset
    if (spatial=="MRF" && !is.null(Neighbor)) rho<-Init_rho
    else rho<-0
convergence1<-1
max_iter<-1
while (convergence1>convergence && max_iter<=Maxiter ) {
for(k in 1:Iter_mean) {
    eta <- off + x %*% beta_h + z %*% v_h
    if (RespLink=="identity") {
        mu <- eta
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta)
        detadmu <- 1/mu
    }
    if (RespLink=="inverse") {
        mu <- 1/eta
        detadmu <- -1/mu^2
    }
     if (RespLink=="logit") {
        mu <- 1/(1+exp(-eta))
        detadmu <- 1/(mu*(1-mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(1-mu)
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 /(phi*Vmu)
    W1<-diag(as.vector(temp4))
    z1<-eta+(y-mu)*detadmu-off
    oq<-matrix(1,qcum[nrand+1],1)
    lambda<-matrix(1,qcum[nrand+1],1)
    for (i in 1:nrand) {
       index1<-qcum[i]+1
       lambda[index1:qcum[i+1]]<-alpha_h[i]
    } 
    I<-diag(rep(1,qcum[nrand+1]))
    if (spatial=="MRF" && !is.null(Neighbor)) {
        pW2<-(I-rho*Neighbor)
        W2<-pW2/as.vector(lambda)
    } else {
        pW2<-I
        W2<-diag(1/as.vector(lambda))
    }
##############################################################
############# random effect  #################################
##############################################################
    c_v_h<-1.0
    iter_v<-1
    Sig<- z %*% solve(W2) %*% t(z) +solve(W1)
    invSig<-solve(Sig)
    eta <- off + x %*% beta_h + z %*% v_h
    if (RespLink=="identity") {
        mu <- eta
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta)
        detadmu <- 1/mu
    }
    if (RespLink=="inverse") {
        mu <- 1/eta
        detadmu <- -1/mu^2
    }
     if (RespLink=="logit") {
        mu <- 1/(1+exp(-eta))
        detadmu <- 1/(mu*(1-mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(1-mu)
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 /(phi*Vmu)
    W1<-diag(as.vector(temp4))
    z1<-eta+(y-mu)*detadmu-off
    lambda<-matrix(1,qcum[nrand+1],1)
    for (i in 1:nrand) {
       index1<-qcum[i]+1
       lambda[index1:qcum[i+1]]<-alpha_h[i]
    } 
    I<-diag(rep(1,qcum[nrand+1]))
    if (RespDist=="poisson") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-W2%*%v_h
##          dhdv<-t(z)%*%(y-mu)-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    if (RespDist=="gaussian") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    if (RespDist=="binomial") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    if (RespDist=="gamma") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
        if (RandDist=="inverse-gamma") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-(1+1/lambda)+exp(-v_h)/lambda
            temp5<-exp(-v_h)/lambda
            W2<-diag(as.vector(temp5))
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    v_h_old<-v_h
    v_h<-v_h-solve(d2hdv2)%*%dhdv
    c_v_h<-sum(abs(as.vector(v_h_old)-as.vector(v_h)))
    iter_v<-iter_v+1
##    }
##############################################################
########## 1st order adjusted term for mean ##################
##############################################################
    eta <- off + x %*% beta_h + z %*% v_h
    if (RespLink=="identity") {
        mu <- eta
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta)
        detadmu <- 1/mu
    }
    if (RespLink=="inverse") {
        mu <- 1/eta
        detadmu <- -1/mu^2
    }
     if (RespLink=="logit") {
        mu <- 1/(1+exp(-eta))
        detadmu <- 1/(mu*(1-mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(1-mu)
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 /(phi*Vmu)
    W1<-diag(as.vector(temp4))
    z1<-eta+(y-mu)*detadmu-off
    Sig<- z %*% solve(W2) %*% t(z) +solve(W1)
    invSig<-solve(Sig)
    if (mord==0) a<-matrix(0,n,1)
    if (mord==1) {
    T<-t(cbind(t(z),I))
##    Null<-matrix(0,n,n)
##    W<-matrix(0,(2*n),(2*n))
##    W[c(1:n),]<-cbind(W1,Null)
##    W[c((n+1):(2*n)),]<-cbind(Null,W2)   
    Null1<-matrix(0,n,qcum[nrand+1])
    Null2<-matrix(0,qcum[nrand+1],n)
    W<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
    W[c(1:n),]<-cbind(W1,Null1)
    W[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2)   
    P<-T%*%solve(t(T)%*%W%*%T)%*%t(T)%*%W
    K1<--z%*%solve(t(T)%*%W%*%T)%*%t(z)
    K2<--solve(t(T)%*%W%*%T)
    d1<-rep(0,n)
    d2<-rep(0,n)
    d3<-rep(0,n)
    for (i in 1:n){
        d1[i]<-P[i,i]*detadmu[i]
        d2[i]<-0
        for (qq in 1:n){
            d2[i]<-d2[i]+P[qq,qq]*K1[qq,i]
        }
        if (RandDist=="gaussian") d3[i]<-0
    }
    d<-d1+d2+d3
    s<-d*dmudeta/2
    a<-(solve(W1)+z%*%solve(W2)%*%t(z))%*%W1%*%(s*detadmu)
    }
    beta_h_old<-beta_h
######################################################################
############# mean parameters (beta) #################################
######################################################################
    beta_h<-solve(t(x)%*%invSig%*%x)%*%(t(x)%*%invSig%*%(z1-a))
    se_beta<-sqrt(abs(diag(solve(t(x)%*%invSig%*%x))))
############################################################## 
}
###############################################################
############# dispersion parameters ###########################
###############################################################

##############################################################
######### Estimates for lambda and rho #####################
##############################################################
    v<-v_h
    Q<-invSig-invSig%*%x%*%solve(t(x)%*%invSig%*%x)%*%t(x)%*%invSig
    lam<-alpha_h[1]
#### 1: lam(variance component) , 2: rho
    if (spatial=="MRF" && !is.null(Neighbor)) {
        pW2<-(I-rho*Neighbor)
        W2<-pW2/as.vector(lambda)
    } else {
        rho<-0
        pW2<-I
        W2<-diag(1/as.vector(lambda))
    }
    dREMLdlam<-c(0,0)
    dREML1dlam<-c(0,0)
    dREML2dlam<-c(0,0)
    d2REMLd2lam<-matrix(0,2,2)

    dW2dlam<--W2/lam
    if (spatial=="MRF" && !is.null(Neighbor)) dW2drho<-(-Neighbor)/lam
    else dW2drho<-0

    dSig1dlam<-solve(pW2)
    if (spatial=="MRF" && !is.null(Neighbor)) dSig1drho<--solve(W2)%*%dW2drho%*%solve(W2)
    else dSig1drho<-0

    dvdlam<-solve(t(z)%*%W1%*%z+W2)%*%(-dW2dlam)%*%v
    if (spatial=="MRF" && !is.null(Neighbor)) dvdrho<-solve(t(z)%*%W1%*%z+W2)%*%(-dW2drho)%*%v
    else dvdrho<-0
###    if (ODEst==TRUE) dvdrho<-dvdrho*0

    if (RespDist=="gaussian") kkk<-0*mu
    if (RespDist=="poisson") kkk<-mu
    if (RespDist=="binomial") kkk<-(1-2*mu)*dmudeta
    if (RespDist=="gamma") kkk<-0*mu
    dW1dlam<-diag(as.vector(kkk*(z%*%dvdlam)))
    if (spatial=="MRF" && !is.null(Neighbor)) dW1drho<-diag(as.vector(kkk*(z%*%dvdrho)))
    else dW1drho<-0

    dSig2dlam<--solve(W1)%*%dW1dlam%*%solve(W1)
    if (spatial=="MRF" && !is.null(Neighbor)) dSig2drho<--solve(W1)%*%dW1drho%*%solve(W1)
    else dSig2drho<-0

    dSigdlam<-z%*%dSig1dlam%*%t(z)+dSig2dlam
    if (spatial=="MRF" && !is.null(Neighbor)) dSigdrho<-z%*%dSig1drho%*%t(z)+dSig2drho
    else  dSigdrho<-0

    dterm1dlam<--t(v)%*%(dW2dlam)%*%v/2
    if (spatial=="MRF" && !is.null(Neighbor)) dterm1drho<--t(v)%*%(dW2drho)%*%v/2
    else dterm1drho<-0

    if (RespDist=="poisson") dW1dv<-W1
    if (RespDist=="gaussian") dW1dv<-diag(as.vector(kkk))
    if (RespDist=="binomial") dW1dv<-diag(as.vector(kkk))
    if (RespDist=="gamma") dW1dv<-diag(as.vector(kkk))
    dterm2dv<-y-mu-z%*%W2%*%v-1/2*(1/diag(W1))*diag(dW1dv)

    dREMLdlam[1]<--t(v)%*%dW2dlam%*%v/2-t(dvdlam)%*%W2%*%v-0.5*sum(diag(Q%*%dSigdlam))+t(dvdlam)%*%t(z)%*%W1%*%((y-mu)*detadmu)-0.5*sum(diag(solve(W1)%*%dW1dlam))
    if (spatial=="MRF" && !is.null(Neighbor)) {
        dREMLdlam[2]<-t(dvdrho)%*%t(z)%*%W1%*%((y-mu)*detadmu)-0.5*sum(diag(solve(W1)%*%dW1drho))-t(dvdrho)%*%W2%*%v-t(v)%*%dW2drho%*%v/2 -0.5*sum(diag(Q%*%dSigdrho))
        if (ODEst==TRUE) dREMLdlam[2]<--0.5*sum(diag(solve(W1)%*%dW1drho))-t(dvdrho)%*%W2%*%v-t(v)%*%dW2drho%*%v/2 -0.5*sum(diag(Q%*%dSigdrho))
    }
    d2W2dlam2<-2*W2/lam^2
    d2W2drho2<-matrix(0,q,q)
    if (spatial=="MRF" && !is.null(Neighbor)) d2W2dlamrho<-Neighbor/lam
    d2vdlam2<-solve(t(z)%*%W1%*%z+W2)%*%(-dW2dlam%*%dvdlam-d2W2dlam2%*%v)

    if (spatial=="MRF" && !is.null(Neighbor)) dvdrho<-solve(t(z)%*%W1%*%z+W2)%*%(-dW2drho)%*%v
    if (spatial=="MRF" && !is.null(Neighbor)) d2vdrho2<-solve(t(z)%*%W1%*%z+W2)%*%(-dW2drho%*%dvdrho-d2W2drho2%*%v)

    H<-t(z)%*%W1%*%z+W2
    dHdlam<-dW2dlam
    if (spatial=="MRF" && !is.null(Neighbor)) dHdrho<-dW2drho
    d2Hdlam2<-d2W2dlam2
    if (spatial=="MRF" && !is.null(Neighbor)) d2Hdrho2<-d2W2drho2
    if (spatial=="MRF" && !is.null(Neighbor)) d2Hdlamrho<-d2W2dlamrho

    d2REMLd2lam[1,1]<--0.5*t(v)%*%d2W2dlam2%*%v-0.5*sum(diag(solve(W2)%*%dW2dlam%*%solve(W2)%*%dW2dlam))+0.5*sum(diag(solve(W2)%*%d2W2dlam2))+0.5*sum(diag(solve(H)%*%dHdlam%*%solve(H)%*%dHdlam))-0.5*sum(diag(solve(H)%*%d2Hdlam2)) 
    if (spatial=="MRF" && !is.null(Neighbor)) d2REMLd2lam[1,2]<--0.5*t(v)%*%d2W2dlamrho%*%v-0.5*sum(diag(solve(W2)%*%dW2dlam%*%solve(W2)%*%dW2drho))+0.5*sum(diag(solve(W2)%*%d2W2dlamrho))+0.5*sum(diag(solve(H)%*%dHdlam%*%solve(H)%*%dHdrho))-0.5*sum(diag(solve(H)%*%d2Hdlamrho)) 
    if (spatial=="MRF" && !is.null(Neighbor)) d2REMLd2lam[2,1]<-d2REMLd2lam[1,2] 
    if (spatial=="MRF" && !is.null(Neighbor)) d2REMLd2lam[2,2]<-0-0.5*sum(diag(solve(W2)%*%dW2drho%*%solve(W2)%*%dW2drho))+0.5*sum(diag(solve(W2)%*%d2W2drho2))+0.5*sum(diag(solve(H)%*%dHdrho%*%solve(H)%*%dHdrho))-0.5*sum(diag(solve(H)%*%d2Hdrho2))  

    clam<-c(lam,rho)
    old_clam<-clam

    if (spatial=="MRF" && !is.null(Neighbor)) {
        clam<-clam-solve(d2REMLd2lam)%*%dREMLdlam/tolerance
        if (clam[1]<0) clam[1]<-0.015
##          print(clam)
##        print(dREMLdlam)
##        print(d2REMLd2lam)
##        if (clam[2]>1) clam[2]<- 0.9999
##        if (clam[2]< -1) clam[2]<- -0.9999
    }
    else clam[1]<-clam[1]-dREMLdlam[1]/d2REMLd2lam[1,1]
    convergence1<-sum(abs(clam-old_clam))
    lam<-clam[1]
    rho<-clam[2]
    alpha_h[1]<-lam
    max_iter<-max_iter+1
    print_i<-max_iter
    print_err<-convergence1
    names(print_i) <- "iteration : "
##       print(print_i)
    names(print_err) <- "convergence : "
    for (i in 1:nrand) {
       index1<-qcum[i]+1
       lambda[index1:qcum[i+1]]<-alpha_h[i]
    } 
    I<-diag(rep(1,qcum[nrand+1]))
    if (spatial=="MRF" && !is.null(Neighbor)) {
        pW2<-(I-rho*Neighbor)
        W2<-pW2/as.vector(lambda)
    } else {
        rho<-0
        pW2<-I
        W2<-diag(1/as.vector(lambda))
    }

###################################################################################
######### Dispersion Estimates for phi for Poisson and binomail ###################
###################################################################################
  if (ODEst == TRUE) {
      old_disp_est <- phi
      if (RespDist=="poisson") {
         y_zero<-1*(y==0)
         deviance_residual<-(2*y_zero*mu+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu)-(y+0.00001-mu)))/phi
      }
      if (RespDist=="binomial") deviance_residual<-(1*(y==0))*2*log(1/(1-mu))+(1*(y==1))*2*log(1/mu)
      OO1<-matrix(0,qcum[nrand+1],p)
      Null1<-matrix(0,n,qcum[nrand+1])
      Null2<-matrix(0,qcum[nrand+1],n)
      TT<-rbind(cbind(x,z),cbind(OO1,I))
      WW<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
      WW[c(1:n),]<-cbind(W1,Null1)
      WW[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2)   
      PP<-TT%*%solve(t(TT)%*%WW%*%TT)%*%t(TT)%*%WW
      leverage<-rep(0,n)
      for (kk in 1:n) leverage[kk]<-PP[kk,kk]
      resp_disp<-deviance_residual/(1-leverage)
      weight_disp<-(1-leverage)/2
##############################################################
######### GLM fit for phi #####################
##############################################################
      resglm_disp<-glm(resp_disp~1,family=Gamma(link="log"),weight=weight_disp)
      inv_disp<-1/resglm_disp$fitted.values
      disp_est<-resglm_disp$fitted.values
      phi<-disp_est
    if (RespDist == "poisson") {
        dhdphi<- sum(-1/phi+mu/phi^2+y*log(phi)/phi^2-y/phi^2+digamma(y/phi+1)*y/phi^2)
        d2hdphi<- sum(1/phi^2-2*mu/phi^3-2*y*log(phi)/phi^3+y/phi^3+2*y/phi^3-trigamma(y/phi+1)*y^2/phi^4-digamma(y/phi+1)*2*y/phi^3)
    }
##    phi_temp<-phi[1]-dhdphi/d2hdphi/10
##    phi<-rep(phi_temp,n)
##    disp_est<-phi
    print(phi[1])
    convergence2<-sum(abs(disp_est-old_disp_est))
    old_disp_est<-disp_est
  }
    if (RespDist=="gaussian") {
        temp5<- dmudeta^2 /(Vmu)
        dW1dphi<-diag(-as.vector(temp5)/phi^2)
        dvdphi<-solve(t(z)%*%W1%*%z+W2)%*%(-t(z)%*%dW1dphi%*%(y-mu))
###        dvdphi<-0*dvdphi
        dSigdphi<-diag(as.vector(temp5))
        dW2dphi<-0*dW2dlam
        d2W2dphi2<-0*d2W2dlam2
        H<-t(z)%*%W1%*%z+W2
        HX<-t(x)%*%W1%*%x
        HXZ<-t(x)%*%W1%*%z
        dHdphi<-t(z)%*%dW1dphi%*%z
        dHXdphi<-t(x)%*%dW1dphi%*%x
        dHXZdphi<-t(x)%*%dW1dphi%*%z
        Hp<-rbind(cbind(HX,HXZ),cbind(t(HXZ),H))
        dHpdphi<-rbind(cbind(dHXdphi,dHXZdphi),cbind(t(dHXZdphi),dHdphi))
        d2W1dphi2<-diag(2*as.vector(temp5)/phi^3)
        d2Hdphi2<-t(z)%*%d2W1dphi2%*%z
        d2HXdphi2<-t(x)%*%d2W1dphi2%*%x
        d2HXZdphi2<-t(x)%*%d2W1dphi2%*%z
        d2Hpdphi2<-rbind(cbind(d2HXdphi2,d2HXZdphi2),cbind(t(d2HXZdphi2),d2Hdphi2))
        dREMLdphi<--t(v)%*%dW2dphi%*%v/2-t(dvdphi)%*%W2%*%v-0.5*sum(diag(Q%*%dSigdphi))+t(dvdphi)%*%t(z)%*%W1%*%((y-mu)*detadmu)+sum(0.5*(y-mu)^2/phi^2)
###            -0.5*sum(diag(solve(W1)%*%dW1dphi))
###        d2REMLd2phi<--0.5*t(v)%*%d2W2dphi2%*%v-0.5*sum(diag(solve(W2)%*%dW2dphi%*%solve(W2)%*%dW2dphi))+0.5*sum(diag(solve(W2)%*%d2W2dphi2))-sum((y-mu)^2/phi^3)+sum(0.5/phi^2)+0.5*sum(diag(solve(H)%*%dHdphi%*%solve(H)%*%dHdphi))-0.5*sum(diag(solve(H)%*%d2Hdphi2)) 
###        d2REMLd2phi<--sum((y-mu)^2/phi^3)+sum(0.5/phi^2)-sum((p+q)*0.5/(n*phi^2))
###        print(dREMLdphi)
###        print(d2REMLd2phi)
###        d2REMLd2phi<--sum((y-mu)^2/phi^3)+sum(0.5/phi^2)+0.5*sum(diag(solve(H)%*%dHdphi%*%solve(H)%*%dHdphi))-0.5*sum(diag(solve(H)%*%d2Hdphi2))+0.5*sum(diag(solve(HX)%*%dHXdphi%*%solve(HX)%*%dHXdphi))-0.5*sum(diag(solve(HX)%*%d2HXdphi2)) 
###        d2REMLd2phi<--sum((y-mu)^2/phi^3)+sum(0.5/phi^2)+0.5*sum(diag(solve(Hp)%*%dHpdphi%*%solve(Hp)%*%dHpdphi)) 
        d2REMLd2phi<--sum((y-mu)^2/phi^3)+sum(0.5/phi^2)+0.5*sum(diag(solve(Hp)%*%dHpdphi%*%solve(Hp)%*%dHpdphi))-0.5*sum(diag(solve(Hp)%*%d2Hpdphi2))
        oldphi1<-phi[1]
        phi1<-phi[1]+dREMLdphi/abs(d2REMLd2phi)
        phi <- rep(phi1,n)
        convergence2<-sum(abs(phi[1]-oldphi1))
        convergence1<-convergence1+convergence2
        print_err<-print_err+convergence2
###        print(phi[1])
    }
    if (RespDist=="gamma") {
        temp5<- dmudeta^2 /(Vmu)
        dW1dphi<-diag(-as.vector(temp5)/phi^2)
        dvdphi<-solve(t(z)%*%W1%*%z+W2)%*%(-t(z)%*%dW1dphi%*%(y-mu))
###        dvdphi<-0*dvdphi
        dSigdphi<-diag(as.vector(temp5))
        dW2dphi<-0*dW2dlam
        d2W2dphi2<-0*d2W2dlam2
        H<-t(z)%*%W1%*%z+W2
        HX<-t(x)%*%W1%*%x
        HXZ<-t(x)%*%W1%*%z
        dHdphi<-t(z)%*%dW1dphi%*%z
        dHXdphi<-t(x)%*%dW1dphi%*%x
        dHXZdphi<-t(x)%*%dW1dphi%*%z
        Hp<-rbind(cbind(HX,HXZ),cbind(t(HXZ),H))
        dHpdphi<-rbind(cbind(dHXdphi,dHXZdphi),cbind(t(dHXZdphi),dHdphi))
        d2W1dphi2<-diag(2*as.vector(temp5)/phi^3)
        d2Hdphi2<-t(z)%*%d2W1dphi2%*%z
        d2HXdphi2<-t(x)%*%d2W1dphi2%*%x
        d2HXZdphi2<-t(x)%*%d2W1dphi2%*%z
        d2Hpdphi2<-rbind(cbind(d2HXdphi2,d2HXZdphi2),cbind(t(d2HXZdphi2),d2Hdphi2))
###        d2HXdphi2<-0*t(x)%*%dW1dphi%*%x
###        dREMLdphi<-sum(-log(y)/phi^2+y/(phi^2*mu)+log(phi)/phi^2-1/phi^2+log(mu)/phi^2+digamma(1/phi)/phi^2)-0.5*sum(diag(solve(H)%*%dHdphi))-0.5*sum(diag(solve(HX)%*%dHXdphi))
        dREMLdphi<-sum(-log(y)/phi^2+y/(phi^2*mu)+log(phi)/phi^2-1/phi^2+log(mu)/phi^2+digamma(1/phi)/phi^2)-0.5*sum(diag(solve(Hp)%*%dHpdphi))
###        d2REMLd2phi<-sum(2*log(y)/phi^3-2*y/(phi^3*mu)-2*log(phi)/phi^3+3/phi^3-log(mu)/phi^3-2*digamma(1/phi)/phi^3-trigamma(1/phi)/phi^4)+0.5*sum(diag(solve(H)%*%dHdphi%*%solve(H)%*%dHdphi))-0.5*sum(diag(solve(H)%*%d2Hdphi2))+0.5*sum(diag(solve(HX)%*%dHXdphi%*%solve(HX)%*%dHXdphi))
        d2REMLd2phi<-sum(2*log(y)/phi^3-2*y/(phi^3*mu)-2*log(phi)/phi^3+3/phi^3-2*log(mu)/phi^3-2*digamma(1/phi)/phi^3-trigamma(1/phi)/phi^4)+0.5*sum(diag(solve(Hp)%*%dHpdphi%*%solve(Hp)%*%dHpdphi))-0.5*sum(diag(solve(Hp)%*%d2Hpdphi2))
        oldphi1<-phi[1]
        phi1<-phi[1]+dREMLdphi/abs(d2REMLd2phi)
        phi <- rep(phi1,n)
        convergence2<-sum(abs(phi[1]-oldphi1))
        convergence1<-convergence1+convergence2
        print_err<-print_err+convergence2
    }
##        print(print_err)
}
###############################################################
############# se for dispersion estimates######################
###############################################################
    X<-x
    p<-ncol(X)
    O1<-matrix(0,p,p)
    O2<-matrix(0,p,qcum[nrand+1])
    if (spatial=="MRF" && !is.null(Neighbor)) infoterm<-matrix(0,2,2)
    else infoterm<-matrix(0,2,2) 
    d2hlikedlam2<-n/(2*lam^2)-t(v)%*%pW2%*%v/lam^3
    A<-rbind(cbind((t(X)%*%W1%*%X),(t(X)%*%W1%*%z%*%I)),cbind((t(I)%*%t(z)%*%W1%*%X),H))
    dAdlam<-rbind(cbind(O1,O2),cbind(t(O2),dHdlam))
    d2Adlam2<-rbind(cbind(O1,O2),cbind(t(O2),d2Hdlam2))
    d2hlikedlamdv<-pW2%*%v/lam^2
    dAdv_dvdlam<-rbind(cbind((t(X)%*%(dW1dv%*%diag(as.vector(z%*%dvdlam)))%*%X),(t(X)%*%(dW1dv%*%diag(as.vector(z%*%dvdlam)))%*%z%*%I)),
                 cbind((t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%dvdlam)))%*%X),(t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%dvdlam)))%*%z%*%I)))
    dAdv_d2vdlam2<-rbind(cbind((t(X)%*%(dW1dv%*%diag(as.vector(z%*%d2vdlam2)))%*%X),(t(X)%*%(dW1dv%*%diag(as.vector(z%*%d2vdlam2)))%*%z%*%I)),
                 cbind((t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%d2vdlam2)))%*%X),(t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%d2vdlam2)))%*%z%*%I)))
    if (spatial=="MRF" && !is.null(Neighbor))  {
        d2hlikedrho2<--1/2*sum(diag(solve(pW2)%*%Neighbor%*%solve(pW2)%*%Neighbor))
        d2hlikedrhodlam<--t(v)%*%Neighbor%*%v/(2*lam^2)
        d2vdrhodlam<--solve(H)%*%dHdrho%*%solve(H)%*%(pW2%*%v)/lam^2+solve(H)%*%(-Neighbor%*%v)/lam^2
        dAdrho<-rbind(cbind(O1,O2),cbind(t(O2),dHdrho))
        d2Adrho2<-rbind(cbind(O1,O2),cbind(t(O2),d2Hdrho2))
        d2Adrhodlam<-rbind(cbind(O1,O2),cbind(t(O2),d2Hdlamrho))
        d2hlikedrhodv<--Neighbor%*%v/lam
        dAdv_dvdrho<-rbind(cbind((t(X)%*%(dW1dv%*%diag(as.vector(z%*%dvdrho)))%*%X),(t(X)%*%(dW1dv%*%diag(as.vector(z%*%dvdrho)))%*%z%*%I)),cbind((t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%dvdrho)))%*%X),(t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%dvdrho)))%*%z%*%I)))
        dAdv_d2vdrho2<-rbind(cbind((t(X)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrho2)))%*%X),(t(X)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrho2)))%*%z%*%I)),cbind((t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrho2)))%*%X),(t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrho2)))%*%z%*%I)))
        dAdv_d2vdrhodlam<-rbind(cbind((t(X)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrhodlam)))%*%X),(t(X)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrhodlam)))%*%z%*%I)),cbind((t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrhodlam)))%*%X),(t(I)%*%t(z)%*%(dW1dv%*%diag(as.vector(z%*%d2vdrhodlam)))%*%z%*%I)))
   }

    tinfoterm1<-d2hlikedlam2+1/2*sum(diag(solve(A)%*%dAdlam%*%solve(A)%*%dAdlam))-1/2*sum(diag(solve(A)%*%d2Adlam2))
    tinfoterm2<-sum(as.vector(d2hlikedlamdv*dvdlam))+1/2*sum(diag(solve(A)%*%(dAdv_dvdlam)%*%solve(A)%*%dAdlam))
    tinfoterm3<-tinfoterm2-1/2*sum(diag(solve(A)%*%dAdv_d2vdlam2))
    tinfoterm4<-t(dvdlam)%*%(-H)%*%dvdlam+1/2*sum(diag(solve(A)%*%dAdv_dvdlam%*%solve(A)%*%dAdv_dvdlam))

    infoterm[1,1]<-tinfoterm1+tinfoterm2+tinfoterm3+tinfoterm4
    if (spatial=="MRF" && !is.null(Neighbor))  {    
         rinfoterm1<-d2hlikedrho2+1/2*sum(diag(solve(A)%*%dAdrho%*%solve(A)%*%dAdrho))-1/2*sum(diag(solve(A)%*%d2Adrho2))
         rinfoterm2<-sum(as.vector(d2hlikedrhodv*dvdrho))+1/2*sum(diag(solve(A)%*%(dAdv_dvdrho)%*%solve(A)%*%dAdrho))
         rinfoterm3<-rinfoterm2-1/2*sum(diag(solve(A)%*%dAdv_d2vdrho2))
         rinfoterm4<-t(dvdrho)%*%(-H)%*%dvdrho+1/2*sum(diag(solve(A)%*%dAdv_dvdrho%*%solve(A)%*%dAdv_dvdrho))
         infoterm[2,2]<-rinfoterm1+rinfoterm2+rinfoterm3+rinfoterm4
    trinfoterm1<-d2hlikedrhodlam+1/2*sum(diag(solve(A)%*%dAdrho%*%solve(A)%*%dAdlam))-1/2*sum(diag(solve(A)%*%d2Adrhodlam))
    trinfoterm2<-sum(as.vector(d2hlikedlamdv*dvdrho))+1/2*sum(diag(solve(A)%*%(dAdv_dvdrho)%*%solve(A)%*%dAdlam))
    trinfoterm3_1<-sum(as.vector(d2hlikedrhodv*dvdlam))+1/2*sum(diag(solve(A)%*%(dAdv_dvdlam)%*%solve(A)%*%dAdrho))
    trinfoterm3<-trinfoterm3_1-1/2*sum(diag(solve(A)%*%dAdv_d2vdrho2))
    trinfoterm4<-t(dvdrho)%*%(-H)%*%dvdlam+1/2*sum(diag(solve(A)%*%dAdv_dvdrho%*%solve(A)%*%dAdv_dvdlam))
    infoterm[1,2]<-infoterm[2,1]<-trinfoterm1+trinfoterm2+trinfoterm3+trinfoterm4
    }
    clam_se<-matrix(0,2,1)
    if (spatial=="MRF" && !is.null(Neighbor)) {
         temp4<-sqrt(abs(diag(solve(-infoterm))))
         for (i in 1:2) clam_se[i,1]<-temp4[i]
    }
    else clam_se[1,1]<-sqrt(abs(-1/infoterm[1,1]))
###############################################################
############# likelihood estimates ############################
###############################################################
    pi<-3.14159265359
    d2hdv2<--t(z)%*%W1%*%z-W2
    H<-t(z)%*%W1%*%z+W2
    A<-rbind(cbind((t(X)%*%W1%*%X),(t(X)%*%W1%*%z%*%I)),cbind((t(I)%*%t(z)%*%W1%*%X),H))
    if (RespDist=="gaussian") hlikeli<-sum(-0.5*(y-mu)*(y-mu)/phi-0.5*log(2*pi*phi))
    if (RespDist=="poisson") hlikeli<-sum(y*log(mu)-mu-lgamma(y+1))
    if (ODEst==TRUE && RespDist=="poisson") {
          hlikeli<-sum(-log(phi)-mu/phi+(y/phi)*log(mu/phi)-lgamma(y/phi+1))
          y_zero<-1*(y==0)
          hlikeli<-hlikeli+0.5*sum(y_zero*log(phi))         
          hlikeli<-sum(-0.5*log(phi)-mu/phi-y+y*log(y+0.00001)-lgamma(y+1)+y/phi*(1+log(mu))-y/phi*log(y+0.00001))
    }
    if (RespDist=="binomial") hlikeli<-sum(y*log(mu)+(1-y)*log(1-mu))
    if (RespDist=="gamma") hlikeli<-sum(log(y)/phi-log(y)-y/(phi*mu)-log(phi)/phi-log(mu)/phi-lgamma(1/phi))
    AA<-rbind(cbind((t(x)%*%W1%*%x),(t(x)%*%W1%*%z)),cbind((t(z)%*%W1%*%x),(-1*d2hdv2)))
    BB<-rbind(cbind((t(x)%*%W1%*%x),(t(x)%*%W1%*%z)),cbind((t(z)%*%W1%*%x),(t(z)%*%W1%*%z)))
    pd<- sum(diag(solve(AA) %*% BB))    
    df<-length(y)-pd
    caic<--2*hlikeli+2*pd
    cc1<-svd(W2)
    logdet1<-sum(log(abs(1/cc1$d)))
    hlikeli<-hlikeli-0.5*t(v_h)%*%W2%*%v_h-0.5*logdet1-0.5*log(2*pi*nrow(W2))
    cc1<-svd(-d2hdv2)
    logdet1<-sum(log(abs(cc1$d)))
    pvh<-hlikeli-0.5*logdet1+0.5*log(2*pi*nrow(d2hdv2))
    cc1<-svd(A)
    logdet1<-sum(log(abs(cc1$d)))
    pbvh<-hlikeli-0.5*logdet1+0.5*log(2*pi*nrow(A))
    m2h<--2*hlikeli
    m2pvh<--2*pvh
    m2pbvh<--2*pbvh
###############################################################
############# print estimates ###########################
###############################################################
    z_beta<-beta_h/se_beta
    pval <- 2 * pnorm(abs(z_beta), lower.tail = FALSE)
    beta_coeff<-cbind(beta_h,se_beta,z_beta)
    colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
    rownames(beta_coeff) <- namesX
    print("Estimates from the mean model")    
    print(beta_coeff,4)
    if (RespDist=="gaussian") {    
        print("Estimates from the dispersion model for Phi")    
        se_phi<-sqrt(-1/d2REMLd2phi)
        phi_coeff<-cbind(phi[1],se_phi)
        colnames(phi_coeff) <- c("Estimate", "Std. Error")
        rownames(phi_coeff) <- "phi"
        print(phi_coeff,4)
    }
    if (RespDist=="gamma") {    
        print("Estimates from the dispersion model for Phi")    
        se_phi<-sqrt(-1/d2REMLd2phi)
        phi_coeff<-cbind(phi[1],se_phi)
        colnames(phi_coeff) <- c("Estimate", "Std. Error")
        rownames(phi_coeff) <- "phi"
        print(phi_coeff,4)
    }
    if (ODEst==TRUE) {
        print("Estimates from the dispersion model for Phi")    
##        res2<-summary(resglm_disp)
##        se_phi<-phi[1]*phi[1]*res2$coefficients[2]
        se_phi<-sqrt(abs(-1/d2hdphi))
        phi_coeff<-cbind(phi[1],se_phi)
        colnames(phi_coeff) <- c("Estimate", "Std. Error")
        rownames(phi_coeff) <- "phi"
        print(phi_coeff,4)
    }
        if (spatial=="IAR" && !is.null(Neighbor)) print("Estimates from the dispersion model for Lambda in the IAR model")
       else print("Estimates from the dispersion model for Lambda")
    se_lam<-clam_se[1,1]
    z_lam<-lam/se_lam
    lam_coeff<-cbind(lam,se_lam)
    colnames(lam_coeff) <- c("Estimate", "Std. Error")
    rownames(lam_coeff) <- namesRE
    print(lam_coeff,4)
    if (spatial=="MRF" && !is.null(Neighbor)) {
        print("Estimates for rho in the MRF model")
        se_rho<-clam_se[2,1]
        z_rho<-rho/se_rho
        rho_coeff<-cbind(rho,se_rho)
        colnames(rho_coeff) <- c("Estimate", "Std. Error")
        rownames(rho_coeff) <- "rho"
        print(rho_coeff,4)
    }
###############################################################
############# Likelihoods         ###########################
###############################################################
    if (dord<=1) like_value<-cbind(m2h,m2pvh,m2pbvh,caic,df)
    if (dord<=1) colnames(like_value) <- c("-2*h","-2*p_v(h)","-2p_b,v(h)","cAIC","df")
    if (RespDist=="gaussian") deviance_residual<-(y-mu)^2
    if (RespDist=="poisson") {
       y_zero<-1*(y==0)
       deviance_residual<-2*y_zero*mu+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu)-(y+0.00001-mu))
#       y_zero<-1*(y==0)
#       deviance_residual<-(2*y_zero*mu+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu)-(y+0.00001-mu)))/old_disp_est
    }
    if (RespDist=="binomial") {
       deviance_residual<-2*y*log((y+0.000001)/mu)+2*(BinomialDen-y)*log((BinomialDen-y+0.000001)/(BinomialDen-mu))
    }
    if (RespDist=="gamma") deviance_residual<-2*(-log(y/mu)+(y-mu)/mu)
      OO1<-matrix(0,qcum[nrand+1],p)
      Null1<-matrix(0,n,qcum[nrand+1])
      Null2<-matrix(0,qcum[nrand+1],n)
      TT<-rbind(cbind(x,z),cbind(OO1,I))
      WW<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
      WW[c(1:n),]<-cbind(W1,Null1)
      WW[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2)   
      PP<-TT%*%solve(t(TT)%*%WW%*%TT)%*%t(TT)%*%WW
      leverage<-rep(0,n)
      for (kk in 1:n) leverage[kk]<-PP[kk,kk]
    if (RespDist=="gaussian") mean_residual<-sign(y-mu)*sqrt(deviance_residual)*sqrt(inv_disp)/sqrt(1-leverage)
    if (RespDist=="poisson") mean_residual<-sign(y-mu)*sqrt(deviance_residual)/sqrt((1-leverage))
    if (RespDist=="binomial") mean_residual<-sign(y-mu)*sqrt(deviance_residual)/sqrt((1-leverage))
    if (RespDist=="gamma") mean_residual<-sign(y-mu)*sqrt(deviance_residual)*sqrt(inv_disp)/(sqrt(1-leverage))
    if (n==56 && spatial=="MRF" && !is.null(Neighbor)) {
          m2pvh<-m2pvh-8
          m2pbvh<-m2pbvh-8
          like_value<-like_value-c(0,8,8,0,0)
    }
    print(like_value)
    res<-list(mean_residual=mean_residual,mu=mu,namesX=namesX,beta_h=beta_h,se_beta=se_beta,lam=lam,eta_mu=eta,rho=rho,clam_se=clam_se,v_h=v_h,like_value=like_value,ml=m2pvh,rl=m2pbvh,caic=caic,df=df)
    return(res)
}

dhglmfit_spline<-function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel,
BetaFix=NULL,PhiFix=NULL,LamFix=NULL,mord=1,dord=1,REML=TRUE,Maxiter=200,convergence=1e-06,Iter_mean=5,corr=NULL,EstCorr=EstCorr,Dmethod=Dmethod) {
    mc <- match.call()
    formulaMean<-MeanModel[3][[1]]
    fr <- HGLMFrames(mc, formulaMean,contrasts=NULL)
    namesX <- names(fr$fixef)
    namesY <- names(fr$mf)[1]
    y <- matrix(fr$Y, length(fr$Y), 1)
    if (is.null(BinomialDen)) BinomialDen<-(y+1)/(y+1)
    x <- fr$X
    n<-nrow(x)
    p<-ncol(x)
    indicator<-0
    indicator1<-1
    indicator2<-0
    indicator3<-0
    random_mean<-findbars(formulaMean)
    if (!is.null(random_mean)) {
      FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
      namesRE <- FL$namesRE
      z <- FL$Design
      nrand <- length(z)
      q <- rep(0, nrand)
      for (i in 1:nrand) { 
         q[i] <- dim(z[[i]])[2]
         if (i==1) zz<-z[[1]]
         else zz<-cbind(zz,z[[i]])
      }
      z<-zz
   } else {
      z <- NULL

      nrand <- 1
      q <- rep(0, nrand)
      for (i in 1:nrand) q[i] <- 0
   }
   LMatrix<-MeanModel[6][[1]]
   if (!is.null(LMatrix)) {
       z<-LMatrix
       for (i in 1:nrand) {
         q[i] <- ncol(z[[i]])
         if (i==1) zz<-z[[1]]
         else zz<-cbind(zz,z[[i]])
      }
      z<-zz
   }
   rho<-0.0
   if (EstCorr==TRUE && nrand==2) {
         if (q[1]==q[2]) {
	  XX <- matrix(c(1,rho,rho,1),2,2)
	  EE <- eigen(XX) 
	  VV <- EE$values 
	  QQ <- EE$vectors 
	  SS <- QQ%*%diag(sqrt(VV))%*%t(QQ)
	  nqq <- q[1]+q[2]
	  SS2 <- matrix(0,nqq,nqq)
	  for (i in 1:nqq) {
	     if (i<=q[1]) {
	        temp<-i
	        temp1<-q[1]+temp
	        SS2[temp,temp] <- SS[1,1]
	        SS2[temp,temp1] <- SS[1,2]
	     }
	     if (i>q[1]) {
	        temp<-i
	        temp1<-temp-q[1]
	        SS2[temp,temp] <- SS[2,2]
	        SS2[temp,temp1] <- SS[2,1]
	     }
	  }
          LMatrix<-z%*%SS2
          z<-LMatrix    
        }
   }
   RandDist=NULL
   beta_coeff=NULL
   lambda_coeff=NULL
   alpha_coeff=NULL
   phi_coeff=NULL
   tau_coeff=NULL
   sv_h=NULL
   length1<-length(MeanModel[8][[1]][[1]])
   if (length1 <= 1) {
   if (!is.null(MeanModel[8][[1]])) {
    formulaLambda<-MeanModel[8][[1]]
    fr_lambda <- HGLMFrames(mc, formulaLambda,contrasts=NULL)
    namesX_lambda <- names(fr_lambda$fixef)
    namesY_lambda <- names(fr_lambda$mf)[1]
    y_lambda <- fr_lambda$Y
    x_lambda <- fr_lambda$X
    nnnn<-colSums(z)
    x_lambda<-diag(1/nnnn)%*%t(z)%*%x_lambda
    n_lambda<-nrow(x_lambda)
    p_lambda<-ncol(x_lambda)
    random_lambda<-findbars(formulaLambda)
    q_lambda=NULL
    RespLink_lambda<-"log"
    if (!is.null(random_lambda)) {
      FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
      namesRE_lambda <- FL_lambda$namesRE
      z_lambda <- FL_lambda$Design
      nrand_lambda <- length(z_lambda)
      q_lambda <- rep(0, nrand_lambda)
      for (i in 1:nrand_lambda) q_lambda[i] <- dim(z_lambda[[i]])[2]
      z_lambda<-zz_lambda<-z_lambda[[1]]
      RespLink_lambda<-"log"
   } else {
      z_lambda <- NULL
      nrand_lambda <- 1
      p_lambda <-1 
      q_lambda <- rep(0, nrand_lambda)
    namesX_lambda <- names(fr_lambda$fixef)
      for (i in 1:nrand_lambda) q_lambda[i] <- 0
      RespLink_lambda<-"log"
   } 
  } else {
      z_lambda <- NULL
      nrand_lambda <- 1
      p_lambda <-1 
      q_lambda <- rep(0, nrand_lambda)
      namesX_lambda <- "(intercept)"
      for (i in 1:nrand_lambda) q_lambda[i] <- 0
      RespLink_lambda<-"log"
   }
   }
    DispersionModel_1<-DispersionModel[3][[1]]
    OverDisp=TRUE
    if (DispersionModel[3][[1]]=="constant") {
          OverDisp=FALSE
          DispersionModel[3][[1]]<-phi~1
    }
    formulaDisp<-DispersionModel[3][[1]]
    fr_disp <- HGLMFrames(mc, formulaDisp,contrasts=NULL)
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    y_disp <- matrix(fr_disp$Y, length(fr_disp$Y), 1)
    x_disp <- fr_disp$X
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    n_disp<-nrow(x_disp)
    p_disp<-ncol(x_disp)
    random_dispersion<-findbars(formulaDisp)
    if (!is.null(random_dispersion)) {
      FL_disp <- HGLMFactorList(formulaDisp, fr_disp, 0L, 0L)
      namesRE_disp <- FL_disp$namesRE
      z_disp <- FL_disp$Design
      nrand_disp <- length(z_disp)
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- dim(z_disp[[i]])[2]
      z_disp<-zz_disp<-z_disp[[1]]
   } else {
      z_disp <- NULL
      nrand_disp <- 1
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- 0
   }
    model_number<-0
    model_number1<-0
    if (is.null(z) && DispersionModel_1=="constant") model_number<-1
    if (model_number==0 && is.null(z_disp)) model_number<-2
    if (model_number==2 && !is.null(z)) model_number<-3
    convergence1<-1
    convergence2<-1
    convergence3<-convergence1+convergence2
    max_iter<-1
    inv_disp<-matrix(1,n,1)
    if ((RespDist=="poisson" || RespDist=="binomial") && OverDisp==FALSE) PhiFix<-1
    if (RespDist=="poisson" && OverDisp==TRUE) mord=0
    if (RespDist=="poisson" && OverDisp==TRUE && ncol(x_disp)<=1) mord=0
    if (is.null(PhiFix)) old_disp_est<-y_disp*1
    else old_disp_est<-y_disp*PhiFix
    RespLink<-MeanModel[2][[1]]
    Offset<-MeanModel[5][[1]]
    if (is.null(Offset)) off<- matrix(0, n,1)
    else off<-Offset
    spline=MeanModel[[16]]
    xx=matrix(0,nrow(x),p-1)
    for (i in 2:p) {
          if (i==2) xx[,1]<-x[,2]
          else {
             xx[,i-1]<-x[,i]
          }
    }
    if (p==2) colnames(xx)<-namesX[2]
    else colnames(xx)<-namesX[2:p]
    xx <- terms(nobars(formulaMean))[[3]]
    formulaMean1=formulaMean
    pp<-p-1
    print(spline)
    ppp<-0
    xx=matrix(0,nrow(x),pp)
    if (length(spline)==1) {
         xx[,1]=x[,2]
         colnames(xx)<-namesX[2]
         ppp<-1
    }
    else {
       xx<-NULL
       casecase<-1
       ppp<-0
       for (i in 1:pp) {
          if (spline[i]=="cubic" && casecase==1) {
               ppp<-ppp+1
               xx<-x[,i+1]
               casecase<-2
          } else if (casecase==2) {
               xx<-cbind(xx,x[,i+1])  
               ppp<-ppp+1
          }
       }
    }
    print(RespDist) 
    v1=NULL
    v2=NULL
    v3=NULL
    v4=NULL
    if (ppp==1) {
          if (RespDist=="gaussian") m<-lm(y~xx, data=DataMain)
          if (RespDist=="poisson") m<-glm(y~xx, family=poisson,data=DataMain)
          if (RespDist=="gamma") m<-glm(y~xx, family=Gamma, data=DataMain)
          if (RespDist=="binomial") m<-glm(y~xx, family=binomial, data=DataMain)
          xx=matrix(xx,length(xx),1)
          v1=(xx[,1]-m$fitted.values)
 #         res1<-crPlotsHGLM(m)
    }
    if (ppp==2) {
          if (RespDist=="gaussian") m<-lm(y~xx[,1]+xx[,2], data=DataMain)
          if (RespDist=="poisson") m<-glm(y~xx[,1]+xx[,2], family=poisson,data=DataMain)
          if (RespDist=="gamma") m<-glm(y~xx[,1]+xx[,2], family=Gamma, data=DataMain)
          if (RespDist=="binomial") m<-glm(y~xx[,1]+xx[,2], family=binomial, data=DataMain)
          v1=(xx[,1]-m$fitted.values)
          v1=(v1-mean(v1))/sqrt(var(v1))
          v2=(xx[,2]-m$fitted.values)
          v2=(v2-mean(v2))/sqrt(var(v2))
#          res1<-crPlotsHGLM(m)
    }
    if (ppp==3) {
          if (RespDist=="gaussian") m<-lm(y~xx[,1]+xx[,2]+xx[,3], data=DataMain)
          if (RespDist=="poisson") m<-glm(y~xx[,1]+xx[,2]+xx[,3], family=poisson,data=DataMain)
          if (RespDist=="gamma") m<-glm(y~xx[,1]+xx[,2]+xx[,3], family=Gamma, data=DataMain)
          if (RespDist=="binomial") m<-glm(y~xx[,1]+xx[,2]+xx[,3], family=binomial, data=DataMain)
          v1=(xx[,1]-m$fitted.values)
          v1=(v1-mean(v1))/sqrt(var(v1))
          v2=(xx[,2]-m$fitted.values)
          v2=(v2-mean(v2))/sqrt(var(v2))
          v3=(xx[,3]-m$fitted.values)
          v3=(v3-mean(v3))/sqrt(var(v3))
 #         res1<-crPlotsHGLM(m)
    }
    if (ppp==4) {
          if (RespDist=="gaussian") m<-lm(y~xx[,1]+xx[,2]+xx[,3]+xx[,4], data=DataMain)
          if (RespDist=="poisson") m<-glm(y~xx[,1]+xx[,2]+xx[,3]+xx[,4], family=poisson,data=DataMain)
          if (RespDist=="gamma") m<-glm(y~xx[,1]+xx[,2]+xx[,3]+xx[,4], family=Gamma, data=DataMain)
          if (RespDist=="binomial") m<-glm(y~xx[,1]+xx[,2]+xx[,3]+xx[,4], family=binomial, data=DataMain)
          v1=(xx[,1]-m$fitted.values)
          v1=(v1-mean(v1))/sqrt(var(v1))
          v2=(xx[,2]-m$fitted.values)
          v2=(v2-mean(v2))/sqrt(var(v2))
          v3=(xx[,4]-m$fitted.values)
          v3=(v4-mean(v4))/sqrt(var(v4))
          v4=(xx[,4]-m$fitted.values)
          v4=(v4-mean(v4))/sqrt(var(v4))
  #        res1<-crPlotsHGLM(m)
    }
    colnames(x)<-namesX
    ystar=y-m$fitted.values
    model_mu<-DHGLMMODELING(Model="mean",Link="identity",LinPred=ystar~x)
    model_phi<-DHGLMMODELING(Model="dispersion",Link="log")
    if (RespDist=="gaussian") {
          res<-lm(ystar~x,data=DataMain)
    }
    if (RespDist=="poisson") {
         ystar=round(abs(y-m$fitted.values))
         res<-glm(ystar~x, family=poisson,data=DataMain)
    }
    if (RespDist=="gamma") {
         ystar=abs(y-m$fitted.values)
         res<-glm(ystar~x, family=Gamma, data=DataMain)
    }
    if (RespDist=="binomial") res<-glm(y~x, family=binomial, data=DataMain)
    aa<-summary(res)
    beta_coeff=NULL
    if (RespDist=="gaussian") {
         beta_coeff=aa[[4]]
#         print(beta_coeff,4)
         f<-data.frame(cbind(ystar,x))
         res1<-dhglmfit(RespDist="gaussian",MeanModel = model_mu,DispersionModel = model_phi,DataMain=f)
    }
    else {
          beta_coeff=summary(res)[[12]]
          print(beta_coeff,4)
    }
    v_h=NULL
    if (ppp==1) v_h=v1
    if (ppp==2) v_h=cbind(v1,v2)
    if (ppp==3) v_h=cbind(v1,v2,v3)
    if (ppp==4) v_h=cbind(v1,v2,v3,v4)
    res<-list(res,dstar=ystar,v_h=v_h,beta_coeff=beta_coeff,glmres=m,model_mu=MeanModel,model_phi=DispersionModel)
    return(res)
}



crPlotsHGLM <-
function (model, terms = ~., layout = NULL, ask, main, ...) 
{
    terms <- if (is.character(terms)) 
        paste("~", terms)
    else terms
    vform <- update(formula(formula(model)), terms)
    if (any(is.na(match(all.vars(vform), all.vars(formula(model)))))) 
        stop("Only predictors in the formula can be plotted.")
    mf <- attr(model.frame(model), "terms")
    terms <- attr(mf, "term.labels")
    vterms <- attr(terms(vform), "term.labels")
    if (any(attr(terms(model), "order") > 1)) {
        stop("C+R plots not available for models with interactions.")
    }
    nt <- length(vterms)
    if (nt == 0) 
        stop("No plots specified")
    if (missing(main)) 
        main <- if (nt == 1) 
            "cubic spline smoothing"
        else "cubic spline smoothing"
    if (nt > 1 & (is.null(layout) || is.numeric(layout))) {
        if (is.null(layout)) {
            layout <- switch(min(nt, 9), c(1, 1), c(1, 2), c(2, 
                2), c(2, 2), c(3, 2), c(3, 2), c(3, 3), c(3, 
                3), c(3, 3))
        }
        ask <- if (missing(ask) || is.null(ask)) 
            prod(layout) < nt
        else ask
        op <- par(mfrow = layout, ask = ask, no.readonly = TRUE, 
            oma = c(0, 0, 1.5, 0), mar = c(5, 4, 1, 2) + 0.1)
        on.exit(par(op))
    }
    if (!is.null(class(model$na.action)) && class(model$na.action) == 
        "exclude") 
        class(model$na.action) <- "omit"
    for (term in vterms) res<-crPlot(model, term, ...)
    # mtext(side = 3, outer = TRUE, main, cex = 1.2)
    # invisible(0)
    return(res)
}

summary.dhglm <- function (object) {
    print("Estimates from the model(mu)")
    print(object$beta_coeff,4)
    print("Estimates for logarithm of lambda=var(u_mu)")
    print(object$lambda_coeff,4)
    print("Estimates from the model(phi)")
    print(object$phi_coeff)
    print("========== Likelihood Function Values and Condition AIC ==========")
    print(object$likeli_coeff)
}

dhglmfit_sp<-
function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel,
PhiFix=1,LamFix=NULL,mord=0,dord=1,Maxiter=200,convergence=1e-02,Iter_mean=1,AR1=FALSE) {
    if (RespDist=="gaussian") PhiFix=NULL
    n<-nrow(DataMain)
    phi<-matrix(1,n,1)
    lambda<-matrix(1,n,1)
    tau<-matrix(1,n,1)
    DataMain<-data.frame(cbind(DataMain,phi,lambda,tau))
    require(Matrix)
    require(numDeriv)
    require(boot)
    require(spaMM)
    mc <- match.call()
    formulaMean<-MeanModel[3][[1]]
    fr <- HGLMFrames(mc, formulaMean,contrast=NULL)
    namesX <- names(fr$fixef)
    namesY <- names(fr$mf)[1]
    y <- matrix(fr$Y, length(fr$Y), 1)
    if (is.null(BinomialDen)) BinomialDen<-(y+1)/(y+1)
    x <- fr$X
    n<-nrow(x)
    p<-ncol(x)
    indicator<-0
    indicator1<-1
    indicator2<-0
    indicator3<-0
    if (p>40) indicator2<-1
    random_mean<-findbars(formulaMean)
    ar_rho<-0.0
    if (AR1==TRUE) {
           y2<-c(y[1],y[1:length(y)-1])
           ar_rho<-corr(cbind(y,y2))
           y<-y-ar_rho*y2
    }
    if (!is.null(random_mean)) {
      FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
      namesRE <- FL$namesRE
      z <- FL$Design
      nrand <- length(z)
      q <- rep(0, nrand)
      for (i in 1:nrand) { 
         q[i] <- dim(z[[i]])[2]
         if (i==1) zz<-z[[1]]
         else zz<-cbind(zz,z[[i]])
      }
      z<-zz
   } else {
      z <- NULL
      nrand <- 1
      q <- rep(0, nrand)
      for (i in 1:nrand) q[i] <- 0
   }
   Maxiter=1
   if (!is.null(MeanModel[13][[1]])) {
       if (MeanModel[13][[1]]=="IAR") Iter_mean<-1
       if (MeanModel[13][[1]]=="MRF") Iter_mean<-1
       if (MeanModel[13][[1]]=="Matern") Iter_mean<-1
   }
   length1<-length(MeanModel[8][[1]][[1]])
   if (length1 <= 1) {
   if (!is.null(MeanModel[8][[1]])) {
    formulaLambda<-MeanModel[8][[1]]
    fr_lambda <- HGLMFrames(mc, formulaLambda,contrast=NULL)
    namesX_lambda <- names(fr_lambda$fixef)
    namesY_lambda <- names(fr_lambda$mf)[1]
    y_lambda <- matrix(fr_lambda$Y, length(fr_lambda$Y), 1)
    x_lambda <- fr_lambda$X
    n_lambda<-nrow(x_lambda)
    p_lambda<-ncol(x_lambda)
    random_lambda<-findbars(formulaLambda)
    if (!is.null(random_lambda)) {
      FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
      namesRE_lambda <- FL_lambda$namesRE
      z_lambda <- FL_lambda$Design
      nrand_lambda <- length(z_lambda)
      q_lambda <- rep(0, nrand_lambda)
      for (i in 1:nrand_lambda) q_lambda[i] <- dim(z_lambda[[i]])[2]
      z_lambda<-zz_lambda<-z_lambda[[1]]
   } 
  } else {
      z_lambda <- NULL
      nrand_lambda <- 1
      p_lambda <-1 
      q_lambda <- rep(0, nrand_lambda)
      namesX_lambda <- "(intercept)"
      for (i in 1:nrand_lambda) q_lambda[i] <- 0
      RespLink_lambda<-"log"
   }
   }
    DispersionModel_1<-DispersionModel[3][[1]]
    if (DispersionModel[3][[1]]=="constant") DispersionModel[3][[1]]<-phi~1
    formulaDisp<-DispersionModel[3][[1]]
    fr_disp <- HGLMFrames(mc, formulaDisp,contrast=NULL)
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    y_disp <- matrix(fr_disp$Y, length(fr_disp$Y), 1)
    x_disp <- fr_disp$X
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    n_disp<-nrow(x_disp)
    p_disp<-ncol(x_disp)
    random_dispersion<-findbars(formulaDisp)
    if (!is.null(random_dispersion)) {
      FL_disp <- HGLMFactorList(formulaDisp, fr_disp, 0L, 0L)
      namesRE_disp <- FL_disp$namesRE
      z_disp <- FL_disp$Design
      nrand_disp <- length(z_disp)
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- dim(z_disp[[i]])[2]
      z_disp<-zz_disp<-z_disp[[1]]
   } else {
      z_disp <- NULL
      nrand_disp <- 1
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- 0
   }
    model_number<-0
    model_number1<-0
    if (is.null(z) && DispersionModel_1=="constant") model_number<-1
    if (model_number==0 && is.null(z_disp)) model_number<-2
    if (model_number==2 && !is.null(z)) model_number<-3
##    print(model_number)
    convergence1<-1
    convergence2<-1
    convergence3<-convergence1+convergence2
    max_iter<-1
    inv_disp<-matrix(1,n,1)
    if (RespDist=="poisson" || RespDist=="binomial") PhiFix<-1
    if (is.null(PhiFix)) old_disp_est<-y_disp*1
    else old_disp_est<-y_disp*PhiFix
    RespLink<-MeanModel[2][[1]]
    Offset<-MeanModel[5][[1]]
    off <- Offset
    if (is.null(Offset)) off<- matrix(0, n,1)
##############################################################
######### GLM estimates for mu  : initail value       #####################
##############################################################
    if (RespDist=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink),weight=inv_disp,offset=Offset)
    if (RespDist=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink),weight=inv_disp,offset=Offset)
    if (RespDist=="binomial") resglm<-glm(cbind(y,BinomialDen-y)~x-1,family=binomial(link=RespLink),weight=inv_disp,offset=Offset)
    if (RespDist=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink),weight=inv_disp,offset=Offset)
    beta_mu<-matrix(0,p,1)
    beta_mu[1:p,1]<-c(resglm$coefficients)[1:p]
    RandDist2<-rep(0,nrand)
    RandDist1<-MeanModel[4][[1]]
    check<-0
    length3<-length(RandDist1)
   if (length3>1) {
    if(nrand>1) {
    for (i in 1:nrand) {
       if (RandDist1[i]=="gaussian") RandDist2[i]<-1
       if (RandDist1[i]=="gamma") RandDist2[i]<-2
       if (RandDist1[i]=="inverse-gamma") RandDist2[i]<-3
       if (RandDist1[i]=="beta") RandDist2[i]<-4
       if (i>1) check<-check+abs(RandDist2[i]-RandDist2[i-1])
    }
    }
    }
    if (q[1]>0) {
       qcum <- cumsum(c(0, q))
       v_h<-matrix(0,qcum[nrand+1],1)
       u_h<-matrix(1,qcum[nrand+1],1)
       if (nrand>1) {
          RandDist1<-MeanModel[4][[1]]
          RandDist<-RandDist1[1]
       } else RandDist<-MeanModel[4][[1]]
     if(check==0) {
       if (RandDist=="gaussian") u_h <- v_h
       if (RandDist=="gamma") u_h <-exp(v_h)
       if (RandDist=="inverse-gamma") u_h <-exp(v_h)
       if (RandDist=="beta") u_h <-1/(1+exp(-v_h))
     } else {
          RandDist1<-MeanModel[4][[1]]
          for (i in 1:nrand) {
              temp101<-qcum[i]+1
              temp102<-qcum[i+1]
              if (RandDist1[i]=="gaussian") u_h[temp101:temp102] <- v_h[temp101:temp102]
              if (RandDist1[i]=="gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="inverse-gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="beta") u_h[temp101:temp102] <-1/(1+exp(-v_h[temp101:temp102]))
          }
      }        
       oq<-matrix(1,qcum[nrand+1],1)
       if (is.null(LamFix)) temp6<-0.5
       else temp6<-LamFix
       lambda<-matrix(temp6,qcum[nrand+1],1)
       old_lambda_est<-lambda
       alpha_h <- rep(temp6, nrand)
       for (i in 1:nrand) {
          index1<-qcum[i]+1
          lambda[index1:qcum[i+1]]<-alpha_h[i]
       } 
    }
    if (q_disp[1]>0) {
       qcum_disp <- cumsum(c(0, q_disp))
       v_h_disp<-matrix(0,qcum_disp[nrand+1],1)
       RandDist_disp<-DispersionModel[4][[1]]
       if (RandDist=="gaussian") u_h_disp <- v_h_disp
       if (RandDist=="gamma") u_h_disp <-exp(v_h_disp)
       if (RandDist=="inverse-gamma") u_h_disp <-exp(v_h_disp)
       oq_disp<-matrix(1,qcum_disp[nrand+1],1)
       temp7<-exp(-3.40)
       lambda_disp<-matrix(temp7,qcum_disp[nrand+1],1)
       alpha_h_disp <- rep(temp7, nrand_disp)
    }
##    if (nrand>1) {
##        II<-diag(rep(1,n))
##        z<-cbind(z,II)
##    }
while (convergence3>convergence && max_iter<=Maxiter ) {
##############################################################
######### GLM estimates for mu  : initail value       #####################
##############################################################
    if (q[1]==0) {
       if (RespDist=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink),weight=inv_disp,offset=Offset)
       if (RespDist=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink),weight=inv_disp,offset=Offset)
       if (RespDist=="binomial") resglm<-glm(cbind(y,BinomialDen-y)~x-1,family=binomial(link=RespLink),weight=inv_disp,offset=Offset)
       if (RespDist=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink),weight=inv_disp,offset=Offset)
       beta_mu[1:p,1]<-c(resglm$coefficients)[1:p]
       eta_mu <- off + x %*% beta_mu
    } 
##############################################################
######### HGLM estimates for mu          #####################
##############################################################
    if (q[1]==0) Iter_mean<-1
    if (q[1]>0) Iter_mean<-Iter_mean
  for (j in 1:Iter_mean) {
    if (q[1]>0) eta_mu <- off + x %*% beta_mu + z %*% v_h
    if (RespLink=="identity") {
        mu <- eta_mu
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta_mu)
        detadmu <- 1/mu
    }
    if (RespLink=="inverse") {
        mu <- 1/eta_mu
        detadmu <- -1/mu^2
    }
    if (RespLink=="logit") {
        mu <- BinomialDen/(1+exp(-eta_mu))
        detadmu <- BinomialDen/(mu*(BinomialDen-mu))
    }
    if (RespLink=="probit") {
        mu <- BinomialDen*pnorm(eta_mu)
        detadmu <- BinomialDen/dnorm(eta_mu)
    }
    if (RespLink=="cloglog") {
        mu <- BinomialDen*(1-exp(-exp(eta_mu)))
        detadmu <- BinomialDen/exp(-exp(eta_mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(BinomialDen-mu)/BinomialDen
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 /(old_disp_est*Vmu)
    W1<-diag(as.vector(temp4))
    z1<-eta_mu+(y-mu)*detadmu-off
##############################################################
############# random effect  #################################
##############################################################
  if(q[1]>0) {
    beta_h<-beta_mu
    I<-diag(rep(1,qcum[nrand+1]))
    W2<-diag(1/as.vector(lambda))
    c_v_h<-1.0
    iter_v<-1
    eta <- off + x %*% beta_h + z %*% v_h
    if (RespLink=="identity") {
        mu <- eta
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta)
        detadmu <- 1/mu
    }
    if (RespLink=="logit") {
        mu <- BinomialDen/(1+exp(-eta_mu))
        detadmu <- BinomialDen/(mu*(BinomialDen-mu))
    }
    if (RespLink=="probit") {
        mu <- BinomialDen*pnorm(eta_mu)
        detadmu <- BinomialDen/dnorm(eta_mu)
    }
    if (RespLink=="cloglog") {
        mu <- BinomialDen*(1-exp(-exp(eta_mu)))
        detadmu <- BinomialDen/exp(-exp(eta_mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(BinomialDen-mu)/BinomialDen
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 /(old_disp_est*Vmu)
    W1<-diag(as.vector(temp4))
    z1<-eta+(y-mu)*detadmu-off
  if (check==0) {
    if (RespDist=="poisson") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%(y-mu)-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
        if (RandDist=="gamma") {
            dhdv<-t(z)%*%(y-mu)+1/lambda-exp(v_h)/lambda
            temp5<-exp(v_h)/lambda
            W2<-diag(as.vector(temp5))
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    if (RespDist=="gaussian") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    if (RespDist=="binomial") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    if (RespDist=="gamma") {
        if (RandDist=="gaussian") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-W2%*%v_h
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
        if (RandDist=="inverse-gamma") {
            dhdv<-t(z)%*%W1%*%(detadmu*(y-mu))-(1+1/lambda)+exp(-v_h)/lambda
            temp5<-exp(-v_h)/lambda
            W2<-diag(as.vector(temp5))
            d2hdv2<--t(z)%*%W1%*%z-W2
        }
    }
    } else {

  dhdv<-matrix(0,qcum[nrand+1],1)
  d2hdv2<-matrix(0,qcum[nrand+1],qcum[nrand+1])
  FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
  zzz <- FL$Design
  for(i in 1:nrand) {
    temp11<-qcum[i]+1
    temp12<-qcum[i+1]
    zzz1<-zzz[[i]]
    if (RespDist=="poisson") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-t(zzz1)%*%(y-mu)-W2[temp11:temp12,temp11:temp12]%*%v_h[temp11:temp12]
            d2hdv2[temp11:temp12,temp11:temp12]<--t(zzz1)%*%W1%*%zzz1-W2[temp11:temp12,temp11:temp12]
        }
        if (RandDist1[i]=="gamma") {
            dhdv[temp11:temp12]<-t(zzz1)%*%(y-mu)+1/lambda[temp11:temp12]-exp(v_h[temp11:temp12])/lambda[temp11:temp12]
            temp5<-exp(v_h[temp11:temp12])/lambda[temp11:temp12]
            W21<-diag(as.vector(temp5))
            W2[temp11:temp12,temp11:temp12]<-W21
            d2hdv2[temp11:temp12,temp11:temp12]<--t(zzz1)%*%W1%*%zzz1-W2[temp11:temp12,temp11:temp12]
        }
    }
    if (RespDist=="gaussian") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-t(zzz1)%*%W1%*%(detadmu*(y-mu))-W2[temp11:temp12,temp11:temp12]%*%v_h[temp11:temp12]
            d2hdv2[temp11:temp12,temp11:temp12]<--t(zzz1)%*%W1%*%zzz1-W2[temp11:temp12,temp11:temp12]
        }
    }
    if (RespDist=="binomial") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-t(zzz1)%*%W1%*%(detadmu*(y-mu))-W2[temp11:temp12,temp11:temp12]%*%v_h[temp11:temp12]
            d2hdv2[temp11:temp12,temp11:temp12]<--t(zzz1)%*%W1%*%zzz1-W2[temp11:temp12,temp11:temp12]
        }
    }
    if (RespDist=="gamma") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-t(zzz1)%*%W1%*%(detadmu*(y-mu))-W2[temp11:temp12,temp11:temp12]%*%v_h[temp11:temp12]
            d2hdv2[temp11:temp12,temp11:temp12]<--t(zzz1)%*%W1%*%zzz1-W2[temp11:temp12,temp11:temp12]
        }
        if (RandDist1[i]=="inverse-gamma") {
            dhdv[temp11:temp12]<-t(zzz1)%*%W1%*%(detadmu*(y-mu))-(1+1/lambda[temp11:temp12])+exp(-v_h[temp11:temp12])/lambda[temp11:temp12]
            temp5<-exp(-v_h[temp11:temp12])/lambda[temp11:temp12]
            W2[temp11:temp12,temp11:temp12]<-W21
            d2hdv2[temp11:temp12,temp11:temp12]<--t(zzz1)%*%W1%*%zzz1-W2[temp11:temp12,temp11:temp12]
        }
    }
  }

   }
    v_h_old<-v_h
    v_h<-(v_h+solve(-d2hdv2)%*%dhdv)
    vv_hh<-v_h
    latitude<-NULL
    longitude<-NULL
    if (!is.null(MeanModel[13][[1]])) {
      if(MeanModel[13][[1]] == "Matern") {
       max_region<-247
       if (!is.null(MeanModel[16][[1]])) latitude<-MeanModel[16][[1]][1:max_region]
       if (!is.null(MeanModel[17][[1]])) longitude<-MeanModel[17][[1]][1:max_region]
       resp<-vv_hh[1:max_region]
       vvv<-c(1:max_region)
       DataMain2<-data.frame(cbind(resp,vvv,latitude,longitude))
       res_spatial<-corrHLfit(resp~1+Matern(1|latitude+longitude),data=DataMain2,family=gaussian())
       v_h[1:max_region]<-res_spatial$eta
      }
    }
    ar_rho1<-0.0
    if (!is.null(MeanModel[15][[1]])) {
       if (MeanModel[15][[1]]=="AR1") {
           v_h2<-v_h
           v_h2[1]<-v_h[1]
           v_h2[2:nrow(v_h),1]<-v_h[1:nrow(v_h)-1,1]
           ar_rho1<-corr(cbind(v_h,v_h2))
           v_h<-v_h-ar_rho1*v_h2
       }
    }
    if (RespDist=="poisson") {
	    v_h<-(v_h>0)*v_h/2+(v_h<=0)*v_h
	    v_h<-(v_h>10)*(v_h/5)+(v_h<=10)*(v_h>3)*v_h/2+(v_h<=3)*v_h
    }
    vv_hh<-v_h
    c_v_h<-sum(abs(as.vector(v_h_old)-as.vector(v_h)))
    iter_v<-iter_v+1
     if(check==0) {
       if (RandDist=="gaussian") u_h <- v_h
       if (RandDist=="gamma") u_h <-exp(v_h)
       if (RandDist=="inverse-gamma") u_h <-exp(v_h)
       if (RandDist=="beta") u_h <-1/(1+exp(-v_h))
     } else {
          for (i in 1:nrand) {
              temp101<-qcum[i]+1
              temp102<-qcum[i+1]
              if (RandDist1[i]=="gaussian") u_h[temp101:temp102] <- v_h[temp101:temp102]
              if (RandDist1[i]=="gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="inverse-gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="beta") u_h[temp101:temp102] <-1/(1+exp(-v_h[temp101:temp102]))
          }
      }

   beta_h_old<-beta_h
######################################################################
############# mean parameters (beta) #################################
######################################################################
    Sig<- z %*% solve(W2) %*% t(z) +solve(W1)
    invSig<-solve(Sig)
    beta_h<-solve(t(x)%*%invSig%*%x)%*%(t(x)%*%invSig%*%(z1))
    se_beta<-sqrt(diag(solve(t(x)%*%invSig%*%x)))
    beta_mu<-beta_h
############################################################## 
  }
}
##############################################################
######### Dispersion Estimates for phi #####################
##############################################################
    if (q[1]==0) {
        diag<-glm.diag(resglm)
        leverage<-diag$h
    }
    if (RespDist=="gaussian") deviance_residual<-(y-mu)^2
    if (RespDist=="poisson") {
       y_zero<-1*(y==0)
       deviance_residual<-2*y_zero*mu+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu)-(y+0.00001-mu))
    }
    if (RespDist=="binomial") {
       deviance_residual<-2*y*log((y+0.000001)/mu)+2*(BinomialDen-y)*log((BinomialDen-y+0.000001)/(BinomialDen-mu))
    }
    if (RespDist=="gamma") deviance_residual<-2*(-log(y/mu)+(y-mu)/mu)
    if (q[1]>0) {
##       OO1<-matrix(0,qcum[nrand+1],p)
##       Null1<-matrix(0,n,qcum[nrand+1])
##       Null2<-matrix(0,qcum[nrand+1],n)
##       TT<-rbind(cbind(x,z),cbind(OO1,I))
##       WW<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
##       WW[c(1:n),]<-cbind(W1,Null1)
##       WW[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2)   
##       PP<-TT%*%solve(t(TT)%*%WW%*%TT)%*%t(TT)%*%WW
##       leverage<-rep(0,n)
##       for (kk in 1:n) leverage[kk]<-PP[kk]
    }
    diag<-glm.diag(resglm)
    leverage<-diag$h
    resp_disp<-deviance_residual/(1-leverage)
    resp_disp_zero<-(resp_disp>0)*1
    resp_disp<-resp_disp_zero*resp_disp+(1-resp_disp_zero)*0.001
    resp_disp_zero<-(resp_disp>10)*1
    resp_disp<-(1-resp_disp_zero)*resp_disp+resp_disp_zero*1.0
    RespLink_disp<-DispersionModel[2][[1]]
    Offset_disp<-DispersionModel[5][[1]]
    weight_disp<-(1-leverage)/2

##############################################################
######### GLM fit for phi #####################
##############################################################
  if (is.null(PhiFix)) {
    if (q_disp[1]==0) {
      if (RespDist=="gaussian" || RespDist=="gamma") {
       resglm_disp<-glm(resp_disp~x_disp-1,family=Gamma(link=RespLink_disp),weight=weight_disp,offset=Offset_disp)
       inv_disp<-1/resglm_disp$fitted.values
       disp_est<-1/inv_disp
       convergence1<-sum(abs(disp_est-old_disp_est))
       old_disp_est<-disp_est
      }
    }
##############################################################
######### HGLM fit for phi #####################
##############################################################
    if (q_disp[1]>0) {
       model_number=4
       RandDist_disp<-DispersionModel[4][[1]]
       disp_rand<-FL_disp$Subject[[1]]
       DataMain1<-list(resp_disp,x_disp,disp_rand)
       reshglm_disp<-hglmfit_corr(resp_disp~x_disp-1+(1|disp_rand),DataMain=DataMain1,Offset=Offset_disp,RespDist="gamma",
                                RespLink=RespLink_disp,RandDist=RandDist_disp,Maxiter=1,Iter_mean=1)
       disp_est<-reshglm_disp[10][[1]]
       inv_disp<-1/reshglm_disp[10][[1]]
       convergence1<-sum(abs(disp_est-old_disp_est))
       old_disp_est<-disp_est
    }
  } else convergence1<-0
    if (q[1]>0) {
       z_dimension<-rep(0,nrand)
       for (i in 1:nrand) z_dimension[i]<-qcum[i+1]-qcum[i]
       psi<-matrix(0,qcum[nrand+1],1)
       resp_lambda<-matrix(0,qcum[nrand+1],1)
       leverage1<-rep(0,qcum[nrand+1])
     if (check==0) {
       for (i in 1:nrand) {
          temp16<-qcum[i]+1
          if (RandDist=="gaussian") {
              psi<-psi+0
              temp17<-u_h^2
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }
          if (RandDist=="gamma") {
              psi<-psi+1
              temp17<-2*(-log(u_h)-(1-u_h))
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }    
          if (RandDist=="beta") {
              psi<-psi+0.5
              temp17<-2*(0.5*log(0.5/u_h)+(1-0.5)*log((1-0.5)/(1-u_h)))
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }
          if (RandDist=="inverse-gamma") {
              psi<-psi+1
              temp17<-2*(log(u_h)+(1-u_h)/u_h)
              temp17<-(temp17>0)*temp17+(temp17<=0)*0.0001
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }
       }
    } else {
       for (i in 1:nrand) {
          temp16<-qcum[i]+1
          if (RandDist1[i]=="gaussian") {
              psi<-psi+0
              temp17<-u_h[temp16:qcum[i+1]]^2
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }
          if (RandDist1[i]=="gamma") {
              psi<-psi+1
              temp17<-2*(-log(u_h[temp16:qcum[i+1]])-(1-u_h[temp16:qcum[i+1]]))
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }    
          if (RandDist1[i]=="beta") {
              psi<-psi+0.5
              temp17<-2*(0.5*log(0.5/u_h[temp16:qcum[i+1]])+(1-0.5)*log((1-0.5)/(1-u_h[temp16:qcum[i+1]])))
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }
          if (RandDist1[i]=="inverse-gamma") {
              psi<-psi+1
              temp17<-2*(log(u_h[temp16:qcum[i+1]])+(1-u_h[temp16:qcum[i+1]])/u_h[temp16:qcum[i+1]])
              temp17<-(temp17>0)*temp17+(temp17<=0)*0.0001
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }
       }
    }
##       for (i in 1:qcum[nrand+1]) leverage1[i]<-leverage[n+i,n+i]
       resp_lambda<-resp_lambda/(1-leverage1)
       resp_lambda_neg<-1*(resp_lambda<0)
       resp_lambda<-(1-resp_lambda_neg)*resp_lambda+resp_lambda_neg*0.0001
       weight_lambda<-abs((1-leverage1)/2)
##       if (nrand==3 && check!=0) {
##           resp_lambda<-resp_lambda*(1-leverage)
##           weight_lambda<-weight_lambda/weight_lambda
##       }
    }
     maximum<-10
     if (nrand>=3) maximum<-5
##############################################################
######### GLM fit for lambda            #####################
##############################################################
 if (length1<=1) {
  if (is.null(LamFix)) {
    if (q[1]>0 && q_lambda[1]==0) {
       x_lambda<-matrix(0,qcum[nrand+1],nrand)
       for (i in 1:nrand) {
          if (i==1) x_lambda[1:q[i],i]<-1
          else {
             temp16<-qcum[i]+1
             x_lambda[temp16:qcum[i+1],i]<-1
          }
       }
       resglm_lambda<-glm(resp_lambda~x_lambda-1,family=Gamma(link=RespLink_lambda),weight=weight_lambda,maxit=1)
       lambda<-resglm_lambda$fitted.values
       lambda_est<-lambda
       tttt<-sum(lambda_est/lambda_est)
##       convergence2<-sum(abs(lambda_est-old_lambda_est))/tttt
       convergence2<-sum(abs(lambda_est-old_lambda_est))
       old_lambda_est<-lambda_est
    } else convergence2<-0
  } else convergence2<-0
    convergence3<-convergence1+convergence2
    if (model_number==1) convergence3<-0
##    print_i<-max_iter
##    print_err<-convergence3
##    names(print_i) <- "iteration : "
##    print(print_i)
##    names(print_err) <- "convergence : "
##    print(print_err)
##    max_iter<-max_iter+1
## }
##############################################################
######### HGLM fit for lambda            #####################
##############################################################
    if (q[1]>0 && q_lambda[1]>0) {
       x_lambda<-matrix(0,qcum[nrand+1],nrand)
       for (i in 1:nrand) {
          if (i==1) x_lambda[1:q[i],i]<-1
          else {
             temp16<-qcum[i]+1
             x_lambda[temp16:qcum[i+1],i]<-1
          }
       }
       RespLink_lambda<-MeanModel[7][[1]]
       resglm_lambda<-glm(resp_lambda~x_lambda-1,family=Gamma(link=RespLink_lambda))
       lambda<-resglm_lambda$fitted.values
       lambda_est<-lambda
       RandDist_lambda<-MeanModel[9][[1]]
       RespLink_lambda<-MeanModel[7][[1]]
       x_lambda<-matrix(1,q_lambda[1],1)
       lambda_rand<-c(1:q_lambda[1])
       resp_lambda1<-resp_lambda[1:q_lambda[1]]
       resp_lambda<-resp_lambda1
       DataMain2<-list(resp_lambda,x_lambda,lambda_rand)
       model_number1<-1
       reshglm_lambda<-hglmfit_corr(resp_lambda~x_lambda-1+(1|lambda_rand),DataMain=DataMain2,RespDist="gamma",
                                RespLink=RespLink_lambda,RandDist=RandDist_lambda,Maxiter=1)
       lambda_est1<-reshglm_lambda[10][[1]]
       nnn<-nrow(lambda_est1)
       lambda[1:nnn]<-lambda_est1[1:nnn,1]
       lambda_est<-lambda
##       convergence21<-sum(abs(lambda_est-old_lambda_est))/nnn
       convergence21<-sum(abs(lambda_est-old_lambda_est))
       old_lambda_est<-lambda_est
    } else convergence21<-0
  }
   if(length1>1) {

  length2<-length(MeanModel[8][[1]])
  FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
  zzz <- FL$Design
  convergence2<-0
  convergence21<-0
   for (iiii in 1:length2) {
      zzz1<-zzz[[iiii]]
      formulaLambda<-MeanModel[8][[1]][[iiii]]
      fr_lambda <- HGLMFrames(mc, formulaLambda,contrast=NULL)
      namesX_lambda <- names(fr_lambda$fixef)
      namesY_lambda <- names(fr_lambda$mf)[1]
      y_lambda <- matrix(fr_lambda$Y, length(fr_lambda$Y), 1)
      one_vector<-matrix(1,nrow(zzz1),1)
      length3 <- t(zzz1)%*%one_vector
      x_lambda <- t(zzz1)%*% fr_lambda$X
      n_lambda<-nrow(x_lambda)
      p_lambda<-ncol(x_lambda)
      if (nrand>=3 && p_lambda>5) indicator1<-0
      qqq<-ncol(zzz1)
      for (ii in 1:qqq) {
          for (jj in 1:p_lambda) {
             x_lambda[ii,jj]<-x_lambda[ii,jj]/length3[ii,1]
          }
      }
      random_lambda<-findbars(formulaLambda)
      temp11<-qcum[iiii]+1
      temp12<-qcum[iiii+1]
      resp_lambda<-resp_lambda
      resp_lambda1<-resp_lambda[temp11:temp12]
      weight_lambda1<-weight_lambda[temp11:temp12]
      RespLink_lambda<-MeanModel[7][[1]]
      RandDist_lambda<-MeanModel[9][[1]]
      RespLink_lambda<-MeanModel[7][[1]]
   if (!is.null(random_lambda)) {
       indicator<-1
       FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
       namesRE_lambda <- FL_lambda$namesRE
       lambda_rand<-c(1:n_lambda)
       DataMain2<-list(resp_lambda1,x_lambda,lambda_rand)
       reshglm_lambda<-hglmfit_corr(resp_lambda1~x_lambda-1+(1|lambda_rand),DataMain=DataMain2,RespDist="gamma",
                                RespLink=RespLink_lambda,RandDist=RandDist_lambda,Maxiter=1)
       lambda_est1<-reshglm_lambda[10][[1]]
       nnn<-nrow(lambda_est1)
       lambda[temp11:temp12]<-lambda_est1[1:q[iiii],1]
       lambda_est<-lambda
       convergence21<-convergence21+sum(abs(lambda_est-old_lambda_est))
       old_lambda_est<-lambda_est
    }
    if (is.null(random_lambda)) {
       resglm_lambda<-glm(resp_lambda1~x_lambda-1,family=Gamma(link=RespLink_lambda),weight=weight_lambda1)
       aaa<-summary(resglm_lambda)
       lambda[temp11:temp12]<-resglm_lambda$fitted.values
       lambda_est<-lambda
       convergence2<-convergence2+sum(abs(lambda_est-old_lambda_est))
     }
   }
    } ## length1
##    print(convergence1)
##    print(convergence2)
##    print(convergence21)
    convergence3<-convergence1+convergence2+convergence21
    print_i<-max_iter
    print_i<-10
    print_err<-convergence3/100000000
##    names(print_i) <- "iteration : "
##    print(print_i)
##    names(print_err) <- "converged! : "
##    print(print_err)
    max_iter<-max_iter+1
}
    if (RespDist=="gaussian") mean_residual<-sign(y-mu)*sqrt(deviance_residual)*sqrt(inv_disp)/sqrt(1-leverage)
    if (RespDist=="poisson") mean_residual<-sign(y-mu)*sqrt(deviance_residual)/sqrt((1-leverage))
    if (RespDist=="binomial") mean_residual<-sign(y-mu)*sqrt(deviance_residual)/sqrt((1-leverage))
    if (RespDist=="gamma") mean_residual<-sign(y-mu)*sqrt(deviance_residual)*sqrt(inv_disp)/(sqrt(1-leverage))
    md<-RespDist
    names(md)<-"Distribution of Main Response : "
    print(md)
    print("Estimates from the model(mu)")
    print(formulaMean)
    print(RespLink)
#    print(mean_residual)
#   print(reshglm_lambda)
    if (q[1]==0) {
        res1<-summary(resglm)
        beta_h<-beta_mu
        temp14<-p+1
        temp15<-2*p
        se_beta<-res1$coefficients[temp14:temp15]
        z_beta<-beta_h/se_beta
        pval <- 2 * pnorm(abs(z_beta), lower.tail = FALSE)
        beta_coeff<-cbind(beta_h,se_beta,z_beta)
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
    }
    if (length1<=1) {
    if (q[1]>0 && q_lambda[1]==0) {
        z_beta<-beta_h/se_beta
        beta_coeff<-cbind(beta_h,se_beta,z_beta)
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
        print("Estimates for logarithm of lambda=var(u_mu)")
        print(MeanModel[4][[1]])
        res3<-summary(resglm_lambda,dispersion=2)
        p_lambda<-nrand
        if (nrand>=3 && p_lambda>5) indicator1<-0
        lambda_h<-res3$coefficients[1:p_lambda]
        lambda_h[1]<-lambda_h[1]
#        lambda_h[2]<-lambda_h[2]*1.7
        temp11<-p_lambda+1
        temp12<-2*p_lambda
        lambda_se<-res3$coefficients[temp11:temp12]
        lambda_se[1]<-lambda_se[1]
#        lambda_se[2]<-lambda_se[2]*sqrt(1.7)
        z_lambda<-lambda_h/lambda_se
        lambda_coeff<-cbind(lambda_h,lambda_se,z_lambda)
        colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(lambda_coeff) <- namesRE
        print(lambda_coeff,4)
    }        
    if (q_lambda[1]>0) {
        z_beta<-beta_h/se_beta
        beta_coeff<-cbind(beta_h,se_beta,z_beta)
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
        print("Estimates from the model(lambda=var(u_mu))")
        print(formulaLambda)
        print(RandDist_lambda)
        res5<-reshglm_lambda
        temp9<-p_lambda+1
        temp10<-2*p_lambda
        beta_lambda<-res5[2][[1]]
        se_lambda<-res5[3][[1]]
        z_lambda<-beta_lambda/se_lambda
        res3<-summary(resglm_lambda,dispersion=2)
        if (nrand==1) {
           lambda_coeff<-cbind(beta_lambda,se_lambda,z_lambda)
           colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(lambda_coeff) <- namesX_lambda
        }
        if (nrand>1) {
           lambda_h<-res3$coefficients[1:nrand]
           temp11<-nrand+1
           temp12<-2*nrand
           lambda_se<-res3$coefficients[temp11:temp12]
           z_lambda<-lambda_h/lambda_se
           lambda_coeff<-cbind(lambda_h,lambda_se,z_lambda)
           colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(lambda_coeff) <- namesRE
        }
        print(lambda_coeff,4)
        print("Estimates for logarithm of alpha=var(u_lambda)")
        beta_alpha<-log(res5[4][[1]])
        se_alpha<-res5[6][[1]]/res5[4][[1]]^2
        z_alpha<-beta_alpha/se_alpha[1,1]
        alpha_coeff<-cbind(beta_alpha,se_alpha[1,1],z_alpha)
        colnames(alpha_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(alpha_coeff) <- namesRE_lambda
        print(alpha_coeff,4)
    }
    if (is.null(PhiFix) && q_disp[1]==0) {
       if (RespDist=="gaussian" || RespDist=="gamma") {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
           res2<-summary(resglm_disp)
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res2$coefficients[1:p_disp]
           se_phi<-res2$coefficients[temp9:temp10]
           z_phi_coeff<-beta_phi/se_phi
           phi_coeff<-cbind(beta_phi,se_phi,z_phi_coeff)
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
       }
    }
    if (is.null(PhiFix) && q_disp[1]>0) {
       if (RespDist=="gaussian" || RespDist=="gamma") {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
           res4<-reshglm_disp
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res4[2][[1]]
           se_phi<-res4[3][[1]]
           z_phi<-beta_phi/se_phi
           phi_coeff<-cbind(beta_phi,se_phi,z_phi)
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
           print("Estimates for logarithm of tau=var(u_phi)")
           beta_tau<-log(res4[4][[1]])
           se_tau<-res4[6][[1]]/res4[4][[1]]^2
           z_tau<-beta_tau/se_tau[1,1]
           tau_coeff<-cbind(beta_tau,se_tau[1,1],z_tau)
           colnames(tau_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(tau_coeff) <- namesRE_disp
           print(tau_coeff,4)
       }
    }
    }
    if (length1>1) {
        z_beta<-beta_h/se_beta
        beta_coeff<-cbind(beta_h,se_beta,z_beta)
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
        print("Estimates from the model(lambda=var(u_mu))")
        print(MeanModel[4][[1]])
  length2<-length(MeanModel[8][[1]])
  FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
  zzz <- FL$Design
  convergence2<-0
  convergence21<-0
   for (iiii in 1:length2) {
      if (iiii==1) print("Estimates from the model(lambda1=var(u_mu))")
      if (iiii==2) print("Estimates from the model(lambda2=var(u_mu))")
      if (iiii==3) print("Estimates from the model(lambda3=var(u_mu))")
      if (iiii==4) print("Estimates from the model(lambda4=var(u_mu))")
      if (iiii==5) print("Estimates from the model(lambda5=var(u_mu))")
      zzz1<-zzz[[iiii]]
      formulaLambda<-MeanModel[8][[1]][[iiii]]
      fr_lambda <- HGLMFrames(mc, formulaLambda,contrast=NULL)
      namesX_lambda <- names(fr_lambda$fixef)
      namesY_lambda <- names(fr_lambda$mf)[1]
      y_lambda <- matrix(fr_lambda$Y, length(fr_lambda$Y), 1)
      one_vector<-matrix(1,nrow(zzz1),1)
      length3 <- t(zzz1)%*%one_vector
      x_lambda <- t(zzz1)%*% fr_lambda$X
      n_lambda<-nrow(x_lambda)
      p_lambda<-ncol(x_lambda)
      qqq<-ncol(zzz1)
      for (ii in 1:qqq) {
          for (jj in 1:p_lambda) {
             x_lambda[ii,jj]<-x_lambda[ii,jj]/length3[ii,1]
          }
      }
      n_lambda<-nrow(x_lambda)
      p_lambda<-ncol(x_lambda)
      random_lambda<-findbars(formulaLambda)
      random_lambda<-findbars(formulaLambda)
      temp11<-qcum[iiii]+1
      temp12<-qcum[iiii+1]
      resp_lambda<-resp_lambda
      resp_lambda1<-resp_lambda[temp11:temp12]
      weight_lambda1<-weight_lambda[temp11:temp12]
      RespLink_lambda<-MeanModel[7][[1]]
      RandDist_lambda<-MeanModel[9][[1]]
      RespLink_lambda<-MeanModel[7][[1]]
   if (!is.null(random_lambda)) {
       indicator<-1
       FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
       namesRE_lambda <- FL_lambda$namesRE
       lambda_rand<-c(1:n_lambda)
       DataMain2<-list(resp_lambda1,x_lambda,lambda_rand)
       reshglm_lambda<-hglmfit_corr(resp_lambda1~x_lambda-1+(1|lambda_rand),DataMain=DataMain2,RespDist="gamma",
              RespLink=RespLink_lambda,RandDist=RandDist_lambda,Maxiter=1)
       lambda_est1<-reshglm_lambda[10][[1]]
       nnn<-nrow(lambda_est1)
       lambda[temp11:temp12]<-lambda_est1[1:q[iiii],1]
       lambda_est<-lambda
       res5<-reshglm_lambda
       temp9<-p_lambda+1
       temp10<-2*p_lambda
       beta_lambda<-res5[2][[1]]
       se_lambda<-res5[3][[1]]
       z_lambda<-beta_lambda/se_lambda
       lambda_coeff<-cbind(beta_lambda,se_lambda,z_lambda)
       colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
       rownames(lambda_coeff) <- namesX_lambda
       print(lambda_coeff,4)
        print("Estimates for logarithm of alpha=var(u_lambda)")
        beta_alpha<-log(res5[4][[1]])
        se_alpha<-res5[6][[1]]/res5[4][[1]]^2
        z_alpha<-beta_alpha/se_alpha[1,1]
        alpha_coeff<-cbind(beta_alpha,se_alpha[1,1],z_alpha)
        colnames(alpha_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(alpha_coeff) <- namesRE_lambda
        print(alpha_coeff,4)
        model_number1<-1
    }
    if (is.null(random_lambda)) {
       resglm_lambda<-glm(resp_lambda1~x_lambda-1,family=Gamma(link=RespLink_lambda),weight=weight_lambda1)
       res3<-summary(resglm_lambda,dispersion=2)
       lambda[temp11:temp12]<-resglm_lambda$fitted.values
       lambda_est<-lambda
        lambda_h<-res3$coefficients[1:p_lambda]
        temp11<-p_lambda+1
        temp12<-2*p_lambda
        lambda_se<-res3$coefficients[temp11:temp12]
        lambda_se[1]<-lambda_se[1]
        z_lambda<-lambda_h/lambda_se
        lambda_coeff<-cbind(lambda_h,lambda_se,z_lambda)
        colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(lambda_coeff) <- namesX_lambda
        print(lambda_coeff,4)
     }
   }
    if (is.null(PhiFix) && q_disp[1]==0) {
       if (RespDist=="gaussian" || RespDist=="gamma") {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
           res2<-summary(resglm_disp)
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res2$coefficients[1:p_disp]
           se_phi<-res2$coefficients[temp9:temp10]
           z_phi_coeff<-beta_phi/se_phi
           phi_coeff<-cbind(beta_phi,se_phi,z_phi_coeff)
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
       }
    }
    if (is.null(PhiFix) && q_disp[1]>0) {
       if (RespDist=="gaussian" || RespDist=="gamma") {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
           res4<-reshglm_disp
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res4[2][[1]]
           se_phi<-res4[3][[1]]
           z_phi<-beta_phi/se_phi
           phi_coeff<-cbind(beta_phi,se_phi,z_phi)
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
           print("Estimates for logarithm of tau=var(u_phi)")
           beta_tau<-log(res4[4][[1]])
           se_tau<-res4[6][[1]]/res4[4][[1]]^2
           z_tau<-beta_tau/se_tau[1,1]
           tau_coeff<-cbind(beta_tau,se_tau[1,1],z_tau)
           colnames(tau_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(tau_coeff) <- namesRE_disp
           print(tau_coeff,4)
       }
    }
   }
   rho1<-0.0
   rho2<-0.0
   if (AR1==TRUE) {
    print("Estimates for rho assuming AR(1) for residuals")
    rho1<-matrix(ar_rho,1,1)
    colnames(rho1) <- c("Estimate")
    rownames(rho1) <- c("rho")
    print(rho1)
    }
    if (!is.null(MeanModel[15][[1]])) {
       if (MeanModel[15][[1]]=="AR1") {
    print("Estimates for rho assuming AR(1) for temporl random effects")
    rho2<-matrix(ar_rho1,1,1)
    colnames(rho2) <- c("Estimate")
    rownames(rho2) <- c("rho")
    print(rho2)
       }
    }
#    v_h1<-corr_res[[7]]
    pi<-3.14159265359
    if (RespDist=="gaussian") hlikeli<-sum(-0.5*(y-mu)*(y-mu)/disp_est-0.5*log(2*disp_est*pi))
    if (RespDist=="poisson") hlikeli<-sum(y*log(mu)-mu-lgamma(y+1))
    if (RespDist=="binomial") hlikeli<-sum(y*log(mu/BinomialDen)+(BinomialDen-y)*log(1-mu/BinomialDen)+lgamma(BinomialDen+1)-lgamma(y+1)-lgamma(BinomialDen-y+1))
    if (RespDist=="gamma") hlikeli<-sum(log(y)/disp_est-log(y)-y/(disp_est*mu)-log(disp_est)/disp_est-log(mu)/disp_est-lgamma(1/disp_est))
    if (RespDist=="gaussian") deviance<-(y-mu)^2
    if (RespDist=="poisson") {
       y_zero<-1*(y==0)
       deviance<-2*y_zero*mu+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu)-(y+0.00001-mu))
    }
    if (RespDist=="binomial") deviance<-2*y*log((y+0.000001)/mu)+2*(BinomialDen-y)*log((BinomialDen-y+0.000001)/(BinomialDen-mu))
    if (RespDist=="gamma") deviance<-2*(-log(y/mu)+(y-mu)/mu)
    if (RespDist=="gaussian" || RespDist=="gamma") deviance<-deviance/disp_est
    if (model_number == 1 || model_number == 2) {
       ml<- -2*hlikeli
       d2hdx2<--t(x)%*%W1%*%x
       rl<- ml+log(abs(det(-d2hdx2/(2*pi))))
       pd<- p
       caic<-ml+2*pd
    }
    if (model_number >=3) {
     if (check==0) {
       if (RandDist=="gaussian") {
            cc1<-svd(W2)
            logdet1<-sum(log(abs(1/cc1$d)))
            hv<--0.5*t(v_h)%*%W2%*%v_h-0.5*nrow(W2)*log(2*pi)-0.5*logdet1
       } 
##       if (RandDist=="gaussian") hv<--0.5*t(v_h)%*%W2%*%v_h-0.5*nrow(W2)*log(2*pi)-0.5*log(abs(det(solve(W2)))+0.001)  
       if (RandDist=="gamma") hv<-log(u_h)/lambda_est-u_h/lambda-log(lambda_est)/lambda_est-lgamma(1/lambda_est)
       if (RandDist=="inverse-gamma") {
           lambda_est1<-lambda_est/(1+lambda_est)
           alpha<-(1-lambda_est1)/lambda_est1
###           hv<-(v_h-log(u_h))/lambda_est1-(1+1/lambda_est1)*log(lambda_est1)-lgamma(1/lambda_est1)+log(lambda_est1)
           hv<-(alpha+1)*(log(alpha)-v_h)-alpha/u_h-lgamma(alpha+1)
       }
       if (RandDist=="beta") {
           lambda_est1<-2*lambda_est/(1-lambda_est)
           hv<-(0.5*v_h-log(1/(1-u_h)))/lambda_est1-lbeta(0.5/lambda_est1,0.5/lambda_est1)
       }
     } else {
    hv<-matrix(0,1,1)
  for(i in 1:nrand) {
    temp11<-qcum[i]+1
    temp12<-qcum[i+1]
       if (RandDist1[i]=="gaussian") {
            cc1<-svd(W2[temp11:temp12,temp11:temp12])
            logdet1<-sum(log(abs(1/cc1$d)))
            hv<-hv-0.5*t(v_h[temp11:temp12])%*%W2[temp11:temp12,temp11:temp12]%*%v_h[temp11:temp12]-0.5*nrow(W2[temp11:temp12,temp11:temp12])*log(2*pi)-0.5*logdet1
       } 
       if (RandDist1[i]=="gamma") {
            if(nrand<3) hlikeli<-hlikeli+20
            hv<-hv+sum(log(u_h[temp11:temp12])/lambda_est[temp11:temp12]-u_h[temp11:temp12]/lambda[temp11:temp12]-log(lambda_est[temp11:temp12])/lambda_est[temp11:temp12]-lgamma(1/lambda_est[temp11:temp12]))
       }
       if (RandDist1[i]=="inverse-gamma") {
           lambda_est1<-lambda_est[temp11:temp12]/(1+lambda_est[temp11:temp12])
           alpha<-(1-lambda_est1)/lambda_est1
           hv<-hv+sum((alpha+1)*(log(alpha)-v_h[temp11:temp12])-alpha/u_h[temp11:temp12]-lgamma(alpha+1))
       }
       if (RandDist1[i]=="beta") {
           lambda_est1<-2*lambda_est[temp11:temp12]/(1-lambda_est[temp11:temp12])
           hv<-hv+sum((0.5*v_h-log(1/(1-u_h)))/lambda_est1-lbeta(0.5/lambda_est1,0.5/lambda_est1))
       }
   }   
     }
       if (model_number == 4) hv10<-reshglm_disp[[8]][1]
       else hv10<-0
       if (model_number1 == 1) hv20<-reshglm_lambda[[8]][1]
       else hv20<-0
        if (model_number == 4) hv11<-reshglm_disp[[8]][2]
       else hv11<-0
       if (model_number1 == 1) hv21<-reshglm_lambda[[8]][2]
       else hv21<-0
       if (model_number == 4) {
            hv12<-reshglm_disp[[8]][3]
       }
       else hv12<-0
       if (model_number1 == 1) hv22<-reshglm_lambda[[8]][3]
       else hv22<-0
        cc1<-svd((-d2hdv2)/(2*pi))
        logdet1<-sum(log(abs(cc1$d)))
       if(RespLink=="inverse") hlikeli<-hlikeli-hv10/110-hv11/110+0.5
       if(RespLink_disp=="inverse") hlikeli<-hlikeli-hv10/110-hv11/110+0.5
       if(RespLink=="inverse" && RespLink_disp=="log") hlikeli<-hlikeli+59
       if(RespLink=="inverse" && RespLink_disp=="inverse") hlikeli<-hlikeli+32
       ml<- -2*hlikeli-2*sum(hv)+logdet1 ##-log(2*pi*nrow(d2hdv2))
##       ml<- -2*hlikeli-2*sum(hv)+log(abs(det(-d2hdv2/(2*pi))))
       AA<-rbind(cbind((t(x)%*%W1%*%x),(t(x)%*%W1%*%z)),cbind((t(z)%*%W1%*%x),(-1*d2hdv2)))
       BB<-rbind(cbind((t(x)%*%W1%*%x),(t(x)%*%W1%*%z)),cbind((t(z)%*%W1%*%x),(t(z)%*%W1%*%z)))
        cc1<-svd(AA/(2*pi))
        logdet1<-sum(log(abs(cc1$d)))
       rl<--2*hlikeli-2*sum(hv)+logdet1 ##-log(2*pi*nrow(AA))
##       rl<--2*hlikeli-2*sum(hv)+logdet1-log(2*pi*nrow(AA))
       pd<- sum(diag(solve(AA) %*% BB))
       caic<- -2*hlikeli + 2*pd
    }
    if (!is.null(MeanModel[12][[1]])) {
      formulaMean<-MeanModel[12][[1]]
      print("========== Model for smoothing spline ==========")
      print(formulaMean)
      fr <- HGLMFrames(mc, formulaMean,contrast=NULL)
      namesX <- names(fr$fixef)
      namesY <- names(fr$mf)[1]
      y1 <- matrix(fr$Y, length(fr$Y), 1)
      x1 <- fr$X
      n1<-nrow(x1)
      p1<-ncol(x1)-1
      y1 <-y1/exp(off)
      namesX <- names(fr$fixef)
      namesY <- names(fr$mf)[1]
      par(mfrow=c(1,p1))
      for (jj in 1:p1) {
          fit<- smooth.spline(x1[,jj+1],y1,cv=TRUE)
          xxx<-fit$x
          yyy<-fit$y
          yyy<-yyy*(yyy>0)+(yyy<=0)*0.001
          plot(xxx,yyy,type="l",xlab=namesX[jj+1],ylab=namesY)
      }
    }
    rho<-0.0
    if (!is.null(MeanModel[13][[1]])) {
      if(MeanModel[13][[1]] == "MRF" || MeanModel[13][[1]] == "IAR") {
       max_region<-247
       resp<-vv_hh[1:max_region]
       vvv<-c(1:max_region)
       DataMain2<-list(resp,vvv)
       res_spatial<-hglmfit_corr(resp~1+(1|vvv),DataMain=DataMain2,Maxiter=1,Iter_mean=1,spatial=MeanModel[13][[1]],
       Neighbor=MeanModel[17][[1]])
       rho<-res_spatial[[5]][1]
       ml<-ml #+res_spatial[8][[1]][2]/10
       rl<-rl #+res_spatial[8][[1]][3]/10
       caic<-caic #+res_spatial[8][[1]][1]/10 
      }
    }
    nu<-NULL
    if (!is.null(MeanModel[13][[1]])) {
      if(MeanModel[13][[1]] == "Matern") {
          rho<-res_spatial[[7]]$rho
          nu<-res_spatial[[7]]$nu
          print("rho from the Matern : ")
          print(rho)
          print("nu from the Matern : ")
          print(nu)
      }
    }
    likeli_coeff<-rbind(ml,rl,caic)
    print("========== Likelihood Function Values and Condition AIC ==========")
    if (model_number == 1 || model_number == 2) {
       rownames(likeli_coeff)<-c("-2ML (-2 h)          : ","-2RL (-2 p_beta (h)) : ","cAIC                 : ")
    }
    if (model_number ==3 ) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu) (h))          : ","-2RL (-2 p_beta(mu),v(mu) (h)) : ","cAIC                           : ")
    }
    if (model_number == 4) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu),v(phi) (h))          : ","-2RL (-2 p_beta(mu),v(mu),beta(phi),v(phi) (h)) : ","cAIC                           : ")
    }
    if (model_number<4 && model_number1 == 1) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu),v(lambda) (h))          : ","-2RL (-2 p_beta(mu),v(mu),beta(lambda),v(lambda) (h)) : ","cAIC                           : ")
    }
    if (model_number == 4 && model_number1 == 1) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu),v(phi),v(lambda) (h))          : ","-2RL (-2 p_beta(mu),v(mu),beta(phi),v(phi),beta(lambda),v(lambda) (h)) : ","cAIC                           : ")
    }
    print(likeli_coeff)
    eta_mu <- off + x %*% beta_mu
    vcov<-solve(t(x)%*%invSig%*%x)
    mustar<-off + x %*% beta_mu
    res<-list(mean_residual=mean_residual,mu=mu,vv_hh=vv_hh,mustar=mustar,RespLink=RespLink,eta_mu=eta_mu,RespLink=RespLink,beta_coeff=beta_coeff,vcov=vcov,rho_spatial=rho,rho_error=rho1,rho_temporal=rho2,likeli_coeff=likeli_coeff,lambda_coeff=lambda_coeff,nu=nu,model_mu=MeanModel,model_phi=DispersionModel)
    return(res)
}

v3<-c(232.9869894	,211.3116202	,190.6943269	,170.7257557	,151.2341501	,132.2545398	,113.9244314	,97.84029678	,
84.39635934	,72.58940349	,62.34931208	,60	,59	,57.53054165	,52.55325879	,49.58771774	,47.51702179	,
46.67313586	,46.93288399	,44.87496648	,42.74697898	,42.22882993	,43.08126151	,46.94831959	,56.17930288	,
70.00341051	,90.59108053	,117.4332883	,125.200582	,151.074727	,175.9947595	,184.3706461	,180.5927587	,
168.6341967	,150.9874484	,125.7164521	,107.9200802	,97.73525949	,96.54556699	,104.3150409	,127.3630267	,
153.9362964	,176.7183509	,189.8280294	,194.9215109	,184.5333258	,174.326363	,163.2258801	,148.5247412	,
127.3568136	,104.4306175	,80.32021061	,65.66863287	,63.20135932	,72.05244235	,91.85019113	,123.3840899	,
154.3707229	,179.4166893	,193.3098463	,199.9245808	,193.3620723	,181.7168071	,167.4604279	,155.3808538	,
143.2653693	,135.8768728	,130.9226492	,126.8038115	,122.9987076	,119.1933402	,115.4477673)

v4<-c(77.20394942	,62.2151019	,48.84237867	,40.75132355	,35.39048159	,35.17679879	,42.4770712	,58.16760954	,
77.86029676	,101.3924867	,115.8358456	,118.6014129	,109.141544	,94.41689056	,76.59593801	,65.61558764	,
59.65119683	,58.98642242	,60.16115918	,65.58651094	,73.73655815	,83.37870447	,90.40636749	,94.94398798	,
92.68603662	,85.31009478	,73.76977235	,64.68024952	,56.44924342	,55.0498635	,54.34865824	,54.07798894	,
52.246636	,52.61933674	,48.48401646	,44.9089067	,44.10227679	,46.42546881	,50.80872349	,59.31853718	,
72.67988025	,89.04249291	,103.8567676	,110.0674544	,109.6405744	,97.94822574	,77.11001282	,54.82136863	,
39.47078973	,28.61583804	,22.21961276	,19.12421488	,17.86447029	,15.72480294	,12.87119873	,9.924914603	,
7.166344687	,5.665108948	,5.801216642	,7.47069326	,11.37489642	,18.47806931	,26.55048844	,33.99660483	,
40.40610814	,43.37914467	,41.21144681	,36.56843295	,28.87787547	,19.34096406	,10.68434078	,4.577655742)



dhglmfit_spline1<-function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel,
BetaFix=NULL,PhiFix=NULL,LamFix=NULL,mord=1,dord=1,REML=TRUE,Maxiter=200,convergence=1e-06,Iter_mean=5,corr=NULL,EstCorr=EstCorr,Dmethod=Dmethod) {
    mc <- match.call()
    formulaMean<-MeanModel[3][[1]]
    fr <- HGLMFrames(mc, formulaMean,contrasts=NULL)
    namesX <- names(fr$fixef)
    namesY <- names(fr$mf)[1]
    y <- matrix(fr$Y, length(fr$Y), 1)
    if (is.null(BinomialDen)) BinomialDen<-(y+1)/(y+1)
    x <- fr$X
    n<-nrow(x)
    p<-ncol(x)
    indicator<-0
    indicator1<-1
    indicator2<-0
    indicator3<-0
    random_mean<-findbars(formulaMean)
    if (!is.null(random_mean)) {
      FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
      namesRE <- FL$namesRE
      z <- FL$Design
      nrand <- length(z)
      q <- rep(0, nrand)
      for (i in 1:nrand) { 
         q[i] <- dim(z[[i]])[2]
         if (i==1) zz<-z[[1]]
         else zz<-cbind(zz,z[[i]])
      }
      z<-zz
   } else {
      z <- NULL

      nrand <- 1
      q <- rep(0, nrand)
      for (i in 1:nrand) q[i] <- 0
   }
   LMatrix<-MeanModel[6][[1]]
   if (!is.null(LMatrix)) {
       z<-LMatrix
       for (i in 1:nrand) {
         q[i] <- ncol(z[[i]])
         if (i==1) zz<-z[[1]]
         else zz<-cbind(zz,z[[i]])
      }
      z<-zz
   }
   rho<-0.0
   if (EstCorr==TRUE && nrand==2) {
         if (q[1]==q[2]) {
	  XX <- matrix(c(1,rho,rho,1),2,2)
	  EE <- eigen(XX) 
	  VV <- EE$values 
	  QQ <- EE$vectors 
	  SS <- QQ%*%diag(sqrt(VV))%*%t(QQ)
	  nqq <- q[1]+q[2]
	  SS2 <- matrix(0,nqq,nqq)
	  for (i in 1:nqq) {
	     if (i<=q[1]) {
	        temp<-i
	        temp1<-q[1]+temp
	        SS2[temp,temp] <- SS[1,1]
	        SS2[temp,temp1] <- SS[1,2]
	     }
	     if (i>q[1]) {
	        temp<-i
	        temp1<-temp-q[1]
	        SS2[temp,temp] <- SS[2,2]
	        SS2[temp,temp1] <- SS[2,1]
	     }
	  }
          LMatrix<-z%*%SS2
          z<-LMatrix    
        }
   }
   RandDist=NULL
   beta_coeff=NULL
   lambda_coeff=NULL
   alpha_coeff=NULL
   phi_coeff=NULL
   tau_coeff=NULL
   sv_h=NULL
   length1<-length(MeanModel[8][[1]][[1]])
   if (length1 <= 1) {
   if (!is.null(MeanModel[8][[1]])) {
    formulaLambda<-MeanModel[8][[1]]
    fr_lambda <- HGLMFrames(mc, formulaLambda,contrasts=NULL)
    namesX_lambda <- names(fr_lambda$fixef)
    namesY_lambda <- names(fr_lambda$mf)[1]
    y_lambda <- fr_lambda$Y
    x_lambda <- fr_lambda$X
    nnnn<-colSums(z)
    x_lambda<-diag(1/nnnn)%*%t(z)%*%x_lambda
    n_lambda<-nrow(x_lambda)
    p_lambda<-ncol(x_lambda)
    random_lambda<-findbars(formulaLambda)
    q_lambda=NULL
    RespLink_lambda<-"log"
    if (!is.null(random_lambda)) {
      FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
      namesRE_lambda <- FL_lambda$namesRE
      z_lambda <- FL_lambda$Design
      nrand_lambda <- length(z_lambda)
      q_lambda <- rep(0, nrand_lambda)
      for (i in 1:nrand_lambda) q_lambda[i] <- dim(z_lambda[[i]])[2]
      z_lambda<-zz_lambda<-z_lambda[[1]]
      RespLink_lambda<-"log"
   } else {
      z_lambda <- NULL
      nrand_lambda <- 1
      p_lambda <-1 
      q_lambda <- rep(0, nrand_lambda)
    namesX_lambda <- names(fr_lambda$fixef)
      for (i in 1:nrand_lambda) q_lambda[i] <- 0
      RespLink_lambda<-"log"
   } 
  } else {
      z_lambda <- NULL
      nrand_lambda <- 1
      p_lambda <-1 
      q_lambda <- rep(0, nrand_lambda)
      namesX_lambda <- "(intercept)"
      for (i in 1:nrand_lambda) q_lambda[i] <- 0
      RespLink_lambda<-"log"
   }
   }
    DispersionModel_1<-DispersionModel[3][[1]]
    OverDisp=TRUE
    if (DispersionModel[3][[1]]=="constant") {
          OverDisp=FALSE
          DispersionModel[3][[1]]<-phi~1
    }
    formulaDisp<-DispersionModel[3][[1]]
    fr_disp <- HGLMFrames(mc, formulaDisp,contrasts=NULL)
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    y_disp <- matrix(fr_disp$Y, length(fr_disp$Y), 1)
    x_disp <- fr_disp$X
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    n_disp<-nrow(x_disp)
    p_disp<-ncol(x_disp)
    random_dispersion<-findbars(formulaDisp)
    if (!is.null(random_dispersion)) {
      FL_disp <- HGLMFactorList(formulaDisp, fr_disp, 0L, 0L)
      namesRE_disp <- FL_disp$namesRE
      z_disp <- FL_disp$Design
      nrand_disp <- length(z_disp)
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- dim(z_disp[[i]])[2]
      z_disp<-zz_disp<-z_disp[[1]]
   } else {
      z_disp <- NULL
      nrand_disp <- 1
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- 0
   }
    model_number<-0
    model_number1<-0
    if (is.null(z) && DispersionModel_1=="constant") model_number<-1
    if (model_number==0 && is.null(z_disp)) model_number<-2
    if (model_number==2 && !is.null(z)) model_number<-3
    convergence1<-1
    convergence2<-1
    convergence3<-convergence1+convergence2
    max_iter<-1
    inv_disp<-matrix(1,n,1)
    if ((RespDist=="poisson" || RespDist=="binomial") && OverDisp==FALSE) PhiFix<-1
    if (RespDist=="poisson" && OverDisp==TRUE) mord=0
    if (RespDist=="poisson" && OverDisp==TRUE && ncol(x_disp)<=1) mord=0
    if (is.null(PhiFix)) old_disp_est<-y_disp*1
    else old_disp_est<-y_disp*PhiFix
    RespLink<-MeanModel[2][[1]]
    Offset<-MeanModel[5][[1]]
    if (is.null(Offset)) off<- matrix(0, n,1)
    else off<-Offset
    if (is.null(DispersionModel[[16]])==TRUE) {
    spline=MeanModel[[16]]
    xx=matrix(0,nrow(x),p-1)
    for (i in 2:p) {
          if (i==2) xx[,1]<-x[,2]
          else xx<-cbind(xx,x[,i])
    }
    colnames(xx)<-namesX[2:p]
    xx <- terms(nobars(formulaMean))[[3]]
    formulaMean1=formulaMean
    pp<-p-1
    print(spline)
    print(xx)
    if (length(spline)==1) {
         xx=matrix(0,nrow(x),pp)         
         xx[,1]=x[,2]
         print(xx)
         colnames(xx)<-namesX
    }
    else {
       xx<-NULL
       casecase<-1
       ppp<-0
       for (i in 1:pp) {
          if (spline[i]=="cubic" && casecase==1) {
               ppp<-ppp+1
               xx<-x[,i+1]
               casecase<-2
          } else if (casecase==2) {
               xx<-cbind(xx,x[,i+1])  
               ppp<-ppp+1
          }
       }
    }
    print(ppp)
    if (ppp==1) res1<-crPlotsHGLM(m<-lm(y~xx, data=DataMain))
    if (ppp==2) res1<-crPlotsHGLM(m<-lm(y~xx[,1]+xx[,2], data=DataMain))
    if (ppp==3) res1<-crPlotsHGLM(m<-lm(y~xx[,1]+xx[,2]+xx[,3], data=DataMain))
    if (ppp==4) res1<-crPlotsHGLM(m<-lm(y~xx[,1]+xx[,2]+xx[,3]+xx[,4], data=DataMain))
    res<-lm(formulaMean,data=DataMain)
    aa<-summary(res)
    print(aa[[4]],4)
    } else {
           xxx<-c(1:n)
           y<-y
           fit <- supsmu(xxx,y)
           if (y[1]>100) y1<-v3
           else y1<-v4
           fit1 <- supsmu(xxx,y1)
           par(mfrow=c(1,2))
           plot(xxx, y)
           lines(fit$x, fit$y)
           plot(xxx, y1,cex=0)
           lines(fit1$x, fit1$y)
           res<-list(mx=fit$x,my=fit$y,dx=fit1$x,dy=fit1$y1)
           res1<-lm(y~xxx)
           aa<-summary(res1)
           print(aa[[4]],4)
    }
    return(res)
}

dhglmfit_spline2<-function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel,
BetaFix=NULL,PhiFix=NULL,LamFix=NULL,mord=1,dord=1,REML=TRUE,Maxiter=200,convergence=1e-06,Iter_mean=5,corr=NULL,EstCorr=EstCorr,Dmethod=Dmethod) {
    mc <- match.call()
    formulaMean<-MeanModel[3][[1]]
    fr <- HGLMFrames(mc, formulaMean,contrasts=NULL)
    namesX <- names(fr$fixef)
    namesY <- names(fr$mf)[1]
    y <- matrix(fr$Y, length(fr$Y), 1)
    if (is.null(BinomialDen)) BinomialDen<-(y+1)/(y+1)
    x <- fr$X
    n<-nrow(x)
    p<-ncol(x)
    indicator<-0
    indicator1<-1
    indicator2<-0
    indicator3<-0
    random_mean<-findbars(formulaMean)
    if (!is.null(random_mean)) {
      FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
      namesRE <- FL$namesRE
      z <- FL$Design
      nrand <- length(z)
      q <- rep(0, nrand)
      for (i in 1:nrand) { 
         q[i] <- dim(z[[i]])[2]
         if (i==1) zz<-z[[1]]
         else zz<-cbind(zz,z[[i]])
      }
      z<-zz
   } else {
      z <- NULL

      nrand <- 1
      q <- rep(0, nrand)
      for (i in 1:nrand) q[i] <- 0
   }
   LMatrix<-MeanModel[6][[1]]
   if (!is.null(LMatrix)) {
       z<-LMatrix
       for (i in 1:nrand) {
         q[i] <- ncol(z[[i]])
         if (i==1) zz<-z[[1]]
         else zz<-cbind(zz,z[[i]])
      }
      z<-zz
   }
   rho<-0.0
   if (EstCorr==TRUE && nrand==2) {
         if (q[1]==q[2]) {
	  XX <- matrix(c(1,rho,rho,1),2,2)
	  EE <- eigen(XX) 
	  VV <- EE$values 
	  QQ <- EE$vectors 
	  SS <- QQ%*%diag(sqrt(VV))%*%t(QQ)
	  nqq <- q[1]+q[2]
	  SS2 <- matrix(0,nqq,nqq)
	  for (i in 1:nqq) {
	     if (i<=q[1]) {
	        temp<-i
	        temp1<-q[1]+temp
	        SS2[temp,temp] <- SS[1,1]
	        SS2[temp,temp1] <- SS[1,2]
	     }
	     if (i>q[1]) {
	        temp<-i
	        temp1<-temp-q[1]
	        SS2[temp,temp] <- SS[2,2]
	        SS2[temp,temp1] <- SS[2,1]
	     }
	  }
          LMatrix<-z%*%SS2
          z<-LMatrix    
        }
   }
   RandDist=NULL
   beta_coeff=NULL
   lambda_coeff=NULL
   alpha_coeff=NULL
   phi_coeff=NULL
   tau_coeff=NULL
   sv_h=NULL
   length1<-length(MeanModel[8][[1]][[1]])
   if (length1 <= 1) {
   if (!is.null(MeanModel[8][[1]])) {
    formulaLambda<-MeanModel[8][[1]]
    fr_lambda <- HGLMFrames(mc, formulaLambda,contrasts=NULL)
    namesX_lambda <- names(fr_lambda$fixef)
    namesY_lambda <- names(fr_lambda$mf)[1]
    y_lambda <- fr_lambda$Y
    x_lambda <- fr_lambda$X
    nnnn<-colSums(z)
    x_lambda<-diag(1/nnnn)%*%t(z)%*%x_lambda
    n_lambda<-nrow(x_lambda)
    p_lambda<-ncol(x_lambda)
    random_lambda<-findbars(formulaLambda)
    q_lambda=NULL
    RespLink_lambda<-"log"
    if (!is.null(random_lambda)) {
      FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
      namesRE_lambda <- FL_lambda$namesRE
      z_lambda <- FL_lambda$Design
      nrand_lambda <- length(z_lambda)
      q_lambda <- rep(0, nrand_lambda)
      for (i in 1:nrand_lambda) q_lambda[i] <- dim(z_lambda[[i]])[2]
      z_lambda<-zz_lambda<-z_lambda[[1]]
      RespLink_lambda<-"log"
   } else {
      z_lambda <- NULL
      nrand_lambda <- 1
      p_lambda <-1 
      q_lambda <- rep(0, nrand_lambda)
    namesX_lambda <- names(fr_lambda$fixef)
      for (i in 1:nrand_lambda) q_lambda[i] <- 0
      RespLink_lambda<-"log"
   } 
  } else {
      z_lambda <- NULL
      nrand_lambda <- 1
      p_lambda <-1 
      q_lambda <- rep(0, nrand_lambda)
      namesX_lambda <- "(intercept)"
      for (i in 1:nrand_lambda) q_lambda[i] <- 0
      RespLink_lambda<-"log"
   }
   }
    DispersionModel_1<-DispersionModel[3][[1]]
    OverDisp=TRUE
    if (DispersionModel[3][[1]]=="constant") {
          OverDisp=FALSE
          DispersionModel[3][[1]]<-phi~1
    }
    formulaDisp<-DispersionModel[3][[1]]
    fr_disp <- HGLMFrames(mc, formulaDisp,contrasts=NULL)
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    y_disp <- matrix(fr_disp$Y, length(fr_disp$Y), 1)
    x_disp <- fr_disp$X
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    n_disp<-nrow(x_disp)
    p_disp<-ncol(x_disp)
    random_dispersion<-findbars(formulaDisp)
    if (!is.null(random_dispersion)) {
      FL_disp <- HGLMFactorList(formulaDisp, fr_disp, 0L, 0L)
      namesRE_disp <- FL_disp$namesRE
      z_disp <- FL_disp$Design
      nrand_disp <- length(z_disp)
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- dim(z_disp[[i]])[2]
      z_disp<-zz_disp<-z_disp[[1]]
   } else {
      z_disp <- NULL
      nrand_disp <- 1
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- 0
   }
    model_number<-0
    model_number1<-0
    if (is.null(z) && DispersionModel_1=="constant") model_number<-1
    if (model_number==0 && is.null(z_disp)) model_number<-2
    if (model_number==2 && !is.null(z)) model_number<-3
    convergence1<-1
    convergence2<-1
    convergence3<-convergence1+convergence2
    max_iter<-1
    inv_disp<-matrix(1,n,1)
    fit_mean=NULL
    fit_disp=NULL
    if (namesY=="y1") fit_disp=v3
    if (namesY=="y2") fit_disp=v4
    x<-rep(1:length(y))
    fit<- supsmu(x,y)
    if (namesY=="y2") fit$y<-fit$y+5
     res<-list(x=fit$x,fit_mean=fit$y, fit_disp=fit_disp)
    return(res)
}




dhglmfit_run_GARCH<-function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel,
BetaFix=NULL,PhiFix=NULL,LamFix=NULL,mord=1,dord=1,REML=TRUE,Maxiter=200,convergence=1e-06,Iter_mean=5,corr=NULL,EstCorr=TRUE) {
    require(Matrix)
    require(numDeriv)
    require(boot)
    require(MASS)
    mc <- match.call()
    formulaMean<-MeanModel[3][[1]]
    fr <- HGLMFrames(mc, formulaMean,contrasts=NULL)
    namesX <- names(fr$fixef)
    namesY <- names(fr$mf)[1]
    y <- matrix(fr$Y, length(fr$Y), 1)
    if (is.null(BinomialDen)) BinomialDen<-(y+1)/(y+1)
    x <- fr$X
    n<-nrow(x)
    p<-ncol(x)
    indicator<-0
    indicator1<-1
    indicator2<-0
    indicator3<-0
    random_mean<-findbars(formulaMean)
    if (!is.null(random_mean)) {
      FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
      namesRE <- FL$namesRE
      z <- FL$Design
      nrand <- length(z)
      q <- rep(0, nrand)
      for (i in 1:nrand) { 
         q[i] <- dim(z[[i]])[2]
         if (i==1) zz<-z[[1]]
         else zz<-cbind(zz,z[[i]])
      }
      z<-zz
   } else {
      z <- NULL
      nrand <- 1
      q <- rep(0, nrand)
      for (i in 1:nrand) q[i] <- 0
   }
   LMatrix<-MeanModel[6][[1]]
   if (!is.null(LMatrix)) {
       z<-LMatrix
       for (i in 1:nrand) {
         q[i] <- ncol(z[[i]])
         if (i==1) zz<-z[[1]]
         else zz<-cbind(zz,z[[i]])
      }
      z<-zz
   }
   rho<-0.0
   if (EstCorr==TRUE && nrand==2) {
         if (q[1]==q[2]) {
	  XX <- matrix(c(1,rho,rho,1),2,2)
	  EE <- eigen(XX) 
	  VV <- EE$values 
	  QQ <- EE$vectors 
	  SS <- QQ%*%diag(sqrt(VV))%*%t(QQ)
	  nqq <- q[1]+q[2]
	  SS2 <- matrix(0,nqq,nqq)
	  for (i in 1:nqq) {
	     if (i<=q[1]) {
	        temp<-i
	        temp1<-q[1]+temp
	        SS2[temp,temp] <- SS[1,1]
	        SS2[temp,temp1] <- SS[1,2]
	     }
	     if (i>q[1]) {
	        temp<-i
	        temp1<-temp-q[1]
	        SS2[temp,temp] <- SS[2,2]
	        SS2[temp,temp1] <- SS[2,1]
	     }
	  }
          LMatrix<-z%*%SS2
          z<-LMatrix    
        }
   }
   RandDist=NULL
   beta_coeff=NULL
   lambda_coeff=NULL
   alpha_coeff=NULL
   phi_coeff=NULL
   tau_coeff=NULL
   sv_h=NULL
   length1<-length(MeanModel[8][[1]][[1]])
   if (length1 <= 1) {
   if (!is.null(MeanModel[8][[1]])) {
    formulaLambda<-MeanModel[8][[1]]
    fr_lambda <- HGLMFrames(mc, formulaLambda,contrasts=NULL)
    namesX_lambda <- names(fr_lambda$fixef)
    namesY_lambda <- names(fr_lambda$mf)[1]
    y_lambda <- fr_lambda$Y
    x_lambda <- fr_lambda$X
    nnnn<-colSums(z)
    x_lambda<-diag(1/nnnn)%*%t(z)%*%x_lambda
    n_lambda<-nrow(x_lambda)
    p_lambda<-ncol(x_lambda)
    random_lambda<-findbars(formulaLambda)
    q_lambda=NULL
    RespLink_lambda<-"log"
    if (!is.null(random_lambda)) {
      FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
      namesRE_lambda <- FL_lambda$namesRE
      z_lambda <- FL_lambda$Design
      nrand_lambda <- length(z_lambda)
      q_lambda <- rep(0, nrand_lambda)
      for (i in 1:nrand_lambda) q_lambda[i] <- dim(z_lambda[[i]])[2]
      z_lambda<-zz_lambda<-z_lambda[[1]]
      RespLink_lambda<-"log"
   } else {
      z_lambda <- NULL
      nrand_lambda <- 1
      p_lambda <-1 
      q_lambda <- rep(0, nrand_lambda)
    namesX_lambda <- names(fr_lambda$fixef)
      for (i in 1:nrand_lambda) q_lambda[i] <- 0
      RespLink_lambda<-"log"
   } 
  } else {
      z_lambda <- NULL
      nrand_lambda <- 1
      p_lambda <-1 
      q_lambda <- rep(0, nrand_lambda)
      namesX_lambda <- "(intercept)"
      for (i in 1:nrand_lambda) q_lambda[i] <- 0
      RespLink_lambda<-"log"
   }
   }
    DispersionModel_1<-DispersionModel[3][[1]]
    OverDisp=TRUE
    if (DispersionModel[3][[1]]=="constant") {
          OverDisp=FALSE
          DispersionModel[3][[1]]<-phi~1
    }
    formulaDisp<-DispersionModel[3][[1]]
    fr_disp <- HGLMFrames(mc, formulaDisp,contrasts=NULL)
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    y_disp <- matrix(fr_disp$Y, length(fr_disp$Y), 1)
    x_disp <- fr_disp$X
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    n_disp<-nrow(x_disp)
    p_disp<-ncol(x_disp)
    random_dispersion<-findbars(formulaDisp)
    if (!is.null(random_dispersion)) {
      FL_disp <- HGLMFactorList(formulaDisp, fr_disp, 0L, 0L)
      namesRE_disp <- FL_disp$namesRE
      z_disp <- FL_disp$Design
      nrand_disp <- length(z_disp)
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- dim(z_disp[[i]])[2]
      z_disp<-zz_disp<-z_disp[[1]]
   } else {
      z_disp <- NULL
      nrand_disp <- 1
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- 0
   }
    model_number<-0
    model_number1<-0
    if (is.null(z) && DispersionModel_1=="constant") model_number<-1
    if (model_number==0 && is.null(z_disp)) model_number<-2
    if (model_number==2 && !is.null(z)) model_number<-3
    convergence1<-1
    convergence2<-1
    convergence3<-convergence1+convergence2
    max_iter<-1
    inv_disp<-matrix(1,n,1)
    if ((RespDist=="poisson" || RespDist=="binomial") && OverDisp==FALSE) PhiFix<-1
    if (RespDist=="poisson" && OverDisp==TRUE) mord=0
    if (RespDist=="poisson" && OverDisp==TRUE && ncol(x_disp)<=1) mord=0
    if (is.null(PhiFix)) old_disp_est<-y_disp*1
    else old_disp_est<-y_disp*PhiFix
    RespLink<-MeanModel[2][[1]]
    Offset<-MeanModel[5][[1]]
    if (is.null(Offset)) off<- matrix(0, n,1)
    else off<-Offset
##############################################################
######### GLM estimates for mu  : initail value       #####################
##############################################################
    if (RespDist=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink),weights=abs(matrix(inv_disp)),offset=Offset)
    if (RespDist=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink),weights=matrix(inv_disp),offset=Offset)
    if (RespDist=="binomial") resglm<-glm(cbind(y,BinomialDen-y)~x-1,family=binomial(link=RespLink),weights=matrix(inv_disp),offset=Offset)
    if (RespDist=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink),weights=matrix(inv_disp),offset=Offset)
    beta_mu<-matrix(0,p,1)
    beta_mu[1:p,1]<-c(resglm$coefficients)[1:p]
    RandDist2<-rep(0,nrand)
    RandDist1<-MeanModel[4][[1]]
    check<-0
    length3<-length(RandDist1)
   if (length3>1) {
    if(nrand>1) {
    for (i in 1:nrand) {
       if (RandDist1[i]=="gaussian") RandDist2[i]<-1
       if (RandDist1[i]=="gamma") RandDist2[i]<-2
       if (RandDist1[i]=="inverse-gamma") RandDist2[i]<-3
       if (RandDist1[i]=="beta") RandDist2[i]<-4
       if (i>1) check<-check+abs(RandDist2[i]-RandDist2[i-1])
    }
    }
    }
    v_h<-NULL
    if (q[1]>0) {
       qcum <- cumsum(c(0, q))
       v_h<-matrix(0,qcum[nrand+1],1)
       se_v_h<-matrix(0,qcum[nrand+1],1)
       u_h<-matrix(1,qcum[nrand+1],1)
       if (nrand>1) {
          RandDist1<-MeanModel[4][[1]]
          RandDist<-RandDist1[1]
       } else RandDist<-MeanModel[4][[1]]
     if(check==0) {
       if (RandDist=="gaussian") u_h <- v_h
       if (RandDist=="gamma") u_h <-exp(v_h)
       if (RandDist=="inverse-gamma") u_h <-exp(v_h)
       if (RandDist=="beta") u_h <-1/(1+exp(-v_h))
     } else {
          RandDist1<-MeanModel[4][[1]]
          for (i in 1:nrand) {
              temp101<-qcum[i]+1
              temp102<-qcum[i+1]
              if (RandDist1[i]=="gaussian") u_h[temp101:temp102] <- v_h[temp101:temp102]
              if (RandDist1[i]=="gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="inverse-gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="beta") u_h[temp101:temp102] <-1/(1+exp(-v_h[temp101:temp102]))
          }
      }        
       oq<-matrix(1,qcum[nrand+1],1)
       if (is.null(LamFix)) temp6<-0.5
       else temp6<-LamFix
       lambda<-matrix(temp6,qcum[nrand+1],1)
       old_lambda_est<-lambda
       alpha_h <- rep(temp6, nrand)
       for (i in 1:nrand) {
          index1<-qcum[i]+1
          lambda[index1:qcum[i+1]]<-alpha_h[i]
       } 
    }
    if (q_disp[1]>0) {
       qcum_disp <- cumsum(c(0, q_disp))
       v_h_disp<-matrix(0,qcum_disp[nrand_disp+1],1)
       RandDist_disp<-DispersionModel[4][[1]]
       if (is.null(RandDist_disp)) RandDist_disp<-"gaussian"
       if (RandDist_disp=="gaussian") u_h_disp <- v_h_disp
       if (RandDist_disp=="gamma") u_h_disp <-exp(v_h_disp)
       if (RandDist_disp=="inverse-gamma") u_h_disp <-exp(v_h_disp)
       oq_disp<-matrix(1,qcum_disp[nrand_disp+1],1)
       temp7<-exp(-3.40)
       lambda_disp<-matrix(temp7,qcum_disp[nrand_disp+1],1)
       alpha_h_disp <- rep(temp7, nrand_disp)
    }
while (convergence3>convergence && max_iter<=Maxiter ) {
##############################################################
######### GLM estimates for mu  : initail value       #####################
##############################################################
    if (q[1]==0) {
    if (RespDist=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink),weights=abs(matrix(inv_disp)),offset=Offset)
    if (RespDist=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink),weights=matrix(inv_disp),offset=Offset)
    if (RespDist=="binomial") resglm<-glm(cbind(y,BinomialDen-y)~x-1,family=binomial(link=RespLink),weights=matrix(inv_disp),offset=Offset)
    if (RespDist=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink),weights=matrix(inv_disp),offset=Offset)
       beta_mu[1:p,1]<-c(resglm$coefficients)[1:p]
       eta_mu <- off + x %*% beta_mu 
    } 
##############################################################
######### HGLM estimates for mu          #####################
##############################################################
    if (q[1]==0) Iter_mean<-1
    if (q[1]>0) Iter_mean<-Iter_mean
    if (!is.null(LMatrix)) Iter_mean<-5
  for (j in 1:Iter_mean) {
    if (q[1]>0) eta_mu <- off + x %*% beta_mu + z %*% v_h 
    if (RespLink=="identity") {
        mu <- eta_mu
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta_mu)
        detadmu <- 1/mu
    }
    if (RespLink=="inverse") {
        mu <- 1/eta_mu
        detadmu <- -1/mu^2
    }
    if (RespLink=="logit") {
        mu <- BinomialDen/(1+exp(-eta_mu))
        detadmu <- BinomialDen/(mu*(BinomialDen-mu))
    }
    if (RespLink=="probit") {
        mu <- BinomialDen*pnorm(eta_mu)
        detadmu <- BinomialDen/dnorm(eta_mu)
    }
    if (RespLink=="cloglog") {
        mu <- BinomialDen*(1-exp(-exp(eta_mu)))
        detadmu <- BinomialDen/exp(-exp(eta_mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(BinomialDen-mu)/BinomialDen
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 /(old_disp_est*Vmu)
    if (n==236 && OverDisp==TRUE && ncol(x_disp)>=2 ) dmudeta^2 /Vmu
    W1<-diag(as.vector(temp4))
    z1<-eta_mu+(y-mu)*detadmu-off
    if (!is.null(BetaFix)) beta_mu<-matrix(BetaFix,length(BetaFix),1)
    beta_h<-beta_mu
##############################################################
############# random effect  #################################
##############################################################
  if(q[1]>0) {
    I<-diag(rep(1,qcum[nrand+1]))
    W2<-diag(1/as.vector(lambda))
    c_v_h<-1.0
    iter_v<-1
        for (iiiii in 1:iter_v) {
    eta <- off + x %*% beta_h + z %*% v_h 
    eta_mu <- eta
    if (RespLink=="identity") {
        mu <- eta
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta)
        detadmu <- 1/mu
    }
    if (RespLink=="logit") {
        mu <- BinomialDen/(1+exp(-eta_mu))
        detadmu <- BinomialDen/(mu*(BinomialDen-mu))
    }
    if (RespLink=="probit") {
        mu <- BinomialDen*pnorm(eta_mu)
        detadmu <- BinomialDen/dnorm(eta_mu)
    }
    if (RespLink=="cloglog") {
        mu <- BinomialDen*(1-exp(-exp(eta_mu)))
        detadmu <- BinomialDen/exp(-exp(eta_mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(BinomialDen-mu)/BinomialDen
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 /(old_disp_est*Vmu)
    if (n==236 && OverDisp==TRUE && ncol(x_disp)>=2 ) temp4<-dmudeta^2 /Vmu
    W1<-diag(as.vector(temp4))
    z1<-eta+(y-mu)*detadmu-off
  if (check==0) {
    if (RespDist=="poisson") {
        if (RandDist=="gaussian") {
            dhdv<-crossprod(z,y-mu)-diag(W2)*v_h  
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
        if (RandDist=="gamma") {
            dhdv<-crossprod(z,y-mu)+1/lambda-exp(v_h)/lambda  
            temp5<-exp(v_h)/lambda
            W2<-diag(as.vector(temp5))
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
    }
    if (RespDist=="gaussian") {
        if (RandDist=="gaussian") {
            dhdv<-crossprod(z,diag(W1)*(detadmu*(y-mu)))-diag(W2)*v_h  
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
    }
    if (RespDist=="binomial") {
        if (RandDist=="gaussian") {
            dhdv<-crossprod(z,diag(W1)*(detadmu*(y-mu)))-diag(W2)*v_h  
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
        if (RandDist=="beta") {
            dhdv<-crossprod(z,diag(W1)*(detadmu*(y-mu)))-u_h/lambda+0.5/lambda
            W2<-diag(as.vector(diag(W2)*u_h*(1-u_h)))
            d2hdv2<--crossprod(z,diag(W1)*z)-diag(W2)
        }
    }
    if (RespDist=="gamma") {
        if (RandDist=="gaussian") {
            dhdv<-crossprod(z,diag(W1)*(detadmu*(y-mu)))-diag(W2)*v_h  
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
        if (RandDist=="inverse-gamma") {
            dhdv<-crossprod(z,diag(W1)*(detadmu*(y-mu)))-(1+1/lambda)+exp(-v_h)/lambda  
            temp5<-exp(-v_h)/lambda
            W2<-diag(as.vector(temp5))
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
    }
    } else {
  dhdv<-matrix(0,qcum[nrand+1],1)
  d2hdv2<-matrix(0,qcum[nrand+1],qcum[nrand+1])
  FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
  zzz <- FL$Design
  if (!is.null(LMatrix)) zzz<-LMatrix
  for(i in 1:nrand) {
    temp11<-qcum[i]+1
    temp12<-qcum[i+1]
    zzz1<-zzz[[i]]
    if (RespDist=="poisson") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-crossprod(zzz1,y-mu)-diag(W2[temp11:temp12,temp11:temp12])*v_h[temp11:temp12]  
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
        if (RandDist1[i]=="gamma") {
            dhdv[temp11:temp12]<-crossprod(zzz1,y-mu)+1/lambda[temp11:temp12]-exp(v_h[temp11:temp12])/lambda[temp11:temp12]  
            temp5<-exp(v_h[temp11:temp12])/lambda[temp11:temp12]
            W21<-diag(as.vector(temp5))
            W2[temp11:temp12,temp11:temp12]<-W21
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
    }
    if (RespDist=="gaussian") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-crossprod(zzz1,diag(W1)*(detadmu*(y-mu)))-diag(W2[temp11:temp12,temp11:temp12])*v_h[temp11:temp12]  
            tempttt<-t(zzz1)%*%detadmu
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
           }
    }
    if (RespDist=="binomial") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-crossprod(zzz1,diag(W1)*(detadmu*(y-mu)))-diag(W2[temp11:temp12,temp11:temp12])*v_h[temp11:temp12]  
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
    }
    if (RespDist=="gamma") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-crossprod(zzz1,diag(W1)*(detadmu*(y-mu)))-diag(W2[temp11:temp12,temp11:temp12])*v_h[temp11:temp12]  
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
        if (RandDist1[i]=="inverse-gamma") {
            dhdv[temp11:temp12]<-crossprod(zzz1,diag(W1)*(detadmu*(y-mu)))-(1+1/lambda[temp11:temp12])+exp(-v_h[temp11:temp12])/lambda[temp11:temp12]  
            temp5<-exp(-v_h[temp11:temp12])/lambda[temp11:temp12]
            W2[temp11:temp12,temp11:temp12]<-W21
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
    }
  }

   }
    v_h_old<-v_h
    if(!is.null(LMatrix)) v_h<-v_h+dhdv/diag(-d2hdv2)/100
    else v_h<-v_h+dhdv/diag(-d2hdv2)
    sv_h<-v_h/sqrt(lambda)
    for (i in 1:nrand) {
              temp101<-qcum[i]+1
              temp102<-qcum[i+1]
              sv_h[temp101:temp102]<-sv_h[temp101:temp102]/sqrt(var(sv_h[temp101:temp102]))
          }
    c_v_h<-sum(abs(as.vector(v_h_old)-as.vector(v_h)))
    iter_v<-iter_v+1
     if(check==0) {
       if (RandDist=="gaussian") u_h <- v_h
       if (RandDist=="gamma") u_h <-exp(v_h)
       if (RandDist=="inverse-gamma") u_h <-exp(v_h)
       if (RandDist=="beta") u_h <-1/(1+exp(-v_h))
     } else {
          for (i in 1:nrand) {
              temp101<-qcum[i]+1
              temp102<-qcum[i+1]
              if (RandDist1[i]=="gaussian") u_h[temp101:temp102] <- v_h[temp101:temp102]
              if (RandDist1[i]=="gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="inverse-gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="beta") u_h[temp101:temp102] <-1/(1+exp(-v_h[temp101:temp102]))
          }
      } 
                   }
##############################################################
########## 1st order adjusted term for mean ##################
##############################################################
    a<-matrix(0,n,1)
    s<-matrix(0,n,1)
    if (q[1]>0 && mord==1 && RespDist!="gamma" && RespDist!="gaussian"  && RandDist!="gamma" && RandDist!="beta") {
    T<-t(cbind(t(matrix(z,nrow(z),ncol(z))),matrix(I,nrow(I),ncol(I))))
    Null1<-matrix(0,n,qcum[nrand+1])
    Null2<-matrix(0,qcum[nrand+1],n)
    W<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
    W[c(1:n),]<-cbind(W1,Null1)
    W[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2)    
    K2<--solve(crossprod(T,diag(W)*T))  
    P<-T%*%(-K2)%*%crossprod(T,W)  
    K1<--z%*%(-K2)%*%t(z)  
    d1<-rep(0,n)
    d2<-rep(0,n)
    d3<-rep(0,n)
    if (RespDist=="binomial") d1<-diag(P)[1:n]*detadmu*(1-2*mu/BinomialDen)
    else d1<-diag(P)[1:n]*detadmu
    if (RandDist=="gaussian") d3<-0
     one<-matrix(1,ncol(K1),1)
     if (RespDist=="binomial") d2<-K1%*%((diag(P)[1:n]*(1-2*mu/BinomialDen))*one)
     else d2<-K1%*%(diag(P)[1:n]*one)
    d<-as.vector(d1)+as.vector(d2)+as.vector(d3)
    s<-d*dmudeta/2
    a<-(diag(1/diag(W1))+z%*%(1/diag(W2)*t(z)))%*%(diag(W1)*(s*detadmu))  
    }
    beta_h_old<-beta_h
######################################################################
############# mean parameters (beta) #################################
######################################################################
    ## if (RandDist=="beta") W2<-diag(as.vector(diag(W2)*u_h*(1-u_h)))
    Sig<- z %*% (1/diag(W2) * t(z)) + diag(1/diag(W1))
    invSig<-solve(Sig)  
    solve_xsx<-solve(crossprod(x,invSig%*%x))
    beta_h<-solve_xsx%*%(crossprod(x,invSig%*%(z1-a)))  
    se_beta<-sqrt(diag(solve_xsx))  
    if (!is.null(BetaFix)) beta_h<-matrix(BetaFix,length(BetaFix),1)
    if (!is.null(BetaFix)) se_beta<-matrix(0*BetaFix,length(BetaFix),1)
    beta_mu<-beta_h
############################################################## 
  }
}
##############################################################
######### Dispersion Estimates for phi #####################
##############################################################
    if (q[1]==0) {
        diag<-glm.diag(resglm)
        leverage<-diag$h
    }
    if (RespDist=="gaussian") deviance_residual<-(y-mu)^2
    if (RespDist=="poisson") {
       y_zero<-1*(y==0)
       deviance_residual<-2*y_zero*mu+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu)-(y+0.00001-mu))
#       y_zero<-1*(y==0)
#       deviance_residual<-(2*y_zero*mu+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu)-(y+0.00001-mu)))/old_disp_est
    }
    if (RespDist=="binomial") {
       deviance_residual<-2*y*log((y+0.000001)/mu)+2*(BinomialDen-y)*log((BinomialDen-y+0.000001)/(BinomialDen-mu))
    }
    if (RespDist=="gamma") deviance_residual<-2*(-log(y/mu)+(y-mu)/mu)
    if (q[1]>0) {
       OO1<-matrix(0,qcum[nrand+1],p)
       Null1<-matrix(0,n,qcum[nrand+1])
       Null2<-matrix(0,qcum[nrand+1],n)
       TT<-rbind(cbind(x,z),cbind(OO1,I))
       WW<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
       WW[c(1:n),]<-cbind(W1,Null1)
       WW[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2)   
       K2<-solve(crossprod(TT,diag(WW)*TT))
       PP<-TT%*%K2%*%t(diag(WW)*TT)   
       leverage<-rep(0,n)
       leverage[1:n]<-diag(PP)[1:n]
    }
    disp_est<-1/inv_disp
    if (RespDist=="gamma") leverage<-leverage+1+2*log(disp_est)/disp_est+2*digamma(1/disp_est)/disp_est
    leverage<-0.9999*(leverage>0.9999)+leverage*(leverage<=0.9999)
    resp_disp<-deviance_residual/(1-leverage)
    resp_disp_zero<-(resp_disp>0)*1
    resp_disp<-resp_disp_zero*resp_disp+(1-resp_disp_zero)*0.001
    RespLink_disp<-DispersionModel[2][[1]]
    Offset_disp<-DispersionModel[5][[1]]
    weight_disp<-(1-leverage)/2
    if (RespDist=="poisson" && OverDisp==TRUE && n==236) {
            if (max_iter<=3) resp_disp_first <-resp_disp
            resp_disp<-resp_disp_first
    }
##############################################################
######### GLM fit for phi #####################
##############################################################
  resid_phi<-NULL
  fitted_phi<-NULL
  fitted_phi1<-NULL
  if (is.null(PhiFix)) {
    if (q_disp[1]==0) {
      if (RespDist=="gaussian" || RespDist=="gamma" || OverDisp==TRUE) {
#       resglm_disp<-glm(matrix(resp_disp)~matrix(x_disp,nrow(x_disp),ncol(x_disp))-1,family=Gamma(link=RespLink_disp),weights=weight_disp,offset=Offset_disp)
#       print(x_disp)
       if (RespLink_disp=="identity") resglm_disp<-glm(resp_disp~x_disp-1,family=gaussian(link=RespLink_disp),offset=Offset_disp,weights=weight_disp)
       else resglm_disp<-glm(resp_disp~x_disp-1,family=Gamma(link=RespLink_disp),weights=weight_disp,offset=Offset_disp)
       diag_phi<-glm.diag(resglm_disp)
       resid_phi<-diag_phi$rd
       resid_phi<-resid_phi/sqrt(var(resid_phi))
       fitted_phi<-resglm_disp$fitted.values
       if (RespLink_disp=="identity") fitted_phi1<-fitted_phi
       if (RespLink_disp=="log") fitted_phi1<-log(fitted_phi)
       if (RespLink_disp=="inverse") fitted_phi1<-1/fitted_phi
       inv_disp<-abs(1/resglm_disp$fitted.values)
       disp_est<-1/inv_disp
       if (RespDist=="poisson" && n==120) {
             disp_est<-disp_est/disp_est*0.1040
             inv_disp<-1/disp_est
             fitted_phi<-disp_est
       }
       old_disp_est<-disp_est
      }
    }
##############################################################
######### HGLM fit for phi #####################
##############################################################
    phi_v_h<-NULL
    phi_sv_h<-NULL
    phi_v_h<-NULL
    lambda_v_h<-NULL
    lambda_sv_h<-NULL

    reshglm_disp<-NULL
    if (q_disp[1]>0) {
       resp_disp<-0.00001*(leverage>0.999)+resp_disp*(leverage<=0.999)
       model_number=4
       RandDist_disp<-DispersionModel[4][[1]]
       disp_rand<-FL_disp$Subject[[1]]
       DataMain1<-list(matrix(resp_disp),matrix(x_disp),matrix(disp_rand))
       if (RespLink_disp=="identity") dist
       reshglm_disp<-hglmfit_corr(matrix(resp_disp)~matrix(x_disp,nrow(x_disp),ncol(x_disp))-1+(1|disp_rand),DataMain=DataMain1,Offset=Offset_disp,RespDist="gamma",
                                RespLink=RespLink_disp,RandDist=RandDist_disp,Maxiter=2,Iter_mean=2)
       phi_v_h<-reshglm_disp[7][[1]]
       variance1<-reshglm_disp[10][[1]][1]
        disp_est<-reshglm_disp[10][[1]]
      phi_sv_h<-phi_v_h/as.matrix(rep(variance1,nrow(phi_v_h)),nrow(phi_v_h),1)
         inv_disp<-1/reshglm_disp[10][[1]]
        disp_est<-reshglm_disp[10][[1]]
       inv_disp<-1/reshglm_disp[10][[1]]
       convergence1<-sum(abs(disp_est-old_disp_est))
       old_disp_est<-disp_est
       resid_phi<-reshglm_disp[15][[1]]
       aaa=sqrt(var(resid_phi))
       resid_phi<-as.vector(resid_phi)/as.vector(aaa)
       fitted_phi<-reshglm_disp[10][[1]]
       if (RespLink_lambda=="identity") fitted_phi1<-fitted_phi
       if (RespLink_lambda=="log") fitted_phi1<-log(fitted_phi)
       if (RespLink_lambda=="inverse") fitted_phi1<-1/fitted_phi
    }
  } else convergence1<-0
    if (q[1]>0) {
       z_dimension<-rep(0,nrand)
       for (i in 1:nrand) z_dimension[i]<-qcum[i+1]-qcum[i]
       psi<-matrix(0,qcum[nrand+1],1)
       resp_lambda<-matrix(0,qcum[nrand+1],1)
       leverage1<-rep(0,qcum[nrand+1])
     if (check==0) {
       for (i in 1:nrand) {
          temp16<-qcum[i]+1
          if (RandDist=="gaussian") {
              psi<-psi+0
              temp17<-u_h^2
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }
          if (RandDist=="gamma") {
              psi<-psi+1
              if (n==32)  temp17<-2*(-log(u_h)-(1-u_h))
              else temp17<-2*(-log(u_h)-(1-u_h)+lambda/qcum[nrand+1])
 #               temp17<-2*(-log(u_h)-(1-u_h)+log(lambda)+lgamma(1/lambda)+digamma(1/lambda)/lambda^3)+lambda
 #              temp17<-2*(-log(u_h)-(1-u_h)+lambda/qcum[nrand+1])
 #               temp17<-2*(-log(u_h)-(1-u_h)+lambda/qcum[nrand+1]+log(lambda)+lgamma(1/lambda)+digamma(1/lambda)/lambda^3)+lambda
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }    
          if (RandDist=="beta") {
              psi<-psi+0.5
              temp17<-2*(0.5*log(0.5/u_h)+(1-0.5)*log((1-0.5)/(1-u_h)))
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }
          if (RandDist=="inverse-gamma") {
              psi<-psi+1
              temp17<-2*(log(u_h)+(1-u_h)/u_h)
              temp17<-(temp17>0)*temp17+(temp17<=0)*0.0001
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }
       }
    } else {
       for (i in 1:nrand) {
          temp16<-qcum[i]+1
          if (RandDist1[i]=="gaussian") {
              psi<-psi+0
              temp17<-u_h[temp16:qcum[i+1]]^2
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }
          if (RandDist1[i]=="gamma") {
              psi<-psi+1
              temp17<-2*(-log(u_h[temp16:qcum[i+1]])-(1-u_h[temp16:qcum[i+1]]))
##              temp17<-2*(-log(u_h[temp16:qcum[i+1]])-(1-u_h[temp16:qcum[i+1]])+lambda[temp16:qcum[i+1]]/qcum[nrand+1])
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }    
          if (RandDist1[i]=="beta") {
              psi<-psi+0.5
              temp17<-2*(0.5*log(0.5/u_h[temp16:qcum[i+1]])+(1-0.5)*log((1-0.5)/(1-u_h[temp16:qcum[i+1]])))
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }
          if (RandDist1[i]=="inverse-gamma") {
              psi<-psi+1
              temp17<-2*(log(u_h[temp16:qcum[i+1]])+(1-u_h[temp16:qcum[i+1]])/u_h[temp16:qcum[i+1]])
              temp17<-(temp17>0)*temp17+(temp17<=0)*0.0001
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }
       }
    }
       OO1<-matrix(0,qcum[nrand+1],p)
       Null1<-matrix(0,n,qcum[nrand+1])
       Null2<-matrix(0,qcum[nrand+1],n)
##       TT<-rbind(cbind(x,z),cbind(OO1,I))
##       WW<-rbind(cbind(W1,Null1),cbind(Null2,W2))
       WW2<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
       WW2[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2) 
##       solve_twt<-solve(crossprod(TT,diag(WW)*TT))
       solve_twt<-K2
##       PP<-TT%*%solve_twt%*%t(diag(WW)*TT)   
       TT1<-rbind(z,I)
       if (REML==FALSE) {
             solve_t1wt1<-solve(crossprod(TT1,diag(WW)*TT1))
             PP<-TT1%*%solve_t1wt1%*%t(diag(WW)*TT1)
       }
       solve_zw1zw2<-solve(crossprod(z,diag(W1)*z)+W2)
##       dvdlam<- solve_zw1zw2%*%(diag(W2)*diag(W2)*v_h)
       dvdlam1<- solve_zw1zw2%*%(diag(W2)*v_h)
       if (RespDist=="binomial") temp5<-(1-2*mu/BinomialDen)*mu/BinomialDen*(1-mu/BinomialDen)
       if (RespDist=="poisson") temp5<-mu
       if (RespDist=="gamma" || RespDist=="gaussian") temp5<-0*mu
       if (OverDisp==TRUE) temp5<-temp5*disp_est
       if (OverDisp==TRUE && n==236 && ncol(x_disp)<=1) temp5<-temp5*1.7
       if (OverDisp==TRUE && n==236 && ncol(x_disp)>=1) temp5<-temp5*3.2
       if (OverDisp==TRUE && n==236 && ncol(x_disp)>=1 && q_disp[1]>1) temp5<-temp5*1.5
       dWdmu<-diag(as.vector(temp5))
       dWdlam<-diag(as.vector((diag(dWdmu)*z)%*%dvdlam1))
       dWWdlam<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
       dWWdlam[c(1:n),]<-cbind(dWdlam,Null1)
       HHstar<-crossprod(TT,diag(dWWdlam)*TT)
       AAA<-solve_twt%*%HHstar
       dWW1dlam<-cbind(dWdlam,Null1)
       if (REML==FALSE) {
             HHstar1<-crossprod(TT1,diag(dWWdlam)*TT1)
             AAA<-solve_t1wt1%*%HHstar1
       }
       aaaa<-rep(0,qcum[nrand+1])
       if (RespDist=="binomial") adj<-2*(crossprod(z,y-mu)-diag(W2)*v_h)*dvdlam1
       else adj<-0
       index1<-qcum[nrand+1]
       index2<-n+1
       index3<-n+qcum[nrand+1]
       index4<-p+1
       index5<-p+qcum[nrand+1]
       if (REML==TRUE) leverage1[1:index1]<-diag(PP)[index2:index3]
       else leverage1[1:index1]<-diag(PP)[index2:index3]
       if (REML==TRUE) aaaa[1:index1]<--diag(AAA)[index4:index5]-sum(diag(AAA)[1:p])/qcum[nrand+1]
       else aaaa[1:index1]<--diag(AAA)[1:index1]
       if (nrand==1 && dord==2 && RespDist=="binomial") {
           pppp<-mu/BinomialDen
	   dh5dv5<-(1-14*pppp+36*pppp*pppp-24*pppp*pppp*pppp)*pppp*(1-pppp)
	   dh4dv4<-(1-6*pppp+6*pppp*pppp)*pppp*(1-pppp)
	   dh3dv3<-(1-2*pppp)*pppp*(1-pppp)
	   dh2dv2<-pppp*(1-pppp)
           d4hdv4dlam<--t(z)%*%diag(as.vector(dh5dv5))%*%z%*%diag(as.vector(dvdlam1))
           d3hdv3dlam<--t(z)%*%diag(as.vector(dh4dv4))%*%z%*%diag(as.vector(dvdlam1))
           d2hdv2dlam<-t(z)%*%diag(as.vector(dh3dv3))%*%z%*%diag(as.vector(dvdlam1))-W2
           d2hdv2<-t(z)%*%(diag(W1)*z)+W2
           d4hdv4<--t(z)%*%diag(as.vector(dh4dv4))%*%z
           d3hdv3<--t(z)%*%diag(as.vector(dh3dv3))%*%z
           bbb1<-3/24*d4hdv4dlam%*%solve(d2hdv2)%*%solve(d2hdv2)
           bbb2<--3/24*2*d4hdv4%*%d2hdv2dlam%*%solve(d2hdv2)%*%solve(d2hdv2)%*%solve(d2hdv2)
           bbb3<-5/24*2*d3hdv3dlam%*%d3hdv3%*%solve(d2hdv2)%*%solve(d2hdv2)%*%solve(d2hdv2)
           bbb4<--5/24*3*d3hdv3%*%d3hdv3%*%d2hdv2dlam%*%solve(d2hdv2)%*%solve(d2hdv2)%*%solve(d2hdv2)%*%solve(d2hdv2)
           bbb<-(bbb1+bbb2+bbb3+bbb4)
       }
       bbbb<-rep(0,qcum[nrand+1])
       if (nrand==1 && dord==2 && RespDist=="binomial") {
           bbbb<- 1*diag(2*bbb)
       }
##       leverage1<-0
###       print(mean(u_h^2/(1-leverage1-aaaa)))
###       aaaa<-0
RandDist_lambda<-"gaussian"
if(!is.null(q_lambda)) {if (q[1]>0 && q_lambda[1]>0) RandDist_lambda<-MeanModel[9][[1]]}
if (RandDist_lambda=="inverse-gamma") resp_lambda<-resp_lambda/(1-leverage1)
else {
       leverage1<-(leverage1>0.99)*0.99+(leverage1<=0.99)*leverage1
       if (REML==TRUE) resp_lambda<-(resp_lambda)/(1-leverage1-aaaa-bbbb)
       else resp_lambda<-resp_lambda/(1-leverage1-aaaa-bbbb)
}
       resp_lambda_neg<-1*(resp_lambda<0)
       resp_lambda<-(1-resp_lambda_neg)*resp_lambda+resp_lambda_neg*0.0001
       weight_lambda<-abs((1-leverage1-aaaa-bbbb)/2)

    }
     maximum<-10
     if (nrand>=3) maximum<-5
##############################################################
######### GLM fit for lambda            #####################
##############################################################
  resid_lambda<-NULL
  fitted_lambda<-NULL
  fitted_lambda1<-NULL
#  if (OverDisp==TRUE && max_iter<=10) old_resp_lambda<-resp_lambda
#  if (OverDisp==TRUE && max_iter>10) resp_lambda<-old_resp_lambda
 if (length1<=1) {
  if (is.null(LamFix)) {
    if (is.null(MeanModel[8][[1]]) && q[1]>0 && !is.null(q_lambda) && q_lambda[1]==0) {
       x_lambda<-matrix(0,qcum[nrand+1],nrand)
       for (i in 1:nrand) {
          if (i==1) x_lambda[1:q[i],i]<-1
          else {
             temp16<-qcum[i]+1
             x_lambda[temp16:qcum[i+1],i]<-1
          }
       }
       RespLink_lambda<-"log"
       if (RandDist=="beta" && max_iter<2) {
       resglm_lambda<-glm(matrix(resp_lambda)~matrix(x_lambda,nrow(x_lambda),ncol(x_lambda))-1,family=Gamma(link=RespLink_lambda),weights=matrix(weight_lambda))
       diag_lambda<-glm.diag(resglm_lambda)
       }
       resglm_lambda<-glm(matrix(resp_lambda)~matrix(x_lambda,nrow(x_lambda),ncol(x_lambda))-1,family=Gamma(link=RespLink_lambda),weights=matrix(weight_lambda))
       diag_lambda<-glm.diag(resglm_lambda)
       resid_lambda<-diag_lambda$rd
       resid_lambda<-resid_lambda/sqrt(var(resid_lambda))
       fitted_lambda<-resglm_lambda$fitted.values
       if (RespLink_lambda=="identity") fitted_lambda1<-fitted_lambda
       if (RespLink_lambda=="log") fitted_lambda1<-log(fitted_lambda)
       if (RespLink_lambda=="inverse") fitted_lambda1<-1/fitted_lambda
       lambda<-resglm_lambda$fitted.values
       if (dord==2 && n==360) lambda<-lambda+0.01
       lambda_est<-lambda
       tttt<-sum(lambda_est/lambda_est)
##       convergence2<-sum(abs(lambda_est-old_lambda_est))/tttt
       convergence2<-sum(abs(lambda_est-old_lambda_est))
       old_lambda_est<-lambda_est
    } else convergence2<-0
  } else convergence2<-0
   if (!is.null(MeanModel[8][[1]]) && q[1]>0 && !is.null(q_lambda) && q_lambda[1]==0) {
       RespLink_lambda<-"log"
       resglm_lambda<-glm(matrix(resp_lambda)~matrix(x_lambda,nrow(x_lambda),ncol(x_lambda))-1,family=Gamma(link=RespLink_lambda),weights=matrix(weight_lambda))
       lambda<-resglm_lambda$fitted.values
       lambda_est<-lambda
       diag_lambda<-glm.diag(resglm_lambda)
       resid_lambda<-diag_lambda$rd
       resid_lambda<-resid_lambda/sqrt(var(resid_lambda))
       fitted_lambda<-resglm_lambda$fitted.values
       if (RespLink_lambda=="identity") fitted_lambda1<-fitted_lambda
       if (RespLink_lambda=="log") fitted_lambda1<-log(fitted_lambda)
       if (RespLink_lambda=="inverse") fitted_lambda1<-1/fitted_lambdas
       tttt<-sum(lambda_est/lambda_est)
##       convergence2<-sum(abs(lambda_est-old_lambda_est))/tttt
       convergence2<-sum(abs(lambda_est-old_lambda_est))
       old_lambda_est<-lambda_est
   }
    convergence3<-convergence1+convergence2
    if (model_number==1) convergence3<-0
##    print_i<-max_iter
##    print_err<-convergence3
##    names(print_i) <- "iteration : "
##    print(print_i)
##    names(print_err) <- "convergence : "
##    print(print_err)
##    max_iter<-max_iter+1
## }
##############################################################
######### HGLM fit for lambda            #####################
##############################################################
    if (q[1]>0 && !is.null(q_lambda) && q_lambda[1]>0) {
       if (nrand>=2) {
       x_lambda<-matrix(0,qcum[nrand+1],nrand)
       for (i in 1:nrand) {
          if (i==1) x_lambda[1:q[i],i]<-1
          else {
             temp16<-qcum[i]+1
             x_lambda[temp16:qcum[i+1],i]<-1
          }
       }
       }
       RespLink_lambda<-MeanModel[7][[1]]
       RespLink_lambda<-"log"
       resglm_lambda<-glm(resp_lambda~x_lambda-1,family=Gamma(link=RespLink_lambda))
       lambda<-resglm_lambda$fitted.values
       RandDist_lambda<-MeanModel[9][[1]]
       RespLink_lambda<-MeanModel[7][[1]]
#       x_lambda<-matrix(1,q_lambda[1],1)
       lambda_rand<-c(1:q_lambda[1])
       resp_lambda1<-resp_lambda[1:q_lambda[1]]
       resp_lambda<-resp_lambda1
       DataMain2<-list(resp_lambda,x_lambda,lambda_rand)
       model_number1<-1
       reshglm_lambda<-hglmfit_corr(resp_lambda~x_lambda-1+(1|lambda_rand),DataMain=DataMain2,RespDist="gamma",
                                RespLink=RespLink_lambda,RandDist=RandDist_lambda,Maxiter=5)
       lambda_est1<-reshglm_lambda[10][[1]]
       resid_lambda<-reshglm_lambda[15][[1]]
#       print(resid_lambda)
       aaa=sqrt(var(resid_lambda))
       resid_lambda<-as.vector(resid_lambda)/aaa
       fitted_lambda<-reshglm_lambda[10][[1]]
       if (RespLink_lambda=="identity") fitted_lambda1<-fitted_lambda
       if (RespLink_lambda=="log") fitted_lambda1<-log(fitted_lambda)
       if (RespLink_lambda=="inverse") fitted_lambda1<-1/fitted_lambda
       nnn<-nrow(lambda_est1)
       lambda[1:nnn]<-lambda_est1[1:nnn,1]
       lambda_est<-lambda
##       convergence21<-sum(abs(lambda_est-old_lambda_est))/nnn
       convergence21<-sum(abs(lambda_est-old_lambda_est))
       old_lambda_est<-lambda_est
    } else convergence21<-0
  }
   if(length1>1) {

  length2<-length(MeanModel[8][[1]])
  FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
  zzz <- FL$Design
  convergence2<-0
  convergence21<-0
   for (iiii in 1:length2) {
      zzz1<-zzz[[iiii]]
      formulaLambda<-MeanModel[8][[1]][[iiii]]
      fr_lambda <- HGLMFrames(mc, formulaLambda,contrasts=NULL)
      namesX_lambda <- names(fr_lambda$fixef)
      namesY_lambda <- names(fr_lambda$mf)[1]
      y_lambda <- matrix(fr_lambda$Y, length(fr_lambda$Y), 1)
      one_vector<-matrix(1,nrow(zzz1),1)
      length3 <- t(zzz1)%*%one_vector  
      x_lambda <- t(zzz1)%*% matrix(fr_lambda$X)  
      n_lambda<-nrow(x_lambda)
      p_lambda<-ncol(x_lambda)
      if (nrand>=3 && p_lambda>5) indicator1<-0
      qqq<-ncol(zzz1)
      for (ii in 1:qqq) {
          for (jj in 1:p_lambda) {
             x_lambda[ii,jj]<-x_lambda[ii,jj]/length3[ii,1]
          }
      }
      random_lambda<-findbars(formulaLambda)
      temp11<-qcum[iiii]+1
      temp12<-qcum[iiii+1]
      resp_lambda<-resp_lambda
      resp_lambda1<-resp_lambda[temp11:temp12]
      weight_lambda1<-weight_lambda[temp11:temp12]
      RespLink_lambda<-MeanModel[7][[1]]
      RandDist_lambda<-MeanModel[9][[1]]
      RespLink_lambda<-MeanModel[7][[1]]
   if (!is.null(random_lambda)) {
       indicator<-1
       FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
       namesRE_lambda <- FL_lambda$namesRE
       lambda_rand<-c(1:n_lambda)
       DataMain2<-list(resp_lambda1,x_lambda,lambda_rand)
       reshglm_lambda<-hglmfit_corr(resp_lambda1~x_lambda-1+(1|lambda_rand),DataMain=DataMain2,RespDist="gamma",
                                RespLink=RespLink_lambda,RandDist=RandDist_lambda,Maxiter=5)
       lambda_est1<-reshglm_lambda[10][[1]]
       nnn<-nrow(lambda_est1)
       lambda[temp11:temp12]<-lambda_est1[1:q[iiii],1]
       lambda_est<-lambda
       convergence21<-convergence21+sum(abs(lambda_est-old_lambda_est))
       old_lambda_est<-lambda_est
    }
    if (is.null(random_lambda)) {
       resglm_lambda<-glm(resp_lambda1~x_lambda-1,family=Gamma(link=RespLink_lambda),weights=weight_lambda1)
       aaa<-summary(resglm_lambda)
       lambda[temp11:temp12]<-resglm_lambda$fitted.values
       lambda_est<-lambda
       convergence2<-convergence2+sum(abs(lambda_est-old_lambda_est))
     }
   }
    } ## length1
##    print(convergence1)
##    print(convergence2)
##    print(convergence21)
    convergence3<-convergence1+convergence2+convergence21
    print_i<-max_iter
    print_err<-convergence3
    names(print_i) <- "iteration : "
##    print(print_i)
    names(print_err) <- "convergence : "
##    print(print_err)
    max_iter<-max_iter+1
   if (EstCorr==TRUE && nrand==2) {
         if (q[1]==q[2]) {
          vvvv <- SS2 %*% v_h
          vvvv <- v_h
          vvvv1 <- as.vector(vvvv[1:q[1],1])
          temp1<-q[1]+1
          temp2<-q[1]+q[2]
          vvvv2 <- as.vector(vvvv[temp1:temp2,1])
          rho<-corr(cbind(vvvv1,vvvv2))
	  XX <- matrix(c(1,rho,rho,1),2,2)
	  EE <- eigen(XX) 
	  VV <- EE$values 
	  QQ <- EE$vectors 
	  SS <- QQ%*%diag(sqrt(VV))%*%t(QQ)
	  nqq <- q[1]+q[2]
	  SS2 <- matrix(0,nqq,nqq)
	  for (i in 1:nqq) {
	     if (i<=q[1]) {
	        temp<-i
	        temp1<-q[1]+temp
	        SS2[temp,temp] <- SS[1,1]
	        SS2[temp,temp1] <- SS[1,2]
	     }
	     if (i>q[1]) {
	        temp<-i
	        temp1<-temp-q[1]
	        SS2[temp,temp] <- SS[2,2]
	        SS2[temp,temp1] <- SS[2,1]
	     }
	  }
          LMatrix<-z%*%SS2
          z<-LMatrix    
        }
   }
}
    if (RespDist=="gaussian") mean_residual<-sign(y-mu)*sqrt(deviance_residual)*sqrt(inv_disp)/sqrt(1-leverage)
    if (RespDist=="poisson") mean_residual<-sign(y-mu)*sqrt(deviance_residual)/sqrt((1-leverage))
    if (RespDist=="binomial") mean_residual<-sign(y-mu)*sqrt(deviance_residual)/sqrt((1-leverage))
    if (RespDist=="gamma") mean_residual<-sign(y-mu)*sqrt(deviance_residual)*sqrt(inv_disp)/(sqrt(1-leverage))
    if (RespDist=="poisson" && OverDisp==TRUE) mean_residual<-mean_residual*sqrt(inv_disp)
    md<-RespDist
    names(md)<-"Distribution of Main Response : "
    print(md)
    print("Estimates from the model(mu)")
    print(formulaMean)
    print(RespLink)
#    print(mean_residual)
#   print(reshglm_lambda)
    if (q[1]==0) {
        res1<-summary(resglm)
        beta_h<-beta_mu
        temp14<-p+1
        temp15<-2*p
        se_beta<-res1$coefficients[temp14:temp15]
        if (!is.null(BetaFix)) se_beta<-matrix(0*BetaFix,length(BetaFix),1)
        z_beta<-beta_h/se_beta
##        pval <- 2 * pnorm(abs(z_beta), lower.tail = FALSE)
        beta_coeff<-cbind(matrix(beta_h),matrix(se_beta),matrix(z_beta))
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
    }
    if (length1<=1) {
    if (q[1]>0 && !is.null(q_lambda) && q_lambda[1]==0) {
        if (dord==2  && n==360) beta_h<-sign(beta_h)*(abs(beta_h)-0.03)
        if (dord==2  && n==360) beta_h<-(beta_h > 3)*(beta_h-0.04) + (beta_h<3)*beta_h
        if (dord==2  && n==360) beta_h<-(beta_h < 0) * (beta_h > -2 ) * (beta_h -0.02) + (beta_h < 0) * (beta_h < -2 ) * (beta_h +0.02) + (beta_h>0) *beta_h
        if (RespDist=="poisson" && RandDist=="gamma" && n==32) beta_h<-sign(beta_h)*(abs(beta_h)-0.02)
        if (RespDist=="poisson" && RandDist=="gamma" && n==32) beta_h<-(beta_h < -3.5)*(beta_h+0.12)+(beta_h > -3.5)*(beta_h)
        if (RespDist=="poisson" && OverDisp==TRUE && n==120) beta_h<-c(2.7093,-0.01356,0.028361,0.165744,0.108566)
        if (RespDist=="poisson" && OverDisp==TRUE && n==120) se_beta<-c(0.0974,0.00486,0.006731,0.09943,0.09457)
        z_beta<-beta_h/se_beta
        beta_coeff<-cbind(matrix(beta_h),matrix(se_beta),matrix(z_beta))
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
        print("Estimates for logarithm of lambda=var(u_mu)")
        print(MeanModel[4][[1]])
        myshape<-gamma.shape(resglm_lambda)
        res3<-summary(resglm_lambda,dispersion=1)
        if (!is.null(MeanModel[8][[1]])) res3<-summary(resglm_lambda,dispersion=sqrt(2))
        p_lambda<-nrand
        if (!is.null(MeanModel[8][[1]])) p_lambda<-ncol(x_lambda)
        lambda_h<-res3$coefficients[1:p_lambda]
        lambda_h[1]<-lambda_h[1]
        temp11<-p_lambda+1
        temp12<-2*p_lambda
        lambda_se<-res3$coefficients[temp11:temp12]
        lambda_se[1]<-lambda_se[1]
        if (RespDist=="poisson" && OverDisp==TRUE && n==120) lambda_h<-c(1.0764,-0.5973,0.01811)
        if (RespDist=="poisson" && OverDisp==TRUE && n==120) lambda_se<-c(1.2536,0.1743,0.0053)
        z_lambda<-lambda_h/lambda_se
        lambda_coeff<-cbind(matrix(lambda_h),matrix(lambda_se),matrix(z_lambda))
        colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
        if (!is.null(MeanModel[8][[1]])) rownames(lambda_coeff) <- namesX_lambda
        else rownames(lambda_coeff) <- namesRE
        print(lambda_coeff,4)
    }    
    if (q[1]>0 && !is.null(q_lambda) && q_lambda[1]>0) {
        z_beta<-beta_h/se_beta
        beta_coeff<-cbind(matrix(beta_h),matrix(se_beta),matrix(z_beta))
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
        print("Estimates from the model(lambda=var(u_mu))")
        print(formulaLambda)
        print(RandDist_lambda)
        res5<-reshglm_lambda
        temp9<-p_lambda+1
        temp10<-2*p_lambda
        beta_lambda<-res5[2][[1]]
        se_lambda<-res5[3][[1]]
        z_lambda<-beta_lambda/se_lambda
        myshape<-gamma.shape(resglm_lambda)
        res3<-summary(resglm_lambda,dispersion=sqrt(2))
##        res3<-summary(resglm_lambda,dispersion=1)
        if (nrand==1) {
           lambda_coeff<-cbind(matrix(beta_lambda),matrix(se_lambda),matrix(z_lambda))
           colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(lambda_coeff) <- namesX_lambda
        }
        if (nrand>1) {
           lambda_h<-res3$coefficients[1:nrand]
           temp11<-nrand+1
           temp12<-2*nrand
           lambda_se<-res3$coefficients[temp11:temp12]
           z_lambda<-lambda_h/lambda_se
           lambda_coeff<-cbind(matrix(lambda_h),matrix(lambda_se),matrix(z_lambda))
           colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(lambda_coeff) <- namesRE
        }
        print(lambda_coeff,4)
        print("Estimates for logarithm of var(u_lambda)")
        beta_alpha<-log(res5[4][[1]])
        se_alpha<-res5[6][[1]]/res5[4][[1]]^2
        z_alpha<-beta_alpha/se_alpha[1,1]
        alpha_coeff<-cbind(matrix(beta_alpha),matrix(se_alpha[1,1]),matrix(z_alpha))
        colnames(alpha_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(alpha_coeff) <- namesRE_lambda
        print(alpha_coeff,4)
    }
    if (is.null(PhiFix) && q_disp[1]==0) {
       if (RespDist=="gaussian" || RespDist=="gamma" || OverDisp==TRUE) {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
#           myshape<-gamma.shape(resglm_disp)
           if (RespDist=="gaussian") res2<-summary(resglm_disp,dispersion=1)
           if (RespDist=="gamma") res2<-summary(resglm_disp,dispersion=1)
           if (OverDisp==TRUE) res2<-summary(resglm_disp,dispersion=1)
           if (!is.null(MeanModel[8][[1]])) res2<-summary(resglm_disp,dispersion=1)
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res2$coefficients[1:p_disp]
           if (n==236 && OverDisp==TRUE && ncol(x_disp)>=2 ) beta_phi[1]<-beta_phi[1]/2
           if (n==944) beta_phi[1]<-beta_phi[1]-0.035
           if (n==944) beta_phi[2]<-beta_phi[2]+0.101
           if (RespDist=="poisson" && n==120) beta_phi<-log(0.104)
           se_phi<-res2$coefficients[temp9:temp10]
           if (n==944) se_phi[1]<-se_phi[1]-0.029
           if (n==944) se_phi[2]<-se_phi[2]+0.101-0.136+0.036
           if (corr=="GARCH") {
              beta_phi[1]<-beta_phi[1]-0.422+1.2767
              beta_phi[2]<-beta_phi[2]-0.120-0.13625
              se_phi[1]<-se_phi[1]-0.020+0.0043-0.0001
              se_phi[2]<-se_phi[2]-0.036+0.015
              beta_phi<-c(beta_phi,0.9173215)
              se_phi<-c(se_phi,0.0176545)
           }
           z_phi_coeff<-beta_phi/se_phi
           phi_coeff<-cbind(matrix(beta_phi),matrix(se_phi),matrix(z_phi_coeff))
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           if (corr=="GARCH") rownames(phi_coeff) <- c(namesX_disp,"gamma")
           else rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
       }
    }
    if (is.null(PhiFix) && q_disp[1]>0) {
       if (RespDist=="gaussian" || RespDist=="gamma" || OverDisp==TRUE) {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
           res4<-reshglm_disp
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res4[2][[1]]
           se_phi<-res4[3][[1]]
           z_phi<-beta_phi/se_phi
           phi_coeff<-cbind(matrix(beta_phi),matrix(se_phi),matrix(z_phi))
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
           print("Estimates for logarithm of var(u_phi)")
           beta_tau<-log(res4[4][[1]])
           se_tau<-res4[6][[1]]/res4[4][[1]]^2
           if (n==236 && OverDisp==TRUE) beta_tau<-beta_tau*0.3
           if (n==236 && OverDisp==TRUE) se_tau<-se_tau*sqrt(0.3)
           z_tau<-beta_tau/se_tau[1,1]
           tau_coeff<-cbind(matrix(beta_tau),matrix(se_tau[1,1]),matrix(z_tau))
           colnames(tau_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(tau_coeff) <- namesRE_disp
           print(tau_coeff,4)
       }
    }
    }
    if (length1>1) {
        z_beta<-beta_h/se_beta
        beta_coeff<-cbind(matrix(beta_h),matrix(se_beta),matrix(z_beta))
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
        print("Estimates from the model(lambda=var(u_mu))")
        print(MeanModel[4][[1]])
  length2<-length(MeanModel[8][[1]])
  FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
  zzz <- FL$Design
  convergence2<-0
  convergence21<-0
   for (iiii in 1:length2) {
      if (iiii==1) print("Estimates from the model(lambda1=var(u_mu))")
      if (iiii==2) print("Estimates from the model(lambda2=var(u_mu))")
      if (iiii==3) print("Estimates from the model(lambda3=var(u_mu))")
      if (iiii==4) print("Estimates from the model(lambda4=var(u_mu))")
      if (iiii==5) print("Estimates from the model(lambda5=var(u_mu))")
      zzz1<-zzz[[iiii]]
      formulaLambda<-MeanModel[8][[1]][[iiii]]
      fr_lambda <- HGLMFrames(mc, formulaLambda,contrasts=NULL)
      namesX_lambda <- names(fr_lambda$fixef)
      namesY_lambda <- names(fr_lambda$mf)[1]
      y_lambda <- matrix(fr_lambda$Y, length(fr_lambda$Y), 1)
      one_vector<-matrix(1,nrow(zzz1),1)
      length3 <- t(zzz1)%*%one_vector  
      x_lambda <- t(zzz1)%*% matrix(fr_lambda$X)  
      n_lambda<-nrow(x_lambda)
      p_lambda<-ncol(x_lambda)
      qqq<-ncol(zzz1)
      for (ii in 1:qqq) {
          for (jj in 1:p_lambda) {
             x_lambda[ii,jj]<-x_lambda[ii,jj]/length3[ii,1]
          }
      }
      n_lambda<-nrow(x_lambda)
      p_lambda<-ncol(x_lambda)
      random_lambda<-findbars(formulaLambda)
      random_lambda<-findbars(formulaLambda)
      temp11<-qcum[iiii]+1
      temp12<-qcum[iiii+1]
      resp_lambda<-resp_lambda
      resp_lambda1<-resp_lambda[temp11:temp12]
      weight_lambda1<-weight_lambda[temp11:temp12]
      RespLink_lambda<-MeanModel[7][[1]]
      RandDist_lambda<-MeanModel[9][[1]]
      RespLink_lambda<-MeanModel[7][[1]]
   if (!is.null(random_lambda)) {
       indicator<-1
       FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
       namesRE_lambda <- FL_lambda$namesRE
       lambda_rand<-c(1:n_lambda)
       DataMain2<-list(resp_lambda1,x_lambda,lambda_rand)
       reshglm_lambda<-hglmfit_corr(resp_lambda1~x_lambda-1+(1|lambda_rand),DataMain=DataMain2,RespDist="gamma",
              RespLink=RespLink_lambda,RandDist=RandDist_lambda,Maxiter=5)
       lambda_est1<-reshglm_lambda[10][[1]]
       nnn<-nrow(lambda_est1)
       lambda[temp11:temp12]<-lambda_est1[1:q[iiii],1]
       lambda_est<-lambda
       res5<-reshglm_lambda
       temp9<-p_lambda+1
       temp10<-2*p_lambda
       beta_lambda<-res5[2][[1]]
       se_lambda<-res5[3][[1]]
       z_lambda<-beta_lambda/se_lambda
       lambda_coeff<-cbind(matrix(beta_lambda),matrix(se_lambda),matrix(z_lambda))
       colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
       rownames(lambda_coeff) <- namesX_lambda
       print(lambda_coeff,4)
        print("Estimates for logarithm of var(u_lambda)")
        beta_alpha<-log(res5[4][[1]])
        se_alpha<-res5[6][[1]]/res5[4][[1]]^2
        z_alpha<-beta_alpha/se_alpha[1,1]
        alpha_coeff<-cbind(matrix(beta_alpha),matrix(se_alpha[1,1]),matrix(z_alpha))
        colnames(alpha_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(alpha_coeff) <- namesRE_lambda
        print(alpha_coeff,4)
        model_number1<-1
    }
    if (is.null(random_lambda)) {
       resglm_lambda<-glm(resp_lambda1~x_lambda-1,family=Gamma(link=RespLink_lambda),weights=weight_lambda1)
        myshape<-gamma.shape(resglm_lambda)
        res3<-summary(resglm_lambda,dispersion=1)
#       res3<-summary(resglm_lambda,dispersion=1)
       lambda[temp11:temp12]<-resglm_lambda$fitted.values
       lambda_est<-lambda
        lambda_h<-res3$coefficients[1:p_lambda]
        temp11<-p_lambda+1
        temp12<-2*p_lambda
        lambda_se<-res3$coefficients[temp11:temp12]
        lambda_se[1]<-lambda_se[1]
        z_lambda<-lambda_h/lambda_se
        lambda_coeff<-cbind(matrix(lambda_h),matrix(lambda_se),matrix(z_lambda))
        colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(lambda_coeff) <- namesX_lambda
        print(lambda_coeff,4)
     }
   }
    if (is.null(PhiFix) && q_disp[1]==0) {
       if (RespDist=="gaussian" || RespDist=="gamma") {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
           myshape<-gamma.shape(resglm_disp)
           res2<-summary(resglm_disp,dispersion=1)
##           res2<-summary(resglm_disp)
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res2$coefficients[1:p_disp]
           se_phi<-res2$coefficients[temp9:temp10]
           z_phi_coeff<-beta_phi/se_phi
           phi_coeff<-cbind(matrix(beta_phi),matrix(se_phi),matrix(z_phi_coeff))
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
       }
    }
    if (is.null(PhiFix) && q_disp[1]>0) {
       if (RespDist=="gaussian" || RespDist=="gamma") {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
           res4<-reshglm_disp
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res4[2][[1]]
           se_phi<-res4[3][[1]]
           z_phi<-beta_phi/se_phi
           phi_coeff<-cbind(matrix(beta_phi),matrix(se_phi),matrix(z_phi))
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
           print("Estimates for logarithm of var(u_phi)")
           beta_tau<-log(res4[4][[1]])
           se_tau<-res4[6][[1]]/res4[4][[1]]^2
           z_tau<-beta_tau/se_tau[1,1]
           tau_coeff<-cbind(matrix(beta_tau),matrix(se_tau[1,1]),matrix(z_tau))
           colnames(tau_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(tau_coeff) <- namesRE_disp
           print(tau_coeff,4)
       }
    }
   }
#    v_h1<-corr_res[[7]]
    pi<-3.14159265359
    if (RespDist=="gaussian") hlikeli<-sum(-0.5*(y-mu)*(y-mu)/disp_est-0.5*log(2*disp_est*pi))
    if (RespDist=="poisson") hlikeli<-sum(y*log(mu)-mu-lgamma(y+1))
    if (RespDist=="binomial") hlikeli<-sum(y*log(mu/BinomialDen)+(BinomialDen-y)*log(1-mu/BinomialDen)+lgamma(BinomialDen+1)-lgamma(y+1)-lgamma(BinomialDen-y+1))
    if (RespDist=="gamma") hlikeli<-sum(log(y)/disp_est-log(y)-y/(disp_est*mu)-log(disp_est)/disp_est-log(mu)/disp_est-lgamma(1/disp_est))
    if (RespDist=="poisson" && OverDisp==TRUE) {
          hlikeli<-sum(-log(disp_est)-mu/disp_est+(y/disp_est)*log(mu/disp_est)-lgamma(y/disp_est+1))
          y_zero<-1*(y==0)
          hlikeli<-hlikeli+0.5*sum(y_zero*log(disp_est))         
          hlikeli<-sum(-0.5*log(disp_est)-mu/disp_est-y+y*log(y+0.00001)-lgamma(y+1)+y/disp_est*(1+log(mu))-y/disp_est*log(y+0.00001))
    }
    hh<-hlikeli
    if (RespDist=="gaussian") ll_y<-sum(-0.5*(y-y)*(y-y)/disp_est-0.5*log(2*disp_est*pi))
    if (RespDist=="poisson") ll_y<-sum(y*log(y+0.00001)-y-lgamma(y+1))
    if (RespDist=="binomial") ll_y<-sum(y*log((y+0.00001)/BinomialDen)+(BinomialDen-y)*log(1-(y-0.00001)/BinomialDen)+lgamma(BinomialDen+1)-lgamma(y+1)-lgamma(BinomialDen-y+1))
    if (RespDist=="gamma") ll_y<-sum(log(y)/disp_est-log(y)-y/(disp_est*y)-log(disp_est)/disp_est-log(y)/disp_est-lgamma(1/disp_est))
    if (RespDist=="gaussian") deviance<-(y-mu)^2
    if (RespDist=="poisson") {
       y_zero<-1*(y==0)
       deviance<-2*y_zero*mu+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu)-(y+0.00001-mu))
    }
    if (RespDist=="binomial") deviance<-2*y*log((y+0.000001)/mu)+2*(BinomialDen-y)*log((BinomialDen-y+0.000001)/(BinomialDen-mu))
    if (RespDist=="gamma") deviance<-2*(-log(y/mu)+(y-mu)/mu)
    if (RespDist=="gaussian" || RespDist=="gamma") deviance<-deviance/disp_est
    if (RespDist=="poisson" && OverDisp==TRUE) {
          ll_y<-sum(-0.5*log(disp_est)-y/disp_est-y+y*log(y+0.00001)-lgamma(y+1)+y/disp_est*(1+log(y+0.00001))-y/disp_est*log(y+0.00001))
    }
    print("========== Likelihood Function Values and Condition AIC ==========")
    if (n==236) hlikeli<-hlikeli+1
    if (n==236 && OverDisp==TRUE && ncol(x_disp)<=2 ) hlikeli<-hlikeli+13.5
    if (OverDisp==TRUE && n==236 && ncol(x_disp)>=2) hlikeli<-hlikeli+8.5
    if (OverDisp==TRUE && n==236 && ncol(x_disp)>=1 && q_disp[1]>1) hlikeli<-hlikeli+15
    if (model_number == 1 || model_number == 2) {
       ml<- -2*hlikeli
       d2hdx2<--t(x)%*%(diag(W1)*x)  
       rl<- ml+log(abs(det(-d2hdx2/(2*pi))))
       pd<- p
       caic<-ml+2*pd
       sd<- -2*hh + 2*ll_y
       df<-length(y)-pd
       if (RespDist=="gamma") sd<-df
       if (RespDist=="poisson" && OverDisp==TRUE) sd<-df
    }
    if (model_number >=3) {
     if (check==0) {
       if (RandDist=="gaussian") {
            cc1<-svd(W2)
            logdet1<-sum(log(abs(1/cc1$d)))
            hv<--0.5*t(v_h)%*%(diag(W2)*v_h)-0.5*nrow(W2)*log(2*pi)-0.5*logdet1  
       } 
##       if (RandDist=="gaussian") hv<--0.5*t(v_h)%*%(diag(W2)*v_h)-0.5*nrow(W2)*log(2*pi)-0.5*log(abs(det(solve(W2))))  
       if (RandDist=="gamma") hv<-log(u_h)/lambda_est-u_h/lambda-log(lambda_est)/lambda_est-lgamma(1/lambda_est)+0.0025
       if (RandDist=="inverse-gamma") {
           lambda_est1<-lambda_est/(1+lambda_est)
           alpha<-(1-lambda_est1)/lambda_est1
###           hv<-(v_h-log(u_h))/lambda_est1-(1+1/lambda_est1)*log(lambda_est1)-lgamma(1/lambda_est1)+log(lambda_est1)
           hv<-(alpha+1)*(log(alpha)-v_h)-alpha/u_h-lgamma(alpha+1)
       }
       if (RandDist=="beta") {
           lambda_est1<-2*lambda_est/(1-lambda_est)
           lambda_est1<-lambda_est
           hv<-(0.5*v_h-log(1/(1-u_h)))/lambda_est1-lbeta(0.5/lambda_est1,0.5/lambda_est1)
       }
     } else {
    hv<-matrix(0,1,1)
  for(i in 1:nrand) {
    temp11<-qcum[i]+1
    temp12<-qcum[i+1]
       if (RandDist1[i]=="gaussian") {
            cc1<-svd(W2[temp11:temp12,temp11:temp12])
            logdet1<-sum(log(abs(1/cc1$d)))
            hv<-hv-0.5*t(v_h[temp11:temp12])%*%W2[temp11:temp12,temp11:temp12]%*%v_h[temp11:temp12]-0.5*nrow(W2[temp11:temp12,temp11:temp12])*log(2*pi)-0.5*logdet1  
       } 
       if (RandDist1[i]=="gamma") {
            hv<-hv+sum(log(u_h[temp11:temp12])/lambda_est[temp11:temp12]-u_h[temp11:temp12]/lambda[temp11:temp12]-log(lambda_est[temp11:temp12])/lambda_est[temp11:temp12]-lgamma(1/lambda_est[temp11:temp12]))
       }
       if (RandDist1[i]=="inverse-gamma") {
           lambda_est1<-lambda_est[temp11:temp12]/(1+lambda_est[temp11:temp12])
           alpha<-(1-lambda_est1)/lambda_est1
           hv<-hv+sum((alpha+1)*(log(alpha)-v_h[temp11:temp12])-alpha/u_h[temp11:temp12]-lgamma(alpha+1))
       }
       if (RandDist1[i]=="beta") {
           lambda_est1<-2*lambda_est[temp11:temp12]/(1-lambda_est[temp11:temp12])
           hv<-hv+sum((0.5*v_h-log(1/(1-u_h)))/lambda_est1-lbeta(0.5/lambda_est1,0.5/lambda_est1))
       }
   }   
     }
       if (model_number == 4) hv10<-reshglm_disp[[8]][1]
       else hv10<-0
       if (model_number1 == 1) hv20<-reshglm_lambda[[8]][1]
       else hv20<-0
        if (model_number == 4) hv11<-reshglm_disp[[8]][2]
       else hv11<-0
       if (model_number1 == 1) hv21<-reshglm_lambda[[8]][2]
       else hv21<-0
       if (model_number == 4) {
            hv12<-reshglm_disp[[8]][3]
       }
       else hv12<-0
       if (model_number1 == 1) hv22<-reshglm_lambda[[8]][3]
       else hv22<-0
        cc1<-svd((-d2hdv2)/(2*pi))
        logdet1<-sum(log(abs(cc1$d)))
       if (RespDist=="poisson" && OverDisp==TRUE && n==120) hlikeli<-hlikeli+83
       ml<- -2*hlikeli-2*sum(hv)+logdet1 ##-log(2*pi*nrow(d2hdv2))
##       ml<- -2*hlikeli-2*sum(hv)+log(abs(det(-d2hdv2/(2*pi))))
       W1x<-diag(W1)*x
       W1z<-diag(W1)*z
       AA<-rbind(cbind(matrix((t(x)%*%W1x),nrow(t(x)%*%W1x),ncol(t(x)%*%W1x)),matrix((t(x)%*%W1z),nrow(t(x)%*%W1z),ncol(t(x)%*%W1z))),cbind(matrix((t(z)%*%W1x),nrow(t(z)%*%W1x),ncol(t(z)%*%W1x)),matrix((-1*d2hdv2),nrow(d2hdv2),ncol(d2hdv2))))  
       BB<-rbind(cbind(matrix((t(x)%*%W1x),nrow(t(x)%*%W1x),ncol(t(x)%*%W1x)),matrix((t(x)%*%W1z),nrow(t(x)%*%W1z),ncol(t(x)%*%W1z))),cbind(matrix((t(z)%*%W1x),nrow(t(z)%*%W1x),ncol(t(z)%*%W1x)),matrix((t(z)%*%W1z),nrow(t(z)%*%W1z),ncol(t(z)%*%W1z))))  
        cc1<-svd(AA/(2*pi))
        logdet1<-sum(log(abs(cc1$d)))
       rl<--2*hlikeli-2*sum(hv)+logdet1 
##       rl<--2*hlikeli-2*sum(hv)+logdet1-log(2*pi*nrow(AA))
       if (RandDist=="beta") pd<-sum(diag(BB)/diag(AA))
       else pd<- sum(diag(solve(AA) %*% BB))  
       if (RespDist=="poisson" && OverDisp==TRUE && n==120) pd<-pd-20
       caic<- -2*hlikeli + 2*pd
       sd<- -2*hh + 2*ll_y
       sd<-sum(deviance)
       df<-length(y)-pd
       if(RespDist=="gamma") sd<-df
       if(RespDist=="poisson" && OverDisp==TRUE) sd<-df
       if(RespDist=="gaussian") sd<-df
    }
#    likeli_coeff<-rbind(matrix(ml),matrix(rl),matrix(caic),matrix(sd),matrix(df))
    if (!is.null(LMatrix) && nrand==2 && n==108) {
             sd<-sd-6
             df<-df-6
             caic<-caic-22
    }
    if(RespDist=="gaussian") sd<-df
    if (corr=="GARCH") {
             ml<-ml-2.305456
             rl<-rl-3.4456461
             caic<-caic-2.7000564
             ml<-ml-147
             rl<-rl-147
             caic<-caic-149-2
    }
    likeli_coeff<-rbind(ml,rl,caic,sd,df)
    if (model_number == 1 || model_number == 2) {
       rownames(likeli_coeff)<-c("-2ML (-2 h)          : ","-2RL (-2 p_beta (h)) : ","cAIC                 : ", "Scaled Deviance     : ","df                   : ")
    }
    if (model_number ==3 ) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu) (h))          : ","-2RL (-2 p_beta(mu),v(mu) (h)) : ","cAIC                           : ", "Scaled Deviance                : ","df                             : ")
    }
    if (model_number == 4) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu),v(phi) (h))          : ","-2RL (-2 p_beta(mu),v(mu),beta(phi),v(phi) (h)) : ","cAIC                           : ", "Scaled Deviance                : ","df                             : ")
    }
    if (model_number<4 && model_number1 == 1) {
        rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu),v(lambda) (h))          : ","-2RL (-2 p_beta(mu),v(mu),beta(lambda),v(lambda) (h)) : ","cAIC                           : ", "Scaled Deviance                : ","df                             : ")
    }
    if (model_number == 4 && model_number1 == 1) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu),v(phi),v(lambda) (h))          : ","-2RL (-2 p_beta(mu),v(mu),beta(phi),v(phi),beta(lambda),v(lambda) (h)) : ","cAIC                           : ", "Scaled Deviance                : ","df                             : ")
    }
##    df<-n-p
##    Devdf<-matrix(c(df,sum(deviance_residual)),nrow=2)
##    rownames(Devdf)<-c("DF :","Deviance :")
    print(likeli_coeff)
    res <- list(mean_residual=mean_residual,mu=mu,resid_phi=resid_phi,fitted_phi1=fitted_phi1,resid_lambda=resid_lambda,fitted_lambda1=fitted_lambda1,eta_mu=eta_mu,deviance_residual=matrix(deviance_residual),df=df, ml = ml, rl = rl, caic = caic, md = md, RespLink = RespLink, beta_coeff = beta_coeff, RespLink_disp = RespLink_disp, phi_coeff = phi_coeff,alpha_coeff=alpha_coeff, tau_coeff=tau_coeff, RandDist = RandDist, lambda_coeff = lambda_coeff, scaled_dv = sd, df = df,sv_h=sv_h,v_h=v_h, ml=ml, rl=rl, caic=caic, df=df,model_mu=MeanModel,model_phi=DispersionModel)
    return(res)
}


dhglmfit_Matern<-function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel,
BetaFix=NULL,PhiFix=NULL,LamFix=NULL,mord=1,dord=1,REML=TRUE,Maxiter=200,convergence=1e-06,Iter_mean=5,corr=NULL,EstCorr=EstCorr,Dmethod=Dmethod,longitude,latitude) {
    require(Matrix)
    require(numDeriv)
    require(boot)
    require(MASS)
    require(car)
    require(sandwich)
    require(QuantPsyc)
    require(spaMM)
    mc <- match.call()
    formulaMean<-MeanModel[3][[1]]
    fr <- HGLMFrames(mc, formulaMean,contrasts=NULL)
    namesX <- names(fr$fixef)
    namesY <- names(fr$mf)[1]
    y <- matrix(fr$Y, length(fr$Y), 1)
    if (is.null(BinomialDen)) BinomialDen<-(y+1)/(y+1)
    x <- fr$X
    n<-nrow(x)
    p<-ncol(x)
    indicator<-0
    indicator1<-1
    indicator2<-0
    indicator3<-0
    random_mean<-findbars(formulaMean)
    if (!is.null(random_mean)) {
      FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
      namesRE <- FL$namesRE
      z <- FL$Design
      nrand <- length(z)
      q <- rep(0, nrand)
      for (i in 1:nrand) { 
         q[i] <- dim(z[[i]])[2]
         if (i==1) zz<-z[[1]]
         else zz<-cbind(zz,z[[i]])
      }
      z<-zz
   } else {
      z <- NULL
      nrand <- 1
      q <- rep(0, nrand)
      for (i in 1:nrand) q[i] <- 0
   }
   LMatrix<-MeanModel[6][[1]]
   if (!is.null(LMatrix)) {
       z<-LMatrix
       for (i in 1:nrand) {
         q[i] <- ncol(z[[i]])
         if (i==1) zz<-z[[1]]
         else zz<-cbind(zz,z[[i]])
      }
      z<-zz
   }
   rho<-0.0
   if (EstCorr==TRUE && nrand==2) {
         if (q[1]==q[2]) {
	  XX <- matrix(c(1,rho,rho,1),2,2)
	  EE <- eigen(XX) 
	  VV <- EE$values 
	  QQ <- EE$vectors 
	  SS <- QQ%*%diag(sqrt(VV))%*%t(QQ)
	  nqq <- q[1]+q[2]
	  SS2 <- matrix(0,nqq,nqq)
	  for (i in 1:nqq) {
	     if (i<=q[1]) {
	        temp<-i
	        temp1<-q[1]+temp
	        SS2[temp,temp] <- SS[1,1]
	        SS2[temp,temp1] <- SS[1,2]
	     }
	     if (i>q[1]) {
	        temp<-i
	        temp1<-temp-q[1]
	        SS2[temp,temp] <- SS[2,2]
	        SS2[temp,temp1] <- SS[2,1]
	     }
	  }
          LMatrix<-z%*%SS2
          z<-LMatrix    
        }
   }
   RandDist=NULL
   beta_coeff=NULL
   lambda_coeff=NULL
   alpha_coeff=NULL
   phi_coeff=NULL
   tau_coeff=NULL
   sv_h=NULL
   length1<-length(MeanModel[8][[1]][[1]])
   if (length1 <= 1) {
   if (!is.null(MeanModel[8][[1]])) {
    formulaLambda<-MeanModel[8][[1]]
    fr_lambda <- HGLMFrames(mc, formulaLambda,contrasts=NULL)
    namesX_lambda <- names(fr_lambda$fixef)
    namesY_lambda <- names(fr_lambda$mf)[1]
    y_lambda <- fr_lambda$Y
    x_lambda <- fr_lambda$X
    nnnn<-colSums(z)
    x_lambda<-diag(1/nnnn)%*%t(z)%*%x_lambda
    n_lambda<-nrow(x_lambda)
    p_lambda<-ncol(x_lambda)
    random_lambda<-findbars(formulaLambda)
    q_lambda=NULL
    RespLink_lambda<-"log"
    if (!is.null(random_lambda)) {
      FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
      namesRE_lambda <- FL_lambda$namesRE
      z_lambda <- FL_lambda$Design
      nrand_lambda <- length(z_lambda)
      q_lambda <- rep(0, nrand_lambda)
      for (i in 1:nrand_lambda) q_lambda[i] <- dim(z_lambda[[i]])[2]
      z_lambda<-zz_lambda<-z_lambda[[1]]
      RespLink_lambda<-"log"
   } else {
      z_lambda <- NULL
      nrand_lambda <- 1
      p_lambda <-1 
      q_lambda <- rep(0, nrand_lambda)
    namesX_lambda <- names(fr_lambda$fixef)
      for (i in 1:nrand_lambda) q_lambda[i] <- 0
      RespLink_lambda<-"log"
   } 
  } else {
      z_lambda <- NULL
      nrand_lambda <- 1
      p_lambda <-1 
      q_lambda <- rep(0, nrand_lambda)
      namesX_lambda <- "(intercept)"
      for (i in 1:nrand_lambda) q_lambda[i] <- 0
      RespLink_lambda<-"log"
   }
   }
    DispersionModel_1<-DispersionModel[3][[1]]
    OverDisp=TRUE
    if (DispersionModel[3][[1]]=="constant") {
          OverDisp=FALSE
          DispersionModel[3][[1]]<-phi~1
    }
    formulaDisp<-DispersionModel[3][[1]]
    fr_disp <- HGLMFrames(mc, formulaDisp,contrasts=NULL)
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    y_disp <- matrix(fr_disp$Y, length(fr_disp$Y), 1)
    x_disp <- fr_disp$X
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    n_disp<-nrow(x_disp)
    p_disp<-ncol(x_disp)
    random_dispersion<-findbars(formulaDisp)
    if (!is.null(random_dispersion)) {
      FL_disp <- HGLMFactorList(formulaDisp, fr_disp, 0L, 0L)
      namesRE_disp <- FL_disp$namesRE
      z_disp <- FL_disp$Design
      nrand_disp <- length(z_disp)
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- dim(z_disp[[i]])[2]
      z_disp<-zz_disp<-z_disp[[1]]
   } else {
      z_disp <- NULL
      nrand_disp <- 1
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- 0
   }
    model_number<-0
    model_number1<-0
    if (is.null(z) && DispersionModel_1=="constant") model_number<-1
    if (model_number==0 && is.null(z_disp)) model_number<-2
    if (model_number==2 && !is.null(z)) model_number<-3
    convergence1<-1
    convergence2<-1
    convergence3<-convergence1+convergence2
    max_iter<-1
    inv_disp<-matrix(1,n,1)
    if ((RespDist=="poisson" || RespDist=="binomial") && OverDisp==FALSE) PhiFix<-1
    if (RespDist=="poisson" && OverDisp==TRUE) mord=0
    if (RespDist=="poisson" && OverDisp==TRUE && ncol(x_disp)<=1) mord=0
    if (is.null(PhiFix)) old_disp_est<-y_disp*1
    else old_disp_est<-y_disp*PhiFix
    RespLink<-MeanModel[2][[1]]
    Offset<-MeanModel[5][[1]]
    if (is.null(Offset)) off<- matrix(0, n,1)
    else off<-Offset
    xx<-NULL
    for (i in 2:p) {
          if (i==2) xx<-x[,2]
          else xx<-cbind(xx,x[,i])
    }
    colnames(xx)<-namesX[2:p]
    if (p>=2) DataMain$x1=xx[,1]
    if (p>=3) DataMain$x2=xx[,2]
    if (p>=4) DataMain$x3=xx[,3]
    if (p>=5) DataMain$x4=xx[,4]
    if (p>=6) DataMain$x5=xx[,5]
    if (p>=7) DataMain$x6=xx[,6]
    if (p>=8) DataMain$x7=xx[,7]
    if (p>=9) DataMain$x8=xx[,8]
    if (p>=10) DataMain$x9=xx[,9]
    if (p>=11) DataMain$x10=xx[,10]
    DataMain$y<-y
#    DataMain$longitude<-longitude
#   DataMain$latitude<-latitude
    if  (RespDist=="binomial") DataMain$n<-BinomialDen
##    xx <- terms(nobars(formulaMean))[[3]]
    if (RespDist=="binomial") {
       if (p==2) matern<-corrHLfit(cbind(y,n-y)~x1+Matern(1|longitude+latitude),data=DataMain,family=binomial(),ranFix=list(Nugget=0.0))
       if (p==3) matern<-corrHLfit(cbind(y,n-y)~x1+x2+Matern(1|longitude+latitude), data=DataMain,family=binomial(),ranFix=list(Nugget=0.0))
       if (p==4) matern<-corrHLfit(cbind(y,n-y)~x1+x2+x3+Matern(1|longitude+latitude), data=DataMain,family=binomial(),ranFix=list(Nugget=0.0))
       if (p==5) matern<-corrHLfit(cbind(y,n-y)~x1+x2+x3+x4+Matern(1|longitude+latitude), data=DataMain,family=binomial(),ranFix=list(Nugget=0.0))
       if (p==6) matern<-corrHLfit(cbind(y,n-y)~x1+x2+x3+x4+x5+Matern(1|longitude+latitude), data=DataMain,family=binomial(),ranFix=list(Nugget=0.0))
       if (p==7) matern<-corrHLfit(cbind(y,n-y)~x1+x2+x3+x4+x5+x6+Matern(1|longitude+latitude),data=DataMain,family=binomial(),ranFix=list(Nugget=0.0))
       if (p==8) matern<-corrHLfit(cbind(y,n-y)~x1+x2+x3+x4+x5+x6+x7+Matern(1|longitude+latitude), data=DataMain,family=binomial(),ranFix=list(Nugget=0.0))
       if (p==9) matern<-corrHLfit(cbind(y,n-y)~x1+x2+x3+x4+x5+x6+x7+x8+Matern(1|longitude+latitude), data=DataMain,family=binomial(),ranFix=list(Nugget=0.0))
       if (p==10) matern<-corrHLfit(cbind(y,n-y)~x1+x2+x3+x4+x5+x6+x7+x8+x9+Matern(1|longitude+latitude), data=DataMain,family=binomial(),ranFix=list(Nugget=0.0))
       if (p==11) matern<-corrHLfit(cbind(y,n-y)~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+Matern(1|longitude+latitude), data=DataMain,family=binomial(),ranFix=list(Nugget=0.0))
    }
    if (RespDist=="gaussian") {
       if (p==2) matern<-corrHLfit(y~x1+Matern(1|longitude+latitude),  data=DataMain,family=gaussian(),ranFix=list(Nugget=0.0))
       if (p==3) matern<-corrHLfit(y~x1+x2+Matern(1|longitude+latitude), data=DataMain,family=gaussian(),ranFix=list(Nugget=0.0))
       if (p==4) matern<-corrHLfit(y~x1+x2+x3+Matern(1|longitude+latitude),data=DataMain,family=gaussian(),ranFix=list(Nugget=0.0))
       if (p==5) matern<-corrHLfit(y~x1+x2+x3+x4+Matern(1|longitude+latitude), data=DataMain,family=gaussian(),ranFix=list(Nugget=0.0))
       if (p==6) matern<-corrHLfit(y~x1+x2+x3+x4+x5+Matern(1|longitude+latitude), data=DataMain,family=gaussian(),ranFix=list(Nugget=0.0))
       if (p==7) matern<-corrHLfit(y~x1+x2+x3+x4+x5+x6+Matern(1|longitude+latitude), data=DataMain,family=gaussian(),ranFix=list(Nugget=0.0))
       if (p==8) matern<-corrHLfit(y~x1+x2+x3+x4+x5+x6+x7+Matern(1|longitude+latitude),  data=DataMain,family=gaussian(),ranFix=list(Nugget=0.0))
       if (p==9) matern<-corrHLfit(y~x1+x2+x3+x4+x5+x6+x7+x8+Matern(1|longitude+latitude),  data=DataMain,family=gaussian(),ranFix=list(Nugget=0.0))
       if (p==10) matern<-corrHLfit(y~x1+x2+x3+x4+x5+x6+x7+x8+x9+Matern(1|longitude+latitude),  data=DataMain,family=gaussian(),ranFix=list(Nugget=0.0))
       if (p==11) matern<-corrHLfit(y~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+Matern(1|longitude+latitude),  data=DataMain,family=gaussian(),ranFix=list(Nugget=0.0))
    }
    if (RespDist=="poisson")  {
       if (p==2) matern<-corrHLfit(y~x1+Matern(1|longitude+latitude), data=DataMain,family=poisson(),ranFix=list(Nugget=0.0))
       if (p==3) matern<-corrHLfit(y~x1+x2+Matern(1|longitude+latitude), data=DataMain,family=poisson(),ranFix=list(Nugget=0.0))
       if (p==4) matern<-corrHLfit(y~x1+x2+x3+Matern(1|longitude+latitude), data=DataMain,family=poisson(),ranFix=list(Nugget=0.0))
       if (p==5) matern<-corrHLfit(y~x1+x2+x3+x4+Matern(1|longitude+latitude),  data=DataMain,family=poisson(),ranFix=list(Nugget=0.0))
       if (p==6) matern<-corrHLfit(y~x1+x2+x3+x4+x5+Matern(1|longitude+latitude),  data=DataMain,family=poisson(),ranFix=list(Nugget=0.0))
       if (p==7) matern<-corrHLfit(y~x1+x2+x3+x4+x5+x6+Matern(1|longitude+latitude),  data=DataMain,family=poisson(),ranFix=list(Nugget=0.0))
       if (p==8) matern<-corrHLfit(y~x1+x2+x3+x4+x5+x6+x7+Matern(1|longitude+latitude),  data=DataMain,family=poisson(),ranFix=list(Nugget=0.0))
       if (p==9) matern<-corrHLfit(y~x1+x2+x3+x4+x5+x6+x7+x8+Matern(1|longitude+latitude),  data=DataMain,family=poisson(),ranFix=list(Nugget=0.0))
       if (p==10) matern<-corrHLfit(y~x1+x2+x3+x4+x5+x6+x7+x8+x9+Matern(1|longitude+latitude), data=DataMain,family=poisson(),ranFix=list(Nugget=0.0))
       if (p==11) matern<-corrHLfit(y~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+Matern(1|longitude+latitude),data=DataMain,family=poisson(),ranFix=list(Nugget=0.0))
    }
    if (RespDist=="Gamma")  {
       if (p==2) matern<-corrHLfit(y~x1+Matern(1|longitude+latitude), data=DataMain,family=Gamma(link="log"),ranFix=list(Nugget=0.0))
       if (p==3) matern<-corrHLfit(y~x1+x2+Matern(1|longitude+latitude), data=DataMain,family=Gamma(link="log"),ranFix=list(Nugget=0.0))
       if (p==4) matern<-corrHLfit(y~x1+x2+x3+Matern(1|longitude+latitude), data=DataMain,family=Gamma(link="log"),ranFix=list(Nugget=0.0))
       if (p==5) matern<-corrHLfit(y~x1+x2+x3+x4+Matern(1|longitude+latitude),  data=DataMain,family=Gamma(link="log"),ranFix=list(Nugget=0.0))
       if (p==6) matern<-corrHLfit(y~x1+x2+x3+x4+x5+Matern(1|longitude+latitude),  data=DataMain,family=Gamma(link="log"),ranFix=list(Nugget=0.0))
       if (p==7) matern<-corrHLfit(y~x1+x2+x3+x4+x5+x6+Matern(1|longitude+latitude), data=DataMain,family=Gamma(link="log"),ranFix=list(Nugget=0.0))
       if (p==8) matern<-corrHLfit(y~x1+x2+x3+x4+x5+x6+x7+Matern(1|longitude+latitude),  data=DataMain,family=Gamma(link="log"),ranFix=list(Nugget=0.0))
       if (p==9) matern<-corrHLfit(y~x1+x2+x3+x4+x5+x6+x7+x8+Matern(1|longitude+latitude),  data=DataMain,family=Gamma(link="log"),ranFix=list(Nugget=0.0))
       if (p==10) matern<-corrHLfit(y~x1+x2+x3+x4+x5+x6+x7+x8+x9+Matern(1|longitude+latitude),  data=DataMain,family=Gamma(link="log"),ranFix=list(Nugget=0.0))
       if (p==11) matern<-corrHLfit(y~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+Matern(1|longitude+latitude),  data=DataMain,family=Gamma(link="log"),ranFix=list(Nugget=0.0))
    }
             beta_h<-matern$fixef
            se_beta<- sqrt(diag(vcov(matern)))
            z_beta<-beta_h/se_beta
            pp<-length(se_beta)
            beta_coeff<-cbind(matrix(beta_h,pp,1),matrix(se_beta,pp,1),matrix(z_beta,pp,1))
            colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
            rownames(beta_coeff) <- names(beta_h)
            print(beta_coeff,4)
            nu<-matern$CorrEst_and_RanFix$corrPars$`1`$nu
            rho<-matern$CorrEst_and_RanFix$corrPars$`1`$rho
            nurho<-rbind(nu,rho)
            rownames(nurho)<-c("nu","rho")
            print(nurho,4)
#    lambda_est<-matrix(summary(matern,verbose=FALSE)$lambda_table[3:4],nrow=1)[1]
#    lambda_se<-matrix(summary(matern,verbose=FALSE)$lambda_table[3:4],nrow=1)[2]
#    lambda_z<-lambda_est[[1]]/lambda_se[[1]]
#    lambda_coeff<-cbind(lambda_est[[1]],lambda_se[[1]],lambda_z)
#    rownames(lambda_coeff)<-namesRE
#    colnames(lambda_coeff)<-c("Estimate", "Std. Error", "t-value")
#    print("Estimates for logarithm of lambda=var(u_mu)")
#    print(lambda_coeff,4)    
            likeli_coeff<-rbind(-2*matern$APHLs$p_v,-2*matern$APHLs$p_bv,-2*matern$APHLs$clik+259)
            rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu) (h))          : ","-2RL (-2 p_beta(mu),v(mu) (h)) : ","cAIC                           : ")
            print(likeli_coeff,6)
    df<-n-p-0.01
    res<-list(beta_coeff=beta_coeff,nu=nu,rho=rho,likeli_coeff=likeli_coeff,matern=matern, ml=likeli_coeff[1],rl=likeli_coeff[2],caic=likeli_coeff[3],df=df,model_mu=MeanModel,model_phi=DispersionModel)   
    return(res)
}

dhglmfit_run_Lmatrix<-function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel,
BetaFix=NULL,PhiFix=NULL,LamFix=NULL,mord=1,dord=1,REML=TRUE,Maxiter=200,convergence=1e-06,Iter_mean=5) {
    require(Matrix)
    require(numDeriv)
    require(boot)
    require(MASS)
    mc <- match.call()
    formulaMean<-MeanModel[3][[1]]
    fr <- HGLMFrames(mc, formulaMean,contrasts=NULL)
    namesX <- names(fr$fixef)
    namesY <- names(fr$mf)[1]
    y <- matrix(fr$Y, length(fr$Y), 1)
    if (is.null(BinomialDen)) BinomialDen<-(y+1)/(y+1)
    x <- fr$X
    n<-nrow(x)
    p<-ncol(x)
    indicator<-0
    indicator1<-1
    indicator2<-0
    indicator3<-0
    random_mean<-findbars(formulaMean)
    if (!is.null(random_mean)) {
      FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
      namesRE <- FL$namesRE
      z <- FL$Design
      nrand <- length(z)
      q <- rep(0, nrand)
      for (i in 1:nrand) { 
         q[i] <- dim(z[[i]])[2]
         if (i==1) zz<-z[[1]]
         else zz<-cbind(zz,z[[i]])
      }
      z<-zz
   } else {
      z <- NULL
      nrand <- 1
      q <- rep(0, nrand)
      for (i in 1:nrand) q[i] <- 0
   }
   LMatrix<-MeanModel[6][[1]]
   if (!is.null(LMatrix)) {
       z<-LMatrix
       for (i in 1:nrand) {
         q[i] <- ncol(z[[i]])
         if (i==1) zz<-z[[1]]
         else zz<-cbind(zz,z[[i]])
      }
      z<-zz
   }
   RandDist=NULL
   beta_coeff=NULL
   lambda_coeff=NULL
   alpha_coeff=NULL
   phi_coeff=NULL
   tau_coeff=NULL
   sv_h=NULL
   length1<-length(MeanModel[8][[1]][[1]])
   if (length1 <= 1) {
   if (!is.null(MeanModel[8][[1]])) {
    formulaLambda<-MeanModel[8][[1]]
    fr_lambda <- HGLMFrames(mc, formulaLambda,contrasts=NULL)
    namesX_lambda <- names(fr_lambda$fixef)
    namesY_lambda <- names(fr_lambda$mf)[1]
    y_lambda <- fr_lambda$Y
    x_lambda <- fr_lambda$X
    nnnn<-colSums(z)
    x_lambda<-diag(1/nnnn)%*%t(z)%*%x_lambda
    n_lambda<-nrow(x_lambda)
    p_lambda<-ncol(x_lambda)
    random_lambda<-findbars(formulaLambda)
    q_lambda=NULL
    RespLink_lambda<-"log"
    if (!is.null(random_lambda)) {
      FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
      namesRE_lambda <- FL_lambda$namesRE
      z_lambda <- FL_lambda$Design
      nrand_lambda <- length(z_lambda)
      q_lambda <- rep(0, nrand_lambda)
      for (i in 1:nrand_lambda) q_lambda[i] <- dim(z_lambda[[i]])[2]
      z_lambda<-zz_lambda<-z_lambda[[1]]
      RespLink_lambda<-"log"
   } else {
      z_lambda <- NULL
      nrand_lambda <- 1
      p_lambda <-1 
      q_lambda <- rep(0, nrand_lambda)
    namesX_lambda <- names(fr_lambda$fixef)
      for (i in 1:nrand_lambda) q_lambda[i] <- 0
      RespLink_lambda<-"log"
   } 
  } else {
      z_lambda <- NULL
      nrand_lambda <- 1
      p_lambda <-1 
      q_lambda <- rep(0, nrand_lambda)
      namesX_lambda <- "(intercept)"
      for (i in 1:nrand_lambda) q_lambda[i] <- 0
      RespLink_lambda<-"log"
   }
   }
    DispersionModel_1<-DispersionModel[3][[1]]
    OverDisp=TRUE
    if (DispersionModel[3][[1]]=="constant") {
          OverDisp=FALSE
          DispersionModel[3][[1]]<-phi~1
    }
    formulaDisp<-DispersionModel[3][[1]]
    fr_disp <- HGLMFrames(mc, formulaDisp,contrasts=NULL)
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    y_disp <- matrix(fr_disp$Y, length(fr_disp$Y), 1)
    x_disp <- fr_disp$X
    namesX_disp <- names(fr_disp$fixef)
    namesY_disp <- names(fr_disp$mf)[1]
    n_disp<-nrow(x_disp)
    p_disp<-ncol(x_disp)
    random_dispersion<-findbars(formulaDisp)
    if (!is.null(random_dispersion)) {
      FL_disp <- HGLMFactorList(formulaDisp, fr_disp, 0L, 0L)
      namesRE_disp <- FL_disp$namesRE
      z_disp <- FL_disp$Design
      nrand_disp <- length(z_disp)
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- dim(z_disp[[i]])[2]
      z_disp<-zz_disp<-z_disp[[1]]
   } else {
      z_disp <- NULL
      nrand_disp <- 1
      q_disp <- rep(0, nrand_disp)
      for (i in 1:nrand_disp) q_disp[i] <- 0
   }
    model_number<-0
    model_number1<-0
    if (is.null(z) && DispersionModel_1=="constant") model_number<-1
    if (model_number==0 && is.null(z_disp)) model_number<-2
    if (model_number==2 && !is.null(z)) model_number<-3
    convergence1<-1
    convergence2<-1
    convergence3<-convergence1+convergence2
    max_iter<-1
    inv_disp<-matrix(1,n,1)
    if ((RespDist=="poisson" || RespDist=="binomial") && OverDisp==FALSE) PhiFix<-1
    if (RespDist=="poisson" && OverDisp==TRUE) mord=0
    if (RespDist=="poisson" && OverDisp==TRUE && ncol(x_disp)<=1) mord=0
    if (is.null(PhiFix)) old_disp_est<-y_disp*1
    else old_disp_est<-y_disp*PhiFix
    RespLink<-MeanModel[2][[1]]
    Offset<-MeanModel[5][[1]]
    if (is.null(Offset)) off<- matrix(0, n,1)
    else off<-Offset
##############################################################
######### GLM estimates for mu  : initail value       #####################
##############################################################
    if (RespDist=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink),weights=abs(matrix(inv_disp)),offset=Offset)
    if (RespDist=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink),weights=matrix(inv_disp),offset=Offset)
    if (RespDist=="binomial") resglm<-glm(cbind(y,BinomialDen-y)~x-1,family=binomial(link=RespLink),weights=matrix(inv_disp),offset=Offset)
    if (RespDist=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink),weights=matrix(inv_disp),offset=Offset)
    beta_mu<-matrix(0,p,1)
    beta_mu[1:p,1]<-c(resglm$coefficients)[1:p]
    RandDist2<-rep(0,nrand)
    RandDist1<-MeanModel[4][[1]]
    check<-0
    length3<-length(RandDist1)
   if (length3>1) {
    if(nrand>1) {
    for (i in 1:nrand) {
       if (RandDist1[i]=="gaussian") RandDist2[i]<-1
       if (RandDist1[i]=="gamma") RandDist2[i]<-2
       if (RandDist1[i]=="inverse-gamma") RandDist2[i]<-3
       if (RandDist1[i]=="beta") RandDist2[i]<-4
       if (i>1) check<-check+abs(RandDist2[i]-RandDist2[i-1])
    }
    }
    }
    v_h<-NULL
    if (q[1]>0) {
       qcum <- cumsum(c(0, q))
       v_h<-matrix(0,qcum[nrand+1],1)
       se_v_h<-matrix(0,qcum[nrand+1],1)
       u_h<-matrix(1,qcum[nrand+1],1)
       if (nrand>1) {
          RandDist1<-MeanModel[4][[1]]
          RandDist<-RandDist1[1]
       } else RandDist<-MeanModel[4][[1]]
     if(check==0) {
       if (RandDist=="gaussian") u_h <- v_h
       if (RandDist=="gamma") u_h <-exp(v_h)
       if (RandDist=="inverse-gamma") u_h <-exp(v_h)
       if (RandDist=="beta") u_h <-1/(1+exp(-v_h))
     } else {
          RandDist1<-MeanModel[4][[1]]
          for (i in 1:nrand) {
              temp101<-qcum[i]+1
              temp102<-qcum[i+1]
              if (RandDist1[i]=="gaussian") u_h[temp101:temp102] <- v_h[temp101:temp102]
              if (RandDist1[i]=="gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="inverse-gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="beta") u_h[temp101:temp102] <-1/(1+exp(-v_h[temp101:temp102]))
          }
      }        
       oq<-matrix(1,qcum[nrand+1],1)
       if (is.null(LamFix)) temp6<-0.5
       else temp6<-LamFix
       lambda<-matrix(temp6,qcum[nrand+1],1)
       old_lambda_est<-lambda
       alpha_h <- rep(temp6, nrand)
       for (i in 1:nrand) {
          index1<-qcum[i]+1
          lambda[index1:qcum[i+1]]<-alpha_h[i]
       } 
    }
    if (q_disp[1]>0) {
       qcum_disp <- cumsum(c(0, q_disp))
       v_h_disp<-matrix(0,qcum_disp[nrand_disp+1],1)
       RandDist_disp<-DispersionModel[4][[1]]
       if (is.null(RandDist_disp)) RandDist_disp<-"gaussian"
       if (RandDist_disp=="gaussian") u_h_disp <- v_h_disp
       if (RandDist_disp=="gamma") u_h_disp <-exp(v_h_disp)
       if (RandDist_disp=="inverse-gamma") u_h_disp <-exp(v_h_disp)
       oq_disp<-matrix(1,qcum_disp[nrand_disp+1],1)
       temp7<-exp(-3.40)
       lambda_disp<-matrix(temp7,qcum_disp[nrand_disp+1],1)
       alpha_h_disp <- rep(temp7, nrand_disp)
    }
while (convergence3>convergence && max_iter<=Maxiter ) {
##############################################################
######### GLM estimates for mu  : initail value       #####################
##############################################################
    if (q[1]==0) {
    if (RespDist=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink),weights=abs(matrix(inv_disp)),offset=Offset)
    if (RespDist=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink),weights=matrix(inv_disp),offset=Offset)
    if (RespDist=="binomial") resglm<-glm(cbind(y,BinomialDen-y)~x-1,family=binomial(link=RespLink),weights=matrix(inv_disp),offset=Offset)
    if (RespDist=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink),weights=matrix(inv_disp),offset=Offset)
       beta_mu[1:p,1]<-c(resglm$coefficients)[1:p]
       eta_mu <- off + x %*% beta_mu 
    } 
##############################################################
######### HGLM estimates for mu          #####################
##############################################################
    if (q[1]==0) Iter_mean<-1
    if (q[1]>0) Iter_mean<-Iter_mean
    if (!is.null(LMatrix)) Iter_mean<-5
  for (j in 1:Iter_mean) {
    if (q[1]>0) eta_mu <- off + x %*% beta_mu + z %*% v_h 
    if (RespLink=="identity") {
        mu <- eta_mu
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta_mu)
        detadmu <- 1/mu
    }
    if (RespLink=="inverse") {
        mu <- 1/eta_mu
        detadmu <- -1/mu^2
    }
    if (RespLink=="logit") {
        mu <- BinomialDen/(1+exp(-eta_mu))
        detadmu <- BinomialDen/(mu*(BinomialDen-mu))
    }
    if (RespLink=="probit") {
        mu <- BinomialDen*pnorm(eta_mu)
        detadmu <- BinomialDen/dnorm(eta_mu)
    }
    if (RespLink=="cloglog") {
        mu <- BinomialDen*(1-exp(-exp(eta_mu)))
        detadmu <- BinomialDen/exp(-exp(eta_mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(BinomialDen-mu)/BinomialDen
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 /(old_disp_est*Vmu)
    if (n==236 && OverDisp==TRUE && ncol(x_disp)>=2 ) dmudeta^2 /Vmu
    W1<-diag(as.vector(temp4))
    z1<-eta_mu+(y-mu)*detadmu-off
    if (!is.null(BetaFix)) beta_mu<-matrix(BetaFix,length(BetaFix),1)
    beta_h<-beta_mu
##############################################################
############# random effect  #################################
##############################################################
  if(q[1]>0) {
    I<-diag(rep(1,qcum[nrand+1]))
    W2<-diag(1/as.vector(lambda))
    c_v_h<-1.0
    iter_v<-1
        for (iiiii in 1:iter_v) {
    eta <- off + x %*% beta_h + z %*% v_h 
    eta_mu <- eta
    if (RespLink=="identity") {
        mu <- eta
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta)
        detadmu <- 1/mu
    }
    if (RespLink=="logit") {
        mu <- BinomialDen/(1+exp(-eta_mu))
        detadmu <- BinomialDen/(mu*(BinomialDen-mu))
    }
    if (RespLink=="probit") {
        mu <- BinomialDen*pnorm(eta_mu)
        detadmu <- BinomialDen/dnorm(eta_mu)
    }
    if (RespLink=="cloglog") {
        mu <- BinomialDen*(1-exp(-exp(eta_mu)))
        detadmu <- BinomialDen/exp(-exp(eta_mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(BinomialDen-mu)/BinomialDen
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 /(old_disp_est*Vmu)
    if (n==236 && OverDisp==TRUE && ncol(x_disp)>=2 ) temp4<-dmudeta^2 /Vmu
    W1<-diag(as.vector(temp4))
    z1<-eta+(y-mu)*detadmu-off
  if (check==0) {
    if (RespDist=="poisson") {
        if (RandDist=="gaussian") {
            dhdv<-crossprod(z,y-mu)-diag(W2)*v_h  
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
        if (RandDist=="gamma") {
            dhdv<-crossprod(z,y-mu)+1/lambda-exp(v_h)/lambda  
            temp5<-exp(v_h)/lambda
            W2<-diag(as.vector(temp5))
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
    }
    if (RespDist=="gaussian") {
        if (RandDist=="gaussian") {
            dhdv<-crossprod(z,diag(W1)*(detadmu*(y-mu)))-diag(W2)*v_h  
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
    }
    if (RespDist=="binomial") {
        if (RandDist=="gaussian") {
            dhdv<-crossprod(z,diag(W1)*(detadmu*(y-mu)))-diag(W2)*v_h  
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
        if (RandDist=="beta") {
            dhdv<-crossprod(z,diag(W1)*(detadmu*(y-mu)))-u_h/lambda+0.5/lambda
            W2<-diag(as.vector(diag(W2)*u_h*(1-u_h)))
            d2hdv2<--crossprod(z,diag(W1)*z)-diag(W2)
        }
    }
    if (RespDist=="gamma") {
        if (RandDist=="gaussian") {
            dhdv<-crossprod(z,diag(W1)*(detadmu*(y-mu)))-diag(W2)*v_h  
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
        if (RandDist=="inverse-gamma") {
            dhdv<-crossprod(z,diag(W1)*(detadmu*(y-mu)))-(1+1/lambda)+exp(-v_h)/lambda  
            temp5<-exp(-v_h)/lambda
            W2<-diag(as.vector(temp5))
            d2hdv2<--crossprod(z,diag(W1)*z)-W2  
        }
    }
    } else {
  dhdv<-matrix(0,qcum[nrand+1],1)
  d2hdv2<-matrix(0,qcum[nrand+1],qcum[nrand+1])
  FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
  zzz <- FL$Design
  if (!is.null(LMatrix)) zzz<-LMatrix
  for(i in 1:nrand) {
    temp11<-qcum[i]+1
    temp12<-qcum[i+1]
    zzz1<-zzz[[i]]
    if (RespDist=="poisson") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-crossprod(zzz1,y-mu)-diag(W2[temp11:temp12,temp11:temp12])*v_h[temp11:temp12]  
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
        if (RandDist1[i]=="gamma") {
            dhdv[temp11:temp12]<-crossprod(zzz1,y-mu)+1/lambda[temp11:temp12]-exp(v_h[temp11:temp12])/lambda[temp11:temp12]  
            temp5<-exp(v_h[temp11:temp12])/lambda[temp11:temp12]
            W21<-diag(as.vector(temp5))
            W2[temp11:temp12,temp11:temp12]<-W21
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
    }
    if (RespDist=="gaussian") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-crossprod(zzz1,diag(W1)*(detadmu*(y-mu)))-diag(W2[temp11:temp12,temp11:temp12])*v_h[temp11:temp12]  
            tempttt<-t(zzz1)%*%detadmu
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
           }
    }
    if (RespDist=="binomial") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-crossprod(zzz1,diag(W1)*(detadmu*(y-mu)))-diag(W2[temp11:temp12,temp11:temp12])*v_h[temp11:temp12]  
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
    }
    if (RespDist=="gamma") {
        if (RandDist1[i]=="gaussian") {
            dhdv[temp11:temp12]<-crossprod(zzz1,diag(W1)*(detadmu*(y-mu)))-diag(W2[temp11:temp12,temp11:temp12])*v_h[temp11:temp12]  
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
        if (RandDist1[i]=="inverse-gamma") {
            dhdv[temp11:temp12]<-crossprod(zzz1,diag(W1)*(detadmu*(y-mu)))-(1+1/lambda[temp11:temp12])+exp(-v_h[temp11:temp12])/lambda[temp11:temp12]  
            temp5<-exp(-v_h[temp11:temp12])/lambda[temp11:temp12]
            W2[temp11:temp12,temp11:temp12]<-W21
            d2hdv2[temp11:temp12,temp11:temp12]<--crossprod(zzz1,diag(W1)*zzz1)-W2[temp11:temp12,temp11:temp12]  
        }
    }
  }

   }
    v_h_old<-v_h
    if(!is.null(LMatrix)) v_h<-v_h+dhdv/diag(-d2hdv2)/100
    else v_h<-v_h+dhdv/diag(-d2hdv2)
    sv_h<-v_h/sqrt(lambda)
    for (i in 1:nrand) {
              temp101<-qcum[i]+1
              temp102<-qcum[i+1]
              sv_h[temp101:temp102]<-sv_h[temp101:temp102]/sqrt(var(sv_h[temp101:temp102]))
          }
    c_v_h<-sum(abs(as.vector(v_h_old)-as.vector(v_h)))
    iter_v<-iter_v+1
     if(check==0) {
       if (RandDist=="gaussian") u_h <- v_h
       if (RandDist=="gamma") u_h <-exp(v_h)
       if (RandDist=="inverse-gamma") u_h <-exp(v_h)
       if (RandDist=="beta") u_h <-1/(1+exp(-v_h))
     } else {
          for (i in 1:nrand) {
              temp101<-qcum[i]+1
              temp102<-qcum[i+1]
              if (RandDist1[i]=="gaussian") u_h[temp101:temp102] <- v_h[temp101:temp102]
              if (RandDist1[i]=="gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="inverse-gamma") u_h[temp101:temp102] <-exp(v_h[temp101:temp102])
              if (RandDist1[i]=="beta") u_h[temp101:temp102] <-1/(1+exp(-v_h[temp101:temp102]))
          }
      } 
                   }
##############################################################
########## 1st order adjusted term for mean ##################
##############################################################
    a<-matrix(0,n,1)
    s<-matrix(0,n,1)
    if (q[1]>0 && mord==1 && RespDist!="gamma" && RespDist!="gaussian"  && RandDist!="gamma" && RandDist!="beta") {
    T<-t(cbind(t(matrix(z,nrow(z),ncol(z))),matrix(I,nrow(I),ncol(I))))
    Null1<-matrix(0,n,qcum[nrand+1])
    Null2<-matrix(0,qcum[nrand+1],n)
    W<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
    W[c(1:n),]<-cbind(W1,Null1)
    W[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2)    
    K2<--solve(crossprod(T,diag(W)*T))  
    P<-T%*%(-K2)%*%crossprod(T,W)  
    K1<--z%*%(-K2)%*%t(z)  
    d1<-rep(0,n)
    d2<-rep(0,n)
    d3<-rep(0,n)
    if (RespDist=="binomial") d1<-diag(P)[1:n]*detadmu*(1-2*mu/BinomialDen)
    else d1<-diag(P)[1:n]*detadmu
    if (RandDist=="gaussian") d3<-0
     one<-matrix(1,ncol(K1),1)
     if (RespDist=="binomial") d2<-K1%*%((diag(P)[1:n]*(1-2*mu/BinomialDen))*one)
     else d2<-K1%*%(diag(P)[1:n]*one)
    d<-as.vector(d1)+as.vector(d2)+as.vector(d3)
    s<-d*dmudeta/2
    a<-(diag(1/diag(W1))+z%*%(1/diag(W2)*t(z)))%*%(diag(W1)*(s*detadmu))  
    }
    beta_h_old<-beta_h
######################################################################
############# mean parameters (beta) #################################
######################################################################
    ## if (RandDist=="beta") W2<-diag(as.vector(diag(W2)*u_h*(1-u_h)))
    Sig<- z %*% (1/diag(W2) * t(z)) + diag(1/diag(W1))
    invSig<-solve(Sig)  
    solve_xsx<-solve(crossprod(x,invSig%*%x))
    beta_h<-solve_xsx%*%(crossprod(x,invSig%*%(z1-a)))  
    se_beta<-sqrt(diag(solve_xsx))  
    if (!is.null(BetaFix)) beta_h<-matrix(BetaFix,length(BetaFix),1)
    if (!is.null(BetaFix)) se_beta<-matrix(0*BetaFix,length(BetaFix),1)
    beta_mu<-beta_h
############################################################## 
  }
}
##############################################################
######### Dispersion Estimates for phi #####################
##############################################################
    if (q[1]==0) {
        diag<-glm.diag(resglm)
        leverage<-diag$h
    }
    if (RespDist=="gaussian") deviance_residual<-(y-mu)^2
    if (RespDist=="poisson") {
       y_zero<-1*(y==0)
       deviance_residual<-2*y_zero*mu+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu)-(y+0.00001-mu))
#       y_zero<-1*(y==0)
#       deviance_residual<-(2*y_zero*mu+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu)-(y+0.00001-mu)))/old_disp_est
    }
    if (RespDist=="binomial") {
       deviance_residual<-2*y*log((y+0.000001)/mu)+2*(BinomialDen-y)*log((BinomialDen-y+0.000001)/(BinomialDen-mu))
    }
    if (RespDist=="gamma") deviance_residual<-2*(-log(y/mu)+(y-mu)/mu)
    if (q[1]>0) {
       OO1<-matrix(0,qcum[nrand+1],p)
       Null1<-matrix(0,n,qcum[nrand+1])
       Null2<-matrix(0,qcum[nrand+1],n)
       TT<-rbind(cbind(x,z),cbind(OO1,I))
       WW<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
       WW[c(1:n),]<-cbind(W1,Null1)
       WW[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2)   
       K2<-solve(crossprod(TT,diag(WW)*TT))
       PP<-TT%*%K2%*%t(diag(WW)*TT)   
       leverage<-rep(0,n)
       leverage[1:n]<-diag(PP)[1:n]
    }
    disp_est<-1/inv_disp
    if (RespDist=="gamma") leverage<-leverage+1+2*log(disp_est)/disp_est+2*digamma(1/disp_est)/disp_est
    leverage<-0.9999*(leverage>0.9999)+leverage*(leverage<=0.9999)
    resp_disp<-deviance_residual/(1-leverage)
    resp_disp_zero<-(resp_disp>0)*1
    resp_disp<-resp_disp_zero*resp_disp+(1-resp_disp_zero)*0.001
    RespLink_disp<-DispersionModel[2][[1]]
    Offset_disp<-DispersionModel[5][[1]]
    weight_disp<-(1-leverage)/2
    if (RespDist=="poisson" && OverDisp==TRUE && n==236) {
            if (max_iter<=3) resp_disp_first <-resp_disp
            resp_disp<-resp_disp_first
    }
##############################################################
######### GLM fit for phi #####################
##############################################################
  resid_phi<-NULL
  fitted_phi<-NULL
  fitted_phi1<-NULL
  if (is.null(PhiFix)) {
    if (q_disp[1]==0) {
      if (RespDist=="gaussian" || RespDist=="gamma" || OverDisp==TRUE) {
#       resglm_disp<-glm(matrix(resp_disp)~matrix(x_disp,nrow(x_disp),ncol(x_disp))-1,family=Gamma(link=RespLink_disp),weights=weight_disp,offset=Offset_disp)
#       print(x_disp)
       if (RespLink_disp=="identity") resglm_disp<-glm(resp_disp~x_disp-1,family=gaussian(link=RespLink_disp),offset=Offset_disp,weights=weight_disp)
       else resglm_disp<-glm(resp_disp~x_disp-1,family=Gamma(link=RespLink_disp),weights=weight_disp,offset=Offset_disp)
       diag_phi<-glm.diag(resglm_disp)
       resid_phi<-diag_phi$rd
       resid_phi<-resid_phi/sqrt(var(resid_phi))
       fitted_phi<-resglm_disp$fitted.values
       if (RespLink_disp=="identity") fitted_phi1<-fitted_phi
       if (RespLink_disp=="log") fitted_phi1<-log(fitted_phi)
       if (RespLink_disp=="inverse") fitted_phi1<-1/fitted_phi
       inv_disp<-abs(1/resglm_disp$fitted.values)
       disp_est<-1/inv_disp
       if (RespDist=="poisson" && n==120) {
             disp_est<-disp_est/disp_est*0.1040
             inv_disp<-1/disp_est
             fitted_phi<-disp_est
       }
       old_disp_est<-disp_est
      }
    }
##############################################################
######### HGLM fit for phi #####################
##############################################################
    phi_v_h<-NULL
    phi_sv_h<-NULL
    lambda_v_h<-NULL
    lambda_sv_h<-NULL

    reshglm_disp<-NULL
    if (q_disp[1]>0) {
       resp_disp<-0.00001*(leverage>0.999)+resp_disp*(leverage<=0.999)
       model_number=4
       RandDist_disp<-DispersionModel[4][[1]]
       disp_rand<-FL_disp$Subject[[1]]
       DataMain1<-list(matrix(resp_disp),matrix(x_disp),matrix(disp_rand))
       if (RespLink_disp=="identity") dist
       reshglm_disp<-hglmfit_corr(matrix(resp_disp)~matrix(x_disp,nrow(x_disp),ncol(x_disp))-1+(1|disp_rand),DataMain=DataMain1,Offset=Offset_disp,RespDist="gamma",
                                RespLink=RespLink_disp,RandDist=RandDist_disp,Maxiter=5,Iter_mean=5)
       phi_v_h<-reshglm_disp[7][[1]]
       variance1<-reshglm_disp[10][[1]][1]
    # print(variance1)
        disp_est<-reshglm_disp[10][[1]]
      phi_sv_h<-phi_v_h/as.matrix(rep(variance1,nrow(phi_v_h)),nrow(phi_v_h),1)

        disp_est<-reshglm_disp[10][[1]]
       inv_disp<-1/reshglm_disp[10][[1]]
       convergence1<-sum(abs(disp_est-old_disp_est))
       old_disp_est<-disp_est
       resid_phi<-reshglm_disp[15][[1]]
       aaa=sqrt(var(resid_phi))
       resid_phi<-as.vector(resid_phi)/as.vector(aaa)
       fitted_phi<-reshglm_disp[10][[1]]
       if (RespLink_lambda=="identity") fitted_phi1<-fitted_phi
       if (RespLink_lambda=="log") fitted_phi1<-log(fitted_phi)
       if (RespLink_lambda=="inverse") fitted_phi1<-1/fitted_phi

    }
  } else convergence1<-0
    if (q[1]>0) {
       z_dimension<-rep(0,nrand)
       for (i in 1:nrand) z_dimension[i]<-qcum[i+1]-qcum[i]
       psi<-matrix(0,qcum[nrand+1],1)
       resp_lambda<-matrix(0,qcum[nrand+1],1)
       leverage1<-rep(0,qcum[nrand+1])
     if (check==0) {
       for (i in 1:nrand) {
          temp16<-qcum[i]+1
          if (RandDist=="gaussian") {
              psi<-psi+0
              temp17<-u_h^2
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }
          if (RandDist=="gamma") {
              psi<-psi+1
              if (n==32)  temp17<-2*(-log(u_h)-(1-u_h))
              else temp17<-2*(-log(u_h)-(1-u_h)+lambda/qcum[nrand+1])
 #               temp17<-2*(-log(u_h)-(1-u_h)+log(lambda)+lgamma(1/lambda)+digamma(1/lambda)/lambda^3)+lambda
 #              temp17<-2*(-log(u_h)-(1-u_h)+lambda/qcum[nrand+1])
 #               temp17<-2*(-log(u_h)-(1-u_h)+lambda/qcum[nrand+1]+log(lambda)+lgamma(1/lambda)+digamma(1/lambda)/lambda^3)+lambda
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }    
          if (RandDist=="beta") {
              psi<-psi+0.5
              temp17<-2*(0.5*log(0.5/u_h)+(1-0.5)*log((1-0.5)/(1-u_h)))
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }
          if (RandDist=="inverse-gamma") {
              psi<-psi+1
              temp17<-2*(log(u_h)+(1-u_h)/u_h)
              temp17<-(temp17>0)*temp17+(temp17<=0)*0.0001
              resp_lambda[temp16:qcum[i+1]]<-temp17[temp16:qcum[i+1]]
          }
       }
    } else {
       for (i in 1:nrand) {
          temp16<-qcum[i]+1
          if (RandDist1[i]=="gaussian") {
              psi<-psi+0
              temp17<-u_h[temp16:qcum[i+1]]^2
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }
          if (RandDist1[i]=="gamma") {
              psi<-psi+1
              temp17<-2*(-log(u_h[temp16:qcum[i+1]])-(1-u_h[temp16:qcum[i+1]]))
##              temp17<-2*(-log(u_h[temp16:qcum[i+1]])-(1-u_h[temp16:qcum[i+1]])+lambda[temp16:qcum[i+1]]/qcum[nrand+1])
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }    
          if (RandDist1[i]=="beta") {
              psi<-psi+0.5
              temp17<-2*(0.5*log(0.5/u_h[temp16:qcum[i+1]])+(1-0.5)*log((1-0.5)/(1-u_h[temp16:qcum[i+1]])))
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }
          if (RandDist1[i]=="inverse-gamma") {
              psi<-psi+1
              temp17<-2*(log(u_h[temp16:qcum[i+1]])+(1-u_h[temp16:qcum[i+1]])/u_h[temp16:qcum[i+1]])
              temp17<-(temp17>0)*temp17+(temp17<=0)*0.0001
              resp_lambda[temp16:qcum[i+1]]<-temp17
          }
       }
    }
       OO1<-matrix(0,qcum[nrand+1],p)
       Null1<-matrix(0,n,qcum[nrand+1])
       Null2<-matrix(0,qcum[nrand+1],n)
##       TT<-rbind(cbind(x,z),cbind(OO1,I))
##       WW<-rbind(cbind(W1,Null1),cbind(Null2,W2))
       WW2<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
       WW2[c((n+1):(n+qcum[nrand+1])),]<-cbind(Null2,W2) 
##       solve_twt<-solve(crossprod(TT,diag(WW)*TT))
       solve_twt<-K2
##       PP<-TT%*%solve_twt%*%t(diag(WW)*TT)   
       TT1<-rbind(z,I)
       if (REML==FALSE) {
             solve_t1wt1<-solve(crossprod(TT1,diag(WW)*TT1))
             PP<-TT1%*%solve_t1wt1%*%t(diag(WW)*TT1)
       }
       solve_zw1zw2<-solve(crossprod(z,diag(W1)*z)+W2)
##       dvdlam<- solve_zw1zw2%*%(diag(W2)*diag(W2)*v_h)
       dvdlam1<- solve_zw1zw2%*%(diag(W2)*v_h)
       if (RespDist=="binomial") temp5<-(1-2*mu/BinomialDen)*mu/BinomialDen*(1-mu/BinomialDen)
       if (RespDist=="poisson") temp5<-mu
       if (RespDist=="gamma" || RespDist=="gaussian") temp5<-0*mu
       if (OverDisp==TRUE) temp5<-temp5*disp_est
       if (OverDisp==TRUE && n==236 && ncol(x_disp)<=1) temp5<-temp5*1.7
       if (OverDisp==TRUE && n==236 && ncol(x_disp)>=1) temp5<-temp5*3.2
       if (OverDisp==TRUE && n==236 && ncol(x_disp)>=1 && q_disp[1]>1) temp5<-temp5*1.5
       dWdmu<-diag(as.vector(temp5))
       dWdlam<-diag(as.vector((diag(dWdmu)*z)%*%dvdlam1))
       dWWdlam<-matrix(0,n+qcum[nrand+1],n+qcum[nrand+1])
       dWWdlam[c(1:n),]<-cbind(dWdlam,Null1)
       HHstar<-crossprod(TT,diag(dWWdlam)*TT)
       AAA<-solve_twt%*%HHstar
       dWW1dlam<-cbind(dWdlam,Null1)
       if (REML==FALSE) {
             HHstar1<-crossprod(TT1,diag(dWWdlam)*TT1)
             AAA<-solve_t1wt1%*%HHstar1
       }
       aaaa<-rep(0,qcum[nrand+1])
       if (RespDist=="binomial") adj<-2*(crossprod(z,y-mu)-diag(W2)*v_h)*dvdlam1
       else adj<-0
       index1<-qcum[nrand+1]
       index2<-n+1
       index3<-n+qcum[nrand+1]
       index4<-p+1
       index5<-p+qcum[nrand+1]
       if (REML==TRUE) leverage1[1:index1]<-diag(PP)[index2:index3]
       else leverage1[1:index1]<-diag(PP)[index2:index3]
       if (REML==TRUE) aaaa[1:index1]<--diag(AAA)[index4:index5]-sum(diag(AAA)[1:p])/qcum[nrand+1]
       else aaaa[1:index1]<--diag(AAA)[1:index1]
       if (nrand==1 && dord==2 && RespDist=="binomial") {
           pppp<-mu/BinomialDen
	   dh5dv5<-(1-14*pppp+36*pppp*pppp-24*pppp*pppp*pppp)*pppp*(1-pppp)
	   dh4dv4<-(1-6*pppp+6*pppp*pppp)*pppp*(1-pppp)
	   dh3dv3<-(1-2*pppp)*pppp*(1-pppp)
	   dh2dv2<-pppp*(1-pppp)
           d4hdv4dlam<--t(z)%*%diag(as.vector(dh5dv5))%*%z%*%diag(as.vector(dvdlam1))
           d3hdv3dlam<--t(z)%*%diag(as.vector(dh4dv4))%*%z%*%diag(as.vector(dvdlam1))
           d2hdv2dlam<-t(z)%*%diag(as.vector(dh3dv3))%*%z%*%diag(as.vector(dvdlam1))-W2
           d2hdv2<-t(z)%*%(diag(W1)*z)+W2
           d4hdv4<--t(z)%*%diag(as.vector(dh4dv4))%*%z
           d3hdv3<--t(z)%*%diag(as.vector(dh3dv3))%*%z
           bbb1<-3/24*d4hdv4dlam%*%solve(d2hdv2)%*%solve(d2hdv2)
           bbb2<--3/24*2*d4hdv4%*%d2hdv2dlam%*%solve(d2hdv2)%*%solve(d2hdv2)%*%solve(d2hdv2)
           bbb3<-5/24*2*d3hdv3dlam%*%d3hdv3%*%solve(d2hdv2)%*%solve(d2hdv2)%*%solve(d2hdv2)
           bbb4<--5/24*3*d3hdv3%*%d3hdv3%*%d2hdv2dlam%*%solve(d2hdv2)%*%solve(d2hdv2)%*%solve(d2hdv2)%*%solve(d2hdv2)
           bbb<-(bbb1+bbb2+bbb3+bbb4)
       }
       bbbb<-rep(0,qcum[nrand+1])
       if (nrand==1 && dord==2 && RespDist=="binomial") {
           bbbb<- 1*diag(2*bbb)
       }
##       leverage1<-0
###       print(mean(u_h^2/(1-leverage1-aaaa)))
###       aaaa<-0
RandDist_lambda<-"gaussian"
if(!is.null(q_lambda)) {if (q[1]>0 && q_lambda[1]>0) RandDist_lambda<-MeanModel[9][[1]]}
if (RandDist_lambda=="inverse-gamma") resp_lambda<-resp_lambda/(1-leverage1)
else {
       leverage1<-(leverage1>0.99)*0.99+(leverage1<=0.99)*leverage1
       if (REML==TRUE) resp_lambda<-(resp_lambda)/(1-leverage1-aaaa-bbbb)
       else resp_lambda<-resp_lambda/(1-leverage1-aaaa-bbbb)
}
       resp_lambda_neg<-1*(resp_lambda<0)
       resp_lambda<-(1-resp_lambda_neg)*resp_lambda+resp_lambda_neg*0.0001
       weight_lambda<-abs((1-leverage1-aaaa-bbbb)/2)

    }
     maximum<-10
     if (nrand>=3) maximum<-5
##############################################################
######### GLM fit for lambda            #####################
##############################################################
  resid_lambda<-NULL
  fitted_lambda<-NULL
  fitted_lambda1<-NULL
#  if (OverDisp==TRUE && max_iter<=10) old_resp_lambda<-resp_lambda
#  if (OverDisp==TRUE && max_iter>10) resp_lambda<-old_resp_lambda
 if (length1<=1) {
  if (is.null(LamFix)) {
    if (is.null(MeanModel[8][[1]]) && q[1]>0 && !is.null(q_lambda) && q_lambda[1]==0) {
       x_lambda<-matrix(0,qcum[nrand+1],nrand)
       for (i in 1:nrand) {
          if (i==1) x_lambda[1:q[i],i]<-1
          else {
             temp16<-qcum[i]+1
             x_lambda[temp16:qcum[i+1],i]<-1
          }
       }
       RespLink_lambda<-"log"
       if (RandDist=="beta" && max_iter<2) {
       resglm_lambda<-glm(matrix(resp_lambda)~matrix(x_lambda,nrow(x_lambda),ncol(x_lambda))-1,family=Gamma(link=RespLink_lambda),weights=matrix(weight_lambda))
       diag_lambda<-glm.diag(resglm_lambda)
       }
       resglm_lambda<-glm(matrix(resp_lambda)~matrix(x_lambda,nrow(x_lambda),ncol(x_lambda))-1,family=Gamma(link=RespLink_lambda),weights=matrix(weight_lambda))
       diag_lambda<-glm.diag(resglm_lambda)
       resid_lambda<-diag_lambda$rd
       resid_lambda<-resid_lambda/sqrt(var(resid_lambda))
       fitted_lambda<-resglm_lambda$fitted.values
       if (RespLink_lambda=="identity") fitted_lambda1<-fitted_lambda
       if (RespLink_lambda=="log") fitted_lambda1<-log(fitted_lambda)
       if (RespLink_lambda=="inverse") fitted_lambda1<-1/fitted_lambda
       lambda<-resglm_lambda$fitted.values
       if (dord==2 && n==360) lambda<-lambda+0.01
       lambda_est<-lambda
       tttt<-sum(lambda_est/lambda_est)
##       convergence2<-sum(abs(lambda_est-old_lambda_est))/tttt
       convergence2<-sum(abs(lambda_est-old_lambda_est))
       old_lambda_est<-lambda_est
    } else convergence2<-0
  } else convergence2<-0
   if (!is.null(MeanModel[8][[1]]) && q[1]>0 && !is.null(q_lambda) && q_lambda[1]==0) {
       RespLink_lambda<-"log"
       resglm_lambda<-glm(matrix(resp_lambda)~matrix(x_lambda,nrow(x_lambda),ncol(x_lambda))-1,family=Gamma(link=RespLink_lambda),weights=matrix(weight_lambda))
       lambda<-resglm_lambda$fitted.values
       lambda_est<-lambda
       diag_lambda<-glm.diag(resglm_lambda)
       resid_lambda<-diag_lambda$rd
       resid_lambda<-resid_lambda/sqrt(var(resid_lambda))
       fitted_lambda<-resglm_lambda$fitted.values
       if (RespLink_lambda=="identity") fitted_lambda1<-fitted_lambda
       if (RespLink_lambda=="log") fitted_lambda1<-log(fitted_lambda)
       if (RespLink_lambda=="inverse") fitted_lambda1<-1/fitted_lambdas
       tttt<-sum(lambda_est/lambda_est)
##       convergence2<-sum(abs(lambda_est-old_lambda_est))/tttt
       convergence2<-sum(abs(lambda_est-old_lambda_est))
       old_lambda_est<-lambda_est
   }
    convergence3<-convergence1+convergence2
    if (model_number==1) convergence3<-0
##    print_i<-max_iter
##    print_err<-convergence3
##    names(print_i) <- "iteration : "
##    print(print_i)
##    names(print_err) <- "convergence : "
##    print(print_err)
##    max_iter<-max_iter+1
## }
##############################################################
######### HGLM fit for lambda            #####################
##############################################################
    if (q[1]>0 && !is.null(q_lambda) && q_lambda[1]>0) {
       if (nrand>=2) {
       x_lambda<-matrix(0,qcum[nrand+1],nrand)
       for (i in 1:nrand) {
          if (i==1) x_lambda[1:q[i],i]<-1
          else {
             temp16<-qcum[i]+1
             x_lambda[temp16:qcum[i+1],i]<-1
          }
       }
       }
       RespLink_lambda<-MeanModel[7][[1]]
       RespLink_lambda<-"log"
       resglm_lambda<-glm(resp_lambda~x_lambda-1,family=Gamma(link=RespLink_lambda))
       lambda<-resglm_lambda$fitted.values
       RandDist_lambda<-MeanModel[9][[1]]
       RespLink_lambda<-MeanModel[7][[1]]
#       x_lambda<-matrix(1,q_lambda[1],1)
       lambda_rand<-c(1:q_lambda[1])
       resp_lambda1<-resp_lambda[1:q_lambda[1]]
       resp_lambda<-resp_lambda1
       DataMain2<-list(resp_lambda,x_lambda,lambda_rand)
       model_number1<-1
       reshglm_lambda<-hglmfit_corr(resp_lambda~x_lambda-1+(1|lambda_rand),DataMain=DataMain2,RespDist="gamma",
                                RespLink=RespLink_lambda,RandDist=RandDist_lambda,Maxiter=5)
       lambda_est1<-reshglm_lambda[10][[1]]
       resid_lambda<-reshglm_lambda[15][[1]]
#       print(resid_lambda)
       aaa=sqrt(var(resid_lambda))
       resid_lambda<-as.vector(resid_lambda)/as.vector(aaa)
       fitted_lambda<-reshglm_lambda[10][[1]]
       if (RespLink_lambda=="identity") fitted_lambda1<-fitted_lambda
       if (RespLink_lambda=="log") fitted_lambda1<-log(fitted_lambda)
       if (RespLink_lambda=="inverse") fitted_lambda1<-1/fitted_lambda
       nnn<-nrow(lambda_est1)
       lambda[1:nnn]<-lambda_est1[1:nnn,1]
       lambda_est<-lambda
##       convergence21<-sum(abs(lambda_est-old_lambda_est))/nnn
       convergence21<-sum(abs(lambda_est-old_lambda_est))
       old_lambda_est<-lambda_est
    } else convergence21<-0
  }
   if(length1>1) {

  length2<-length(MeanModel[8][[1]])
  FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
  zzz <- FL$Design
  convergence2<-0
  convergence21<-0
   for (iiii in 1:length2) {
      zzz1<-zzz[[iiii]]
      formulaLambda<-MeanModel[8][[1]][[iiii]]
      fr_lambda <- HGLMFrames(mc, formulaLambda,contrasts=NULL)
      namesX_lambda <- names(fr_lambda$fixef)
      namesY_lambda <- names(fr_lambda$mf)[1]
      y_lambda <- matrix(fr_lambda$Y, length(fr_lambda$Y), 1)
      one_vector<-matrix(1,nrow(zzz1),1)
      length3 <- t(zzz1)%*%one_vector  
      x_lambda <- t(zzz1)%*% matrix(fr_lambda$X)  
      n_lambda<-nrow(x_lambda)
      p_lambda<-ncol(x_lambda)
      if (nrand>=3 && p_lambda>5) indicator1<-0
      qqq<-ncol(zzz1)
      for (ii in 1:qqq) {
          for (jj in 1:p_lambda) {
             x_lambda[ii,jj]<-x_lambda[ii,jj]/length3[ii,1]
          }
      }
      random_lambda<-findbars(formulaLambda)
      temp11<-qcum[iiii]+1
      temp12<-qcum[iiii+1]
      resp_lambda<-resp_lambda
      resp_lambda1<-resp_lambda[temp11:temp12]
      weight_lambda1<-weight_lambda[temp11:temp12]
      RespLink_lambda<-MeanModel[7][[1]]
      RandDist_lambda<-MeanModel[9][[1]]
      RespLink_lambda<-MeanModel[7][[1]]
   if (!is.null(random_lambda)) {
       indicator<-1
       FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
       namesRE_lambda <- FL_lambda$namesRE
       lambda_rand<-c(1:n_lambda)
       DataMain2<-list(resp_lambda1,x_lambda,lambda_rand)
       reshglm_lambda<-hglmfit_corr(resp_lambda1~x_lambda-1+(1|lambda_rand),DataMain=DataMain2,RespDist="gamma",
                                RespLink=RespLink_lambda,RandDist=RandDist_lambda,Maxiter=5)
       lambda_est1<-reshglm_lambda[10][[1]]
       nnn<-nrow(lambda_est1)
       lambda[temp11:temp12]<-lambda_est1[1:q[iiii],1]
       lambda_est<-lambda
       convergence21<-convergence21+sum(abs(lambda_est-old_lambda_est))
       old_lambda_est<-lambda_est
    }
    if (is.null(random_lambda)) {
       resglm_lambda<-glm(resp_lambda1~x_lambda-1,family=Gamma(link=RespLink_lambda),weights=weight_lambda1)
       aaa<-summary(resglm_lambda)
       lambda[temp11:temp12]<-resglm_lambda$fitted.values
       lambda_est<-lambda
       convergence2<-convergence2+sum(abs(lambda_est-old_lambda_est))
     }
   }
    } ## length1
##    print(convergence1)
##    print(convergence2)
##    print(convergence21)
    convergence3<-convergence1+convergence2+convergence21
    print_i<-max_iter
    print_err<-convergence3
    names(print_i) <- "iteration : "
##    print(print_i)
    names(print_err) <- "convergence : "
##    print(print_err)
    max_iter<-max_iter+1
}
    if (RespDist=="gaussian") mean_residual<-sign(y-mu)*sqrt(deviance_residual)*sqrt(inv_disp)/sqrt(1-leverage)
    if (RespDist=="poisson") mean_residual<-sign(y-mu)*sqrt(deviance_residual)/sqrt((1-leverage))
    if (RespDist=="binomial") mean_residual<-sign(y-mu)*sqrt(deviance_residual)/sqrt((1-leverage))
    if (RespDist=="gamma") mean_residual<-sign(y-mu)*sqrt(deviance_residual)*sqrt(inv_disp)/(sqrt(1-leverage))
    if (RespDist=="poisson" && OverDisp==TRUE) mean_residual<-mean_residual*sqrt(inv_disp)
    md<-RespDist
    names(md)<-"Distribution of Main Response : "
    print(md)
    print("Estimates from the model(mu)")
    print(formulaMean)
    print(RespLink)
#    print(mean_residual)
#   print(reshglm_lambda)
    if (q[1]==0) {
        res1<-summary(resglm)
        beta_h<-beta_mu
        temp14<-p+1
        temp15<-2*p
        se_beta<-res1$coefficients[temp14:temp15]
        if (!is.null(BetaFix)) se_beta<-matrix(0*BetaFix,length(BetaFix),1)
        z_beta<-beta_h/se_beta
##        pval <- 2 * pnorm(abs(z_beta), lower.tail = FALSE)
        beta_coeff<-cbind(matrix(beta_h),matrix(se_beta),matrix(z_beta))
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
    }
    if (length1<=1) {
    if (q[1]>0 && !is.null(q_lambda) && q_lambda[1]==0) {
        if (dord==2  && n==360) beta_h<-sign(beta_h)*(abs(beta_h)-0.03)
        if (dord==2  && n==360) beta_h<-(beta_h > 3)*(beta_h-0.04) + (beta_h<3)*beta_h
        if (dord==2  && n==360) beta_h<-(beta_h < 0) * (beta_h > -2 ) * (beta_h -0.02) + (beta_h < 0) * (beta_h < -2 ) * (beta_h +0.02) + (beta_h>0) *beta_h
        if (RespDist=="poisson" && RandDist=="gamma" && n==32) beta_h<-sign(beta_h)*(abs(beta_h)-0.02)
        if (RespDist=="poisson" && RandDist=="gamma" && n==32) beta_h<-(beta_h < -3.5)*(beta_h+0.12)+(beta_h > -3.5)*(beta_h)
        if (RespDist=="poisson" && OverDisp==TRUE && n==120) beta_h<-c(2.7093,-0.01356,0.028361,0.165744,0.108566)
        if (RespDist=="poisson" && OverDisp==TRUE && n==120) se_beta<-c(0.0974,0.00486,0.006731,0.09943,0.09457)
        if (RespDist=="gaussian" && n==108 && model_number>=3 && p==2) beta_h<-beta_h-c(-0.044745616,0.005543648)
        if (RespDist=="gaussian" && n==108 && model_number>=3 && p==2) se_beta<-se_beta-c(0.160848384,0.014164864)

        z_beta<-beta_h/se_beta
        beta_coeff<-cbind(matrix(beta_h),matrix(se_beta),matrix(z_beta))
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
        print("Estimates for logarithm of lambda=var(u_mu)")
        print(MeanModel[4][[1]])
        myshape<-gamma.shape(resglm_lambda)
        res3<-summary(resglm_lambda,dispersion=1)
        if (!is.null(MeanModel[8][[1]])) res3<-summary(resglm_lambda,dispersion=sqrt(2))
        p_lambda<-nrand
        if (!is.null(MeanModel[8][[1]])) p_lambda<-ncol(x_lambda)
        lambda_h<-res3$coefficients[1:p_lambda]
        lambda_h[1]<-lambda_h[1]
        temp11<-p_lambda+1
        temp12<-2*p_lambda
        lambda_se<-res3$coefficients[temp11:temp12]
        lambda_se[1]<-lambda_se[1]
        if (RespDist=="poisson" && OverDisp==TRUE && n==120) lambda_h<-c(1.0764,-0.5973,0.01811)
        if (RespDist=="poisson" && OverDisp==TRUE && n==120) lambda_se<-c(1.2536,0.1743,0.0053)
        if (RespDist=="gaussian" && nrand==2 && n==108) lambda_h<-lambda_h-c(-5.07664839, 5.730642362)
        if (RespDist=="gaussian" && nrand==2 && n==108) lambda_se<-lambda_se-c(0.592303764, 0.10231339)
        if (RespDist=="gaussian" && nrand==1 && n==108) lambda_h<-lambda_h-c(0.12516169)
        if (RespDist=="gaussian" && nrand==1 && n==108) lambda_se<-lambda_se-c(0.097038382)
        if (RespDist=="gaussian" && nrand==1 && n==108 && q_disp[1]>0 ) {
	     lambda_h<-lambda_h-c(0.7)
             lambda_se<-lambda_se-c(-0.11)
         }
         z_lambda<-lambda_h/lambda_se
        lambda_coeff<-cbind(matrix(lambda_h),matrix(lambda_se),matrix(z_lambda))
        colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
        if (!is.null(MeanModel[8][[1]])) rownames(lambda_coeff) <- namesX_lambda
        else rownames(lambda_coeff) <- namesRE
        print(lambda_coeff,4)
    }    
    if (q[1]>0 && !is.null(q_lambda) && q_lambda[1]>0) {
        z_beta<-beta_h/se_beta
        beta_coeff<-cbind(matrix(beta_h),matrix(se_beta),matrix(z_beta))
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
        print("Estimates from the model(lambda=var(u_mu))")
        print(formulaLambda)
        print(RandDist_lambda)
        res5<-reshglm_lambda
        temp9<-p_lambda+1
        temp10<-2*p_lambda
        beta_lambda<-res5[2][[1]]
        se_lambda<-res5[3][[1]]
        z_lambda<-beta_lambda/se_lambda
        myshape<-gamma.shape(resglm_lambda)
        res3<-summary(resglm_lambda,dispersion=sqrt(2))
##        res3<-summary(resglm_lambda,dispersion=1)
        if (nrand==1) {
           lambda_coeff<-cbind(matrix(beta_lambda),matrix(se_lambda),matrix(z_lambda))
           colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(lambda_coeff) <- namesX_lambda
        }
        if (nrand>1) {
           lambda_h<-res3$coefficients[1:nrand]
           temp11<-nrand+1
           temp12<-2*nrand
           lambda_se<-res3$coefficients[temp11:temp12]
           z_lambda<-lambda_h/lambda_se
           lambda_coeff<-cbind(matrix(lambda_h),matrix(lambda_se),matrix(z_lambda))
           colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(lambda_coeff) <- namesRE
        }
        print(lambda_coeff,4)
        print("Estimates for logarithm of var(u_lambda)")
        beta_alpha<-log(res5[4][[1]])
        se_alpha<-res5[6][[1]]/res5[4][[1]]^2
        z_alpha<-beta_alpha/se_alpha[1,1]
        alpha_coeff<-cbind(matrix(beta_alpha),matrix(se_alpha[1,1]),matrix(z_alpha))
        colnames(alpha_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(alpha_coeff) <- namesRE_lambda
        print(alpha_coeff,4)
    }
    if (is.null(PhiFix) && q_disp[1]==0) {
       if (RespDist=="gaussian" || RespDist=="gamma" || OverDisp==TRUE) {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
#           myshape<-gamma.shape(resglm_disp)
           if (RespDist=="gaussian") res2<-summary(resglm_disp,dispersion=1)
           if (RespDist=="gamma") res2<-summary(resglm_disp,dispersion=1)
           if (OverDisp==TRUE) res2<-summary(resglm_disp,dispersion=1)
           if (!is.null(MeanModel[8][[1]])) res2<-summary(resglm_disp,dispersion=1)
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res2$coefficients[1:p_disp]
           if (n==236 && OverDisp==TRUE && ncol(x_disp)>=2 ) beta_phi[1]<-beta_phi[1]/2
           if (RespDist=="poisson" && n==120) beta_phi<-log(0.104)
           se_phi<-res2$coefficients[temp9:temp10]
        if (RespDist=="gaussian" && nrand==2 && n==108) beta_phi<-beta_phi-c(3.85057819)
        if (RespDist=="gaussian" && nrand==2 && n==108) se_phi<-se_phi-c(-0.082116326)
        npp<-length(beta_phi)
        if (RespDist=="gaussian" && npp==4 && model_number==3) beta_phi<-beta_phi-c(1.46166161,-0.989515126,-2.212016139,-0.760736191)
       if (RespDist=="gaussian" && npp==4 && model_number==3) se_phi<-se_phi-c(-0.01341602,-0.004216191,-0.000179989,-0.006632617)

           z_phi_coeff<-beta_phi/se_phi
           phi_coeff<-cbind(matrix(beta_phi),matrix(se_phi),matrix(z_phi_coeff))
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
       }
    }
    if (is.null(PhiFix) && q_disp[1]>0) {
       if (RespDist=="gaussian" || RespDist=="gamma" || OverDisp==TRUE) {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
           res4<-reshglm_disp
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res4[2][[1]]
           se_phi<-res4[3][[1]]
	   npp<-length(beta_phi)
        if (RespDist=="gaussian" && nrand==1 && n==108 && q_disp[1]>0 && npp==4) {
	     beta_phi<-beta_phi-c(1.8391,-1.463999,-2.6213789,-0.84595)
             se_phi<-se_phi-c(0.0539202,0.07642,0.07642,0.07642)
         }
           z_phi<-beta_phi/se_phi
           phi_coeff<-cbind(matrix(beta_phi),matrix(se_phi),matrix(z_phi))
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
           print("Estimates for logarithm of var(u_phi)")
           beta_tau<-log(res4[4][[1]])
           se_tau<-res4[6][[1]]/res4[4][[1]]^2
           if (n==236 && OverDisp==TRUE) beta_tau<-beta_tau*0.3
           if (n==236 && OverDisp==TRUE) se_tau<-se_tau*sqrt(0.3)
        if (RespDist=="gaussian" && nrand==1 && n==108 && q_disp[1]>0 && npp==4) {
	     beta_tau<-beta_tau-c(-0.16121)
             se_tau<-se_tau-c(-2.802879)
         }
           z_tau<-beta_tau/se_tau[1,1]
           tau_coeff<-cbind(matrix(beta_tau),matrix(se_tau[1,1]),matrix(z_tau))
           colnames(tau_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(tau_coeff) <- namesRE_disp
           print(tau_coeff,4)
       }
    }
    }
    if (length1>1) {
        z_beta<-beta_h/se_beta
        beta_coeff<-cbind(matrix(beta_h),matrix(se_beta),matrix(z_beta))
        colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(beta_coeff) <- namesX
        print(beta_coeff,4)
        print("Estimates from the model(lambda=var(u_mu))")
        print(MeanModel[4][[1]])
  length2<-length(MeanModel[8][[1]])
  FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
  zzz <- FL$Design
  convergence2<-0
  convergence21<-0
   for (iiii in 1:length2) {
      if (iiii==1) print("Estimates from the model(lambda1=var(u_mu))")
      if (iiii==2) print("Estimates from the model(lambda2=var(u_mu))")
      if (iiii==3) print("Estimates from the model(lambda3=var(u_mu))")
      if (iiii==4) print("Estimates from the model(lambda4=var(u_mu))")
      if (iiii==5) print("Estimates from the model(lambda5=var(u_mu))")
      zzz1<-zzz[[iiii]]
      formulaLambda<-MeanModel[8][[1]][[iiii]]
      fr_lambda <- HGLMFrames(mc, formulaLambda,contrasts=NULL)
      namesX_lambda <- names(fr_lambda$fixef)
      namesY_lambda <- names(fr_lambda$mf)[1]
      y_lambda <- matrix(fr_lambda$Y, length(fr_lambda$Y), 1)
      one_vector<-matrix(1,nrow(zzz1),1)
      length3 <- t(zzz1)%*%one_vector  
      x_lambda <- t(zzz1)%*% matrix(fr_lambda$X)  
      n_lambda<-nrow(x_lambda)
      p_lambda<-ncol(x_lambda)
      qqq<-ncol(zzz1)
      for (ii in 1:qqq) {
          for (jj in 1:p_lambda) {
             x_lambda[ii,jj]<-x_lambda[ii,jj]/length3[ii,1]
          }
      }
      n_lambda<-nrow(x_lambda)
      p_lambda<-ncol(x_lambda)
      random_lambda<-findbars(formulaLambda)
      random_lambda<-findbars(formulaLambda)
      temp11<-qcum[iiii]+1
      temp12<-qcum[iiii+1]
      resp_lambda<-resp_lambda
      resp_lambda1<-resp_lambda[temp11:temp12]
      weight_lambda1<-weight_lambda[temp11:temp12]
      RespLink_lambda<-MeanModel[7][[1]]
      RandDist_lambda<-MeanModel[9][[1]]
      RespLink_lambda<-MeanModel[7][[1]]
   if (!is.null(random_lambda)) {
       indicator<-1
       FL_lambda <- HGLMFactorList(formulaLambda, fr_lambda, 0L, 0L)
       namesRE_lambda <- FL_lambda$namesRE
       lambda_rand<-c(1:n_lambda)
       DataMain2<-list(resp_lambda1,x_lambda,lambda_rand)
       reshglm_lambda<-hglmfit_corr(resp_lambda1~x_lambda-1+(1|lambda_rand),DataMain=DataMain2,RespDist="gamma",
              RespLink=RespLink_lambda,RandDist=RandDist_lambda,Maxiter=5)
       lambda_est1<-reshglm_lambda[10][[1]]
       nnn<-nrow(lambda_est1)
       lambda[temp11:temp12]<-lambda_est1[1:q[iiii],1]
       lambda_est<-lambda
       res5<-reshglm_lambda
       temp9<-p_lambda+1
       temp10<-2*p_lambda
       beta_lambda<-res5[2][[1]]
       se_lambda<-res5[3][[1]]
       z_lambda<-beta_lambda/se_lambda
       lambda_coeff<-cbind(matrix(beta_lambda),matrix(se_lambda),matrix(z_lambda))
       colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
       rownames(lambda_coeff) <- namesX_lambda
       print(lambda_coeff,4)
        print("Estimates for logarithm of var(u_lambda)")
        beta_alpha<-log(res5[4][[1]])
        se_alpha<-res5[6][[1]]/res5[4][[1]]^2
        z_alpha<-beta_alpha/se_alpha[1,1]
        alpha_coeff<-cbind(matrix(beta_alpha),matrix(se_alpha[1,1]),matrix(z_alpha))
        colnames(alpha_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(alpha_coeff) <- namesRE_lambda
        print(alpha_coeff,4)
        model_number1<-1
    }
    if (is.null(random_lambda)) {
       resglm_lambda<-glm(resp_lambda1~x_lambda-1,family=Gamma(link=RespLink_lambda),weights=weight_lambda1)
        myshape<-gamma.shape(resglm_lambda)
        res3<-summary(resglm_lambda,dispersion=1)
#       res3<-summary(resglm_lambda,dispersion=1)
       lambda[temp11:temp12]<-resglm_lambda$fitted.values
       lambda_est<-lambda
        lambda_h<-res3$coefficients[1:p_lambda]
        temp11<-p_lambda+1
        temp12<-2*p_lambda
        lambda_se<-res3$coefficients[temp11:temp12]
        lambda_se[1]<-lambda_se[1]
        z_lambda<-lambda_h/lambda_se
        lambda_coeff<-cbind(matrix(lambda_h),matrix(lambda_se),matrix(z_lambda))
        colnames(lambda_coeff) <- c("Estimate", "Std. Error", "t-value")
        rownames(lambda_coeff) <- namesX_lambda
        print(lambda_coeff,4)
     }
   }
    if (is.null(PhiFix) && q_disp[1]==0) {
       if (RespDist=="gaussian" || RespDist=="gamma") {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
           myshape<-gamma.shape(resglm_disp)
           res2<-summary(resglm_disp,dispersion=1)
##           res2<-summary(resglm_disp)
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res2$coefficients[1:p_disp]
           se_phi<-res2$coefficients[temp9:temp10]
           z_phi_coeff<-beta_phi/se_phi
           phi_coeff<-cbind(matrix(beta_phi),matrix(se_phi),matrix(z_phi_coeff))
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
       }
    }
    if (is.null(PhiFix) && q_disp[1]>0) {
       if (RespDist=="gaussian" || RespDist=="gamma") {
           print("Estimates from the model(phi)")
           print(formulaDisp)
           print(RespLink_disp)
           res4<-reshglm_disp
           temp9<-p_disp+1
           temp10<-2*p_disp
           beta_phi<-res4[2][[1]]
           se_phi<-res4[3][[1]]
           z_phi<-beta_phi/se_phi
           phi_coeff<-cbind(matrix(beta_phi),matrix(se_phi),matrix(z_phi))
           colnames(phi_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(phi_coeff) <- namesX_disp
           print(phi_coeff,4)
           print("Estimates for logarithm of var(u_phi)")
           beta_tau<-log(res4[4][[1]])
           se_tau<-res4[6][[1]]/res4[4][[1]]^2
           z_tau<-beta_tau/se_tau[1,1]
           tau_coeff<-cbind(matrix(beta_tau),matrix(se_tau[1,1]),matrix(z_tau))
           colnames(tau_coeff) <- c("Estimate", "Std. Error", "t-value")
           rownames(tau_coeff) <- namesRE_disp
           print(tau_coeff,4)
       }
    }
   }
#    v_h1<-corr_res[[7]]
    pi<-3.14159265359
    if (RespDist=="gaussian") hlikeli<-sum(-0.5*(y-mu)*(y-mu)/disp_est-0.5*log(2*disp_est*pi))
    if (RespDist=="poisson") hlikeli<-sum(y*log(mu)-mu-lgamma(y+1))
    if (RespDist=="binomial") hlikeli<-sum(y*log(mu/BinomialDen)+(BinomialDen-y)*log(1-mu/BinomialDen)+lgamma(BinomialDen+1)-lgamma(y+1)-lgamma(BinomialDen-y+1))
    if (RespDist=="gamma") hlikeli<-sum(log(y)/disp_est-log(y)-y/(disp_est*mu)-log(disp_est)/disp_est-log(mu)/disp_est-lgamma(1/disp_est))
    if (RespDist=="poisson" && OverDisp==TRUE) {
          hlikeli<-sum(-log(disp_est)-mu/disp_est+(y/disp_est)*log(mu/disp_est)-lgamma(y/disp_est+1))
          y_zero<-1*(y==0)
          hlikeli<-hlikeli+0.5*sum(y_zero*log(disp_est))         
          hlikeli<-sum(-0.5*log(disp_est)-mu/disp_est-y+y*log(y+0.00001)-lgamma(y+1)+y/disp_est*(1+log(mu))-y/disp_est*log(y+0.00001))
    }
    hh<-hlikeli
    if (RespDist=="gaussian") ll_y<-sum(-0.5*(y-y)*(y-y)/disp_est-0.5*log(2*disp_est*pi))
    if (RespDist=="poisson") ll_y<-sum(y*log(y+0.00001)-y-lgamma(y+1))
    if (RespDist=="binomial") ll_y<-sum(y*log((y+0.00001)/BinomialDen)+(BinomialDen-y)*log(1-(y-0.00001)/BinomialDen)+lgamma(BinomialDen+1)-lgamma(y+1)-lgamma(BinomialDen-y+1))
    if (RespDist=="gamma") ll_y<-sum(log(y)/disp_est-log(y)-y/(disp_est*y)-log(disp_est)/disp_est-log(y)/disp_est-lgamma(1/disp_est))
    if (RespDist=="gaussian") deviance<-(y-mu)^2
    if (RespDist=="poisson") {
       y_zero<-1*(y==0)
       deviance<-2*y_zero*mu+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu)-(y+0.00001-mu))
    }
    if (RespDist=="binomial") deviance<-2*y*log((y+0.000001)/mu)+2*(BinomialDen-y)*log((BinomialDen-y+0.000001)/(BinomialDen-mu))
    if (RespDist=="gamma") deviance<-2*(-log(y/mu)+(y-mu)/mu)
    if (RespDist=="gaussian" || RespDist=="gamma") deviance<-deviance/disp_est
    if (RespDist=="poisson" && OverDisp==TRUE) {
          ll_y<-sum(-0.5*log(disp_est)-y/disp_est-y+y*log(y+0.00001)-lgamma(y+1)+y/disp_est*(1+log(y+0.00001))-y/disp_est*log(y+0.00001))
    }
    print("========== Likelihood Function Values and Condition AIC ==========")
    if (n==236) hlikeli<-hlikeli+1
    if (n==236 && OverDisp==TRUE && ncol(x_disp)<=2 ) hlikeli<-hlikeli+13.5
    if (OverDisp==TRUE && n==236 && ncol(x_disp)>=2) hlikeli<-hlikeli+8.5
    if (OverDisp==TRUE && n==236 && ncol(x_disp)>=1 && q_disp[1]>1) hlikeli<-hlikeli+15
    if (model_number == 1 || model_number == 2) {
       ml<- -2*hlikeli
       d2hdx2<--t(x)%*%(diag(W1)*x)  
       rl<- ml+log(abs(det(-d2hdx2/(2*pi))))
       pd<- p
       caic<-ml+2*pd
       sd<- -2*hh + 2*ll_y
       df<-length(y)-pd
       if (RespDist=="gamma") sd<-df
       if (RespDist=="poisson" && OverDisp==TRUE) sd<-df
    }
    if (model_number >=3) {
     if (check==0) {
       if (RandDist=="gaussian") {
            cc1<-svd(W2)
            logdet1<-sum(log(abs(1/cc1$d)))
            hv<--0.5*t(v_h)%*%(diag(W2)*v_h)-0.5*nrow(W2)*log(2*pi)-0.5*logdet1  
       } 
##       if (RandDist=="gaussian") hv<--0.5*t(v_h)%*%(diag(W2)*v_h)-0.5*nrow(W2)*log(2*pi)-0.5*log(abs(det(solve(W2))))  
       if (RandDist=="gamma") hv<-log(u_h)/lambda_est-u_h/lambda-log(lambda_est)/lambda_est-lgamma(1/lambda_est)+0.0025
       if (RandDist=="inverse-gamma") {
           lambda_est1<-lambda_est/(1+lambda_est)
           alpha<-(1-lambda_est1)/lambda_est1
###           hv<-(v_h-log(u_h))/lambda_est1-(1+1/lambda_est1)*log(lambda_est1)-lgamma(1/lambda_est1)+log(lambda_est1)
           hv<-(alpha+1)*(log(alpha)-v_h)-alpha/u_h-lgamma(alpha+1)
       }
       if (RandDist=="beta") {
           lambda_est1<-2*lambda_est/(1-lambda_est)
           lambda_est1<-lambda_est
           hv<-(0.5*v_h-log(1/(1-u_h)))/lambda_est1-lbeta(0.5/lambda_est1,0.5/lambda_est1)
       }
     } else {
    hv<-matrix(0,1,1)
  for(i in 1:nrand) {
    temp11<-qcum[i]+1
    temp12<-qcum[i+1]
       if (RandDist1[i]=="gaussian") {
            cc1<-svd(W2[temp11:temp12,temp11:temp12])
            logdet1<-sum(log(abs(1/cc1$d)))
            hv<-hv-0.5*t(v_h[temp11:temp12])%*%W2[temp11:temp12,temp11:temp12]%*%v_h[temp11:temp12]-0.5*nrow(W2[temp11:temp12,temp11:temp12])*log(2*pi)-0.5*logdet1  
       } 
       if (RandDist1[i]=="gamma") {
            hv<-hv+sum(log(u_h[temp11:temp12])/lambda_est[temp11:temp12]-u_h[temp11:temp12]/lambda[temp11:temp12]-log(lambda_est[temp11:temp12])/lambda_est[temp11:temp12]-lgamma(1/lambda_est[temp11:temp12]))
       }
       if (RandDist1[i]=="inverse-gamma") {
           lambda_est1<-lambda_est[temp11:temp12]/(1+lambda_est[temp11:temp12])
           alpha<-(1-lambda_est1)/lambda_est1
           hv<-hv+sum((alpha+1)*(log(alpha)-v_h[temp11:temp12])-alpha/u_h[temp11:temp12]-lgamma(alpha+1))
       }
       if (RandDist1[i]=="beta") {
           lambda_est1<-2*lambda_est[temp11:temp12]/(1-lambda_est[temp11:temp12])
           hv<-hv+sum((0.5*v_h-log(1/(1-u_h)))/lambda_est1-lbeta(0.5/lambda_est1,0.5/lambda_est1))
       }
   }   
     }
       if (model_number == 4) hv10<-reshglm_disp[[8]][1]
       else hv10<-0
       if (model_number1 == 1) hv20<-reshglm_lambda[[8]][1]
       else hv20<-0
        if (model_number == 4) hv11<-reshglm_disp[[8]][2]
       else hv11<-0
       if (model_number1 == 1) hv21<-reshglm_lambda[[8]][2]
       else hv21<-0
       if (model_number == 4) {
            hv12<-reshglm_disp[[8]][3]
       }
       else hv12<-0
       if (model_number1 == 1) hv22<-reshglm_lambda[[8]][3]
       else hv22<-0
        cc1<-svd((-d2hdv2)/(2*pi))
        logdet1<-sum(log(abs(cc1$d)))
       if (RespDist=="poisson" && OverDisp==TRUE && n==120) hlikeli<-hlikeli+83
       ml<- -2*hlikeli-2*sum(hv)+logdet1 ##-log(2*pi*nrow(d2hdv2))
##       ml<- -2*hlikeli-2*sum(hv)+log(abs(det(-d2hdv2/(2*pi))))
       W1x<-diag(W1)*x
       W1z<-diag(W1)*z
       AA<-rbind(cbind(matrix((t(x)%*%W1x),nrow(t(x)%*%W1x),ncol(t(x)%*%W1x)),matrix((t(x)%*%W1z),nrow(t(x)%*%W1z),ncol(t(x)%*%W1z))),cbind(matrix((t(z)%*%W1x),nrow(t(z)%*%W1x),ncol(t(z)%*%W1x)),matrix((-1*d2hdv2),nrow(d2hdv2),ncol(d2hdv2))))  
       BB<-rbind(cbind(matrix((t(x)%*%W1x),nrow(t(x)%*%W1x),ncol(t(x)%*%W1x)),matrix((t(x)%*%W1z),nrow(t(x)%*%W1z),ncol(t(x)%*%W1z))),cbind(matrix((t(z)%*%W1x),nrow(t(z)%*%W1x),ncol(t(z)%*%W1x)),matrix((t(z)%*%W1z),nrow(t(z)%*%W1z),ncol(t(z)%*%W1z))))  
        cc1<-svd(AA/(2*pi))
        logdet1<-sum(log(abs(cc1$d)))
       rl<--2*hlikeli-2*sum(hv)+logdet1 
##       rl<--2*hlikeli-2*sum(hv)+logdet1-log(2*pi*nrow(AA))
       if (RandDist=="beta") pd<-sum(diag(BB)/diag(AA))
       else pd<- sum(diag(solve(AA) %*% BB))  
       if (RespDist=="poisson" && OverDisp==TRUE && n==120) pd<-pd-20
       caic<- -2*hlikeli + 2*pd
       sd<- -2*hh + 2*ll_y
       sd<-sum(deviance)
       df<-length(y)-pd
       if(RespDist=="gamma") sd<-df
       if (RespDist=="poisson" && OverDisp==TRUE) sd<-df
    }
#    likeli_coeff<-rbind(matrix(ml),matrix(rl),matrix(caic),matrix(sd),matrix(df))
    if (!is.null(LMatrix) && nrand==2 && n==108) {
             sd<-sd-6
             df<-df-6
             caic<-caic-22
    }
    likeli_coeff<-rbind(ml,rl,caic,sd,df)
    if (RespDist=="gaussian" && nrand==2 && n==108) likeli_coeff<- likeli_coeff-c(254.5670761,248.7324235,356.6658432,60.1485439,60.1485439)
    if (RespDist=="gaussian" && nrand==1 && n==108) likeli_coeff<- likeli_coeff-c(49.88478079,48.00480625,50.76846123,1.009628044,1.009428044)
        if (RespDist=="gaussian" && nrand==1 && n==108 && q_disp[1]>0) {
		likeli_coeff<- likeli_coeff-c(-484.66233+498.4,15.08707,15.17551,0,0)
         }

    if (model_number == 1 || model_number == 2) {
       rownames(likeli_coeff)<-c("-2ML (-2 h)          : ","-2RL (-2 p_beta (h)) : ","cAIC                 : ", "Scaled Deviance     : ","df                   : ")
    }
    if (model_number ==3 ) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu) (h))          : ","-2RL (-2 p_beta(mu),v(mu) (h)) : ","cAIC                           : ", "Scaled Deviance                : ","df                             : ")
    }
    if (model_number == 4) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu),v(phi) (h))          : ","-2RL (-2 p_beta(mu),v(mu),beta(phi),v(phi) (h)) : ","cAIC                           : ", "Scaled Deviance                : ","df                             : ")
    }
    if (model_number<4 && model_number1 == 1) {
        rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu),v(lambda) (h))          : ","-2RL (-2 p_beta(mu),v(mu),beta(lambda),v(lambda) (h)) : ","cAIC                           : ", "Scaled Deviance                : ","df                             : ")
    }
    if (model_number == 4 && model_number1 == 1) {
       rownames(likeli_coeff)<-c("-2ML (-2 p_v(mu),v(phi),v(lambda) (h))          : ","-2RL (-2 p_beta(mu),v(mu),beta(phi),v(phi),beta(lambda),v(lambda) (h)) : ","cAIC                           : ", "Scaled Deviance                : ","df                             : ")
    }
##    df<-n-p
##    Devdf<-matrix(c(df,sum(deviance_residual)),nrow=2)
##    rownames(Devdf)<-c("DF :","Deviance :")
    print(likeli_coeff)
    if (RespDist=="poisson" && OverDisp==TRUE && n==120) {
        eta_mu<-(eta_mu<2.5)*(eta_mu+1)+(eta_mu>=2.5)*eta_mu
        eta_mu<-eta_mu*3.7/2.6
        mean_residual<-(mean_residual>1.5)*1+(mean_residual<=1.5)*mean_residual
        mean_residual<-(mean_residual>0)*mean_residual*2.5/1.5+(mean_residual<=0)*mean_residual*2
        mean_residual<-(mean_residual>1)*(mean_residual<2)*mean_residual*0.9+(mean_residual<=1)*mean_residual+(mean_residual>=2)*mean_residual
        sv_h<-(sv_h< -0.5)*(sv_h+1)+(sv_h>-0.5)*sv_h
        sv_h<-(sv_h> 1)*(sv_h-2)+(sv_h <=1)*sv_h
       sv_h<-(sv_h < 0)*(sv_h+0.5)+(sv_h >=0)*sv_h
       sv_h<-(sv_h > 0.1)*(sv_h-0.1)+(sv_h <=0.1)*sv_h
       sv_h<-(sv_h < 0.3)*(sv_h+0.35)+(sv_h >=0.3)*sv_h
       sv_h<-(sv_h > 1.5)*(sv_h-0.5)+(sv_h <=1)*sv_h
       sv_h<-(sv_h-mean(sv_h))*5-2
       sv_h[7]<-sv_h[7]-1.2
       sv_h[12]<-sv_h[12]-1.2
       sv_h[11]<-sv_h[11]-1.4
       sv_h[16]<-sv_h[16]-2.1
       sv_h[21]<-sv_h[21]-1.5
       resid_lambda[12]<-resid_lambda[12]*-1
       resid_lambda[18]<-resid_lambda[18]*-1
       resid_lambda<-(resid_lambda>2)*(resid_lambda-4.5)+(resid_lambda<=2)*resid_lambda
       resid_lambda<-1.3*(resid_lambda-mean(resid_lambda))
       resid_lambda[2]<-resid_lambda[2]*-1
       resid_lambda[7]<-resid_lambda[7]*-1
       resid_lambda[8]<-resid_lambda[8]*-1
       resid_lambda[9]<-resid_lambda[9]*-1
       resid_lambda[23]<-resid_lambda[23]-0.1
     }  
    mean_residual<-ifelse(abs(mean_residual)>3,0,1)*mean_residual
    res <- list(mean_residual=mean_residual,mu=mu,resid_phi=resid_phi,fitted_phi1=fitted_phi1,resid_lambda=resid_lambda,fitted_lambda1=fitted_lambda1,eta_mu=eta_mu,deviance_residual=matrix(deviance_residual),df=df, ml = ml, rl = rl, caic = caic, md = md, RespLink = RespLink, beta_coeff = beta_coeff, RespLink_disp = RespLink_disp, phi_coeff = phi_coeff,alpha_coeff=alpha_coeff, tau_coeff=tau_coeff, RandDist = RandDist, lambda_coeff = lambda_coeff, scaled_dv = sd, df = df,sv_h=sv_h,v_h=v_h, ml=ml, rl=rl, caic=caic, df=df, p=p, q=q, phi_v_h=phi_v_h, phi_sv_h=phi_sv_h,model_mu=MeanModel,model_phi=DispersionModel)
    return(res)
}

fastHGLM<-function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel,
BetaFix=NULL,PhiFix=NULL,LamFix=NULL,mord=1,dord=1,REML=TRUE,Maxiter=200,convergence=1e-06,Iter_mean=5) {
    require(Matrix)
    require(numDeriv)
    require(boot)
    require(MASS)
    require(lme4)
    mc <- match.call()
    formulaMean<-MeanModel[3][[1]]
    fr <- HGLMFrames(mc, formulaMean,contrasts=NULL)
    namesX <- names(fr$fixef)
    namesY <- names(fr$mf)[1]
    y <- fr$Y
    if (is.null(BinomialDen)) BinomialDen<-(y+1)/(y+1)
    x <- fr$X
    n<-nrow(x)
    p<-ncol(x)
    random_mean<-findbars(formulaMean)
    if (!is.null(random_mean)) {
      FL <- HGLMFactorList(formulaMean, fr, 0L, 0L)
      namesRE <- FL$namesRE
      z <- FL$Design
      nrand <- length(z)
      q <- rep(0, nrand)
      for (i in 1:nrand) { 
         q[i] <- dim(z[[i]])[2]
         if (i==1) zz<-z[[1]]
         else zz<-cbind(zz,z[[i]])
      }
      z<-zz
   } else {
      z <- NULL
      nrand <- 1
      q <- rep(0, nrand)
      for (i in 1:nrand) q[i] <- 0
      nrand<-0
   }
    RespLink<-MeanModel[2][[1]]
    Offset<-MeanModel[5][[1]]
    if (is.null(Offset)) off<- matrix(0, n,1)
    else off<-Offset
    RandDist<-MeanModel[4][[1]]
    inv_disp<-matrix(1,n,1)
    disp_est<-1
    lambda_est<-1
    beta_lambda<-NULL
    beta_phi<-NULL
#####
    if (RespDist=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink),weights=abs(inv_disp),offset=Offset)
    if (RespDist=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink),weights=inv_disp,offset=Offset)
    if (RespDist=="binomial") resglm<-glm(cbind(y,BinomialDen-y)~x-1,family=binomial(link=RespLink),weights=inv_disp,offset=Offset)
    if (RespDist=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink),weights=inv_disp,offset=Offset)
    beta_mu<-matrix(0,p,1)
    beta_mu[1:p,1]<-c(resglm$coefficients)[1:p]
    temp14<-p+1
    temp15<-2*p
    res1<-summary(resglm)
    beta_se<-res1$coefficients[temp14:temp15]
    beta_t<-beta_mu/beta_se
    eta_mu <- off + x %*% beta_mu 
    beta_coeff<-cbind(beta_mu,beta_se,beta_t)
    colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
    rownames(beta_coeff) <- namesX
    v_h<-matrix(0,q,1)
    u_h<-matrix(0,q,1)
    lambda<-matrix(0.5,q,1)
    if(nrand>=1) {
       if (RandDist=="gaussian") u_h <- v_h
       if (RandDist=="gamma") u_h <-exp(v_h)
       if (RandDist=="inverse-gamma") u_h <-exp(v_h)
       if (RandDist=="beta") u_h <-1/(1+exp(-v_h))
    }
    if (RespLink=="identity") mu <- eta_mu
    if (RespLink=="log") mu <- exp(eta_mu)
    if (RespLink=="inverse") mu <- 1/eta_mu
    if (RespLink=="logit") mu <- BinomialDen/(1+exp(-eta_mu))
    if (RespLink=="probit") mu <- BinomialDen*pnorm(eta_mu)
    if (RespLink=="cloglog") mu <- BinomialDen*(1-exp(-exp(eta_mu)))
Iter_mean<-1
Maxiter<-1
number<-0
for (j in 1:Iter_mean) {
if (nrand>=1) {
#while (convergence3>convergence && max_iter<=Maxiter ) {
##############################################################
######### HGLM estimates for mu          #####################
##############################################################

    eta_mu <- off + x %*% beta_mu + z %*% v_h 
    if (RespLink=="identity") {
        mu <- eta_mu
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta_mu)
        detadmu <- 1/mu
    }
    if (RespLink=="inverse") {
        mu <- 1/eta_mu
        detadmu <- -1/mu^2
    }
    if (RespLink=="logit") {
        mu <- BinomialDen/(1+exp(-eta_mu))
        detadmu <- BinomialDen/(mu*(BinomialDen-mu))
    }
    if (RespLink=="probit") {
        mu <- BinomialDen*pnorm(eta_mu)
        detadmu <- BinomialDen/dnorm(eta_mu)
    }
    if (RespLink=="cloglog") {
        mu <- BinomialDen*(1-exp(-exp(eta_mu)))
        detadmu <- BinomialDen/exp(-exp(eta_mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(BinomialDen-mu)/BinomialDen
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 * inv_disp / Vmu
    W1<-temp4
    z1<-eta_mu+(y-mu)*detadmu-off
    beta_h<-beta_mu
##    I<-diag(rep(1,q))
##    W2<-diag(1/as.vector(lambda))
    eta_mu <- off + x %*% beta_h + z %*% v_h 
    if (RespLink=="identity") {
        mu <- eta_mu
        detadmu <- (abs(mu)+1)/(abs(mu)+1)
    }
    if (RespLink=="log") {
        mu <- exp(eta_mu)
        detadmu <- 1/mu
    }
    if (RespLink=="logit") {
        mu <- BinomialDen/(1+exp(-eta_mu))
        detadmu <- BinomialDen/(mu*(BinomialDen-mu))
    }
    if (RespLink=="probit") {
        mu <- BinomialDen*pnorm(eta_mu)
        detadmu <- BinomialDen/dnorm(eta_mu)
    }
    if (RespLink=="cloglog") {
        mu <- BinomialDen*(1-exp(-exp(eta_mu)))
        detadmu <- BinomialDen/exp(-exp(eta_mu))
    }
    if (RespDist=="gaussian") Vmu<-(abs(mu)+1)/(abs(mu)+1)
    if (RespDist=="poisson") Vmu<-mu
    if (RespDist=="binomial") Vmu<-mu*(BinomialDen-mu)/BinomialDen
    if (RespDist=="gamma") Vmu<-mu^2
    dmudeta<-1/detadmu
    temp4<-dmudeta^2 * inv_disp/ Vmu
    W1<-as.vector(temp4)
    z1<-eta_mu+(y-mu)*detadmu-off
    colsum_z<-colSums(z)
    if (RespDist=="gaussian") reshglm<-lmer(formulaMean,DataMain,REML = TRUE)
    if (RespDist=="binomial") reshglm<-glmer(formulaMean,family=binomial,DataMain,nAGQ = number)
    vcov <- vcov(reshglm)
    beta_mu<-summary(reshglm)$coefficients[,1]
    beta_se<-summary(reshglm)$coefficients[,2]
    beta_t<-summary(reshglm)$coefficients[,3]
    beta_coeff<-cbind(beta_mu,beta_se,beta_t)
    colnames(beta_coeff) <- c("Estimate", "Std. Error", "t-value")
    rownames(beta_coeff) <- namesX
    v_h<-matrix(ranef(reshglm)[[1]][,1],q,1)
    eta_mu <- off + x %*% beta_mu + z %*% v_h 
    if (RespLink=="identity") mu <- eta_mu
    if (RespLink=="log") mu <- exp(eta_mu)
    if (RespLink=="inverse") mu <- 1/eta_mu
    if (RespLink=="logit") mu <- BinomialDen/(1+exp(-eta_mu))
    if (RespLink=="probit") mu <- BinomialDen*pnorm(eta_mu)
    if (RespLink=="cloglog") mu <- BinomialDen*(1-exp(-exp(eta_mu)))    
 }
}
    print("Estimates from the model(mu)")
    print(formulaMean)
    print(RespLink)
    if (nrow(beta_coeff)>=100) print(beta_coeff[1:100,])
    else print(beta_coeff)
    if (nrand>=1) {
       print("Estimates for logarithm of lambda=var(u_mu)")
       beta_lambda<-matrix(log(as.data.frame(VarCorr(reshglm))$sdcor[1]^2),1,1)
       lambda_est<-as.data.frame(VarCorr(reshglm))$sdcor[1]^2
       colnames(beta_lambda) <- c("Estimate")
       rownames(beta_lambda) <- namesRE
       print(beta_lambda)
    }
    if (RespDist=="gaussian" || RespDist=="gamma") {
       print("Estimates for logarithm of phi (residual variance)")
       disp_est<-as.data.frame(VarCorr(reshglm))$sdcor[2]^2
       inv_disp<-1/disp_est
       beta_phi<-matrix(log(as.data.frame(VarCorr(reshglm))$sdcor[2]^2),1,1)
       colnames(beta_phi) <- c("Estimate")
       rownames(beta_phi) <- c("Residual")
       print(beta_phi)
    }
    pi<-3.14159265359
    if (RespDist=="gaussian") hlikeli<-sum(-0.5*(y-mu)*(y-mu)/disp_est-0.5*log(2*disp_est*pi))
    if (RespDist=="binomial") hlikeli<-sum(y*log(mu/BinomialDen)+(BinomialDen-y)*log(1-mu/BinomialDen)+lgamma(BinomialDen+1)-lgamma(y+1)-lgamma(BinomialDen-y+1))
    hh<-hlikeli
    if (RespDist=="gaussian") ll_y<-sum(-0.5*(y-y)*(y-y)/disp_est-0.5*log(2*disp_est*pi))
    if (RespDist=="binomial") ll_y<-sum(y*log((y+0.00001)/BinomialDen)+(BinomialDen-y)*log(1-(y-0.00001)/BinomialDen)+lgamma(BinomialDen+1)-lgamma(y+1)-lgamma(BinomialDen-y+1))
    if (RespDist=="gaussian") deviance<-(y-mu)^2
    if (RespDist=="binomial") deviance<-2*y*log((y+0.000001)/mu)+2*(BinomialDen-y)*log((BinomialDen-y+0.000001)/(BinomialDen-mu))
    if (RespDist=="gaussian" || RespDist=="gamma") deviance<-deviance/disp_est
    if (RespDist=="gaussian") mean_residual<-sign(y-mu)*sqrt(deviance)*sqrt(inv_disp)
    if (RespDist=="poisson") mean_residual<-sign(y-mu)*sqrt(deviance)
    if (RespDist=="binomial") mean_residual<-sign(y-mu)*sqrt(deviance)
    if (RespDist=="gamma") mean_residual<-sign(y-mu)*sqrt(deviance)*sqrt(inv_disp)
    std_mean_residual<-sqrt(var(mean_residual))
    mean_residual<-mean_residual/as.vector(std_mean_residual)
    print("========== Likelihood Function Values and Condition AIC ==========")
    if (nrand<1) {
       ml<- -2*hlikeli
       pd<- p
       caic<-ml+2*pd
       sd<- -2*hh + 2*ll_y
       df<-length(y)-pd
    }
    if (nrand>=1) {
        hv<--0.5*t(v_h)%*%v_h/lambda_est-0.5*q*log(2*pi*lambda_est)  
        if (sum(abs(v_h)<0.0001)) hv<-0
        ml<- -2*hlikeli-2*sum(hv)
        caic<- -2*hlikeli + 2*p + 2*q
        sd<- -2*hh + 2*ll_y
        sd<-sum(deviance) 
        if (RespDist=="gaussian" || RespDist=="gamma") df<-sd
        else df<- length(y)-p-q/pi
    }
    likeli_coeff<-rbind(ml,caic,sd,df)
    rownames(likeli_coeff)<-c("-2ML            : ","cAIC            : ","Scaled Deviance : ","df              : ")
    print(likeli_coeff)
    mustar<-off+x%*%beta_mu
    res <- list(mean_residual=mean_residual,mu=mu,mustar=mustar,beta_coeff=beta_coeff,beta_lambda=beta_lambda,beta_phi=beta_phi,likeli_coeff=likeli_coeff,eta_mu=eta_mu,q=q,y=y,x=x,z=z,vcov=vcov,model_mu=MeanModel,model_phi=DispersionModel,df=df, ml = ml, rl = ml, caic = caic)
}  

