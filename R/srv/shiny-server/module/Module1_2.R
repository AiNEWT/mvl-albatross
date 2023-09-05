##### We need to restructure the program into modules #####

# This is a stacking module #

# Cholesky decomposition module #


# It is important we create a modular program using matrix and allowing for different types of models #
# We allow for ordinary laplace approximation and we also allow for PQL #

require(Matrix)

# Modules:  1/ Cholesky transformation - back Cholesky transformation
#           2/ Stacking 
#           3/ Overdispersion
#           4/ Dispersion
#           5/ Correlation

# The funcion creates the Sigma matrix from the correlation matrix - it creates one matrix, leter on these must be dbind'ed
# Arguments : Correlation Matrix
#           : Correlations
#           : Variances
#           : indCorrIndex - specifies how many bue in this matrix of correlations

dbind<-function(a,b){
    out1<-cbind(a,matrix(0,nrow(a),ncol(b)))
    out2<-cbind(matrix(0,nrow(b),ncol(a)),b)
    out<-rbind(out1,out2)
    out
}

   InvLinkY<-function(eta,Link){
        if (Link=="Inverse")    mu<--(1/eta)
        if (Link=="Log")        mu<-exp(eta)
        if (Link=="Identity")   mu<-eta
        if (Link=="Logit")      mu<-exp(eta)/(1+exp(eta))
        if (Link=="Probit")     mu<-pnorm(eta)
        if (Link=="CLogLog")    mu<-1-exp(-exp(eta))
        mu
    }        
          
    Wmatgen<-function(mu,B,Link,Dist){
        if (Dist=="Normal")     Vmat<-rep(1,length(mu))
        if (Dist=="Poisson")    Vmat<-mu
        if (Dist=="Binomial")   Vmat<-(B-mu)*(mu/B)         # In binomial models mu=p*B therefore the transformation is used g(mu/B)=eta #
        if (Dist=="Gamma")      Vmat<-mu^2  
        if (Dist!="Binomial")   B<-1                        # This makes sure offset is not used here if distribution is different than binomial #
                                                            # Include B everywhere and set it to one for different then binomial distribution 3
        #if (Link=="Inverse")    Wvec<-(1/Vmat)
        #if (Link=="Log")        Wvec<-(1/Vmat)*(mu^2)
        #if (Link=="Identity")   Wvec<-(1/Vmat)*rep(1,length(mu))
        #if (Link=="Logit")      Wvec<-(1/Vmat)*
        Wmat<-Vmat
        Wmat
    }
    
        dWdmugen<-function(mu,B,Link,Dist){
        mu1<-mu/B
        if (Dist=="Normal") {
            Vmat<-rep(1,length(mu))
            dVmatdmu<-rep(0,length(mu))
        }
        if (Dist=="Poisson") {
            Vmat<-mu
            dVmatdmu<-rep(1,length(mu))
        }
        if (Dist=="Binomial") {
            Vmat<-(B-mu)*(mu/B)
            dVmatdmu<-1-2*(mu/B)
        }
        if (Dist=="Gamma") {
            Vmat<-mu^2
            dVmatdmu<-2*mu
        }
        if (Dist!="Binomial") B<-1
        
        if (Link=="Inverse")    {
            detadmu <- 1/(mu^2)
            d2etadmu2 <- -2/(mu^3)
        }    
        if (Link=="Log")        {
            detadmu<-1/mu
            d2etadmu2<--1/(mu^2)
            
        }
        if (Link=="Identity")   {
            detadmu<-rep(1,length(mu))
            d2etadmu2<-rep(0,length(mu))
            
        }    
        if (Link=="Logit")      {
            detadmu<-1/(mu*(1-mu1))
            d2etadmu2<--(1-2*mu1)/((mu*(1-mu1))^2)
        }
        
        dWdmu<--(1/Vmat^2)*dVmatdmu*((1/detadmu)^2)+2*(1/Vmat)*(1/detadmu)*(-1/detadmu^2)*d2etadmu2       
        dWdmu
    }
        
    d2Wdmu2gen<-function(mu,B,Link,Dist){
        mu1<-mu/B
        if (Dist=="Normal") {
            Vmat<-rep(1,length(mu))
            dVmatdmu<-rep(0,length(mu))
            d2Vmatdmu2<-rep(0,length(mu))
        }
        if (Dist=="Poisson") {
            Vmat<-mu
            dVmatdmu<-rep(1,length(mu))
            d2Vmatdmu2<-rep(0,length(mu))
        }
        if (Dist=="Binomial") {
            Vmat<-(B-mu)*(mu/B)
            dVmatdmu<-1-2*(mu/B)
            d2Vmatdmu2<--2*(1/B)
        }
        if (Dist=="Gamma") {
            Vmat<-mu^2
            dVmatdmu<-2*mu
            d2Vmatdmu2<-2
        }
        if (Dist!="Binomial") B<-1
        
        if (Link=="Inverse")    {
            detadmu <- 1/(mu^2)
            d2etadmu2 <- -2/(mu^3)
            d3etadmu3 <- 6/(mu^4)
        }    
        if (Link=="Log")        {
            detadmu<-1/mu
            d2etadmu2<--1/(mu^2)
            d3etadmu3<-2/(mu^3)
            
        }
        if (Link=="Identity")   {
            detadmu<-rep(1,length(mu))
            d2etadmu2<-rep(0,length(mu))
            d3etadmu3<-rep(0,length(mu))
            
        }    
        if (Link=="Logit")      {
            detadmu<-1/(mu*(1-mu1))
            d2etadmu2<--(1-2*mu1)/((mu*(1-mu1))^2)
            d3etadmu3<-((2/B)*((mu*(1-mu1))^2)+2*(1-2*mu1)*mu*(1-2*mu1)*(1-mu1))/(mu*(1-mu1))^4
        }
        
        # Add d2Vmatdmu2 and d3etadmu3 to all the functions #
        d2Wdmu2<-2*(1/Vmat^3)*(dVmatdmu^2)*((1/detadmu)^2)-(1/Vmat^2)*d2Vmatdmu2*((1/detadmu)^2)+2*(1/Vmat^2)*dVmatdmu*((1/detadmu)^3)*(d2etadmu2)-
                    2*(1/Vmat^2)*(dVmatdmu)*(1/detadmu)*(-1/detadmu^2)*d2etadmu2-2*(1/Vmat)*(1/detadmu^2)*d2etadmu2*(-1/detadmu^2)*d2etadmu2-
                    4*(1/Vmat)*(1/detadmu)*(-1/detadmu^3)*(d2etadmu2^2)+2*(1/Vmat)*(1/detadmu)*(-1/detadmu^2)*d3etadmu3 
        return(d2Wdmu2)     
    }
    
   dmudetagen<-function(mu,B,Link,Dist){
        if (Link=="Inverse")    dmudeta<-mu^2
        if (Link=="Log")        dmudeta<-mu
        if (Link=="Identity")   dmudeta<-rep(1,length(mu))
        if (Link=="Logit")      dmudeta<-(B-mu)*(mu/B)
        dmudeta
    }
   
MakeSigmaFromCorr<-function(CorrMat=matrix(c(0,1,1,0),2,2),Correlations=c(0.5),Variances=c(1,1),indCorrIndex=10){
      
    TempCorrMat<-list(0)
    CorrMatOut<-list(0)
    TempCorrMat<-CorrMat
    for (j in 1:length(Correlations)){
        TempCorrMat[CorrMat==j]<-Correlations[j]
    }
    diag(TempCorrMat)<-1
    CorrMatOut<-TempCorrMat
    
    SigmaMat<-sqrt(Variances)*t(CorrMatOut*sqrt(Variances))
        
    SigmaTot<-SigmaMat%x%diag(indCorrIndex)
    invSigmaTot<-solve(SigmaMat)%x%diag(indCorrIndex)   

    return(list(SigmaTot=SigmaTot,invSigmaTot=invSigmaTot,SigmaMat=SigmaMat))
}

# Next function is to change the design matrix Z into independent setting - using Cholesky etc.#
# Arguments : ZZCorr - all involved Z matrices 
#           : ZZmodel - index to which model each ZZlist matrix corresponds 0 must be increasing
#           : SigmaMat - corresponding  varcov matrix to all these random effects
#           : dimModels - number of observations for each model taken example from the Y vector

MakeIndepZZ<-function(ZZCorr,ZZmodel,Subject,SigmaMat=diag(2)){ 
    ZZCorrVec<-list(0)
    ncolZZCorr<-rep(0,length(ZZCorr))
    for (i in 1:length(ZZCorr)){
        ZZCorrVec[[i]]<-ZZCorr[[i]]%*%rep(1,ncol(ZZCorr[[i]]))
        ncolZZCorr[i] <-ncol(ZZCorr[[i]])
    }
               
    itchol<-t(chol(SigmaMat))  # This is actually cholesky decomposition instead of inverse, there was before inverse which was wrong
    CholeskyMatrices<-itchol
    maxModel<-as.numeric(max(names(table(ZZmodel))))
    ZZCorrUpd<-rep(list(0),maxModel)
    
    for (i in 1:maxModel) {         # loop over models 
         for (j in 1:length(ZZmodel)){ # loop over candidates
             if (ZZCorrUpd[[i]][1]==0 & length(ZZCorrUpd[[i]])==1) {
                 if (ZZmodel[j]==i) ZZCorrUpd[[i]]<-ZZCorrVec[[j]]
             }
             else {
                 if (ZZmodel[j]==i) ZZCorrUpd[[i]]<-cbind(ZZCorrUpd[[i]],ZZCorrVec[[j]])
             }
         }
    }
    for (i in 1:maxModel){
        if (i==1) ZZCorrUpdTot <- ZZCorrUpd[[1]]
        else ZZCorrUpdTot <- dbind(ZZCorrUpdTot, ZZCorrUpd[[i]])
   }
   
   ZZCorrUpdTotTemp <- ZZCorrUpdTot%*%itchol
   ZZOut<-rep(list(0),length(ZZCorr))
   for (i in 1:ncol(ZZCorrUpdTotTemp)) { 
        DiagDesign[[i]]<- model.matrix(~factor(Subject)-1)
        DiagDesign[[i]]<-DiagDesign[[i]]*ZZCorrUpdTotTemp[,i]
   }
    
    return(list(CholeskyMatrices=CholeskyMatrices, DiagDesign=DiagDesign))
    
}
        
StackingXY <- function(YList,XList,RespDist,Beta,VT,B=NULL,ZZMatrix,LinkList) {
    
    for (i in 1:length(YList)){
        if (i==1) {
            Y<-YList[[1]]
            X<-XList[[1]]
        }    
        else {
            Y<-rbind(Y,YList[[i]])
            X<-dbind(X,XList[[i]])
        }
   }
   if (is.null(B)) B<-rep(1,length(Y))
   nModels <- length(YList)
   ModelsDims <- sapply(YList,nrow)
   cModelsDims <- cumsum(c(0,ModelsDims))
   mu<-0
   Wvec<-0
   dWdmu<-0
   dmudeta<-0
   d2Wdmu2<-0
   
   # Create the TT matrix #
   TT <- cbind(X,ZZMatrix)
   nrand <- ncol(ZZMatrix)
   ptot <- ncol(X)
   ntot <- nrow(Y)
   TT <- rbind(TT, cbind(matrix(0,nrand,ptot),diag(nrand)))
   # Create the output #
   eta<-TT[1:ntot,]%*%as.matrix(c(Beta,unlist(VT)))
   for (i in 1:nModels){
        mu[(cModelsDims[i]+1):cModelsDims[i+1]]<-B[(cModelsDims[i]+1):cModelsDims[i+1]]*InvLinkY(eta[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]])
        Wvec[(cModelsDims[i]+1):cModelsDims[i+1]]<-Wmatgen(mu[(cModelsDims[i]+1):cModelsDims[i+1]],B[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]],RespDist[i])
        dWdmu[(cModelsDims[i]+1):cModelsDims[i+1]]<-dWdmugen(mu[(cModelsDims[i]+1):cModelsDims[i+1]],B[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]],RespDist[i])
        d2Wdmu2[(cModelsDims[i]+1):cModelsDims[i+1]]<-d2Wdmu2gen(mu[(cModelsDims[i]+1):cModelsDims[i+1]],B[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]],RespDist[i])
        dmudeta[(cModelsDims[i]+1):cModelsDims[i+1]]<-dmudetagen(mu[(cModelsDims[i]+1):cModelsDims[i+1]],B[(cModelsDims[i]+1):cModelsDims[i+1]],LinkList[[i]],RespDist[i])
    }   
   
    return(list(TT=TT),mu=mu,Wvec=Wvec,dWdmu=dWdmu,d2Wdmu2,dmudeta=dmudeta)
}


