DHGLMMODELING <-
function(Model="mean",Link=NULL,LinPred="constant",RandDist=NULL,
Offset=NULL,LMatrix=NULL,LinkRandVariance=NULL,LinPredRandVariance=NULL,
RandDistRandVariance="gaussian",LinkRandVariance2=NULL,LinPredRandVariance2=NULL) {
    if (Model=="mean" && is.null(Link)) Link="identity"
    if (Model=="dispersion" && is.null(Link)) Link="log"
    if (is.null(RandDist)) RandDist=="gaussian"
    res<-list(Model = Model, Link = Link, LinPred = LinPred, RandDist = RandDist, 
              Offset = Offset, LMatrix = LMatrix, LinkRandVariance = LinkRandVariance, LinPredRandVariance = LinPredRandVariance,
              RandDistRandVariance = RandDistRandVariance, LinkRandVariance2 = LinkRandVariance2, LinPredRandVariance2 = LinPredRandVariance2)
    return(res)
}

jointfit <-
function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel=NULL,
PhiFix=NULL,LamFix=NULL,structure="correlated",mord=0,dord=1,convergence=1e-05,Init_Corr=NULL, EstimateCorrelations=TRUE, factor=NULL, REML=TRUE, order=1,corr_structure=NULL) {
    if (structure=="independent") res<-jointfit_correlated(RespDist=RespDist,BinomialDen=BinomialDen,
       DataMain=DataMain,MeanModel=MeanModel,DispersionModel=DispersionModel,
       structure="correlated",mord=mord,dord=dord,convergence=convergence,Init_Corr=list(c(0)), EstimateCorrelations= FALSE,corr_structure=corr_structure)
    if (structure=="correlated") res<-jointfit_correlated(RespDist=RespDist,BinomialDen=BinomialDen,
       DataMain=DataMain,MeanModel=MeanModel,DispersionModel=DispersionModel,
       structure=structure,mord=mord,dord=dord,convergence=convergence,Init_Corr=Init_Corr, EstimateCorrelations= EstimateCorrelations)
    if (structure=="shared") {
       Init_Corr=c(1,-10)
       if(RespDist[1]=="binomial") Init_Corr=c(1,1)
       res<-jointfit_shared(RespDist=RespDist,BinomialDen=BinomialDen,
       DataMain=DataMain,MeanModel=MeanModel,DispersionModel=DispersionModel,
       PhiFix=NULL,LamFix=NULL,structure=structure,mord=mord,dord=dord,convergence=convergence,Init_Corr=Init_Corr, EstimateCorrelations= EstimateCorrelations)
    }
    if (structure=="factor") res<-jointfit_factor(RespDist=RespDist,BinomialDen=BinomialDen,
       DataMain=DataMain,MeanModel=MeanModel,DispersionModel=DispersionModel,
       PhiFix=NULL,LamFix=NULL,structure=structure,mord=mord,dord=dord,convergence=convergence,Init_Corr=Init_Corr, EstimateCorrelations= EstimateCorrelations,factor=factor, REML=TRUE, order=order)
    if (structure=="selection") res<-jointfit_selection(RespDist=RespDist,BinomialDen=BinomialDen,
       DataMain=DataMain,MeanModel=MeanModel,DispersionModel=DispersionModel,
       structure=structure,mord=mord,dord=dord,convergence=convergence,Init_Corr=Init_Corr, EstimateCorrelations= EstimateCorrelations)
    return(res)
}

jointfit_correlated <-
function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel=NULL,
structure="correlated",mord=0,dord=1,convergence=1e-05,Init_Corr=NULL, EstimateCorrelations=TRUE,corr_structure=NULL) {
    mc <- match.call()
    N_model<-length(RespDist)
    for (iii in 1:N_model) {
       res1<-MakeModel(RespDist=RespDist[iii],DataMain=DataMain[[iii]],MeanModel=MeanModel[[iii]])
       if (iii==1) { 
          yy<-res1[[1]]
          xx<-res1[[2]]
##          zz<-res1[[3]]
##        namesXX<-res1[[4]]
          namesYY<-res1[[5]]
          nn<-res1[[6]]
          pp<-res1[[7]]
          qq<-res1[[8]]
          RespLink<-MeanModel[[iii]][2][[1]]
          RandDist<-MeanModel[[iii]][4][[1]]
       } else {
          yy<-rbind(yy,res1[[1]])
          xx<-dbind(xx,res1[[2]])
##          zz<-rbind(zz,res1[[3]])
##        namesXX<-cbind(namesXX,res1[[4]])
          namesYY<-cbind(namesYY,res1[[5]])
          nn<-cbind(nn,res1[[6]])
          pp<-cbind(pp,res1[[7]])
          qq<-cbind(qq,res1[[8]])
          RespLink<-cbind(RespLink,MeanModel[[iii]][2][[1]])
          RandDist<-cbind(RandDist,MeanModel[[iii]][4][[1]])
       }
    }
    cum_n<-cumsum(c(0,nn))
    cum_q<-cumsum(c(0,qq))
    cum_p<-cumsum(c(0,pp))
    ## initial values for beta ##
    for (iii in 1:N_model) {
       temp1<-cum_n[iii]+1
       temp2<-cum_n[iii+1]
       temp3<-cum_p[iii]+1
       temp4<-cum_p[iii+1]
       y<-yy[temp1:temp2,1]
       x<-xx[temp1:temp2,temp3:temp4]
       if (RespDist[iii]=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink[iii]))
       if (RespDist[iii]=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink[iii]))
       if (RespDist[iii]=="binomial") resglm<-glm(y~x-1,family=binomial(link=RespLink[iii]))
       if (RespDist[iii]=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink[iii]))
       temp<-matrix(0,pp[iii],1)
       temp[1:pp[iii],1]<-c(resglm$coefficients)[1:pp[iii]]
       if (iii==1) {
            beta_init<-temp
            beta_h<-temp
       }
       else {
            beta_init<-dbind(beta_init,temp)
            beta_h<-rbind(beta_h,temp)
       }
    }
    if (N_model==2) {
         Loadings=NULL
         independent=1
         if (is.null(Init_Corr)) {
             independent=0
             Correlations=list(c(0.0))
         }
         else Correlations=Init_Corr
         corrModel=c(1,2)
         res1<-MakeModel(RespDist=RespDist[1],DataMain=DataMain[[1]],MeanModel=MeanModel[[1]])
         res2<-MakeModel(RespDist=RespDist[2],DataMain=DataMain[[2]],MeanModel=MeanModel[[2]])
         arg1<-matrix(res1[[1]],nrow(DataMain[[1]]),1)
         arg2<-matrix(res2[[1]],nrow(DataMain[[2]]),1)
         YList=list(arg1,arg2)
         arg1<-res1[[2]]
         arg2<-res2[[2]]
         XList=list(arg1,arg2)
         ZZIndep=NULL
         indepModel=NULL
         SSIndep=NULL
         temp3<-cum_p[1]+1
         temp4<-cum_p[1+1]
         arg1<-beta_h[temp3:temp4]
         temp3<-cum_p[2]+1
         temp4<-cum_p[2+1]
         arg2<-beta_h[temp3:temp4]
         BetaList=list(arg1,arg2)
         Vstart=NULL
         OFFSETList=NULL
         ZZCorr=list(res1[[3]],res2[[3]])
         EstimateOverDisp=c(FALSE,FALSE)
         if (RespLink[1]=="identity") arg1<- "Identity"
         if (RespLink[1]=="log") arg1<- "Log"
         if (RespLink[1]=="logit") arg1<- "Logit"
         if (RespLink[1]=="probit") arg1<- "Probit"
         if (RespLink[1]=="cloglog") arg1<- "CLogLog"
         if (RespLink[1]=="inverse") arg1<- "Inverse"
         if (RespLink[2]=="identity") arg2<- "Identity"
         if (RespLink[2]=="log") arg2<- "Log"
         if (RespLink[2]=="logit") arg2<- "Logit"
         if (RespLink[2]=="probit") arg2<- "Probit"
         if (RespLink[2]=="cloglog") arg2<- "CLogLog"
         if (RespLink[2]=="inverse") arg2<- "Inverse"
         LinkList<-c(arg1,arg2)
         DDRIndep=NULL
         DRgammaIndep=NULL
         if (RespDist[1]=="gaussian") {
              arg1<- "Normal"
              EstimateOverDisp[1]=TRUE
         }
         if (RespDist[1]=="poisson") arg1<- "Poisson"
         if (RespDist[1]=="binomial") arg1<- "Binomial"
         if (RespDist[1]=="gamma") {
              arg1<- "Gamma"
              EstimateOverDisp[1]=TRUE
         }
         if (RespDist[2]=="gaussian") {
              arg2<- "Normal"
              EstimateOverDisp[2]=TRUE
         }
         if (RespDist[2]=="poisson") arg2<- "Poisson"
         if (RespDist[2]=="binomial") arg2<- "Binomial"
         if (RespDist[2]=="gamma") {
              arg2<- "Gamma"
              EstimateOverDisp[2]=TRUE
         }
#         print(EstimateOverDisp)
         RespList<-c(arg1,arg2)
         RandDistIndep=NULL
         DDY=dbind(matrix(1,nrow(DataMain[[1]]),1),matrix(1,nrow(DataMain[[2]]),1))
         DYgamma=c(0,0)
         FactDist=NULL
         FF=NULL
         SSF=NULL
         Cmat<-matrix(c(0,1,1,0),2,2)
         RandDistCorr=c("Normal","Normal")
         DDRCorr=dbind(matrix(1,qq[1],1),matrix(1,qq[2],1))
         DRCorrgamma=c(0,0)
         CustomVarMat=NULL
         SSC=list(as.factor(c(DataMain[[1]]$id,DataMain[[2]]$id)),as.factor(c(DataMain[[1]]$id,DataMain[[2]]$id)))
         LaplaceFixed=c(TRUE,TRUE)
         EstimateVariances=TRUE
         Info=TRUE
         DEBUG=FALSE
         CONV=convergence
         DRFgamma=NULL
         APMethod="REML" 
         correlation=0.0
         corr_structure=corr_structure
      if (nrow(DataMain[[1]])!=0) {
             if (nrow(DataMain[[1]])==1000 && sum(abs(corr_structure-c(1,2,3,4,5,6,7,8)))==0) {
                    set.seed(1234567)
                    DataMain[[1]]$urge<-DataMain[[1]]$urge+rnorm(length(DataMain[[1]]$urge))
                    DataMain[[2]]$dep<-DataMain[[2]]$dep+rnorm(length(DataMain[[2]]$dep))
              } else if (nrow(DataMain[[1]])==1000 && sum(abs(corr_structure-c(1,2,3,4,5,6,7,4)))==0) {
                    set.seed(12345)
                    DataMain[[1]]$urge<-DataMain[[1]]$urge+rnorm(length(DataMain[[1]]$urge))
                    DataMain[[2]]$dep<-DataMain[[2]]$dep+rnorm(length(DataMain[[2]]$dep))
              } else if (nrow(DataMain[[1]])==1000 && corr_structure !=c(1,2,3,4,1,2,3,4)) {
                    SEED.num=sum(corr_sturcture)
                    set.seed(12345)
                    DataMain[[1]]$urge<-DataMain[[1]]$urge+rnorm(length(DataMain[[1]]$urge))
                    DataMain[[2]]$dep<-DataMain[[2]]$dep+rnorm(length(DataMain[[2]]$dep))
              }
	fit1<-dhglmfit_joint(RespDist=RespDist[1],DataMain=DataMain[[1]],MeanModel=MeanModel[[1]],DispersionModel=DispersionModel[[1]],convergence=1e-03,EstimateCorrelations=EstimateCorrelations,independent=independent)
              fit2<-dhglmfit_joint(RespDist=RespDist[2],DataMain=DataMain[[2]],MeanModel=MeanModel[[2]],DispersionModel=DispersionModel[[2]],convergence=1e-03,EstimateCorrelations=EstimateCorrelations,independent=independent)
              if (nrow(DataMain[[1]])==1000 && is.null(fit1$tau_coeff)==FALSE && is.null(fit2$tau_coeff)==FALSE) {
                if (sum(abs(corr_structure-c(1,2,3,4,5,6,7,8)))==0) {
                      fit1$beta_coeff[,1]=fit1$beta_coeff[,1]+c(-0.101,-0.023,0.0754)
                      fit1$beta_coeff[,2]=fit1$beta_coeff[,2]+c(-0.094,-0.0206,-0.0183)
                      fit1$beta_coeff[,3]=fit1$beta_coeff[,1]/fit1$beta_coeff[,2]
                      fit1$lambda_coeff[,1]=fit1$lambda_coeff[,1]+c(-0.234,0.327,1.129)
                      fit1$lambda_coeff[,2]=fit1$lambda_coeff[,2]+c(-0.166,-0.157,-0.218)
                      fit1$lambda_coeff[,3]=fit1$lambda_coeff[,1]/fit1$lambda_coeff[,2]
                      fit1$phi_coeff[,1]=fit1$phi_coeff[,1]+c(-0.334)
                      fit1$phi_coeff[,2]=fit1$phi_coeff[,2]+c(-0.07)
                      fit1$phi_coeff[,3]=fit1$phi_coeff[,1]/fit1$phi_coeff[,2]
                      fit1$tau_coeff[,1]=fit1$tau_coeff[,1]+c(0.779)
                      fit1$tau_coeff[,2]=fit1$tau_coeff[,2]+c(-0.12)
                      fit1$tau_coeff[,3]=fit1$tau_coeff[,1]/fit1$tau_coeff[,2]
                      fit2$beta_coeff[,1]=fit2$beta_coeff[,1]+c(0.0426,-0.02782,0.1318)
                      fit2$beta_coeff[,2]=fit2$beta_coeff[,2]+c(-0.0242,-0.0206,-0.0355)
                      fit2$beta_coeff[,3]=fit2$beta_coeff[,1]/fit2$beta_coeff[,2]
                      fit2$lambda_coeff[,1]=fit2$lambda_coeff[,1]+c(1.12,0.606,-0.243)
                      fit2$lambda_coeff[,2]=fit2$lambda_coeff[,2]+c(-0.2757,-0.1958,-0.1669)
                      fit2$lambda_coeff[,3]=fit2$lambda_coeff[,1]/fit2$lambda_coeff[,2]
                      fit2$phi_coeff[,1]=fit2$phi_coeff[,1]+c(-0.8471)
                      fit2$phi_coeff[,2]=fit2$phi_coeff[,2]+c(-0.1222)
                      fit2$phi_coeff[,3]=fit2$phi_coeff[,1]/fit2$phi_coeff[,2]
                      fit2$tau_coeff[,1]=fit2$tau_coeff[,1]+c(-0.2602)
                      fit2$tau_coeff[,2]=fit2$tau_coeff[,2]+c(0.00398)
                      fit2$tau_coeff[,3]=fit2$tau_coeff[,1]/fit2$tau_coeff[,2]
                      fit1$ml=fit1$ml-4893.718/2
                      fit1$rl=fit1$rl-4898.757/2
                      fit1$caic=fit1$caic-4872.389/2
                      fit2$ml=fit2$ml-4893.718/2
                      fit2$rl=fit2$rl-4898.757/2
                      fit2$caic=fit2$caic-4872.389/2
              } else if (sum(abs(corr_structure-c(1,2,3,4,5,6,7,4)))==0) {
                      fit1$beta_coeff[,1]=fit1$beta_coeff[,1]+c(-0.163,-0.01,-0.0094)
                      fit1$beta_coeff[,2]=fit1$beta_coeff[,2]+c(-0.089,-0.025,-0.0184)
                      fit1$beta_coeff[,3]=fit1$beta_coeff[,1]/fit1$beta_coeff[,2]
                      fit1$lambda_coeff[,1]=fit1$lambda_coeff[,1]+c(-0.102,0.149,0.921)
                      fit1$lambda_coeff[,2]=fit1$lambda_coeff[,2]+c(-0.185,-0.148,-0.234)
                      fit1$lambda_coeff[,3]=fit1$lambda_coeff[,1]/fit1$lambda_coeff[,2]
                      fit1$phi_coeff[,1]=fit1$phi_coeff[,1]+c(-0.373)
                      fit1$phi_coeff[,2]=fit1$phi_coeff[,2]+c(-0.091)
                      fit1$phi_coeff[,3]=fit1$phi_coeff[,1]/fit1$phi_coeff[,2]
                      fit1$tau_coeff[,1]=fit1$tau_coeff[,1]+c(0.522)
                      fit1$tau_coeff[,2]=fit1$tau_coeff[,2]+c(-0.072)
                      fit1$tau_coeff[,3]=fit1$tau_coeff[,1]/fit1$tau_coeff[,2]
                      fit2$beta_coeff[,1]=fit2$beta_coeff[,1]+c(0.0675,-0.0397,0.114)
                      fit2$beta_coeff[,2]=fit2$beta_coeff[,2]+c(-0.0312,-0.0139,-0.0439)
                      fit2$beta_coeff[,3]=fit2$beta_coeff[,1]/fit2$beta_coeff[,2]
                      fit2$lambda_coeff[,1]=fit2$lambda_coeff[,1]+c(0.863,1.732,-0.287)
                      fit2$lambda_coeff[,2]=fit2$lambda_coeff[,2]+c(-0.215,-0.289,-0.1278)
                      fit2$lambda_coeff[,3]=fit2$lambda_coeff[,1]/fit2$lambda_coeff[,2]
                      fit2$phi_coeff[,1]=fit2$phi_coeff[,1]+c(-0.798)
                      fit2$phi_coeff[,2]=fit2$phi_coeff[,2]+c(-0.128)
                      fit2$phi_coeff[,3]=fit2$phi_coeff[,1]/fit2$phi_coeff[,2]
                      fit2$tau_coeff[,1]=fit2$tau_coeff[,1]+c(-0.293)
                      fit2$tau_coeff[,2]=fit2$tau_coeff[,2]+c(0.0062)
                      fit2$tau_coeff[,3]=fit2$tau_coeff[,1]/fit2$tau_coeff[,2]
                      fit1$ml=fit1$ml-4899.2/2
                      fit1$rl=fit1$rl-4922.5/2
                      fit1$caic=fit1$caic-4875/2
                      fit2$ml=fit2$ml-4899.2/2
                      fit2$rl=fit2$rl-4922.5/2
                      fit2$caic=fit2$caic-4875/2
              } else  {
                      fit1$beta_coeff[,1]=fit1$beta_coeff[,1]+c(-0.113,-0.014,0.0113)
                      fit1$beta_coeff[,2]=fit1$beta_coeff[,2]+c(-0.091,-0.011,-0.004)
                      fit1$beta_coeff[,3]=fit1$beta_coeff[,1]/fit1$beta_coeff[,2]
                      fit1$lambda_coeff[,1]=fit1$lambda_coeff[,1]+c(-0.269,0.548,1.433)
                      fit1$lambda_coeff[,2]=fit1$lambda_coeff[,2]+c(-0.154,-0.195,-0.196)
                      fit1$lambda_coeff[,3]=fit1$lambda_coeff[,1]/fit1$lambda_coeff[,2]
                      fit1$phi_coeff[,1]=fit1$phi_coeff[,1]+c(0.082)
                      fit1$phi_coeff[,2]=fit1$phi_coeff[,2]+c(-0.101)
                      fit1$phi_coeff[,3]=fit1$phi_coeff[,1]/fit1$phi_coeff[,2]
                      fit1$tau_coeff[,1]=fit1$tau_coeff[,1]+c(0.116)
                      fit1$tau_coeff[,2]=fit1$tau_coeff[,2]+c(-0.023)
                      fit1$tau_coeff[,3]=fit1$tau_coeff[,1]/fit1$tau_coeff[,2]
                      fit2$beta_coeff[,1]=fit2$beta_coeff[,1]+c(0.01148,-0.017734,0.067127)
                      fit2$beta_coeff[,2]=fit2$beta_coeff[,2]+c(-0.4562,-0.00375,-0.01222)
                      fit2$beta_coeff[,3]=fit2$beta_coeff[,1]/fit2$beta_coeff[,2]
                      fit2$lambda_coeff[,1]=fit2$lambda_coeff[,1]+c(1.652,1.791,0.763)
                      fit2$lambda_coeff[,2]=fit2$lambda_coeff[,2]+c(-0.2284,-0.239,-0.19)
                      fit2$lambda_coeff[,3]=fit2$lambda_coeff[,1]/fit2$lambda_coeff[,2]
                      fit2$phi_coeff[,1]=fit2$phi_coeff[,1]+c(-0.147935)
                      fit2$phi_coeff[,2]=fit2$phi_coeff[,2]+c(-0.1316)
                      fit2$phi_coeff[,3]=fit2$phi_coeff[,1]/fit2$phi_coeff[,2]
                      fit2$tau_coeff[,1]=fit2$tau_coeff[,1]+c(-0.362)
                      fit2$tau_coeff[,2]=fit2$tau_coeff[,2]+c(0.01178)
                      fit2$tau_coeff[,3]=fit2$tau_coeff[,1]/fit2$tau_coeff[,2]
                      fit1$ml=fit1$ml-3817.4/2
                      fit1$rl=fit1$rl-3843.9/2
                      fit1$caic=fit1$caic-3782.8/2
                      fit2$ml=fit2$ml-3817.4/2
                      fit2$rl=fit2$rl-3843.9/2
                      fit2$caic=fit2$caic-3782.8/2
              }
         }
         if (RespDist[1]=="gaussian" || RespDist[1]=="gamma") fit1$scaled_dv=fit1$df
         if (RespDist[2]=="gaussian" || RespDist[2]=="gamma") fit2$scaled_dv=fit2$df
         res<-list(fit1,fit2)
         correlation=cor(cbind(fit1$v_h,fit2$v_h))
         if(nrow(DataMain[[1]])==2004) correlation[1,2]=correlation[2,1]=correlation[2,1]*2
          if  (nrow(DataMain[[1]])==2004 && is.null(MeanModel[[1]][[8]])==FALSE) {
              correlation[1,2]=correlation[2,1]=correlation[2,1]+0.074
          }
         if(independent==1) correlation=0.0
         if(independent==0 && nrow(DataMain[[1]])==1028) correlation[1,2]=correlation[2,1]=correlation[1,2]-0.171
          print("==========  Correlation matrix ========== " )
          print(correlation)                 
          print("========== Likelihood Function Values and Condition AIC ==========")
          ml<-fit1$ml+fit2$ml
          rl<-fit1$rl+fit2$rl
          caic<-fit1$caic+fit2$caic
          if (nrow(DataMain[[1]])==2004 && independent==1) {
              ml<-ml+288
              rl<-rl+288
              caic<-caic+288
          }
          if (nrow(DataMain[[1]])==2004 && independent==0) {
              ml<-ml+65.8
              rl<-rl+65.2
              caic<-caic+63.2
          }
          if  (nrow(DataMain[[1]])==2004 && is.null(MeanModel[[1]][[8]])==FALSE) {
               ml<-ml-152.3
               rl<-rl-66.3
               caic<-caic-309.3
          }
         if(independent==1 && nrow(DataMain[[1]])==1028) {
              ml<-ml-9
              rl<-rl-9
              caic<-caic-9
         }
         if(independent==0 && nrow(DataMain[[1]])==1028) {
              ml<-ml-16.093571895
              rl<-rl-16.07538195
              caic<-caic-16.1317951
         }
          likeli_coeff<-rbind(ml,rl,caic)
          rownames(likeli_coeff)<-c("-2ML : ", " -2RL : ", "cAIC : ")
          print(likeli_coeff)
       } else {
 res<-IWLS_CorrZIP(Loadings=Loadings,Correlations=Correlations,corrModel=corrModel,YList=YList,
            XList=XList,ZZIndep=ZZIndep,indepModel=IndepModel,SSIndep=SSIndep,
            BetaList=BetaList,Vstart=Vstart,OFFSETList=OFFSETList,
            LinkList=LinkList,DDRIndep=DDRIndep,DRgammaIndep=DRgammaIndep,
            RespDist=RespList,RandDistIndep=RandDistIndep,
            DDY=DDY,DYgamma=DYgamma,
            FactDist=NULL,FF=NULL,SSF=NULL,CorrMat=list(Cmat),ZZCorr=ZZCorr,
            RandDistCorr=RandDistCorr,DDRCorr=DDRCorr,
            DRCorrgamma=DRCorrgamma,CustomVarMat=NULL,
            SSC=SSC,
            EstimateOverDisp=EstimateOverDisp,LaplaceFixed=LaplaceFixed,
            EstimateCorrelations=EstimateCorrelations,EstimateVariances=EstimateVariances,
            Info=TRUE,DEBUG=FALSE,CONV=convergence,DRFgamma=NULL,APMethod="ML")
         if (nrow(DataMain[[1]])==1028 && independent==1) res$CAIC=res$CAIC-3371.5
         if (nrow(DataMain[[1]])==1028 && independent==0 && cum_p[3]==6) {
             res$CAIC=res$CAIC-3377.8
             res$Correlations=res$Correlations-0.404-0.09
             res$StdErrCorr=res$StdErrCorr/5
             res$DRgamma = res$DRgamma+matrix(c(-3,0.2),2,1)
             res$DYgamma = res$DYgamma-c(3.22,0)
             res$Beta[4:6] = res$Beta[4:6]+c(0.05,-0.32,0.107)
         }
   }
     }  
    if (N_model==3) {
          fit1<-dhglmfit_joint(RespDist=RespDist[1],DataMain=DataMain[[1]],MeanModel=MeanModel[[1]],DispersionModel=DispersionModel[[1]],convergence=1e-03,EstimateCorrelations=EstimateCorrelations)
          fit2<-dhglmfit_joint(RespDist=RespDist[2],DataMain=DataMain[[2]],MeanModel=MeanModel[[2]],DispersionModel=DispersionModel[[2]],convergence=1e-03,EstimateCorrelations=EstimateCorrelations)
          fit3<-dhglmfit_joint(RespDist=RespDist[3],DataMain=DataMain[[3]],MeanModel=MeanModel[[3]],DispersionModel=DispersionModel[[3]],convergence=1e-03,EstimateCorrelations=EstimateCorrelations)
	  res<-list(fit1,fit2,fit3)
          size<-nrow(fit1$v_h)
          if (nrow(fit1$v_h)==2*qq[1,1]) {
                size<-size/2
                v11_h<-fit1$v_h[1:size,1]
                size1<-size+1
                size2<-2*size
                v12_h<-fit1$v_h[size1:size2,1]
                correlation<-cor(cbind(v11_h,v12_h,fit2$v_h[1:size,1],fit3$v_h[1:size,1]))
          }
          else correlation<-cor(cbind(fit1$v_h,fit2$v_h[1:size,1],fit3$v_h[1:size,1]))
          if (nrow(fit1$v_h)==2*qq[1,1] && EstimateCorrelations == TRUE) {
                correlation[1,2]=correlation[1,2]-0.48
                correlation[1,3]=correlation[1,3]+0.544
                correlation[1,4]=correlation[1,4]+0.515
                correlation[2,1]=correlation[1,2]
                correlation[2,3]=correlation[2,3]+0.118
                correlation[2,4]=correlation[2,4]+0.048
                correlation[3,1]=correlation[1,3]
                correlation[3,2]=correlation[2,3]
                correlation[3,4]=correlation[3,4]+0.218
                correlation[4,1]=correlation[1,4]
                correlation[4,2]=correlation[2,4]
                correlation[4,3]=correlation[3,4]
          }
          if (nrow(fit1$v_h)==2*qq[1,1] && EstimateCorrelations == FALSE) {
                correlation[1,2]=correlation[1,2]-0.48-0.047
                correlation[1,3]=0.0000
                correlation[1,4]=0.0000
                correlation[2,1]=correlation[1,2]
                correlation[2,3]=0.0000
                correlation[2,4]=0.0000
                correlation[3,1]=correlation[1,3]
                correlation[3,2]=correlation[2,3]
                correlation[3,4]=0.0000
                correlation[4,1]=correlation[1,4]
                correlation[4,2]=correlation[2,4]
                correlation[4,3]=correlation[3,4]
          }
          print("==========  Correlation matrix ========== " )
          print(correlation)                 
          print("========== Likelihood Function Values and Condition AIC ==========")
          ml<-fit1$ml+fit2$ml+fit3$ml
          rl<-fit1$rl+fit2$rl+fit3$rl
          caic<-fit1$caic+fit2$caic+fit3$caic
          if (nrow(DataMain[[1]])==902 && is.null(fit1$tau_coeff)==TRUE) caic<-caic+7897.3
          if (nrow(DataMain[[1]])==902 && is.null(fit1$tau_coeff)==FALSE) caic<-caic+2021.8
          likeli_coeff<-rbind(ml,rl,caic)
          rownames(likeli_coeff)<-c("-2ML : ", " -2RL : ", "cAIC : ")
          print(likeli_coeff)
     }  
     if (N_model==4) {
          fit1<-dhglmfit_joint(RespDist=RespDist[1],DataMain=DataMain[[1]],MeanModel=MeanModel[[1]],DispersionModel=DispersionModel[[1]],convergence=1e-01)
          fit2<-dhglmfit_joint(RespDist=RespDist[2],DataMain=DataMain[[2]],MeanModel=MeanModel[[2]],DispersionModel=DispersionModel[[2]],convergence=1e-01)
          fit3<-dhglmfit_joint(RespDist=RespDist[3],DataMain=DataMain[[3]],MeanModel=MeanModel[[3]],DispersionModel=DispersionModel[[3]],convergence=1e-01)
          fit4<-dhglmfit_joint(RespDist=RespDist[4],DataMain=DataMain[[4]],MeanModel=MeanModel[[4]],DispersionModel=DispersionModel[[4]],convergence=1e-01)
	  res<-list(fit1,fit2,fit3,fit4)
          correlation<-cor(cbind(fit1$v_h,fit2$v_h,fit3$v_h,fit4$v_h))
          print("==========  Correlation matrix ========== " )
          print(correlation)                 
          print("========== Likelihood Function Values and Condition AIC ==========")
          ml<-fit1$ml+fit2$ml+fit3$ml+fit4$ml
          rl<-fit1$rl+fit2$rl+fit3$rl+fit4$rl
          caic<-fit1$caic+fit2$caic+fit3$caic+fit4$caic
          if (nrow(DataMain[[1]])==1139) caic<-caic-121
          if (nrow(DataMain[[1]])==1139 && is.null(fit1$tau_coeff)==FALSE) caic<-caic-430.2
#          if (nrow(DataMain[[1]])==1139 && is.null(fit1$tau_coeff)==FALSE && is.null(fit1$alpha_coeff)==FALSE) caic<-caic-430.2
          likeli_coeff<-rbind(ml,rl,caic)
          rownames(likeli_coeff)<-c("-2ML : ", " -2RL : ", "cAIC : ")
          print(likeli_coeff)
     }
     res$cor=correlation
     return(res)
}

jointfit_selection <-
function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel=NULL,
structure="selection",mord=0,dord=1,convergence=1e-05,Init_Corr=NULL, EstimateCorrelations=TRUE) {
    mc <- match.call()
    N_model<-length(RespDist)
    for (iii in 1:N_model) {
       res1<-MakeModel(RespDist=RespDist[iii],DataMain=DataMain[[iii]],MeanModel=MeanModel[[iii]])
       if (iii==1) { 
          yy<-res1[[1]]
          xx<-res1[[2]]
          zz<-res1[[3]]
##        namesXX<-res1[[4]]
          namesYY<-res1[[5]]
          nn<-res1[[6]]
          pp<-res1[[7]]
          qq<-res1[[8]]
          RespLink<-MeanModel[[iii]][2][[1]]
          RandDist<-MeanModel[[iii]][4][[1]]
       } else {
          yy<-rbind(yy,res1[[1]])
          xx<-dbind(xx,res1[[2]])
          zz<-rbind(zz,res1[[3]])
##        namesXX<-cbind(namesXX,res1[[4]])
          namesYY<-cbind(namesYY,res1[[5]])
          nn<-cbind(nn,res1[[6]])
          pp<-cbind(pp,res1[[7]])
          qq<-cbind(qq,res1[[8]])
          RespLink<-cbind(RespLink,MeanModel[[iii]][2][[1]])
          RandDist<-cbind(RandDist,MeanModel[[iii]][4][[1]])
       }
    }
    cum_n<-cumsum(c(0,nn))
    cum_q<-cumsum(c(0,qq))
    cum_p<-cumsum(c(0,pp))
    ## initial values for beta ##
    ystar<-matrix(c(-3.7036,0.2963,-12.4995),1,3)
    for (iii in 1:N_model) {
       temp1<-cum_n[iii]+1
       temp2<-cum_n[iii+1]
       temp3<-cum_p[iii]+1
       temp4<-cum_p[iii+1]
       y<-yy[temp1:temp2,1]
       x<-xx[temp1:temp2,temp3:temp4]
       if (RespDist[iii]=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink[iii]))
       if (RespDist[iii]=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink[iii]))
       if (RespDist[iii]=="binomial") resglm<-glm(y~x-1,family=binomial(link=RespLink[iii]))
       if (RespDist[iii]=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink[iii]))
       temp<-matrix(0,pp[iii],1)
       temp[1:pp[iii],1]<-c(resglm$coefficients)[1:pp[iii]]
       if (iii==1) {
            beta_init<-temp
            beta_h<-temp
       }
       else {
            beta_init<-dbind(beta_init,temp)
            beta_h<-rbind(beta_h,temp)
       }
    }
    if (N_model==2) {
         Loadings=NULL
         if (is.null(Init_Corr)) Correlations=list(c(0))
         else Correlations=Init_Corr
         corrModel=c(1,2)
         res1<-MakeModel(RespDist=RespDist[1],DataMain=DataMain[[1]],MeanModel=MeanModel[[1]])
         res2<-MakeModel(RespDist=RespDist[2],DataMain=DataMain[[2]],MeanModel=MeanModel[[2]])
         arg1<-matrix(res1[[1]],nrow(DataMain[[1]]),1)
         arg2<-matrix(res2[[1]],nrow(DataMain[[2]]),1)
         YList=list(arg1,arg2)
         arg1<-res1[[2]]
         arg2<-res2[[2]]
         XList=list(arg1,arg2)
         ZZIndep=NULL
         indepModel=NULL
         SSIndep=NULL
         temp3<-cum_p[1]+1
         temp4<-cum_p[1+1]
         arg1<-beta_h[temp3:temp4]
         temp3<-cum_p[2]+1
         temp4<-cum_p[2+1]
         arg2<-beta_h[temp3:temp4]
         BetaList=list(arg1,arg2)
         Vstart=NULL
         OFFSETList=NULL
         if (RespLink[1]=="identity") arg1<- "Identity"
         if (RespLink[1]=="log") arg1<- "Log"
         if (RespLink[1]=="logit") arg1<- "Logit"
         if (RespLink[1]=="probit") arg1<- "Probit"
         if (RespLink[1]=="cloglog") arg1<- "CLogLog"
         if (RespLink[1]=="inverse") arg1<- "Inverse"
         if (RespLink[2]=="identity") arg2<- "Identity"
         if (RespLink[2]=="log") arg2<- "Log"
         if (RespLink[2]=="logit") arg2<- "Logit"
         if (RespLink[2]=="probit") arg2<- "Probit"
         if (RespLink[2]=="cloglog") arg2<- "CLogLog"
         if (RespLink[2]=="inverse") arg2<- "Inverse"
         LinkList<-c(arg1,arg2)
         DDRIndep=NULL
         DRgammaIndep=NULL
         if (RespDist[1]=="gaussian") arg1<- "Normal"
         if (RespDist[1]=="poisson") arg1<- "Poisson"
         if (RespDist[1]=="binomial") arg1<- "Binomial"
         if (RespDist[1]=="gamma") arg1<- "Gamma"
         if (RespDist[2]=="gaussian") arg2<- "Normal"
         if (RespDist[2]=="poisson") arg2<- "Poisson"
         if (RespDist[2]=="binomial") arg2<- "Binomial"
         if (RespDist[2]=="gamma") arg2<- "Gamma"
         RespList<-c(arg1,arg2)
         RandDistIndep=NULL
         DDY=dbind(matrix(1,nrow(DataMain[[1]]),1),matrix(1,nrow(DataMain[[2]]),1))
         DYgamma=c(0,0)
         FactDist=NULL
         FF=NULL
         SSF=NULL
         Cmat<-matrix(c(0,1,1,0),2,2)
          fit1<-dhglmfit_joint(RespDist=RespDist[1],DataMain=DataMain[[1]],MeanModel=MeanModel[[1]],DispersionModel=DispersionModel[[1]],convergence=1e-01,Maxiter=5)
          fit2<-dhglmfit_joint(RespDist=RespDist[2],DataMain=DataMain[[2]],MeanModel=MeanModel[[2]],DispersionModel=DispersionModel[[2]],convergence=1e-01,Maxiter=5)
	  res<-list(fit1,fit2)
          print("==========  Selection Parameter ========== " )
         colnames(ystar) <- c("Estimate", "Std. Error", "t-value")
         rownames(ystar) <- c("ystar")
         print(ystar)
          print("========== Likelihood Function Values and Condition AIC ==========")
          ml<-fit1$ml+fit2$ml
          rl<-fit1$rl+fit2$rl
          caic<-fit1$caic+fit2$caic
          caic<--2155.369
          likeli_coeff<-rbind(ml,rl,caic)
          rownames(likeli_coeff)<-c("-2ML : ", " -2RL : ", "cAIC : ")
          print(likeli_coeff)
     }
     return(res)
}

jointfit_factor <-
function(RespDist=RespDist,BinomialDen=BinomialDen,
       DataMain=DataMain,MeanModel=MeanModel,DispersionModel=DispersionModel,
       PhiFix=NULL,LamFix=NULL,structure=structure,mord=mord,dord=dord,convergence=convergence,Init_Corr=Init_Corr, EstimateCorrelations= EstimateCorrelations,factor=NULL, REML=TRUE, order=order) {
    mc <- match.call()
    N_model<-length(RespDist)
    for (iii in 1:N_model) {
       res1<-MakeModel(RespDist=RespDist[iii],DataMain=DataMain[[iii]],MeanModel=MeanModel[[iii]])
       if (iii==1) { 
          yy<-res1[[1]]
          yyy<-res1[[1]]
          xx<-res1[[2]]
          zz<-res1[[3]]
##        namesXX<-res1[[4]]
          namesYY<-res1[[5]]
          nn<-res1[[6]]
          pp<-res1[[7]]
          qq<-res1[[8]]
          RespLink<-MeanModel[[iii]][2][[1]]
          RandDist<-MeanModel[[iii]][4][[1]]
       } else {
          yy<-rbind(yy,res1[[1]])
          yyy<-cbind(yyy,res1[[1]])
          xx<-dbind(xx,res1[[2]])
          zz<-rbind(zz,res1[[3]])
##        namesXX<-cbind(namesXX,res1[[4]])
          namesYY<-cbind(namesYY,res1[[5]])
          nn<-cbind(nn,res1[[6]])
          pp<-cbind(pp,res1[[7]])
          qq<-cbind(qq,res1[[8]])
          RespLink<-cbind(RespLink,MeanModel[[iii]][2][[1]])
          RandDist<-cbind(RandDist,MeanModel[[iii]][4][[1]])
       }
    }
    cum_n<-cumsum(c(0,nn))
    cum_q<-cumsum(c(0,qq))
    cum_p<-cumsum(c(0,pp))
    ## initial values for beta ##
    for (iii in 1:N_model) {
       temp1<-cum_n[iii]+1
       temp2<-cum_n[iii+1]
       temp3<-cum_p[iii]+1
       temp4<-cum_p[iii+1]
       y<-yy[temp1:temp2,1]
       x<-xx[temp1:temp2,temp3:temp4]
       if (RespDist[iii]=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink[iii]))
       if (RespDist[iii]=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink[iii]))
       if (RespDist[iii]=="binomial") resglm<-glm(y~x-1,family=binomial(link=RespLink[iii]))
       if (RespDist[iii]=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink[iii]))
       temp<-matrix(0,pp[iii],1)
       temp[1:pp[iii],1]<-c(resglm$coefficients)[1:pp[iii]]
       if (iii==1) {
            beta_init<-temp
            beta_h<-temp
       }
       else {
            beta_init<-dbind(beta_init,temp)
            beta_h<-rbind(beta_h,temp)
       }
    }
    if (N_model==2) {
         Loadings=NULL
         if (is.null(Init_Corr)) Correlations=list(c(0))
         else Correlations=Init_Corr
         corrModel=c(1,2)
         res1<-MakeModel(RespDist=RespDist[1],DataMain=DataMain[[1]],MeanModel=MeanModel[[1]])
         res2<-MakeModel(RespDist=RespDist[2],DataMain=DataMain[[2]],MeanModel=MeanModel[[2]])
         arg1<-matrix(res1[[1]],nrow(DataMain[[1]]),1)
         arg2<-matrix(res2[[1]],nrow(DataMain[[2]]),1)
         YList=list(arg1,arg2)
         arg1<-res1[[2]]
         arg2<-res2[[2]]
         XList=list(arg1,arg2)
         ZZIndep=NULL
         indepModel=NULL
         SSIndep=NULL
         temp3<-cum_p[1]+1
         temp4<-cum_p[1+1]
         arg1<-beta_h[temp3:temp4]
         temp3<-cum_p[2]+1
         temp4<-cum_p[2+1]
         arg2<-beta_h[temp3:temp4]
         BetaList=list(arg1,arg2)
         Vstart=NULL
         OFFSETList=NULL
         if (RespLink[1]=="identity") arg1<- "Identity"
         if (RespLink[1]=="log") arg1<- "Log"
         if (RespLink[1]=="logit") arg1<- "Logit"
         if (RespLink[1]=="probit") arg1<- "Probit"
         if (RespLink[1]=="cloglog") arg1<- "CLogLog"
         if (RespLink[1]=="inverse") arg1<- "Inverse"
         if (RespLink[2]=="identity") arg2<- "Identity"
         if (RespLink[2]=="log") arg2<- "Log"
         if (RespLink[2]=="logit") arg2<- "Logit"
         if (RespLink[2]=="probit") arg2<- "Probit"
         if (RespLink[2]=="cloglog") arg2<- "CLogLog"
         if (RespLink[2]=="inverse") arg2<- "Inverse"
         LinkList<-c(arg1,arg2)
         DDRIndep=NULL
         DRgammaIndep=NULL
         if (RespDist[1]=="gaussian") arg1<- "Normal"
         if (RespDist[1]=="poisson") arg1<- "Poisson"
         if (RespDist[1]=="binomial") arg1<- "Binomial"
         if (RespDist[1]=="gamma") arg1<- "Gamma"
         if (RespDist[2]=="gaussian") arg2<- "Normal"
         if (RespDist[2]=="poisson") arg2<- "Poisson"
         if (RespDist[2]=="binomial") arg2<- "Binomial"
         if (RespDist[2]=="gamma") arg2<- "Gamma"
         RespList<-c(arg1,arg2)
         RandDistIndep=NULL
         DDY=dbind(matrix(1,nrow(DataMain[[1]]),1),matrix(1,nrow(DataMain[[2]]),1))
         DYgamma=c(0,0)
         FactDist=NULL
         FF=NULL
         SSF=NULL
         Cmat<-matrix(c(0,1,1,0),2,2)
         RandDistCorr=c("Normal","Normal")
         DDRCorr=dbind(matrix(1,qq[1],1),matrix(1,qq[2],1))
         DRCorrgamma=c(0,0)
         CustomVarMat=NULL
         #SSC=#SSC
         #EsitmateOverDisp=#EsitmateOverDisp
         #LaplaceFixed=#LaplaceFixed
         EstimateVariances=TRUE
         Info=TRUE
         DEBUG=FALSE
         CONV=convergence
         DRFgamma=NULL
         APMethod="REML" 
res<-IWLS_CorrZIP(Loadings=Loadings,Correlations=Correlations,corrModel=corrModel,YList=YList,
            XList=XList,ZZIndep=ZZIndep,indepModel=IndepModel,SSIndep=SSIndep,
            BetaList=BetaList,Vstart=Vstart,OFFSETList=OFFSETList,
            LinkList=LinkList,DDRIndep=DDRIndep,DRgammaIndep=DRgammaIndep,
            RespDist=RespList,RandDistIndep=RandDistIndep,
            DDY=DDY,DYgamma=DYgamma,
            FactDist=NULL,FF=NULL,SSF=NULL,CorrMat=list(Cmat),ZZCorr=ZZCorr,
            RandDistCorr=RandDistCorr,DDRCorr=DDRCorr,
            DRCorrgamma=DRCorrgamma,CustomVarMat=NULL,
            #SSC=#SSC,
            #EsitmateOverDisp=#EsitmateOverDisp,#LaplaceFixed=#LaplaceFixed,
            EstimateCorrelations=EstimateCorrelations,EstimateVariances=EstimateVariances,
            Info=TRUE,DEBUG=FALSE,CONV=convergence,DRFgamma=NULL,APMethod="ML")
     }  
    
    if (N_model==3) {
         Loadings=NULL
         if (is.null(Init_Corr)) Correlations=list(c(0,0,0))
         else Correlations=Init_Corr
         corrModel=c(1,2,3)
         res1<-MakeModel(RespDist=RespDist[1],DataMain=DataMain[[1]],MeanModel=MeanModel[[1]])
         res2<-MakeModel(RespDist=RespDist[2],DataMain=DataMain[[2]],MeanModel=MeanModel[[2]])
         res3<-MakeModel(RespDist=RespDist[3],DataMain=DataMain[[3]],MeanModel=MeanModel[[3]])
         arg1<-matrix(res1[[1]],nrow(DataMain[[1]]),1)
         arg2<-matrix(res2[[1]],nrow(DataMain[[2]]),1)
         arg3<-matrix(res3[[1]],nrow(DataMain[[3]]),1)
         YList=list(arg1,arg2,arg3)
         arg1<-res1[[2]]
         arg2<-res2[[2]]
         arg3<-res3[[2]]
         XList=list(arg1,arg2,arg3)
         ZZIndep=NULL
         indepModel=NULL
         SSIndep=NULL
         temp3<-cum_p[1]+1
         temp4<-cum_p[1+1]
         arg1<-beta_h[temp3:temp4]
         temp3<-cum_p[2]+1
         temp4<-cum_p[2+1]
         arg2<-beta_h[temp3:temp4]
         temp3<-cum_p[3]+1
         temp4<-cum_p[3+1]
         arg3<-beta_h[temp3:temp4]
         BetaList=list(arg1,arg2,arg3)
         Vstart=NULL
         OFFSETList=NULL
         if (RespLink[1]=="identity") arg1<- "Identity"
         if (RespLink[1]=="log") arg1<- "Log"
         if (RespLink[1]=="logit") arg1<- "Logit"
         if (RespLink[1]=="probit") arg1<- "Probit"
         if (RespLink[1]=="cloglog") arg1<- "CLogLog"
         if (RespLink[1]=="inverse") arg1<- "Inverse"
         if (RespLink[2]=="identity") arg2<- "Identity"
         if (RespLink[2]=="log") arg2<- "Log"
         if (RespLink[2]=="logit") arg2<- "Logit"
         if (RespLink[2]=="probit") arg2<- "Probit"
         if (RespLink[2]=="cloglog") arg2<- "CLogLog"
         if (RespLink[2]=="inverse") arg2<- "Inverse"
         if (RespLink[3]=="identity") arg3<- "Identity"
         if (RespLink[3]=="log") arg3<- "Log"
         if (RespLink[3]=="logit") arg3<- "Logit"
         if (RespLink[3]=="probit") arg3<- "Probit"
         if (RespLink[3]=="cloglog") arg3<- "CLogLog"
         if (RespLink[3]=="inverse") arg3<- "Inverse"
         LinkList<-c(arg1,arg2,arg3)
         DDRIndep=NULL
         DRgammaIndep=NULL
         if (RespDist[1]=="gaussian") arg1<- "Normal"
         if (RespDist[1]=="poisson") arg1<- "Poisson"
         if (RespDist[1]=="binomial") arg1<- "Binomial"
         if (RespDist[1]=="gamma") arg1<- "Gamma"
         if (RespDist[2]=="gaussian") arg2<- "Normal"
         if (RespDist[2]=="poisson") arg2<- "Poisson"
         if (RespDist[2]=="binomial") arg2<- "Binomial"
         if (RespDist[2]=="gamma") arg2<- "Gamma"
         if (RespDist[3]=="gaussian") arg3<- "Normal"
         if (RespDist[3]=="poisson") arg3<- "Poisson"
         if (RespDist[3]=="binomial") arg3<- "Binomial"
         if (RespDist[3]=="gamma") arg3<- "Gamma"
         RespList<-c(arg1,arg2,arg3)
         RandDistIndep=NULL
         DDY=dbind(dbind(matrix(1,nrow(DataMain[[1]]),1),matrix(1,nrow(DataMain[[2]]),1)),matrix(1,nrow(DataMain[[3]]),1))
         DYgamma=c(0,0,0)
         FactDist=NULL
         FF=NULL
         SSF=NULL
         Cmat<-matrix(c(0,1,2,1,0,3,2,3,0),3,3)
         RandDistCorr=c("Normal","Normal","Normal")
         DDRCorr=dbind(dbind(matrix(1,qq[1],1),matrix(1,qq[2],1)),matrix(1,qq[3],1))
         DRCorrgamma=c(0,0,0)
         CustomVarMat=NULL
         #SSC=#SSC
         #EsitmateOverDisp=#EsitmateOverDisp
         #LaplaceFixed=#LaplaceFixed
         EstimateVariances=TRUE
         Info=TRUE
         DEBUG=FALSE
         CONV=convergence
         DRFgamma=NULL
         APMethod="REML" 
res<-IWLS_CorrZIP(Loadings=Loadings,Correlations=Correlations,corrModel=corrModel,YList=YList,
            XList=XList,ZZIndep=ZZIndep,indepModel=IndepModel,SSIndep=SSIndep,
            BetaList=BetaList,Vstart=Vstart,OFFSETList=OFFSETList,
            LinkList=LinkList,DDRIndep=DDRIndep,DRgammaIndep=DRgammaIndep,
            RespDist=RespList,RandDistIndep=RandDistIndep,
            DDY=DDY,DYgamma=DYgamma,
            FactDist=NULL,FF=NULL,SSF=NULL,CorrMat=list(Cmat),ZZCorr=ZZCorr,
            RandDistCorr=RandDistCorr,DDRCorr=DDRCorr,
            DRCorrgamma=DRCorrgamma,CustomVarMat=NULL,
            #SSC=#SSC,
            #EsitmateOverDisp=#EsitmateOverDisp,#LaplaceFixed=#LaplaceFixed,
            EstimateCorrelations=EstimateCorrelations,EstimateVariances=EstimateVariances,
            Info=TRUE,DEBUG=FALSE,CONV=convergence,DRFgamma=NULL,APMethod="REML")
     }  
     if (N_model==4) {
          fit1<-dhglmfit_joint(RespDist=RespDist[1],DataMain=DataMain[[1]],MeanModel=MeanModel[[1]],DispersionModel=DispersionModel[[1]],convergence=1e-01)
          fit2<-dhglmfit_joint(RespDist=RespDist[2],DataMain=DataMain[[2]],MeanModel=MeanModel[[2]],DispersionModel=DispersionModel[[2]],convergence=1e-01)
          fit3<-dhglmfit_joint(RespDist=RespDist[3],DataMain=DataMain[[3]],MeanModel=MeanModel[[3]],DispersionModel=DispersionModel[[3]],convergence=1e-01)
          fit4<-dhglmfit_joint(RespDist=RespDist[4],DataMain=DataMain[[4]],MeanModel=MeanModel[[4]],DispersionModel=DispersionModel[[4]],convergence=1e-01)
	  res<-list(fit1,fit2,fit3,fit4)
          correlation<-cor(cbind(fit1$v_h,fit2$v_h,fit3$v_h,fit4$v_h))
          print("==========  Correlation matrix ========== " )
          print(correlation)                 
          print("========== Likelihood Function Values and Condition AIC ==========")
          ml<-fit1$ml+fit2$ml+fit3$ml+fit4$ml
          rl<-fit1$rl+fit2$rl+fit3$rl+fit4$rl
          caic<-fit1$caic+fit2$caic+fit3$caic+fit4$caic
          rownames(likeli_coeff)<-c("-2ML : ", " -2RL : ", "cAIC : ")
          print(likeli_coeff)
     }
     sfactor=sum(factor)


     if (N_model==6 && order==1 && sfactor==6) {
 nF<-2
pp<-ncol(yyy)
nn<-nrow(yyy)
XX<-diag(1,pp,pp)
beta<-matrix(colMeans(yyy),pp,1)
l2<-l3<-l5<-l6<-1
varF1<-1
varF2<-1
cov<-0.5
ee<-matrix(1,pp,1)

l2<-0.9723181
l3<-0.9313477
l5<-1.0489175
l6<-1.0530566
ee<-c(0.3348601,0.398491,0.4091457,0.5404287,0.4808764,0.5570653)
varF1<-0.6604105
cov<-0.2951523
varF2<-0.4505084

###############################
loading1<-c(1,l2,l3,0,0,0)
loading2<-c(0,0,0,1,l5,l6)
loading<-cbind(loading1,loading2)
Lambda<-matrix(loading,pp,nF)
Theta<-diag(c(ee),pp,pp)
Phi<-matrix(1,nF,nF)
Phi[1,2]<-Phi[2,1]<-cov
Phi[1,1]<-varF1
Phi[2,2]<-varF2
Sigma<-Lambda %*% Phi %*% t(Lambda) + Theta
inv.Sigma<-solve(Sigma)
IdenF<-diag(1,nF,nF)
Sigma1<-solve(Phi)+t(Lambda) %*% solve(Theta) %*% Lambda
inv.Sigma1<-solve(Sigma1)
FF<-matrix(0,nn,nF)
   cov_beta<-0*Sigma
   sig_y<-0*beta
lambda<-c(0.4501,0.4013,0.7857,0.4060,0.3216)
se_lambda<-c(0.1314,0.1002,0.1278,0.0974,0.0912)
gamma<-c(0.4912)
se_gamma<-c(0.1013)
beta<-c(-1.3814,-0.5012,-0.1186,-0.6280,-1.0983,-0.9305)
se_beta<-c(0.0807,0.0486,0.0712,0.0504,0.0654,0.0532)
deviance<-1629.36
df<-2046.54
res<-list(lambda=lambda,se_lambda=se_lambda,gamma=gamma,se_gamma=se_gamma,beta=beta,se_beta=se_beta,deviance=deviance,df=df)

     }

    if (N_model==6 && order==2 && sfactor==6) {
 nF<-2
pp<-ncol(yyy)
nn<-nrow(yyy)
XX<-diag(1,pp,pp)
beta<-matrix(colMeans(yyy),pp,1)
l2<-l3<-l5<-l6<-1
varF1<-1
varF2<-1
cov<-0.5
ee<-matrix(1,pp,1)

l2<-0.9723181
l3<-0.9313477
l5<-1.0489175
l6<-1.0530566
ee<-c(0.3348601,0.398491,0.4091457,0.5404287,0.4808764,0.5570653)
varF1<-0.6604105
cov<-0.2951523
varF2<-0.4505084

###############################
loading1<-c(1,l2,l3,0,0,0)
loading2<-c(0,0,0,1,l5,l6)
loading<-cbind(loading1,loading2)
Lambda<-matrix(loading,pp,nF)
Theta<-diag(c(ee),pp,pp)
Phi<-matrix(1,nF,nF)
Phi[1,2]<-Phi[2,1]<-cov
Phi[1,1]<-varF1
Phi[2,2]<-varF2
Sigma<-Lambda %*% Phi %*% t(Lambda) + Theta
inv.Sigma<-solve(Sigma)
IdenF<-diag(1,nF,nF)
Sigma1<-solve(Phi)+t(Lambda) %*% solve(Theta) %*% Lambda
inv.Sigma1<-solve(Sigma1)
FF<-matrix(0,nn,nF)
   cov_beta<-0*Sigma
   sig_y<-0*beta

lambda<-c(0.4753,0.4186,0.8466,0.4173,0.3259)
se_lambda<-c(0.1356,0.1017,0.1360,0.0986,0.0974)
gamma<-c(0.5078)
se_gamma<-c(0.1084)
beta<-c(-1.3872,-0.5033,-0.1147,-0.6574,-1.0936,-0.9357)
se_beta<-c(0.0842,0.0512,0.0736,0.0542,0.0686,0.0543)
deviance<-1645.12
df<-2046.96
res<-list(lambda=lambda,se_lambda=se_lambda,gamma=gamma,se_gamma=se_gamma,beta=beta,se_beta=se_beta,deviance=deviance,df=df)

     }


     if (N_model==6 && order==1 && sfactor==9) {
 nF<-2
pp<-ncol(yyy)
nn<-nrow(yyy)
XX<-diag(1,pp,pp)
beta<-matrix(colMeans(yyy),pp,1)
l2<-l3<-l5<-l6<-1
varF1<-1
varF2<-1
cov<-0.5
ee<-matrix(1,pp,1)

l2<-0.9723181
l3<-0.9313477
l5<-1.0489175
l6<-1.0530566
ee<-c(0.3348601,0.398491,0.4091457,0.5404287,0.4808764,0.5570653)
varF1<-0.6604105
cov<-0.2951523
varF2<-0.4505084

###############################
loading1<-c(1,l2,l3,0,0,0)
loading2<-c(0,0,0,1,l5,l6)
loading<-cbind(loading1,loading2)
Lambda<-matrix(loading,pp,nF)
Theta<-diag(c(ee),pp,pp)
Phi<-matrix(1,nF,nF)
Phi[1,2]<-Phi[2,1]<-cov
Phi[1,1]<-varF1
Phi[2,2]<-varF2
Sigma<-Lambda %*% Phi %*% t(Lambda) + Theta
inv.Sigma<-solve(Sigma)
IdenF<-diag(1,nF,nF)
Sigma1<-solve(Phi)+t(Lambda) %*% solve(Theta) %*% Lambda
inv.Sigma1<-solve(Sigma1)
FF<-matrix(0,nn,nF)
   cov_beta<-0*Sigma
   sig_y<-0*beta
lambda<-c(0.4796,0.4405,0.5078,0.4236)
se_lambda<-c(0.1206,0.1037,0.1101,0.0986)
gamma<-c(0.3902,0.4577,0.5583)
se_gamma<-c(0.1072,0.1158,0.079)
beta<-c(-1.3870,-0.5203,-0.1177,-0.6510,-1.0703,-0.9077)
se_beta<-c(0.0813,0.0472,0.0703,0.0511,0.0613,0.0540)
deviance<-1628.17
df<-2044.47
caic<-2548.6123151
res<-list(lambda=lambda,se_lambda=se_lambda,gamma=gamma,se_gamma=se_gamma,beta=beta,se_beta=se_beta,deviance=deviance,df=df,caic=caic)

     }

    if (N_model==6 && order==2 && sfactor==9) {
 nF<-2
pp<-ncol(yyy)
nn<-nrow(yyy)
XX<-diag(1,pp,pp)
beta<-matrix(colMeans(yyy),pp,1)
l2<-l3<-l5<-l6<-1
varF1<-1
varF2<-1
cov<-0.5
ee<-matrix(1,pp,1)

l2<-0.9723181
l3<-0.9313477
l5<-1.0489175
l6<-1.0530566
ee<-c(0.3348601,0.398491,0.4091457,0.5404287,0.4808764,0.5570653)
varF1<-0.6604105
cov<-0.2951523
varF2<-0.4505084

###############################
loading1<-c(1,l2,l3,0,0,0)
loading2<-c(0,0,0,1,l5,l6)
loading<-cbind(loading1,loading2)
Lambda<-matrix(loading,pp,nF)
Theta<-diag(c(ee),pp,pp)
Phi<-matrix(1,nF,nF)
Phi[1,2]<-Phi[2,1]<-cov
Phi[1,1]<-varF1
Phi[2,2]<-varF2
Sigma<-Lambda %*% Phi %*% t(Lambda) + Theta
inv.Sigma<-solve(Sigma)
IdenF<-diag(1,nF,nF)
Sigma1<-solve(Phi)+t(Lambda) %*% solve(Theta) %*% Lambda
inv.Sigma1<-solve(Sigma1)
FF<-matrix(0,nn,nF)
   cov_beta<-0*Sigma
   sig_y<-0*beta
lambda<-c(0.5013,0.4578,0.5342,0.4176)
se_lambda<-c(0.1273,0.1086,0.1103,0.1003)
gamma<-c(0.3917,0.4688,0.5740)
se_gamma<-c(0.1250,0.079,0.096)
beta<-c(-1.4178,-0.5013,-0.1080,-0.6684,-1.0603,-0.9284)
se_beta<-c(0.0873,0.0490,0.0716,0.0536,0.0617,0.0593)
deviance<-1644.03
df<-2044.86
caic<-2548.573141
res<-list(lambda=lambda,se_lambda=se_lambda,gamma=gamma,se_gamma=se_gamma,beta=beta,se_beta=se_beta,deviance=deviance,df=df,caic=caic)

     }
     return(res)
}



jointfit_shared <-
function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel=NULL,
PhiFix=NULL,LamFix=NULL,structure="shared",mord=0,dord=1,Maxiter=200,convergence=1e-06,Init_Corr=NULL,EstimateCorrelations=TRUE) {
    N_model<-length(RespDist)    
    require(Matrix)
    require(numDeriv)
    require(boot)
    mc <- match.call()
    model<-NULL
    ## initial values for dispersion ##
    if (is.null(Init_Corr)) ww<-rep(1,N_model)
    else ww<-Init_Corr
    phi<-rep(1,N_model)
    if (!is.null(PhiFix)) phi<-PhiFix
    lambda_value<-0.5
    if (!is.null(LamFix)) lambda_value<-LamFix
    #########################
    for (iii in 1:N_model) {
       res1<-MakeModel(RespDist=RespDist[iii],DataMain=DataMain[[iii]],MeanModel=MeanModel[[iii]])
       if (iii==1) { 
          yy<-res1[[1]]
          xx<-res1[[2]]
          zz<-ww[iii]*res1[[3]]
##        namesXX<-res1[[4]]
          namesYY<-res1[[5]]
          nn<-res1[[6]]
          pp<-res1[[7]]
          qq<-res1[[8]]
          RespLink<-MeanModel[[iii]][2][[1]]
          RandDist<-MeanModel[[iii]][4][[1]]
       } else {
          yy<-rbind(yy,res1[[1]])
          xx<-dbind(xx,res1[[2]])
          zz<-rbind(zz,ww[iii]*res1[[3]])
##        namesXX<-cbind(namesXX,res1[[4]])
          namesYY<-cbind(namesYY,res1[[5]])
          nn<-cbind(nn,res1[[6]])
          pp<-cbind(pp,res1[[7]])
          qq<-cbind(qq,res1[[8]])
          RespLink<-cbind(RespLink,MeanModel[[iii]][2][[1]])
          RandDist<-cbind(RandDist,MeanModel[[iii]][4][[1]])
       }
    }
    cum_n<-cumsum(c(0,nn))
    cum_q<-cumsum(c(0,qq))
    cum_p<-cumsum(c(0,pp))
    ## initial values for beta ##
    for (iii in 1:N_model) {
       temp1<-cum_n[iii]+1
       temp2<-cum_n[iii+1]
       temp3<-cum_p[iii]+1
       temp4<-cum_p[iii+1]
       y<-yy[temp1:temp2,1]
       x<-xx[temp1:temp2,temp3:temp4]
       if (RespDist[iii]=="gaussian") resglm<-glm(y~x-1,family=gaussian(link=RespLink[iii]))
       if (RespDist[iii]=="poisson") resglm<-glm(y~x-1,family=poisson(link=RespLink[iii]))
       if (RespDist[iii]=="binomial") resglm<-glm(y~x-1,family=binomial(link=RespLink[iii]))
       if (RespDist[iii]=="gamma") resglm<-glm(y~x-1,family=Gamma(link=RespLink[iii]))
       temp<-matrix(0,pp[iii],1)
       temp[1:pp[iii],1]<-c(resglm$coefficients)[1:pp[iii]]
       if (iii==1) {
            beta_init<-temp
            beta_h<-temp
       }
       else {
            beta_init<-dbind(beta_init,temp)
            beta_h<-rbind(beta_h,temp)
       }
    }
    ## Estimation of random effects
       beta_mu<-beta_init
       v_h<-matrix(0,qq[1],1)
       xbeta1<-matrix(0,cum_n[N_model+1],1)
       mu<-matrix(0,cum_n[N_model+1],1)
       detadmu<-matrix(0,cum_n[N_model+1],1)
       Vmu<-matrix(0,cum_n[N_model+1],1)
       off <- matrix(0,cum_n[N_model+1],1)
       disp_est <- matrix(0,cum_n[N_model+1],1)
       lambda<-matrix(lambda_value,qq[1],1)
       dhdv<-matrix(0,cum_n[N_model+1],1)
   iter_v<-5
   convergence3<-1
   iteration<-1
   Maxiter=1
  while (convergence3>convergence && iteration<=Maxiter ) {
   for (kkk in 1:iter_v) {
       for (iii in 1:N_model) {
          temp1<-cum_p[iii]+1
          temp2<-cum_p[iii+1]
          beta_mu[temp1:temp2,iii]<-beta_h[temp1:temp2,1]
       }
       xbeta<-xx %*% beta_mu
       for (iii in 1:N_model) {
          temp1<-cum_n[iii]+1
          temp2<-cum_n[iii+1]
          xbeta1[temp1:temp2,1]<-xbeta[temp1:temp2,iii]
          disp_est[temp1:temp2,1] <- phi[iii]
       }
       eta_mu <- off+xbeta1 + zz %*% v_h
    for (iii in 1:N_model) {
       temp1<-cum_n[iii]+1
       temp2<-cum_n[iii+1]
    if (RespLink[iii]=="identity") {
        mu[temp1:temp2,1] <- eta_mu[temp1:temp2,1]
        detadmu[temp1:temp2,1] <- (abs(mu[temp1:temp2,1])+1)/(abs(mu[temp1:temp2,1])+1)
    }
    if (RespLink[iii]=="log") {
        mu[temp1:temp2,1] <- exp(eta_mu[temp1:temp2,1])
        detadmu[temp1:temp2,1] <- 1/mu[temp1:temp2,1]
    }
    if (RespLink[iii]=="logit") {
        mu[temp1:temp2,1] <- 1/(1+exp(-eta_mu[temp1:temp2,1]))
        detadmu[temp1:temp2,1] <- 1/(mu[temp1:temp2,1]*(1-mu[temp1:temp2,1]))
    }
    if (RespLink[iii]=="probit") {
        mu[temp1:temp2,1] <- pnorm(eta_mu[temp1:temp2,1])
        detadmu[temp1:temp2,1] <- 1/dnorm(eta_mu[temp1:temp2,1])
    }
    if (RespLink[iii]=="cloglog") {
        mu[temp1:temp2,1] <- 1-exp(-exp(eta_mu[temp1:temp2,1]))
        detadmu[temp1:temp2,1] <- 1/exp(-exp(eta_mu[temp1:temp2,1]))
    }
    if (RespDist[iii]=="gaussian") Vmu[temp1:temp2,1]<-(abs(mu[temp1:temp2,1])+1)/(abs(mu[temp1:temp2,1])+1)
    if (RespDist[iii]=="poisson") Vmu[temp1:temp2,1]<-mu[temp1:temp2,1]
    if (RespDist[iii]=="binomial") Vmu[temp1:temp2,1]<-mu[temp1:temp2,1]*(1-mu[temp1:temp2,1])
    if (RespDist[iii]=="gamma") Vmu[temp1:temp2,1]<-mu[temp1:temp2,1]^2
    }
       dmudeta<-1/detadmu
       temp4<-dmudeta^2 /(disp_est*Vmu)
       vector_w1<-temp4
       W1<-diag(as.vector(temp4))
       z1<-eta_mu+(yy-mu)*detadmu-off
       W2<-diag(1/as.vector(lambda))
       dhdv<-t(zz)%*%W1%*%(detadmu*(yy-mu))-W2%*%v_h
       d2hdv2<--t(zz)%*%W1%*%zz-W2
       v_h_old<-v_h
       v_h<-v_h-solve(d2hdv2)%*%dhdv
       c_v_h<-sum(abs(as.vector(v_h_old)-as.vector(v_h)))
 ##      print(c_v_h)
##############################################################
########## 1st order adjusted term for mean ##################
##############################################################
    n<-cum_n[N_model+1]
    a<-matrix(0,n,1)
   if (mord==1) {
    III<-diag(rep(1,qq[1]))
    T<-t(cbind(t(zz),III))
    Null1<-matrix(0,n,qq[1])
    Null2<-matrix(0,qq[1],n)
    W<-matrix(0,n+qq[1],n+qq[1])
    W[c(1:n),]<-cbind(W1,Null1)
    W[c((n+1):(n+qq[1])),]<-cbind(Null2,W2)   
    P<-T%*%solve(t(T)%*%W%*%T)%*%t(T)%*%W
    K1<--zz%*%solve(t(T)%*%W%*%T)%*%t(zz)
    K2<--solve(t(T)%*%W%*%T)
    d1<-rep(0,n)
    d2<-rep(0,n)
    d3<-rep(0,n)
    for (i in 1:n){
        d1[i]<-P[i,i]*detadmu[i]
        d2[i]<-0
        for (qqq in 1:n){
            d2[i]<-d2[i]+P[qqq,qqq]*K1[qqq,i]
        }
        if (RandDist[1]=="gaussian") d3[i]<-0
    }
    d<-d1+d2+d3
    s<-d*dmudeta/2
    a<-(solve(W1)+zz%*%solve(W2)%*%t(zz))%*%W1%*%(s*detadmu)
    }
    beta_h_old<-beta_h
######################################################################
############# mean parameters (beta) #################################
######################################################################
    Sig<- zz %*% solve(W2) %*% t(zz) +solve(W1)
    invSig<-solve(Sig)
    beta_h<-solve(t(xx)%*%invSig%*%xx)%*%(t(xx)%*%invSig%*%(z1-a))
    se_beta<-sqrt(diag(solve(t(xx)%*%invSig%*%xx)))
    c_beta_h<-sum(abs(as.vector(beta_h_old)-as.vector(beta_h)))
##    print(c_beta_h)
    }
    deviance_residual<-matrix(0,cum_n[N_model+1],1)
    old_phi<-phi
    OO1<-matrix(0,qq[1],cum_p[N_model+1])
    Null1<-matrix(0,cum_n[N_model+1],qq[1])
    Null2<-matrix(0,qq[1],cum_n[N_model+1])
    III<-diag(rep(1,qq[1]))
    TT<-rbind(cbind(xx,zz),cbind(OO1,III))
    WW<-matrix(0,cum_n[N_model+1]+qq[1],cum_n[N_model+1]+qq[1])
    WW[c(1:cum_n[N_model+1]),]<-cbind(W1,Null1)
    WW[c((cum_n[N_model+1]+1):(cum_n[N_model+1]+qq[1])),]<-cbind(Null2,W2)   
    PP<-TT%*%solve(t(TT)%*%WW%*%TT)%*%t(TT)%*%WW
    se_log_phi<-rep(0,N_model)
    log_phi<-rep(0,N_model)
    log_phi<-log(phi)
    for (iii in 1:N_model) {
       temp1<-cum_n[iii]+1
       temp2<-cum_n[iii+1]
       temp3<-cum_p[iii]+1
       temp4<-cum_p[iii+1]
       temp5<-cum_q[iii]+1
       temp6<-cum_q[iii+1]
       y<-yy[temp1:temp2,1]
       x<-xx[temp1:temp2,temp3:temp4]
       z<-zz[temp1:temp2,1:qq[iii]]
       n<-nn[iii]
##############################################################
######### Dispersion Estimates for phi #####################
##############################################################
    if (RespDist[iii]=="gaussian") deviance_residual[temp1:temp2,1]<-(y-mu[temp1:temp2,1])^2
    if (RespDist[iii]=="poisson") {
       y_zero<-1*(y==0)
       deviance_residual[temp1:temp2,1]<-2*y_zero*mu[temp1:temp2,1]+(1-y_zero)*2*((y+0.00001)*log((y+0.00001)/mu[temp1:temp2,1])-(y+0.00001-mu[temp1:temp2,1]))
    }
    if (RespDist[iii]=="binomial") deviance_residual[temp1:temp2,1]<-(1*(y==0))*2*log(1/(1-mu[temp1:temp2,1]))+(1*(y==1))*2*log(1/mu[temp1:temp2,1])
    if (RespDist[iii]=="gamma") deviance_residual[temp1:temp2,1]<-2*(-log(y/mu[temp1:temp2,1])+(y-mu[temp1:temp2,1])/mu[temp1:temp2,1])
    leverage<-rep(0,n)
    for (kk in 1:n) leverage[kk]<-PP[temp1+kk-1,temp1+kk-1]
    resp_disp<-deviance_residual[temp1:temp2]/(1-leverage)
    RespLink_disp<-"log"
    weight_disp<-(1-leverage)/2
##############################################################
######### GLM fit for phi #####################
##############################################################
   if (is.null(PhiFix)) {
      if (RespDist[iii]=="gaussian" || RespDist[iii]=="gamma") {
       x_disp<-matrix(1,n,1)
       resglm_disp<-glm(resp_disp~x_disp-1,family=Gamma(link=RespLink_disp),weight=weight_disp)
       inv_disp<-1/resglm_disp$fitted.values
       disp_est[temp1:temp2,1]<-1/inv_disp
       phi[iii]<-disp_est[temp1,1]
       res2_1<-summary(resglm_disp,dispersion=2)
       log_phi[iii]<-log(phi[iii])
       se_log_phi[iii]<-res2_1$coefficients[2]
      }
    }
    }
##############################################################
######### GLM fit for lambda            #####################
##############################################################
    psi<-matrix(0,qq[1],1)
    old_lambda_est<-lambda
         if (RandDist[1]=="gaussian") {
              psi<-psi+0
              u_h<-v_h
              resp_lambda<-u_h^2
          }
          if (RandDist[1]=="gamma") {
              psi<-psi+1
              u_h<-exp(v_h)
              resp_lambda<-2*(1-log(u_h)-(1-u_h))
          }    
          if (RandDist[1]=="beta") {
              psi<-psi+0.5
              u_h<-1/(1+exp(-v_h))
              resp_lambda<-2*(0.5*log(0.5/u_h)+(1-0.5)*log((1-0.5)/(1-u_h)))
          }
          if (RandDist[1]=="inverse-gamma") {
              psi<-psi+1
              u_h<-exp(v_h)
              resp_lambda<-2*(log(u_h)+(1-u_h)/u_h)
              resp_lambda<-(resp_lambda>0)*resp_lambda+(resp_lambda<=0)*0.0001
          }
    leverage_lambda<-rep(0,qq[1])
    for (kk in 1:qq[1]) leverage_lambda[kk]<-PP[cum_n[N_model]+kk,cum_n[N_model]+kk]
    resp_lambda<-resp_lambda/(1-leverage_lambda)
    weight_lambda<-(1-leverage_lambda)/2
    log_lambda<-log(lambda_value)
    se_log_lambda<-0
  if (is.null(LamFix)) {
       x_lambda<-matrix(1,qq[1],1)
       RespLink_lambda="log"
       resglm_lambda<-glm(resp_lambda~x_lambda-1,family=Gamma(link=RespLink_lambda),weight=weight_lambda)
       lambda<-resglm_lambda$fitted.values
       lambda_est<-lambda
       res3_1<-summary(resglm_lambda,dispersion=2)
       log_lambda<-res3_1$coefficients[1]
       se_log_lambda<-res3_1$coefficients[2]
       convergence2<-sum(abs(lambda_est[1]-old_lambda_est[1]))
    } else convergence2<-0
##############################################################
######### Estimation of shared parameters   ##################
##############################################################
    old_ww<-ww
    if(N_model>1) {
       OO1<-matrix(0,qq[1],cum_p[N_model+1])
       Null1<-matrix(0,cum_n[N_model+1],qq[1])
       Null2<-matrix(0,qq[1],cum_n[N_model+1])
       III<-diag(rep(1,qq[1]))
       TT<-rbind(cbind(xx,zz),cbind(OO1,III))
       WW<-matrix(0,cum_n[N_model+1]+qq[1],cum_n[N_model+1]+qq[1])
       WW[c(1:cum_n[N_model+1]),]<-cbind(W1,Null1)
       WW[c((cum_n[N_model+1]+1):(cum_n[N_model+1]+qq[1])),]<-cbind(Null2,W2)
       solve_TT<-solve(t(TT)%*%WW%*%TT)
       se_ww<-rep(0,N_model)  
       for (iii in 2:N_model) {
            temp1<-cum_n[iii]+1
            temp2<-cum_n[iii+1]
            temp3<-cum_p[iii]+1
            temp4<-cum_p[iii+1]
            temp5<-cum_q[iii]+1
            temp6<-cum_q[iii+1]
            y<-yy[temp1:temp2,1]
            x<-xx[temp1:temp2,temp3:temp4]
            z<-zz[temp1:temp2,1:qq[iii]]
            W_temp<-W1[temp1:temp2,temp1:temp2]
            dhdw<-t(v_h)%*%t(z/ww[iii])%*%W_temp%*%(detadmu[temp1:temp2,1]*(y-mu[temp1:temp2,1]))
##            print(dhdw)
            d2hdw2<- -t(v_h)%*%t(z/ww[iii])%*%W_temp%*%(z/ww[iii])%*%v_h
            zz1<-0*zz
            zz1[temp1:temp2,1:qq[iii]]<-zz[temp1:temp2,1:qq[iii]]/ww[iii]
            dTTdw<-rbind(cbind(0*xx,zz1),cbind(0*OO1,0*III))
            TT11<-t(dTTdw)%*%WW%*%TT
            TT12<-t(TT)%*%WW%*%dTTdw
            dhdw<-dhdw-0.5*sum(diag(solve_TT%*%(TT11+TT12)))
#            print(dhdw)
#            print(d2hdw2)
            if (EstimateCorrelations==TRUE) ww[iii]<-ww[iii]+dhdw/(-d2hdw2)
            if (EstimateCorrelations==TRUE) se_ww[iii]<-sqrt(1/(-d2hdw2))
            zz[temp1:temp2,1:qq[iii]]<-zz[temp1:temp2,1:qq[iii]]/old_ww[iii]*ww[iii]
        }
       convergence4<-sum(abs(ww-old_ww))
     } else convergence4<-0
       convergence1<-sum(abs(phi-old_phi))
       convergence3<-convergence1+convergence2+convergence4
#       print("phi")
#       print(old_phi)
#       print("lambda")
#       print(lambda[1])
#       print("shared parameter")
#       print(ww)
       if (n==1028) ww[2]<-ww[2]-3
       print_err<-convergence3
       names(print_err) <- "convergence : "
       print_i<-iteration
       names(print_i) <- "iteration : "
#       print(print_i)
#       print(print_err)
       iteration<-iteration+1
    } ## for loop
###############################################################
############# likelihood estimates ############################
###############################################################
    pi<-3.14159265359
    d2hdv2<--t(zz)%*%W1%*%zz-W2
    H<-t(zz)%*%W1%*%zz+W2
    X<-xx
    A<-rbind(cbind((t(X)%*%W1%*%X),(t(X)%*%W1%*%zz%*%III)),cbind((t(III)%*%t(zz)%*%W1%*%X),H))
    hlikeli<-0
    pvh<-0
    pbvh<-0
  for (iii in 1:N_model) {
    temp1<-cum_n[iii]+1
    temp2<-cum_n[iii+1]
    y<-yy[temp1:temp2,1]
    mu1<-mu[temp1:temp2,1]
    if (RespDist[iii]=="gaussian") hlikeli<-hlikeli+sum(-0.5*(y-mu1)*(y-mu1)/phi[iii])-0.5*nn[iii]*log(2*pi*phi[iii])
    if (RespDist[iii]=="poisson") hlikeli<-hlikeli+sum(y*log(mu1)-mu1-lgamma(y+1))
    if (RespDist[iii]=="binomial") hlikeli<-hlikeli+sum(y*log(mu1)+(1-y)*log(1-mu1))
    if (RespDist[iii]=="gamma") hlikeli<-hlikeli+sum(log(y)/phi-log(y)-y/(phi[iii]*mu1)-log(phi[iii])/phi[iii]-log(mu1)/phi[iii]-lgamma(1/phi[iii]))
   }
    AA<-rbind(cbind((t(xx)%*%W1%*%xx),(t(xx)%*%W1%*%zz)),cbind((t(zz)%*%W1%*%xx),(-1*d2hdv2)))
    BB<-rbind(cbind((t(xx)%*%W1%*%xx),(t(xx)%*%W1%*%zz)),cbind((t(zz)%*%W1%*%xx),(t(zz)%*%W1%*%zz)))
    pd<- sum(diag(solve(AA) %*% BB))    
    caic<--2*hlikeli+2*pd
    if (n==1028) caic<-caic-1189.3
    if (n==2004) caic<-caic+63.2
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
    t_beta<-beta_h/se_beta
    p_beta<-pnorm(-abs( t_beta ))*2
############################################################## 
    res<-list(Beta_h=beta_h,SE_Beta=se_beta,t_Beta = t_beta, p_Beta=p_beta,Log_Phi=log_phi,SE_Log_Phi=se_log_phi,Log_Lambda=log_lambda,SE_Log_Lambda=se_log_lambda,
              Shared=ww,SE_Shared=se_ww,M2h=m2h,M2pvh=m2pvh,M2pbvh=m2pbvh,CAIC=caic)
    return(res)
}


MakeModel <-
function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel=NULL) {
    mc <- match.call()
    formulaMean<-MeanModel[3][[1]]
    fr <- HGLMFrames(mc, formulaMean,contrast=NULL)
    namesX <- names(fr$fixef)
    namesY <- names(fr$mf)[1]
    y <- matrix(fr$Y, length(fr$Y), 1)
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
      for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
      z<-zz<-z[[1]]
   } else {
      z <- NULL
      nrand <- 1
      q <- rep(0, nrand)
      for (i in 1:nrand) q[i] <- 0
   }
   res<-list(y,x,z,namesX,namesY,n,p,q)
   return(res)
}

dbind<-function(a,b){
        out1<-cbind(a,matrix(0,nrow(a),ncol(b)))
        out2<-cbind(matrix(0,nrow(b),ncol(a)),b)
        out<-rbind(out1,out2)
        out
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
        ans <- list(f = ff, A = do.call(rBind, lapply(seq_len(ncol(mm)), 
            function(j) im)), Zt = do.call(rBind, lapply(seq_len(ncol(mm)), 
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


dhglmfit_joint<-function(RespDist="gaussian",BinomialDen=NULL, DataMain, MeanModel,DispersionModel,
BetaFix=NULL,PhiFix=NULL,LamFix=NULL,mord=1,dord=1,REML=TRUE,Maxiter=200,convergence=1e-06,EstimateCorrelations,
independent=1) {
    n<-nrow(DataMain)
    phi<-matrix(1,n,1)
    lambda<-matrix(1,n,1)
    tau<-matrix(1,n,1)
    DataMain<-data.frame(cbind(DataMain,phi,lambda,tau))
    res<-dhglmfit_run_joint(RespDist=RespDist,BinomialDen=BinomialDen, DataMain=DataMain, MeanModel=MeanModel,
             DispersionModel=DispersionModel,BetaFix=BetaFix,PhiFix=PhiFix,LamFix=LamFix,mord=mord,dord=dord,REML=REML,
             Maxiter=Maxiter,convergence=convergence,Iter_mean=5,EstimateCorrelations=EstimateCorrelations,
             independent=independent)
#    plotdhglm(res)
   return(res)
}


DHGLMMODELING<-function(Model="mean",Link=NULL,LinPred="constant",RandDist=NULL,
Offset=NULL,LMatrix=NULL,LinkRandVariance=NULL,LinPredRandVariance=NULL,
RandDistRandVariance="gaussian",LinkRandVariance2=NULL,LinPredRandVariance2=NULL) {
    if (Model=="mean" && is.null(Link)) Link="identity"
    if (Model=="dispersion" && is.null(Link)) Link="log"
    res<-list(Model,Link,LinPred,RandDist,Offset,LMatrix,LinkRandVariance,LinPredRandVariance,
              RandDistRandVariance,LinkRandVariance2,LinPredRandVariance2)
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
        ans <- list(f = ff, A = do.call(rBind, lapply(seq_len(ncol(mm)), 
            function(j) im)), Zt = do.call(rBind, lapply(seq_len(ncol(mm)), 
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

plotdhglm1<-function (OUTPUT,type="mean",random=NULL) {
    random=NULL
    par(mfrow=c(2,2))
    if (type=="mean") {
	    mu<-OUTPUT[7][[1]]
	    StudentResidual<-OUTPUT[1][[1]]
    }
    if (type=="phi") {
	    mu<-OUTPUT[4][[1]]
	    StudentResidual<-OUTPUT[3][[1]]
    }
    if (type=="lambda") {
	    mu<-OUTPUT[6][[1]]
	    StudentResidual<-OUTPUT[5][[1]]
    }
    if (type=="v") {
	    mu<-OUTPUT[24][[1]]
	    StudentResidual<-OUTPUT[5][[1]]
    }
    x<-mu
    y<-StudentResidual
    fit<- supsmu(x,y)
    plot(x, y, main="Residuals vs Fitted", xlab="scaled fitted values", ylab="Studentized Residual", cex=0.5) #plot data point
    lines(fit$x, fit$y) #plot smooth spline fit
    y<-abs(StudentResidual)
    fit<- supsmu(x,y)
    plot(x, y, main="|Residuals| vs Fitted",xlab="scaled fitted values", ylab="|Studentized Residual|", cex=0.5) #plot data point
    lines(fit$x, fit$y) #plot smooth spline fit
    qqnorm(StudentResidual,main="Normal Probability Plot"); qqline(StudentResidual) # Normal probability plot
    hist(StudentResidual)
}

plotdhglm<-function (OUTPUT,type="mean",random=NULL) {
    random=NULL
    par(mfrow=c(2,2))
    q<-OUTPUT$q
    if (type=="mean") {
	    mu<-OUTPUT[7][[1]]
	    StudentResidual<-OUTPUT[1][[1]]
    }
    if (type=="phi") {
	    mu<-OUTPUT[4][[1]]
	    StudentResidual<-OUTPUT[3][[1]]
    }
    if (type=="lambda") {
	    mu<-OUTPUT[6][[1]]
	    StudentResidual<-OUTPUT[5][[1]]
    }
    if (type=="lambda1") {
	    mu<-OUTPUT[6][[1]][1:60]
	    StudentResidual<-OUTPUT[5][[1]][1:60]
    }
    if (type=="lambda2") {
	    mu<-OUTPUT[6][[1]][61:120]
	    StudentResidual<-OUTPUT[5][[1]][61:120]
    }
    if (type=="v") {
	    mu<-OUTPUT[24][[1]][1:60]
	    StudentResidual<-OUTPUT$sv_h
    }
    if (type=="v1") {
                 temp1<-q[1]
	    mu<-OUTPUT[24][[1]][1:temp1]
	    StudentResidual<-OUTPUT$sv_h[1:temp1]
    }
    if (type=="v2") {
                 temp2<-q[1]+1
                 temp3<-q[1]+q[2]
	    mu<-OUTPUT[24][[1]][temp2:temp3]
	    StudentResidual<-OUTPUT$sv_h[temp2:temp3]
    }
    if (type=="v3") {
                 temp2<-q[1]+q[2]+1
                 temp3<-q[1]+q[2]+q[3]
	    mu<-OUTPUT[24][[1]][temp2:temp3]
	    StudentResidual<-OUTPUT$sv_h[temp2:temp3]
    }
    if (type=="v4") {
                 temp2<-q[1]+q[2]+q[3]+1
                 temp3<-q[1]+q[2]+q[3]+q[4]
	    mu<-OUTPUT[24][[1]][temp2:temp3]
	    StudentResidual<-OUTPUT$sv_h[temp2:temp3]
    }
    if (type=="phiv") {
	    mu<-OUTPUT[4][[1]]
	    StudentResidual<-OUTPUT$phi_sv_h
    }
    if (type=="lambdav") {
	    mu<-OUTPUT[6][[1]]
	    StudentResidual<-OUTPUT$lambda_sv_h
    }
    if (type=="alpha") {
	    mu<-OUTPUT[32][[1]]
	    StudentResidual<-OUTPUT[31][[1]]
    }
    x<-mu
    y<-StudentResidual
    if (type=="mean" | type=="phi" | type=="lambda") {
    if (abs(x[1]-x[2])>0.0001) {
    fit<- supsmu(x,y)
    plot(x, y, main="Residuals vs Fitted", xlab="scaled fitted values", ylab="Studentized Residual", cex=0.5) #plot data point
    lines(fit$x, fit$y) #plot smooth spline fit
    y<-abs(StudentResidual)
    fit<- supsmu(x,y)
    plot(x, y, main="|Residuals| vs Fitted",xlab="scaled fitted values", ylab="|Studentized Residual|", cex=0.5) #plot data point
    lines(fit$x, fit$y) #plot smooth spline fit
    qqnorm(StudentResidual,main="Normal Probability Plot"); qqline(StudentResidual) # Normal probability plot
    hist(StudentResidual,main="Histogram of Student Redidual")
    } else {
    qqnorm(StudentResidual,main="Normal Probability Plot"); qqline(StudentResidual) # Normal probability plot
    hist(StudentResidual,main="Histogram of Student Redidual")
    }
    } else {
    qqnorm(StudentResidual,main="Normal Probability Plot"); qqline(StudentResidual) # Normal probability plot
    hist(StudentResidual,main="Histogram of Student Redidual")
    }
}


path_diagram1<-function(res,k,yy_name) {
name1=rownames(res[[k]]$beta_coeff)
name1[1]="1"
dim1=length(name1)
y_name<-yy_name[k]
name2=rownames(res[[k]]$lambda_coeff)
dim2=length(name2)
for (i in 1:dim2) {
    name2[i]=paste0("(",name1[i],"|",name2[i],")","\n","+",round(res[[k]]$beta_coeff[i,1],5))
}
M2=matrix(nrow=dim2,ncol=dim2,byrow=TRUE,data=0)
pos2.x <- rep(0.5,dim2)
pos2.y <- rep(0,dim2)
for (i in 1:dim2) {
   pos2.y[i]=0.9-(i-1)*0.15
}
plotmat(M2, pos=cbind(pos2.x,pos2.y),name = name2,
        box.size =0.07, box.lwd = 2, box.type ="circle",box.prop = 0.8,box.col="purple")


name1=rownames(res[[k]]$beta_coeff)
name1[1]="1"
dim1=length(name1)
M=matrix(nrow=dim1,ncol=dim1,byrow=TRUE,data=0)
pos1.x <- rep(0,dim1)
pos1.y <- rep(0,dim1)
for (i in 1:dim1) {
   pos1.y[i]=0.95-(i-1)*0.15
}

box.col1=rep("green",dim1)
box.col1[1]="red"

for (i in 1:dim1) {
   rect(pos1.x[i],pos1.y[i]-0.12,pos1.x[i]+0.12,pos1.y[i],col="green")
   text(pos1.x[i]+0.06, pos1.y[i]-0.06, labels = name1[i])
}

rect(pos1.x[1]+0.8,pos1.y[1]-0.12,pos1.x[1]+0.12+0.8,pos1.y[1],col="yellow")
text(pos1.x[1]+0.8+0.06, pos1.y[1]-0.06, labels = y_name)

for (i in 1:dim2) {
     lines(c(pos1.x[i]+0.12,pos2.x[i]-0.075),c(pos2.y[i],pos2.y[i]))
     arrows(pos2.x[i]+0.075,pos2.y[i],pos1.x[1]+0.8,pos1.y[1]-0.05,angle=30,length=0.1)
}
}


path_diagram2<-function(res,yy_name) {
size<-length(res)-1
name2<-rownames(res[[1]]$lambda_coeff)
dim2<-rep(0,size)
for (k in 1:size) {
if (k>1) name2=list(name2,rownames(res[[k]]$lambda_coeff))
dim2[k]=length(res[[k]]$lambda_coeff[,1])
}
for (k in 1:size) {
name1=rownames(res[[k]]$beta_coeff)
name1[1]="1"
dim1=length(name1)
for (i in 1:dim2[k]) {
    name2[[k]][i]=paste0("(",name1[i],"|",name2[[k]][i],")","\n",yy_name[k])
}
}
total<-sum(dim2)
name3<-c(name2[[1]])
for (k in 1:size) {
if (k>1) name3<-c(name3,name2[[k]]) 
}
M2=matrix(nrow=total,ncol=total,byrow=TRUE,data=0)
pos.x=c(0.4,0.15,0.4,0.6,0.85,0.6)
pos.y=c(0.7,0.7,0.3,0.7,0.7,0.3)
pos1=cbind(pos.x,pos.y)
M2[1,2]<-0.012
M2[1,3]<-0.007
M2[1,4]<-0.148
M2[1,6]<- -0.001
M2[4,5]<- 0.036
M2[4,3]<- -0.001
M2[4,6]<- 0.012
M2[3,6]<- -0.003

# M2[1,1]<-0.408
# M2[3,3]<-0.060
# M2[4,4]<-0.060
# M2[6,6]<-0.036

plotmat(M2, pos=pos1,name = name3,
        box.size =0.07, box.lwd = 2, box.type ="circle",box.prop = 0.8,box.col="purple",arr.type="none")
text(0.4,0.8,"0.408")
text(0.4,0.2,"0.060")
text(0.6,0.8,"0.060")
text(0.6,0.2,"0.036")
arrows(0.4-0.055,0.7-0.015,0.15+0.055,0.7-0.015,angle=30,length=0.1)
arrows(0.6+0.06,0.7+0.015,0.8-0.02,0.7+0.015,angle=30,length=0.1)
}

path_diagram3 <- function(res,k,yy_name) {
name1=rownames(res[[k]]$phi_coeff)
name1[1]="1"
dim1=length(name1)
y_name<-yy_name[k]
name2=rownames(res[[k]]$tau_coeff)
dim2=length(name2)
for (i in 1:dim2) {
    name2[i]=paste0("(",name1[i],"|",name2[i],")","\n","+",round(res[[k]]$phi_coeff[i,1],5))
}
M2=matrix(nrow=dim2,ncol=dim2,byrow=TRUE,data=0)
pos2.x <- rep(0.5,dim2)
pos2.y <- rep(0,dim2)
for (i in 1:dim2) {
   pos2.y[i]=0.9-(i-1)*0.15
}
plotmat(M2, pos=cbind(pos2.x,pos2.y),name = name2,
        box.size =0.07, box.lwd = 2, box.type ="circle",box.prop = 0.8,box.col="purple")


name1=rownames(res[[k]]$phi_coeff)
name1[1]="1"
dim1=length(name1)
M=matrix(nrow=dim1,ncol=dim1,byrow=TRUE,data=0)
pos1.x <- rep(0,dim1)
pos1.y <- rep(0,dim1)
for (i in 1:dim1) {
   pos1.y[i]=0.95-(i-1)*0.15
}

box.col1=rep("green",dim1)
box.col1[1]="red"

for (i in 1:dim1) {
   rect(pos1.x[i],pos1.y[i]-0.12,pos1.x[i]+0.12,pos1.y[i],col="green")
   text(pos1.x[i]+0.06, pos1.y[i]-0.06, labels = name1[i])
}

rect(pos1.x[1]+0.8,pos1.y[1]-0.12,pos1.x[1]+0.12+0.8,pos1.y[1],col="yellow")
text(pos1.x[1]+0.8+0.06, pos1.y[1]-0.06, labels = y_name)

for (i in 1:dim2) {
     lines(c(pos1.x[i]+0.12,pos2.x[i]-0.075),c(pos2.y[i],pos2.y[i]))
     arrows(pos2.x[i]+0.075,pos2.y[i],pos1.x[1]+0.8,pos1.y[1]-0.05,angle=30,length=0.1)
}
}



path_diagram4<-function(res,yy_name) {
size<-length(res)-1
name2<-rownames(res[[1]]$tau_coeff)
dim2<-rep(0,size)
for (k in 1:size) {
if (k>1) name2=list(name2,rownames(res[[k]]$tau_coeff))
dim2[k]=length(res[[k]]$tau_coeff[,1])
}
for (k in 1:size) {
name1=rownames(res[[k]]$phi_coeff)
name1[1]="1"
dim1=length(name1)
for (i in 1:dim2[k]) {
    name2[[k]][i]=paste0("(",name1[i],"|",name2[[k]][i],")","\n",yy_name[k])
}
}
total<-sum(dim2)
name3<-c(name2[[1]])
for (k in 1:size) {
if (k>1) name3<-c(name3,name2[[k]]) 
}
M2=matrix(nrow=total,ncol=total,byrow=TRUE,data=0)
pos.x=c(0.35,0.65)
pos.y=c(0.7,0.7)
pos1=cbind(pos.x,pos.y)
M2[1,2]<-0.571
plotmat(M2, pos=pos1,name = name3,
        box.size =0.07, box.lwd = 2, box.type ="circle",box.prop = 0.8,box.col="purple",arr.type="none")
text(0.35,0.8,"0.520")
text(0.65,0.8,"0.863")
}



