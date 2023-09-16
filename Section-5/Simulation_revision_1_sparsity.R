source("generated_true_model_PLS.R")
source("Function_BSS_PLS2.R")
library(mvtnorm)
library(fields)
library(glmnet)
library(Rfast)
library(mixOmics)
sigma.gamma <- 1

n <- 150
p <- 15
q <- 10
nsimu <- 100
NSR <- 0.1

sigma.e <- 6

TABsparsity <- c(3,7,9,11)

tab.cpls <- NULL
tab.spls <- NULL

for (sparsity in TABsparsity){
  
  spls.true <- 0
  cpls.true <- 0


res.MSE.spls <- matrix(0,nrow=nsimu,ncol=p)
res.MSE.lasso <- matrix(0,nrow=nsimu,ncol=1)
res.MSE.BSS <- matrix(0,nrow=nsimu,ncol=p)


min.res.MSE.spls <- matrix(0,nrow=nsimu,ncol=2)
min.MSE.MCC.spls <- matrix(0,nrow=nsimu,ncol=2)
min.MSE.F1.spls <- matrix(0,nrow=nsimu,ncol=2)
min.MSE.sens.spls <- matrix(0,nrow=nsimu,ncol=2)
min.MSE.spe.spls <- matrix(0,nrow=nsimu,ncol=2)

min.res.MSE.bss <- matrix(0,nrow=nsimu,ncol=2)
min.MSE.MCC.bss <- matrix(0,nrow=nsimu,ncol=2)
min.MSE.F1.bss <- matrix(0,nrow=nsimu,ncol=2) 
min.MSE.sens.bss <- matrix(0,nrow=nsimu,ncol=2)
min.MSE.spe.bss <- matrix(0,nrow=nsimu,ncol=2) 


MCC.spls<-matrix(0,nrow=nsimu,ncol=p)
F1.spls<- matrix(0,nrow=nsimu,ncol=p)
sens.spls<-matrix(0,nrow=nsimu,ncol=p)
spe.spls<- matrix(0,nrow=nsimu,ncol=p)

MCC.bss<- matrix(0,nrow=nsimu,ncol=p)
F1.bss<- matrix(0,nrow=nsimu,ncol=p)
sens.bss<- matrix(0,nrow=nsimu,ncol=p)
spe.bss<- matrix(0,nrow=nsimu,ncol=p)






tab.snr <- rep(0,nsimu)
for (kk in 1:nsimu){
  print(kk)
  set.seed(kk)
  simu <- simulation1(n,sigma.gamma,sigma.e,p,q,sparsity)
  #simu <- simu.true.model.paper(n,sigma.gamma,sigma.e,p,q,sparsity)
  ind <- 1:100
  X <- simu$X[ind,]  
  Y <- simu$Y[ind,]  
  X.test <- simu$X[-ind,]  
  Y.test <- simu$Y[-ind,]  
  True.set <- simu$true.set
  
  ################################
  ###### estimate signal to noise
  ###############################
  scale=TRUE
  X.s <- scale(X,center=scale,scale=scale)
  Y.s <- scale(Y,center=scale,scale=scale)
  
  model <- lm(Y.s~X.s)
  #########################
  (Y.s-predict(model))^2 -> RSS
  1-sum(RSS)/sum(Y.s**2) -> R2
  tab.snr[kk] <- R2/(1-R2)
  #########################

  
  ############################################################
  #### Compute the true subset by using an exhausitive search
  ##############################################################
  True_list <- vector(mode = "list", length = p)
  scale <- TRUE
  X.s <- scale(X,center=scale,scale=scale)
  Y.s <- scale(Y,center=scale,scale=scale)
  for(k in 1:p){
    res<-comb_n(1:p, k)
    sol <- apply(res,MARGIN=2,FUN=criterion,m=n,data=X.s,Y=Y.s)
    min.set <- which.min(sol)
    True_list[[k]] <- res[,min.set]
  }
  
  
  
  
  ############
  ## With package mixOmics
  ############
  
  ####################
  ### first component
  ###################
  colnames(X.test) <- paste("X",1:p,sep="")
  spls_list.mix <- vector(mode = "list", length = p)
  MSEP.spls.mix <- NULL
  for(k in 1:p){
    model.spls.mix <- mixOmics::spls(X.s,Y.s,ncomp=1,keepX=c(k),scale = scale)
    res.spls.mix <- which(model.spls.mix$loadings$X!=0)#mixOmics::selectVar(model.spls)$X$name
    spls_list.mix[[k]] <- as.vector(res.spls.mix)
    res.spls.mix <- predict(model.spls.mix,newdata=X.test)
    MSE.spls <- mean((res.spls.mix$predict[,,1] - Y.test)**2)
    MSEP.spls.mix <- c(MSEP.spls.mix,MSE.spls)
    res.MSE.spls[kk,k] <- MSE.spls
    subsetspls <- rep(FALSE,p)
    subsetspls[spls_list.mix[[k]]]<-TRUE
    result.spls <- performance.measure(True.set,subsetspls)
    MCC.spls[kk,k] <- result.spls$MCC
    F1.spls[kk,k] <- result.spls$f1
    sens.spls[kk,k] <- result.spls$sensitivity
    spe.spls[kk,k] <- result.spls$specificity
  }
  
  ### collect result for best model according to MSE
  min.res.MSE.spls[kk,1] <- min(MSEP.spls.mix)
  subsetspls <- rep(FALSE,p)
  subsetspls[spls_list.mix[[which.min(MSEP.spls.mix)]]]<-TRUE
  result.spls <- performance.measure(True.set,subsetspls)
  min.MSE.MCC.spls[kk,1] <- result.spls$MCC
  min.MSE.F1.spls[kk,1] <- result.spls$f1
  min.MSE.sens.spls[kk,1] <- result.spls$sensitivity
  min.MSE.spe.spls[kk,1] <- result.spls$specificity
  
  ###########
  ## BSS
  ###########
  
  ### component 1
  if(T){
    scale <- TRUE
    X.s <- scale(X,center=scale,scale=scale)
    Y.s <- scale(Y,center=scale,scale=scale)
    result.cPLS <- PLS.Lambda.Grid(Nlammax=50,Kmax=p,X=X.s,Y=as.matrix(Y.s),tau=0.5,Niter=1000,alpha=0.01,psy=c(0.9,0.999),epoch=10,tol=0.0001,t0=NULL,collect = 1,Kchoice = p) 
    cpls_list <- result.cPLS$list.BSS
    ta <- rep(0.7,p)
    result.cPLS.2 <- PLS.Lambda.Grid(Nlammax=50,Kmax=p,X=X.s,Y=as.matrix(Y.s),tau=0.5,Niter=1000,alpha=0.01,psy=c(0.9,0.999),epoch=10,tol=0.001,t0=ta,collect = 1,Kchoice=p) 
    result.merge <- merge.model(result.cPLS,result.cPLS.2,X=X.s,Y=as.matrix(Y.s))
    ta <- rep(0.3,p)
    #result.cPLS.3 <- PLS.Lambda.Grid(Nlammax=50,Kmax=p,X=X.s,Y=as.matrix(Y.s),tau=0.5,Niter=1000,alpha=0.01,psy=c(0.9,0.999),epoch=10,tol=0.001,t0=ta,collect = 1,Kchoice=p) 
    #result.merge <- merge.model(result.merge,result.cPLS.3,X=X.s,Y=as.matrix(Y.s))
    
    cpls_list <- result.merge$list.BSS
    
  }
  
  #Bss_list <- vector(mode = "list", length = p)
  MSE.Bss <- NULL
  for(k in 1:p){
    #model.bss <- BSS.PLS.object(cpls_list[[k]],TRUE,X=X.s,Y=Y.s,ncomp=1)
    model.bss <- BSS.PLS.Object(X=X.s,Y=Y.s,ncomp=1,subset = list(cpls_list[[k]]))
    res.bss <- predict.PLS.COMBSS(model.bss,newdata=X.test)
    MSE.Bss.temp <- mean((res.bss$predict[,,1] - Y.test)**2)
    MSE.Bss <- c(MSE.Bss,MSE.Bss.temp)
    
    res.MSE.BSS[kk,k] <- MSE.Bss.temp
    subsetbss <- rep(FALSE,p)
    subsetbss[cpls_list[[k]]]<-TRUE
    result.bss <- performance.measure(True.set,subsetbss)
    MCC.bss[kk,k] <- result.bss$MCC
    F1.bss[kk,k] <- result.bss$f1
    sens.bss[kk,k] <- result.bss$sensitivity
    spe.bss[kk,k] <- result.bss$specificity
    
  }
  
  
  ### collect result for best model according to MSE
  min.res.MSE.bss[kk,1] <- min(MSE.Bss)
  subsetbss <- rep(FALSE,p)
  subsetbss[cpls_list[[which.min(MSE.Bss)]]]<-TRUE
  result.bss <- performance.measure(True.set,subsetbss)
  min.MSE.MCC.bss[kk,1] <- result.bss$MCC
  min.MSE.F1.bss[kk,1] <- result.bss$f1
  min.MSE.sens.bss[kk,1] <- result.bss$sensitivity
  min.MSE.spe.bss[kk,1] <- result.bss$specificity

  
#############
##. Comparison
#############  
spls.true <- spls.true+ compare(spls_list.mix,True_list)
cpls.true <- cpls.true + compare(cpls_list,True_list)
}

result <- list(spls=spls.true ,cpls=cpls.true ,noise=sigma.e,MCC.bss=MCC.bss,
               F1.bss=F1.bss,MSE.Bss=res.MSE.BSS,
               sens.bss=sens.bss,spe.bss=spe.bss, MCC.spls=MCC.spls,
               F1.spls=F1.spls,MSE.spls=res.MSE.spls,
               sens.spls=sens.spls,spe.spls=spe.spls,SNR=mean(tab.snr)) 
               



save(result,file=paste("result_SIMU1_sparsity",sparsity,".Rdata",sep=""))  
}    
