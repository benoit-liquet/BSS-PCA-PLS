source("generated_true_model_PLS.R")
source("Function_BSS_PLS2.R")
library(mvtnorm)
library(fields)
library(Rfast)
library(mixOmics)

sigma.gamma <- 1
TABdimp <- c(50,100,200,500) 
sigma.e <- 5
n <- 150
#p <- 15
q <- 10
nsimu <- 100
NSR <- 0.1

ntruepred <- 10

tab.cpls <- NULL
tab.spls <- NULL

spls.true <- 0
cpls.true <- 0



for (p in TABdimp){
  res.true <- NULL
  res.spls.crit <- NULL
  res.cpls.crit <- NULL  
    
  sparsity <- p - ntruepred
  res.MSE.spls <- matrix(0,nrow=nsimu,ncol=1)
  res.MSE.BSS <- matrix(0,nrow=nsimu,ncol=1)


  MCC.spls<-matrix(0,nrow=nsimu,ncol=1)
  F1.spls<- matrix(0,nrow=nsimu,ncol=1)
  sens.spls<-matrix(0,nrow=nsimu,ncol=1)
  spe.spls<- matrix(0,nrow=nsimu,ncol=1)

  MCC.bss<- matrix(0,nrow=nsimu,ncol=1)
  F1.bss<- matrix(0,nrow=nsimu,ncol=1)
  sens.bss<- matrix(0,nrow=nsimu,ncol=1)
  spe.bss<- matrix(0,nrow=nsimu,ncol=1)


  True_list <- vector(mode = "list", length = 1)
  True_list[[1]] <- (sparsity+1):p
  true.dim <- length(True_list[[1]])



  tab.snr <- rep(0,nsimu)
  for (kk in 1:nsimu){
    print(kk)
    set.seed(kk)
    simu <- simulation1(n,sigma.gamma,sigma.e,p,q,sparsity)
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
    
    #########################
    model <- lm(Y.s~X.s)
    (Y.s-predict(model))^2 -> RSS
    1-sum(RSS)/sum(Y.s**2) -> R2
    tab.snr[kk] <- R2/(1-R2)
    #########################
    
  
  
    scale <- TRUE
    X.s <- scale(X,center=scale,scale=scale)
    Y.s <- scale(Y,center=scale,scale=scale)
    #####################
    ### Result true model
    ntrain <- dim(X)[1]
    true.crit <- criterion(x=True_list[[1]],m=ntrain,data=X.s,Y=Y.s)
    res.true <- c(res.true,true.crit)
    ###################################
  
  
    ############
    ## With package mixOmics
    ############
    
    colnames(X.test) <- paste("X",1:p,sep="")
    spls_list.mix <- vector(mode = "list", length = 1)
    MSEP.spls.mix <- NULL
    k <- 10
    
      model.spls.mix <- mixOmics::spls(X.s,Y.s,ncomp=1,keepX=c(k),scale = scale)
      res.spls.mix <- which(model.spls.mix$loadings$X!=0)#mixOmics::selectVar(model.spls)$X$name
      spls_list.mix[[1]] <- as.vector(res.spls.mix)
      res.spls.mix <- predict(model.spls.mix,newdata=X.test)
      MSE.spls <- mean((res.spls.mix$predict[,,1] - Y.test)**2)
      res.MSE.spls[kk,1] <- MSE.spls
      
      subsetspls <- rep(FALSE,p)
      subsetspls[spls_list.mix[[1]]]<-TRUE
      result.spls <- performance.measure(True.set,subsetspls)
      MCC.spls[kk,1] <- result.spls$MCC
      F1.spls[kk,1] <- result.spls$f1
      sens.spls[kk,1] <- result.spls$sensitivity
      spe.spls[kk,1] <- result.spls$specificity
    
      spls.crit <- criterion(x=spls_list.mix[[1]],m=ntrain,data=X.s,Y=Y.s)
      res.spls.crit <- c(res.spls.crit,spls.crit)
    ###########
    ## BSS
    ###########
    
    ### component 1
      scale <- TRUE
      X.s <- scale(X,center=scale,scale=scale)
      Y.s <- scale(Y,center=scale,scale=scale)
      result.cPLS <- PLS.Lambda.Grid(Nlammax=50,Kmax=10,X=X.s,Y=as.matrix(Y.s),tau=0.5,Niter=1000,alpha=0.01,psy=c(0.9,0.999),epoch=10,tol=0.0001,t0=NULL,collect = 1,Kchoice = 10) 
      cpls_list <- result.cPLS$list.BSS
      ta <- rep(0.7,p)
      result.cPLS.2 <- PLS.Lambda.Grid(Nlammax=50,Kmax=10,X=X.s,Y=as.matrix(Y.s),tau=0.5,Niter=1000,alpha=0.01,psy=c(0.9,0.999),epoch=10,tol=0.001,t0=ta,collect = 1,Kchoice=10) 
      result.merge <- merge.model(result.cPLS,result.cPLS.2,X=X.s,Y=as.matrix(Y.s))
      ta <- rep(0.3,p)
      #result.cPLS.3 <- PLS.Lambda.Grid(Nlammax=50,Kmax=p,X=X.s,Y=as.matrix(Y.s),tau=0.5,Niter=1000,alpha=0.01,psy=c(0.9,0.999),epoch=10,tol=0.001,t0=ta,collect = 1,Kchoice=p) 
      #result.merge <- merge.model(result.merge,result.cPLS.3,X=X.s,Y=as.matrix(Y.s))
      
      cpls_list <- result.merge$list.BSS
      
  
    
    #Bss_list <- vector(mode = "list", length = p)
    MSE.Bss <- NULL
    k <- 10
      #model.bss <- BSS.PLS.object(cpls_list[[k]],TRUE,X=X.s,Y=Y.s,ncomp=1)
      model.bss <- BSS.PLS.Object(X=X.s,Y=Y.s,ncomp=1,subset = list(cpls_list[[k]]))
      res.bss <- predict.PLS.COMBSS(model.bss,newdata=X.test)
      MSE.Bss.temp <- mean((res.bss$predict[,,1] - Y.test)**2)
      res.MSE.BSS[kk,1] <- MSE.Bss.temp
      subsetbss <- rep(FALSE,p)
      subsetbss[cpls_list[[k]]]<-TRUE
      result.bss <- performance.measure(True.set,subsetbss)
      MCC.bss[kk,1] <- result.bss$MCC
      F1.bss[kk,1] <- result.bss$f1
      sens.bss[kk,1] <- result.bss$sensitivity
      spe.bss[kk,1] <- result.bss$specificity
       
      
      cpls_list <- cpls_list[[true.dim]]
      cpls.crit <- criterion(x=cpls_list,m=ntrain,data=X.s,Y=Y.s)
      res.cpls.crit <- c(res.cpls.crit,cpls.crit)
  
  } 

#############
##. Comparison
#############  


  result <- list(tabcrit.cpls=res.cpls.crit,tabcrit.spls=res.spls.crit,tabcrit.true=res.true,noise=sigma.e,MCC.bss=MCC.bss,
                 F1.bss=F1.bss,MSE.Bss=res.MSE.BSS,
                 sens.bss=sens.bss,spe.bss=spe.bss, MCC.spls=MCC.spls,
                 F1.spls=F1.spls,MSE.spls=res.MSE.spls,
                 sens.spls=sens.spls,spe.spls=spe.spls,SNR=mean(tab.snr)) 
                 


  save(result,file=paste("result_SIMU1_high_noise5",p,".Rdata",sep=""))  
  } 
  
