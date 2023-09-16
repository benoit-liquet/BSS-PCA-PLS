source("generated_true_model_PLS.R")
source("Function_BSS_PLS2.R")
library(mvtnorm)
library(fields)
library(glmnet)

sigma.gamma <- 1
TABsigma.e <- c(1.5,3,6) #c(3,4,5,6,7)
sparsity <- 20
n <- 150
p <- 30
q <- 10
nsimu <- 100
NSR <- 0.1

for (sigma.e in TABsigma.e){


res.MSE.spls <- matrix(0,nrow=nsimu,ncol=2)
res.MSE.lasso <- matrix(0,nrow=nsimu,ncol=1)
res.MSE.BSS <- matrix(0,nrow=nsimu,ncol=2)

MCC.spls<-matrix(0,nrow=nsimu,ncol=2)
F1.spls<- matrix(0,nrow=nsimu,ncol=2)
sens.spls<-matrix(0,nrow=nsimu,ncol=2)
spe.spls<- matrix(0,nrow=nsimu,ncol=2)

MCC.bss<- matrix(0,nrow=nsimu,ncol=2)
F1.bss<- matrix(0,nrow=nsimu,ncol=2)
sens.bss<-matrix(0,nrow=nsimu,ncol=2)
spe.bss<- matrix(0,nrow=nsimu,ncol=2)


#MCC.lasso<- matrix(0,nrow=nsimu,ncol=1)
#F1.lasso<- matrix(0,nrow=nsimu,ncol=1)

tab.snr <- rep(0,nsimu)
for (kk in 1:nsimu){
  print(kk)
  set.seed(kk)
  #simu <- simu.true.model.Y.1(n,p,q,NSR)
  #simu <- simu.true.model.2(n,sigma.gamma=1,sigma.e=1.5,p=20,q=1)
  #simu <- simu.true.model.paper(n,sigma.gamma,sigma.e,p,q,sparsity)
  simu <- simu.true.model.paper.2(n,sigma.gamma,sigma.e,p,q,sparsity)
  #simu <- simu.true.model.Chun1(n,p=30,q=10)
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
  snr <- NULL
  for (j in 1:q){
    Y1<-Y[,j] 
    model1 <- lm(Y1~X)
    nn <- dim(X)[1]
    simaX <- var(X)
    beta <- model1$coefficients[-1]
    sima.es2 <- sum(model1$residuals**2)/(nn-p-1)
    snr <- c(snr,(t(beta)%*%simaX%*%beta)/sima.es2)
  }
  model <- lm(Y.s~X.s)
  tab.snr[kk] <- mean(snr)
  #########################
  (Y.s-predict(model))^2 -> RSS
  1-sum(RSS)/sum(Y.s**2) -> R2
  tab.snr[kk] <- R2/(1-R2)
  #########################
  
  ### first component
  colnames(X.test) <- paste("X",1:p,sep="")
  spls_list.mix <- vector(mode = "list", length = p)
  MSEP.spls.mix <- NULL
  for(k in 1:p){
    model.spls.mix <- mixOmics::spls(X,Y,ncomp=1,keepX=c(k))
    res.spls.mix <- which(model.spls.mix$loadings$X!=0)#mixOmics::selectVar(model.spls)$X$name
    spls_list.mix[[k]] <- as.vector(res.spls.mix)
    res.spls.mix <- predict(model.spls.mix,newdata=X.test)
    MSEP.spls.mix <- c(MSEP.spls.mix,mean((res.spls.mix$predict[,,1] - Y.test)**2))
  }
  res.MSE.spls[kk,1] <- min(MSEP.spls.mix)
  subsetspls <- rep(FALSE,p)
  subsetspls[spls_list.mix[[which.min(MSEP.spls.mix)]]]<-TRUE
  keepx1<- which.min(MSEP.spls.mix)
  result.spls <- performance.measure(True.set,subsetspls)
  MCC.spls[kk,1] <- result.spls$MCC
  F1.spls[kk,1] <- result.spls$f1
  sens.spls[kk,1] <- result.spls$sensitivity
  spe.spls[kk,1] <- result.spls$specificity
  ### second component
  
  MSEP.spls.mix <- NULL
  for(k in 1:p){
    #model.spls.mix <- mixOmics::spls(X,Y,ncomp=2,keepX=c(7,k))
    model.spls.mix <- mixOmics::spls(X,Y,ncomp=2,keepX=c(keepx1,k))
    
    res.spls.mix <- c(which(model.spls.mix$loadings$X[,1]!=0),which(model.spls.mix$loadings$X[,2]!=0))#mixOmics::selectVar(model.spls)$X$name
    spls_list.mix[[k]] <- unique(as.vector(res.spls.mix))
    res.spls.mix <- predict(model.spls.mix,newdata=X.test)
    MSEP.spls.mix <- c(MSEP.spls.mix,mean((res.spls.mix$predict[,,2] - Y.test)**2))
  }
  
  
  res.MSE.spls[kk,2] <- min(MSEP.spls.mix)
  subsetspls <- rep(FALSE,p)
  subsetspls[spls_list.mix[[which.min(MSEP.spls.mix)]]]<-TRUE
  
  result.spls <- performance.measure(True.set,subsetspls)
  MCC.spls[kk,2] <- result.spls$MCC
  F1.spls[kk,2] <- result.spls$f1
  sens.spls[kk,2] <- result.spls$sensitivity
  spe.spls[kk,2] <- result.spls$specificity
  
  #colnames(X.test) <- colnames(model.spls.mixo$X)
  #res.spls.mixo <- predict(model.spls.mixo,newdata=X.test)
  #res.spls.mixo$predict
  #MSEP <- NULL
  #for (i in 1:ncomp){
  #  MSEP <- c(MSEP,mean((res.spls.mixo$predict[,,i] - Y.test)**2))
  #}   
  #MSEP





  
  ################
  ## Lasso model
  ################
  if(q==1){
  model.lasso <- glmnet(x=X,y=Y,intercept = TRUE,nlambda=50)
  pred.lasso <- predict(model.lasso,newx=X.test)
  MSE.test.lasso <- apply(pred.lasso,MARGIN=2,FUN=function(x) mean((Y.test-x)**2))
  ind.best <- which.min(MSE.test.lasso)
  best.lambda.lasso <- model.lasso$lambda[ind.best]
  model.lasso$s <- as.vector(model.lasso$beta[,ind.best]!=0)
  #ind <- which(model.lasso$s==TRUE)
  #model.lasso$beta.hat.s <- model.lasso$beta[ind,ind.best] #[model.lasso$s]
  result.lasso <- performance.measure(True.set,model.lasso$s)
  res.MSE.lasso[kk,1] <-min(MSE.test.lasso)
  
  MCC.lasso[kk,1] <- result.lasso$MCC
  F1.lasso[kk,1] <- result.lasso$f1
  
  ###########
  ## BSS
  ###########
  }
  
  ### component 1
  if(T){
    scale <- TRUE
    X.s <- scale(X,center=scale,scale=scale)
    Y.s <- scale(Y,center=scale,scale=scale)
    result.cPLS <- PLS.Lambda.Grid(Nlammax=50,Kmax=p,X=X.s,Y=as.matrix(Y.s),tau=0.5,Niter=1000,alpha=0.01,psy=c(0.9,0.999),epoch=10,tol=0.0001,t0=NULL,collect = 1,Kchoice = p) 
    cpls_list <- result.cPLS$list.BSS
  }
  
  #Bss_list <- vector(mode = "list", length = p)
  MSEP.Bss <- NULL
  for(k in 1:p){
    #model.bss <- BSS.PLS.object(cpls_list[[k]],TRUE,X=X.s,Y=Y.s,ncomp=1)
    model.bss <- BSS.PLS.Object(X=X,Y=Y,ncomp=1,subset = list(cpls_list[[k]]))
    
    res.bss <- predict.PLS.COMBSS(model.bss,newdata=X.test)
    MSEP.Bss <- c(MSEP.Bss,mean((res.bss$predict[,,1] - Y.test)**2))
  }
  
  
  
  res.MSE.BSS[kk,1] <- min(MSEP.Bss)
  subsetbss <- rep(FALSE,p)
  subsetbss[cpls_list[[which.min(MSEP.Bss)]]]<-TRUE
  
  result.bss <- performance.measure(True.set,subsetbss)
  MCC.bss[kk,1] <- result.bss$MCC
  F1.bss[kk,1] <- result.bss$f1
  sens.bss[kk,1] <- result.bss$sensitivity
  spe.bss[kk,1] <- result.bss$specificity
  ### Component 2
  #index.1 <- cpls_list[[10]]
  subset1 <- cpls_list[[which.min(MSEP.Bss)]] # variable for component1
  #print(subset1)
  ####
  Xselect <- X.s[,subset1]
  res.svd <- svd(t(Xselect)%*%Y.s,nu=1,nv=1)
  ###loading
  load.u <- rep(0,p)
  load.u[subset1] <- res.svd$u
  load.v <- res.svd$v
  res.deflat <- deflatation.step(X = X.s, Y = Y.s, load.u, 
                                 load.v, mode = "regression")
  
  if(T){
    result.cPLS.2 <- PLS.Lambda.Grid(Nlammax=p,Kmax=p,X=res.deflat$X.h,Y=as.matrix(res.deflat$Y.h),tau=0.5,Niter=1000,alpha=0.01,psy=c(0.9,0.999),epoch=10,tol=0.0001,t0=NULL,collect = 1,Kchoice = p) 
    cpls_list.2 <- result.cPLS.2$list.BSS
  }
  
  MSEP.Bss.2 <- NULL
  subsetlist <- vector("list",2)
  subsetlist[[1]] <-subset1
  
  for(k in 1:p){
    subsetlist[[2]] <- cpls_list.2[[k]]
    model.bss <- BSS.PLS.Object(X=X,Y=Y,ncomp=2,subset =subsetlist )
    #model.bss <- BSS.PLS.Object(subset1=cpls_list[[which.min(MSEP.Bss)]],subset2=cpls_list.2[[k]],scale=TRUE,X=X.s,Y=Y.s,ncomp=2)
    #model.bss <- BSS.PLS.object.2(subset1=1:4,subset2=cpls_list.2[[k]],scale=TRUE,X=X.s,Y=Y.s,ncomp=2)
    res.bss <- predict.PLS.COMBSS(model.bss,newdata=X.test)
    MSEP.Bss.2 <- c(MSEP.Bss.2,mean((res.bss$predict[,,2] - Y.test)**2))
  }
  concat.subset <- unique(c(cpls_list[[which.min(MSEP.Bss.2)]],cpls_list[[which.min(MSEP.Bss)]]))
  
  #print(cpls_list[[which.min(MSEP.Bss.2)]])
  res.MSE.BSS[kk,2] <- min(MSEP.Bss.2)
  subsetbssf <- rep(FALSE,p)
  subsetbssf[concat.subset]<-TRUE
  
  result.bss <- performance.measure(True.set,subsetbssf)
  MCC.bss[kk,2] <- result.bss$MCC
  F1.bss[kk,2] <- result.bss$f1
  sens.bss[kk,2] <- result.bss$sensitivity
  spe.bss[kk,2] <- result.bss$specificity
  
}

result <- list(noise=sigma.e,MCC.bss=MCC.bss,
               F1.bss=F1.bss,MSE.Bss=res.MSE.BSS,
               sens.bss=sens.bss,spe.bss=spe.bss, MCC.spls=MCC.spls,
               F1.spls=F1.spls,MSE.spls=res.MSE.spls,
               sens.spls=sens.spls,spe.spls=spe.spls,SNR=mean(tab.snr)) 



save(result,file=paste("result_SIMU2new_noise",sigma.e,".Rdata",sep=""))  
}    




