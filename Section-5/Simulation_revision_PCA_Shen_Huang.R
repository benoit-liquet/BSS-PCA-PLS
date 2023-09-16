#####
## Code simulation revision paper: PCA
#####

source("Code_PCA_function.R")
# you need to install this package
#devtools::install_github("chavent/sparsePCA")
library(sparsePCA)


criterion.PCA <- function(x,m,data=X){
  Xselect <- data[,as.vector(x)]
  res.svd <- svd((t(Xselect)%*%Xselect)/m,nu=1,nv=1)
  cov.crit <- -(res.svd$d[1])/m
  return(cov.crit)
}

compare <- function(a,b){
  a.length <- length(a)
  result <- rep(FALSE, a.length)
  for (i in 1:a.length) {
    result[i] <- all(a[[i]]== b[[i]])
  }
  return(result)}


#####
## Simulation setting as Shen and Huan 
####

# True eigenvalues
eigenval <- c(200,100,50,50,6,5,4,3,2,1)
# True loading matrix
z1 <- c(1,1,1,1,0,0,0,0,0.9,0.9)
z2 <- c(0,0,0,0,1,1,1,1,-0.3,0.3)
Ztrue <- cbind(z1/sqrt(sum(z1^2)),z2/sqrt(sum(z2^2))) 
set.seed(1)

p <-10



true.set.1 <- c(1:4,9:10)
true.set.2 <- c(5:10)

nsimu <- 100


tabn <- c(50,100,300)

for (n in tabn){

  spca.true  <- 0
  cpca.true  <- 0
  elastpca.true <- 0
  result.spca<-result.cpca <-result.spca.e<- matrix(0,ncol=2,nrow=nsimu)
  varpca <-varspca <- varbss <- varspca.e <- matrix(0,ncol=2,nrow=nsimu)
  
for (isimu in 1:nsimu){
  print(isimu)
  A <- sparsePCA::simuPCA(n=n, Ztrue, eigenval, seed=isimu)
  X <- scale(A, center=TRUE, scale=FALSE)

  #Total variance with FULL PCA 1 and 2 component
  Totalvar <-sum(eigen(t(X)%*%X)$values)
  PCA.FULL <- mixOmics::pca(X,ncomp=2,scale=FALSE,center=TRUE)
  varpca.FULL1 <- explainedVar(B=PCA.FULL$X, Z=matrix(PCA.FULL$loadings$X[,1],ncol=1), method = "subspVar")
  varpca.FULL2 <- explainedVar(B=PCA.FULL$X, Z=PCA.FULL$loadings$X, method = "subspVar")
  varpca[isimu,1] <- varpca.FULL1/Totalvar
  varpca[isimu,2] <- varpca.FULL2/Totalvar
  
  
  ############################################################
  #### Compute the true subset by using an exhausitive search
  #### First component 
  ##############################################################
  True_list <- vector(mode = "list", length = p)
  for(k in 1:p){
    res<-combn(1:p, k)
    sol <- apply(res,MARGIN=2,FUN=criterion.PCA,m=n,data=X)
    min.set <- which.min(sol)
    True_list[[k]] <- res[,min.set]
  }
  
  
  ############
  ## Result with sparse PCA
  ############
  spca_list <- vector(mode = "list", length = p)
  per.spca <- cpev<-rep(0,10)
  for(k in 1:p){
    model.spca <- mixOmics::spca(X,ncomp=1,keepX=c(k),scale=FALSE,center=TRUE)
    per.spca[k] <- model.spca$prop_expl_var$X
    cpev[k] <- sparsePCA::explainedVar(X,Z=model.spca$loadings$X,method="subspVar")
    res.spca <- which(model.spca$loadings$X!=0)
    spca_list[[k]] <- as.vector(res.spca)
  }
  
  ############
  ## Result with sparse PCA (Zhou)
  ############
  
  spca_list_elasticnet <- vector(mode = "list", length = p)
  per.spca <- cpev<-rep(0,10)
  for(k in 1:p){
    model.spca.zhou <- elasticnet::spca(X,sparse="varnum",K=1,para=k)
    res.spca.zhou <- which(model.spca.zhou$loadings!=0)
    spca_list_elasticnet[[k]] <- as.vector(res.spca.zhou)
  }
  
  
  result.cPCA <- PCA.Lambda.Grid.adapt(Nlammax=50,Kmax=10,X,tau=0.5,Niter=1000,alpha=0.005,psy=c(0.9,0.999),epoch=10,tol=0.0001,t0=NULL,collect=1,Kchoice=NULL)
  cpca_list <- result.cPCA$list.BSS
  
  
  
  
  spca.true <- spca.true+ compare(spca_list,True_list)
  cpca.true <- cpca.true + compare(cpca_list,True_list)
  elastpca.true <- elastpca.true + compare(spca_list_elasticnet,True_list)
  
  
  subset.spca <- spca_list[[6]]
  subset.cpca <- cpca_list[[6]]
  subset.spca.e <- spca_list_elasticnet[[6]]
   
  
  
  
  if(F){
  Xselect <- X[,cpca_list[[4]]]
  res.svd <- svd((t(Xselect)%*%Xselect)/n,nu=1,nv=1)
  vpca.1 <- rep(0,p)
  vpca.1[cpca_list[[4]]] <- res.svd$v
  angleBSS[,1] <- sum(vpca.1*Ztrue[,1])
  }
  
  
  result.spca[isimu,1] <- length(intersect(true.set.1,subset.spca))
  result.cpca[isimu,1] <- length(intersect(true.set.1,subset.cpca))
  result.spca.e[isimu,1] <- length(intersect(true.set.1,subset.spca.e))
  
  ### Second component
  model.spca <- mixOmics::spca(X,ncomp=2,keepX=c(6,6),scale=FALSE,center=TRUE)
  res.spca2 <- as.vector(which(model.spca$loadings$X[,2]!=0))
  result.spca[isimu,2] <- length(intersect(true.set.2,res.spca2))
  
  ### variance explained
  varspca.1 <- explainedVar(B=model.spca$X, Z=matrix(model.spca$loadings$X[,1],ncol=1), method = "subspVar")
  varspca.2 <- explainedVar(B=model.spca$X, Z=model.spca$loadings$X, method = "subspVar")
  varspca[isimu,1] <- varspca.1/Totalvar
  varspca[isimu,2] <- varspca.2/Totalvar
  
  
  ###
  model.spca.zhou <- elasticnet::spca(X,sparse="varnum",K=2,para=c(6,6))
  res.spca.zhou2 <- which(model.spca.zhou$loadings[,2]!=0)
  result.spca.e[isimu,2] <- length(intersect(true.set.2,res.spca.zhou2))
  
  varspcae.1 <- explainedVar(B=X, Z=matrix(model.spca.zhou$loadings[,1],ncol=1), method = "subspVar")
  varspcae.2 <- explainedVar(B=X, Z=model.spca.zhou$loadings[,1:2], method = "subspVar")
  varspca.e[isimu,1] <- varspcae.1/Totalvar
  varspca.e[isimu,2] <- varspcae.2/Totalvar
  
  #BSS second dimension
  Deflat <- Deflation.step(subset.cpca,X,data=X)
  result.cPCA2  <- PCA.Lambda.Grid.adapt(Nlammax=50,Kmax=10,Deflat$Xd,tau=0.5,Niter=1000,alpha=0.005,psy=c(0.9,0.999),epoch=10,tol=0.0001,t0=NULL,Kchoice=10) 
  res.cpca2 <- result.cPCA2$list.BSS[[6]]
  result.cpca[isimu,2] <- length(intersect(true.set.2,res.cpca2))
  varbss.1 <- explainedVar(B=X, Z=matrix(Deflat$loading,ncol=1), method = "subspVar")
  
  ###
  loading.bss.2 <- rep(0,dim(X)[2])
  X.select <- X[,res.cpca2]  
  model <- svd(X.select,nv=1,nu=1)
  loading.bss.2[res.cpca2] <- model$v[,1]
  ###
  load.bss.mat <- cbind(Deflat$loading,loading.bss.2)
  
  varbss.2 <- explainedVar(B=X, Z=load.bss.mat, method = "subspVar")
  varbss[isimu,1] <- varbss.1/Totalvar
  varbss[isimu,2] <- varspca.2/Totalvar
}

BSS.pca <- c(mean(varbss[,1]),mean(result.cpca[,1])/6,1-mean(result.cpca[,1])/6,mean(varbss[,2]),mean(result.cpca[,2])/6,1-mean(result.cpca[,2])/6)
spls.pca <- c(mean(varspca[,1]),mean(result.spca[,1])/6,1-mean(result.spca[,1])/6,mean(varspca[,2]),mean(result.spca[,2])/6,1-mean(result.spca[,2])/6)
zhou.pca <- c(mean(varspca.e[,1]),mean(result.spca.e[,1])/6,1-mean(result.spca.e[,1])/6,mean(varspca.e[,2]),mean(result.spca.e[,2])/6,1-mean(result.spca.e[,2])/6)
full.pca <- c(mean(varpca[,1]),NA,NA,mean(varpca[,2]),NA,NA)
selection <- rbind(full.pca,BSS.pca,spls.pca,zhou.pca)



freq <- cbind(cpca.true,spca.true, elastpca.true)
colnames(freq) <- c("BSS","spca","spca.e")
result <- list(selection=selection,freq=freq)
save(result,file=paste("Result_PCA_",n,".RData",sep=""))

}  
