#####################
## Simulation data
#####################
library(mvtnorm)




### simulation in 2010 paper 
### simulation with univariate Y
### two snr 

simu.true.model.Y.1 <- function(n,p,q,SNR){
H1 <- rnorm(n,sd=5)
H2 <- rnorm(n,sd=5)
H3 <- rnorm(n,sd=5)
H <- cbind(H1,H2,H3)
n0 <- 0
n1 <- (p-q)/2
n2 <- p-q 
n3 <- p
nj0 <- c(n0,n1,n2,n3)
X <- matrix(0,ncol=p,nrow=n)



for (j in 1:3){
  rangei <- (nj0[j]+1):nj0[j+1]
  for (i in rangei){
  X[,i] <- H[,j]+rnorm(n)
  }
}



Y <- 3*H1-4*H2+rnorm(n,sd=SNR)


true.set <- rep(FALSE,p)
true.set[1:(p-q)] <- TRUE
result <- list(X=X,Y=matrix(Y,ncol=1),true.set=true.set)
}








###########
## Simulation data for simulation section 5
##########


simu.true.model.paper <- function(n,sigma.gamma,sigma.e,p,q,sparsity){
  
  theta.x1 <- c(rep(0,sparsity),rep(c(1,-1),(p-sparsity)/2)) #p- sparsity odd number
  theta.y1 <- runif(q,0.5,2) # as in Bioinformatics
  #theta.y1 <- rep(-1,q)  
  
  Sigmax <- matrix(0, nrow = p, ncol = p)
  diag(Sigmax) <- sigma.e ^ 2
  Sigmay <- matrix(0, nrow = q, ncol = q)
  diag(Sigmay) <- sigma.e ^ 2
  gam1 <- rnorm(n,2,sd=sigma.gamma)
  
  
  
  
  X <- matrix(gam1, ncol = 1, byrow = FALSE) %*% matrix(theta.x1,nrow = 1, byrow = TRUE) +
    mvtnorm::rmvnorm(n, mean = rep(0, p), sigma =Sigmax, method = "svd")
             
  Y <- matrix(gam1, ncol = 1, byrow = FALSE) %*% matrix(theta.y1,nrow = 1, byrow = TRUE)
  + mvtnorm::rmvnorm(n, mean = rep(0, q), sigma =Sigmay, method = "svd")
  
  true.set <- rep(FALSE,p)
  true.set[(sparsity+1):p] <- TRUE
  
  result <- list(X=X,Y=Y,true.set=true.set,sparsity=sparsity)
  
  }



#### need to put p=30 and sparsity =20 with give infact 30-20=true zero

simu.true.model.paper.2 <- function(n,sigma.gamma,sigma.e,p,q,sparsity){
  
  theta.x1 <- c(rep(c(1,-1),(p-sparsity)/2),rep(0,sparsity)) #p- sparsity odd number
  theta.x2 <- c(rep(0,sparsity/2),rep(c(1,-1.5),(p-sparsity)/2),rep(0,sparsity/2))
  
  
  theta.y1 <- runif(q,-1,3)
  theta.y2 <- runif(q,-1,3)
  
  
  
  Sigmax <- matrix(0, nrow = p, ncol = p)
  diag(Sigmax) <- sigma.e ^ 2
  Sigmay <- matrix(0, nrow = q, ncol = q)
  diag(Sigmay) <- sigma.e ^ 2
  gam1 <- rnorm(n,2,sd=sigma.gamma)
  
  gam1 <- runif(n,-1,3)
  gam2 <- runif(n,-1,3)
  #gam1 <- rnorm(n)
  #gam2 <- rnorm(n)
  gam <- cbind(gam1,gam2)
  
  S1 <- gam1 + rnorm(n)
  S2 <- gam2 +rnorm(n)
  
  DBtransp <- matrix(runif(2*q,0.5,10),nrow=2,ncol=10,byrow=TRUE)
  
  X <- gam %*% matrix(c(theta.x1, theta.x2),nrow = 2, ncol=p,byrow = TRUE) + mvtnorm::rmvnorm(n, mean = rep(0, p), sigma =
                                                                                       Sigmax, method = "svd")
   Y <- gam %*% DBtransp + mvtnorm::rmvnorm(n, mean = rep(0, q)
                                                                                              , sigma =Sigmay, method = "svd")
  

  true.set <- rep(FALSE,p)
  true.set[1:sparsity] <- TRUE
  
  result <- list(X=X,Y=Y,true.set=true.set,sparsity=sparsity)
  
}




simulation1 <- function(n,sigma.gamma,sigma.e,p,q,sparsity){
   theta.x1 <- c(rep(0,sparsity),rep(c(1,-1),(p-sparsity)/2)) #p- sparsity odd number
  
  theta.y1 <- runif(q,-1,3)
  
  
  
  
  Sigmax <- matrix(0, nrow = p, ncol = p)
  diag(Sigmax) <- sigma.e ^ 2
  Sigmay <- matrix(0, nrow = q, ncol = q)
  diag(Sigmay) <- sigma.e ^ 2
  gam1 <- rnorm(n,2,sd=sigma.gamma)
  
  
  gam <- cbind(gam1)#,gam2)
  
  S1 <- gam1 + rnorm(n)
  
  DBtransp <- matrix(runif(q,0.5,10),nrow=1,ncol=q,byrow=TRUE)
  
  X <- gam %*% matrix(c(theta.x1),nrow = 1, ncol=p,byrow = TRUE) + mvtnorm::rmvnorm(n, mean = rep(0, p), sigma =
                                                                                      Sigmax, method = "svd")
   Y <- gam %*% DBtransp + mvtnorm::rmvnorm(n, mean = rep(0, q)
                                           , sigma =Sigmay, method = "svd")
  
  
  true.set <- rep(TRUE,p)
  true.set[1:sparsity] <- FALSE
  
  result <- list(X=X,Y=Y,true.set=true.set,sparsity=sparsity)
  
}




performance.measure <- function(True.set,selected){
  est.set <- selected
  True.set <- as.factor(True.set)
  if(sum(est.set)==length(True.set)){
    est.set <- as.factor(est.set)
    levels(est.set) <- c(FALSE,TRUE) 
    est.set[1:length(True.set)] <- TRUE
  }else{
    est.set <- as.factor(est.set)
    levels(est.set) <- c(FALSE,TRUE)   
  }
  
  
  confusion_table <- table(True.set,est.set)
  TN <- confusion_table[1,1] #= 'TN'
  FN <- confusion_table[2,1] #= 'FN'
  FP <- confusion_table[1,2] #= 'FP'
  TP <- confusion_table[2,2] #= 'TP'
  accuracy <- (TP + TN) / sum(TP,FP,TN,FN)
  precision <- TP / (TP + FP)
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  if(TP==0){f1_score <-0}else{
    f1_score <- round((2 * precision * sensitivity) / (precision + sensitivity), 4)}
  MCC <- (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  
  
  
  precision <- round(TP / (TP + FP), 4)
  sensitivity <- round(TP / (TP + FN), 4)
  specificity <- round(TN / (TN + FP), 4)
  
  result <- list(MCC=MCC,Accuracy=accuracy,sensitivity=sensitivity,specificity=specificity,f1=f1_score)
  return(result)
}











BSS.PLS.object <- function(subset,scale=TRUE,X,Y,ncomp=1){
  #X.s <- scale(X,center=scale,scale=scale)
  #Y.s <- scale(Y,center=scale,scale=scale)
  X.s <- X
  Y.s <- Y
  p <- dim(X.s)[2]
  n <- dim(X.s)[1]
  mat.c <- matrix(nrow = p, ncol = ncomp)
  #mat.d <- matrix(nrow = q, ncol = ncomp)
  #mat.e <- matrix(nrow = q, ncol = ncomp)
  mat.t <- matrix(nrow = n, ncol = ncomp)
  mat.u <- matrix(nrow = n, ncol = ncomp)
  Xselect <- X.s[,subset]
  res.svd <- svd(t(Xselect)%*%Y.s,nu=1,nv=1)
  
  ###loading
  load.u <- rep(0,p)
  load.u[subset] <- res.svd$u
  load.v <- res.svd$v
  load.u <- matrix(load.u,ncol=1)
  load.v <- matrix(load.v,ncol=1)
  ### variates
  mat.t[, 1] <- X.s[,subset] %*% matrix(load.u[subset],ncol=1)
  mat.u[, 1] <- Y.s %*% load.v
  
  xi.h <- X.s[,subset] %*% matrix(c(res.svd$u), ncol = 1)/((normv(res.svd$u))^2) # useless xi.h is variates x
  #  w.h <- Y %*% matrix(c(res.svd$v), ncol = 1)/((normv(res.svd$v))^2)
  c.h <- t(X.s) %*% matrix(xi.h, ncol = 1)/((normv(xi.h))^2)
  mat.c[, 1] <- c.h
  
  result <- list( X = X.s, Y = Y.s, ncomp = ncomp, mat.c = mat.c, 
                  loadings = list(X = load.u, 
                                  Y = load.v), variates = list(X = mat.t, Y = mat.u))
  # names = list(X = X.names, Y = Y.names, indiv = ind.names), 
  # tol = tol,mat.index=mat.index)
  class(result) = c("sPLS")
  return(invisible(result))
  
}
