#####################
## Simulation data
#####################
library(mvtnorm)

###########
## Simulation data for simulation section 5
##########

simu.true.model <- function(ntrain,ntest,sigma.gamma,sigma.e,p,q,sparsity){
  
  theta.x1 <- c(rep(0,sparsity),rep(c(1,-1),(p-sparsity)/2)) #p- sparsity odd number
  theta.y1 <- runif(q,0.5,2) # as in Bioinformatics
  #theta.y1 <- rep(-1,q)  
  
  Sigmax <- matrix(0, nrow = p, ncol = p)
  diag(Sigmax) <- sigma.e ^ 2
  Sigmay <- matrix(0, nrow = q, ncol = q)
  diag(Sigmay) <- sigma.e ^ 2
  gam1 <- rnorm(ntrain+ntest,2,sd=sigma.gamma)
  X <- matrix(gam1, ncol = 1, byrow = FALSE) %*% matrix(theta.x1,nrow = 1, byrow = TRUE) +
    mvtnorm::rmvnorm(ntrain+ntest, mean = rep(0, p), sigma =Sigmax, method = "svd")
             
  Y <- matrix(gam1, ncol = 1, byrow = FALSE) %*% matrix(theta.y1,nrow = 1, byrow = TRUE)
  + mvtnorm::rmvnorm(ntrain+ntest, mean = rep(0, q), sigma =Sigmay, method = "svd")
  
  X.train <- scale(X[1:ntrain,],center=TRUE,scale=F)
  Y.train <- scale(Y[1:ntrain,],center=TRUE,scale=F)
  
  X.test <- sweep(X[(ntrain+1):ntest,], 2, STATS = colMeans(X[1:ntrain,]))
  Y.test <- sweep(Y[(ntrain+1):ntest,], 2, STATS = colMeans(Y[1:ntrain,]))
  
  result <- list(X.train=X.train,Y.train=Y.train,X.test=X.test,Y.test=Y.test,sparsity=sparsity)
  
  }


