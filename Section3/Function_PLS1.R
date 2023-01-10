####
## R function for PLS1 in order to reproduce surface loss when p=2
####



map.t.to.w <- function (t) sqrt(-(log(1-t)))
map.w.to.t <- function (w) 1-exp(-(w^2))



######
## Loss for PLS1
#####

loss.pls1 <- function(t,X,Y,lambda=0){
  Xt <- X%*%diag(t)
  Mt <- t(Xt)%*%matrix(Y,ncol=1)
  result <- -(1/(dim(X)[1])**2)*(sum(Mt**2))+lambda*sum(t)
  return(result)
}


#############
## PLS1 gradient
#############

PLS1.grad <- function(t,Y,X,lambda=100){
  n <- dim(X)[1]
  matT <- diag(t)
  Xt <- X%*%matT
  Mt <- t(Xt)%*%matrix(Y,ncol=1)
  #norm.Mt <- (sum(Mt**2))^0.5
  Xy <- t(X)%*%matrix(Y,ncol=1)
  grad <- -2*t*(c(Xy)**2)/(n^2)
  w <- map.t.to.w(t) 
  grad <- grad 
  graddelta.2 <- (grad+lambda)*2*w*exp(-w**2)
  res <- list(gradient=graddelta.2,t=t,gradt=grad)
  return(res)
}




gradient.descent.PLS1 <- function(X,Y,alpha,t0,tol=0.00001,Niter=1000,epoch=100,lambda=0){
  p <- dim(X)[2]
  
  if(is.null(t0)) t0 <- rep(0.5,p)
  w <- map.t.to.w(t0)
  t.new <- t0
  counter <- 0
  i <-0
  tab.t <- t0
  while (i<Niter&counter< epoch){
    i <- i+1
    t.old <- t.new  
    
    res.grad <- PLS1.grad(t=c(t.new),Y,X,lambda  = lambda)
    
    w <- w - alpha*res.grad$gradient
    
    t.new <-  map.w.to.t(w)
    tab.t <- rbind(tab.t,c(t.new))
    maxnorm <- max(abs(t.old-t.new))
    
    if(maxnorm < tol){counter <- counter + 1} else counter <-0
  }
  
  res <- list(tab.t=tab.t,Niter=i)
  return(res)
}





#####
## Specific function for landscape loss with p=2
####

plot.landscape.path.BSS <- function(lambda,t0){
  
t1 = seq(0, 1, length.out = 50 )
t2 = seq(0, 1, length.out = 50)

z = matrix(data=NA, nrow=50, ncol=50)
for(i in 1:50)
{
  for(j in 1:50)
  {
    
    z[i,j] = loss.pls1(c(t1[i],t2[j]),X,Y,lambda = lambda)
    
  }
}

#####
## Gradient Descent 
####


alpha <- 0.001
epoch <-10
tol <- 0.00001
Niter <- 10000
if(is.null(t0)) t0 <- c(0.5,0.5)
PLS1 <- gradient.descent.PLS1(X,Y,alpha,t0=t0,tol=tol,Niter=Niter,epoch=10,lambda=lambda)
#plot(PLS1$tab.t[,2],col="red",type="l")
#lines(PLS1$tab.t[,1])

result <- NULL
for(i in 1:dim(PLS1$tab.t)[1]){
  result <- rbind(result,c(c(PLS1$tab.t[i,]),loss.pls1(c(PLS1$tab.t[i,]),X,Y,lambda=lambda)))
}


######
## Plot surface and path of the gradient descent 
#####


persp3d(t1,t2,z, theta=0, phi=0, r=2, shade=0.4, axes=TRUE,scale=TRUE, box=FALSE, 
        
        nticks=5, ticktype="detailed" , col="cyan", xlab="t1", 
        
        ylab="t2", zlab=expression(f[lambda]^PLS1~(t)), main="",zlim=c(min(z)-0.1))
segments3d(rbind(c(0.0, 0.0, min(z)),c(0.0, 0.0, loss.pls1(c(0,0),X,Y,lambda = lambda))),add=TRUE, lwd = 4,col="red") 
segments3d(rbind(c(1, 0, min(z)-0.1),c(1, 0, loss.pls1(c(1,0),X,Y,lambda = lambda))), lwd = 4,add = TRUE,col="red") 
segments3d(rbind(c(1, 1, min(z)-0.1),c(1, 1, loss.pls1(c(1,1),X,Y,lambda = lambda))), lwd = 4,add = TRUE,col="red") 
segments3d(rbind(c(0, 1, min(z)-0.1),c(0, 1, loss.pls1(c(0,1),X,Y,lambda = lambda))), lwd = 4,add = TRUE,col="red") 
points3d(result, lwd = 4,add = TRUE,col="black") 

}

