library(rlist)
library(R.utils)
library(fields)


criterion <- function(x,m,data=X,Y=Y){
  Xselect <- data[,as.vector(x)]
  res.svd <- svd(t(Xselect)%*%Y,nu=1,nv=1)
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


organize_list <- function(obj,p){
  if(dim(obj)[1]!=p) stop("not all model visited")
  list.cPLS <- vector(mode = "list", length = p)
  for (i in 1:p){
    list.cPLS[[i]] <- as.vector(which(obj[i,-c(1:3)]==1))
  }
  return(list.cPLS)
}


####################
## Power method + gradient 
####################
#r1.tilde.t <- function(u,t,A) 2*u*t*(A%*%matrix(u,ncol=1))



r1.tilde.t <- function(u,t,A) 2*c(u)*c(A%*%matrix(c(t*u),ncol=1))

numerical.grad.bis <- function(t,Y,X,h,pos){
  matT <- diag(t)
  Xt <- X%*%matT
  Mt <- t(Xt)%*%Y
  At <- t(Mt)%*%Mt  
  res <- eigen(At)
  maTh <- matT
  maTh[pos,pos] <- matT[pos,pos]+h
  Xt <- X%*%maTh
  Mt <- t(Xt)%*%Y
  At <- t(Mt)%*%Mt  
  resph <- eigen(At)
  maTh[pos,pos] <- matT[pos,pos]-h
  Xt <- X%*%maTh
  Mt <- t(Xt)%*%Y
  At <- Mt%*%t(Mt)  
  resmh <- eigen(At)
  grad <- (resph$values[1]-resmh$values[1])/(2*h)
  result<- list(delta2=res$values[1],grad=grad)#/(dim(X)[1])**2)#/sqrt(sum(t**2)))
  return(result)
}


numerical.grad <- function(t,Y,X,h,pos){
  matT <- diag(t)
  Xt <- X%*%matT
  Mt <- t(Xt)%*%Y
  At <- Mt%*%t(Mt)  
  res <- eigen(At)
  maTh <- matT
  maTh[pos,pos] <- matT[pos,pos]+h
  Xt <- X%*%maTh
  Mt <- t(Xt)%*%Y
  At <- Mt%*%t(Mt)  
  resph <- eigen(At)
  maTh[pos,pos] <- matT[pos,pos]-h
  Xt <- X%*%maTh
  Mt <- t(Xt)%*%Y
  At <- Mt%*%t(Mt)  
  resmh <- eigen(At)
  grad <- (resph$values[1]-resmh$values[1])/(2*h)
  result<- list(delta2=res$values[1],grad=grad)#/(dim(X)[1])**2)#/sqrt(sum(t**2)))
  return(result)
}


PLS1.grad <- function(t,Y,X,lambda=100){
  n <- dim(X)[1]
  matT <- diag(t)
  Xt <- X%*%matT
  Mt <- t(Xt)%*%matrix(Y,ncol=1)
  norm.Mt <- (sum(Mt**2))^0.5
  grad <- t*(Mt**2)/(n*norm.Mt)
  w <- map.t.to.w(t) 
  graddelta.2 <- (grad+lambda)*2*w*exp(-w**2)
  res <- list(grad=graddelta.2,t=t)
  return(res)
}

r1.tilde.t.v <- function(v,t,M){
  Mv <- c(M%*%matrix(v,ncol=1))
  return(2*c(t)*Mv*Mv)
}

power.grad.new <- function(t,Y,X,tol=0.001,Niter=20,seed=10,lambda=100){
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  M <- t(X)%*%Y
  
  #Xt <- t(c(t)%d*%t(X))
  MtranpT2 <- t(c(t**2)%d*%M)
  
  #Mt <- t(Xt)%*%Y
  At <- MtranpT2%*%M
  
  #print(head(At))
  #At <- t(M)%*%diag(c(t**2))%*%M
  #print(head(At))
  set.seed(seed)
  v <- rnorm(dim(Y)[2])
  v <- v/sqrt(sum(v**2))
  norm.sol <- 1
  Niter <- 20
  k <-0
  
  #Bt <- t(t%d*%A)
  Gt <- matrix(0,ncol=p,nrow=q)
  
  temp.Mv <- M%*%matrix(v,ncol=1)
  g1t <- t(M)%*%diag(c(t*c(temp.Mv)))
  
  while((k<Niter) & (norm.sol> tol)) {
    k<- k+1
    
    
    
    temp.vec <- c(At%*%matrix(v,ncol=1))
    temp.norm.Atv <- sqrt(sum(temp.vec**2))
    temp.Mv <- M%*%matrix(v,ncol=1)
    temp <- temp.vec/temp.norm.Atv
    
    Ft <- 2*(t(M)%*%diag(c(t*c(temp.Mv))))/temp.norm.Atv
    dt <- 0.5*temp.vec/(temp.norm.Atv^3)
    gt <- 2*c(t(Gt)%*%(At%*%(matrix(temp.vec,ncol=1))))+2*(t*c(M%*%matrix(v,ncol=1)))*c(M%*%matrix(temp.vec,ncol=1)) 
    Ht <- t(outer(as.vector(gt),as.vector(dt)))
    Gt<- (At%*%Gt)/temp.norm.Atv+Ft-Ht 
    
    norm.sol <- sqrt(sum((v-temp)**2))
    delta1tes <- sum(v*temp.vec)
    v <- temp
  }
  
  #Zt <- t(c(u)%d*%t(Bt))
  #temp.diag <- c(A%*%matrix(t*u,ncol=1))
  temp.vec <- c(At%*%matrix(v,ncol=1))
  temp.norm.Atv <- sqrt(sum(temp.vec**2))
  temp.Mv <- M%*%matrix(v,ncol=1)
  #temp.norm.Atu <- sqrt(sum((At%*%matrix(u,ncol=1))**2))
  Ft <- 2*(t(M)%*%diag(c(t*c(temp.Mv))))/temp.norm.Atv
  dt <- 0.5*temp.vec/(temp.norm.Atv^3)
  gt <- 2*c(t(Gt)%*%(At%*%(matrix(temp.vec,ncol=1))))+2*c(t*c(M%*%matrix(v,ncol=1)))*c(M%*%matrix(temp.vec,ncol=1)) 
  Ht <- t(outer(as.vector(gt),as.vector(dt)))
  Gt<- (At%*%Gt)/temp.norm.Atv+Ft-Ht 
  
  #Gt<- (At%*%Gt)/temp.norm.Atu+Ft-Ht
  #r2 <- c(2*t(Gt)%*%At%*%matrix(u,ncol=1))
  
  r2bis <- c(2*t(Gt)%*%(At%*%matrix(v,ncol=1)))#/2
  r1bis <- r1.tilde.t.v(v,t,M)#/2
  print(dim(Gt))
  print(gt)
  print(Ht)
  print(Ft)
  print(Gt)
  print(dt)
  print(temp.norm.Atv)
  #print(At)
  #print(At%*%matrix(v,ncol=1))
  #print(sqrt(delta1tes)*At%*%matrix(v,ncol=1))
  #print(c(r1bis,r2bis))
  #print(delta1tes)
  w <- map.t.to.w(t) 
  graddelta.2 <- -((1/n^2)*(r1bis+r2bis)-lambda)*2*w*exp(-w**2)
  graddelta.t <- (r1bis+r2bis)
  res <- list(grad=graddelta.2,vtilde=v,delta1=delta1tes,t=t,iter=k,graddelta.t=graddelta.t)
  return(res)
}


myouter <- function(x,y){
  pp <- length(x)
  qq <- length(y)
  h <- matrix(0,ncol=pp,nrow=qq)
  
  for (i in 1: pp){
    h[,i] <- x[i]*y 
  }
  return(h)
}



power.grad.old <- function(t,Y,X,tol=0.001,Niter=20,seed=10,lambda=100){
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  M <- t(X)%*%Y
  
  #Xt <- t(c(t)%d*%t(X))
  MtranpT2 <- t(c(t**2)%d*%M)
  
  #Mt <- t(Xt)%*%Y
  At <- MtranpT2%*%M
  
  #print(head(At))
  #At <- t(M)%*%diag(c(t**2))%*%M
  #print(head(At))
  set.seed(seed)
  v <- rnorm(dim(Y)[2])
  v <- v/sqrt(sum(v**2))
  norm.sol <- 1
  Niter <- 20
  k <-0
  
  #Bt <- t(t%d*%A)
  Gt <- matrix(0,ncol=p,nrow=q)
  Htbis <- matrix(0,ncol=p,nrow=q)
  while((k<Niter) & (norm.sol> tol)) {
    k<- k+1
    #temp.diag <- c(A%*%matrix(t*u,ncol=1))
    #temp.vec <- c(At%*%matrix(u,ncol=1))
    #temp.norm.Atu <- sqrt(sum(temp.vec**2))
    #temp <- temp.vec/temp.norm.Atu
    
    temp.vec <- c(At%*%matrix(v,ncol=1))
    temp.norm.Atv <- sqrt(sum(temp.vec**2))
    temp.Mv <- M%*%matrix(v,ncol=1)
    temp <- temp.vec/temp.norm.Atv
    #print("temp vec")
    #print(c(temp.vec,k))
    Ft <- 2*(t(M)%*%diag(c(t*c(temp.Mv))))/temp.norm.Atv
    #print("Ft times  tempvec")
    #print(Ft%*%matrix(temp.vec,ncol=1))
    dt <- 0.5*temp.vec/(temp.norm.Atv^3)
    gt <- 2*c(t(Gt)%*%(At%*%(matrix(temp.vec,ncol=1))))#+2*(c(t)*c(M%*%matrix(v,ncol=1)))*c(M%*%matrix(temp.vec,ncol=1)) 
    gt <- gt + 4*c(t)*c(temp.Mv)*c(M%*%matrix(temp.vec,ncol=1))
    
    Ht <- t(outer(as.vector(gt),as.vector(dt)))
    #Htbis <- myouter(as.vector(gt),as.vector(dt))
    #print("bullshit")
    #print((At%*%Gt%*%matrix(temp.vec,ncol=1))/temp.norm.Atv)
    Gt<- (At%*%Gt)/temp.norm.Atv+Ft-Ht 
    #print("Ft Ht times  tempvec")
    #print((Ft-Ht)%*%matrix(temp.vec,ncol=1))
    
    # print("dimGt")
    # print(dim(Gt))
    # print("gt")
    # print(gt)
    #print("Ht")
    #print(Ht)
    #print("Htbis")
    #print(Htbis)
    # print("Gt")
    # print(Gt)
    # print("dt")
    # print(dt)
    # print("temp.norm")
    
    norm.sol <- sqrt(sum((v-temp)**2))
    print(norm.sol)
    delta1tes <- sum(v*temp.vec)
    v <- temp
    r2bis <- c(2*t(Gt)%*%(At%*%matrix(v,ncol=1)))#/2
    r1bis <- r1.tilde.t.v(v,t,M)#/2
    #print(c(r1bis,r2bis))
    
  }
  
  #Zt <- t(c(u)%d*%t(Bt))
  #temp.diag <- c(A%*%matrix(t*u,ncol=1))
  temp.vec <- c(At%*%matrix(v,ncol=1))
  temp.norm.Atv <- sqrt(sum(temp.vec**2))
  temp.Mv <- M%*%matrix(v,ncol=1)
  #temp.norm.Atu <- sqrt(sum((At%*%matrix(u,ncol=1))**2))
  Ft <- 2*(t(M)%*%diag(c(t*c(temp.Mv))))/temp.norm.Atv
  dt <- 0.5*temp.vec/(temp.norm.Atv^3)
  gt <- 2*c(t(Gt)%*%(At%*%(matrix(temp.vec,ncol=1))))#+2*c(t*c(M%*%matrix(v,ncol=1)))*c(M%*%matrix(temp.vec,ncol=1)) 
  temp.gt.2 <- 2*c(t*c(M%*%matrix(v,ncol=1)))*c(M%*%matrix(temp.vec,ncol=1)) 
  
  gt <- gt + 4*c(t)*c(temp.Mv)*c(M%*%matrix(temp.vec,ncol=1))
  
  Ht <- t(outer(as.vector(gt),as.vector(dt)))
  Gt<- (At%*%Gt)/temp.norm.Atv+Ft-Ht 
  
  
  
  
  #Gt<- (At%*%Gt)/temp.norm.Atu+Ft-Ht
  #r2 <- c(2*t(Gt)%*%At%*%matrix(u,ncol=1))
  
  r2bis <- c(2*t(Gt)%*%(At%*%matrix(v,ncol=1)))#/2
  r1bis <- r1.tilde.t.v(v,t,M)#/2
  
  #print(v)
  #print(At)
  #print(At%*%matrix(v,ncol=1))
  #print(sqrt(delta1tes)*At%*%matrix(v,ncol=1))
  #print(c(r1bis,r2bis))
  #print(delta1tes)
  w <- map.t.to.w(t) 
  graddelta.2 <- -((1/n^2)*(r1bis+r2bis)-lambda)*2*w*exp(-w**2)
  graddelta.t <- (r1bis+r2bis)
  print(c(r1bis,r2bis,k))
  res <- list(grad=graddelta.2,vtilde=v,delta1=delta1tes,t=t,iter=k,graddelta.t=graddelta.t)
  return(res)
}

power.grad <- function(t,Y,X,tol=0.001,Niter=20,seed=10,lambda=100){
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  M <- t(X)%*%Y
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  M <- t(X)%*%Y
  
  #Xt <- t(c(t)%d*%t(X))
  MtranpT2 <- t(c(t**2)%d*%M)
  
  At <- MtranpT2%*%M
  
  set.seed(seed)
  v <- rnorm(dim(Y)[2])
  v <- v/sqrt(sum(v**2))
  
  norm.sol <- 1
  Niter <- 20
  k <-0
  
 
  while((k<Niter) & (norm.sol> tol)) {
    k<- k+1
    temp.vec <- c(At%*%matrix(v,ncol=1))
    #print(temp.vec)
    temp.norm.Atv <- sqrt(sum(temp.vec**2))
    temp <- temp.vec/temp.norm.Atv
    norm.sol <- sqrt(sum((v-temp)**2))
    
    delta1tes <- sum(v*temp.vec)
    v <- temp
   
  }
  
  r1bis <- r1.tilde.t.v(v,t,M)
  w <- map.t.to.w(t) 
  graddelta.2 <- -((1/n^2)*(r1bis)-lambda)*2*w*exp(-w**2)
  graddelta.t <- r1bis
  
  res <- list(grad=graddelta.2,vtilde=v,delta1=delta1tes,t=t,iter=k,graddelta.t=graddelta.t)
  return(res)
}

power.grad.ut <- function(t,A,Y,X,tol=0.001,Niter=20,seed=10,lambda=100){
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  Xt <- t(c(t)%d*%t(X))
  Mt <- t(Xt)%*%Y
  At <- Mt%*%t(Mt)  
  set.seed(seed)
  u <- rnorm(dim(X)[2])
  u <- u/sqrt(sum(u**2))
  norm.sol <- 1
  Niter <- 20
  k <-0
  
  Bt <- t(t%d*%A)
  Gt <- matrix(0,ncol=p,nrow=p)
  while((k<Niter) & (norm.sol> tol)) {
    k<- k+1
    print(paste("iteration",k))
    temp.diag <- c(A%*%matrix(t*u,ncol=1))
    temp.vec <- c(At%*%matrix(u,ncol=1))
    temp.norm.Atu <- sqrt(sum(temp.vec**2))
    temp <- temp.vec/temp.norm.Atu
    print(c(temp.norm.Atu,k))
    Zt <- t(c(u)%d*%t(Bt))
    
    Ft <- (Zt+diag(temp.diag))/temp.norm.Atu
    dt <- 0.5*temp.vec/(temp.norm.Atu^3)
    print("Ft times  tempvec")
    print(Ft%*%matrix(temp.vec,ncol=1))
    gt <- 2*(t(Gt)%*%(At%*%(matrix(temp.vec,ncol=1))))+2*(t(Zt)%*%matrix(temp.vec,ncol=1)+matrix(temp.diag*temp.vec,ncol=1))
    Ht <- t(outer(as.vector(gt),as.vector(dt)))
    print("Ft - Ht times  tempvec")
    print((Ft-Ht)%*%matrix(temp.vec,ncol=1))
    print("bullshit")
    print((At%*%Gt%*%matrix(temp.vec,ncol=1))/temp.norm.Atu)
    
    Gt<- (At%*%Gt)/temp.norm.Atu+Ft-Ht        
    norm.sol <- sqrt(sum((u-temp)**2))
    delta1tes <- sum(u*temp.vec)
    u <- temp
    r1 <- r1.tilde.t(u,t,A)
    r2 <- c(2*t(Gt)%*%At%*%matrix(u,ncol=1))
    print(c(r1,r2))
    
  }
  
  Zt <- t(c(u)%d*%t(Bt))
  temp.diag <- c(A%*%matrix(t*u,ncol=1))
  temp.vec <- c(At%*%matrix(u,ncol=1))
  temp.norm.Atu <- sqrt(sum(temp.vec**2))
  
  #temp.norm.Atu <- sqrt(sum((At%*%matrix(u,ncol=1))**2))
  Ft <- (Zt+diag(temp.diag))/temp.norm.Atu
  dt <- 0.5*temp.vec/(temp.norm.Atu^3)
  gt <- 2*(t(Gt)%*%(At%*%(matrix(temp.vec,ncol=1))))+2*(t(Zt)%*%matrix(temp.vec,ncol=1)+matrix(temp.diag*temp.vec,ncol=1))
  Ht <- t(outer(as.vector(gt),as.vector(dt)))
  Gt<- (At%*%Gt)/temp.norm.Atu+Ft-Ht 
  
  print("dimGt")
  print(dim(Gt))
  print("gt")
  print(gt)
  print("Ht")
  print(Ht)
  print("Ft")
  print(Ft)
  print("Gt")
  print(Gt)
  print("dt")
  print(dt)
  print("temp.norm")
  print(temp.norm.Atu)
  
  r1 <- r1.tilde.t(u,t,A)
  r2 <- c(2*t(Gt)%*%At%*%matrix(u,ncol=1))
  print(c(r1,r2))
  w <- map.t.to.w(t) 
  graddelta.2 <- -((1/n^2)*(r1+r2)-lambda)*2*w*exp(-w**2)
  graddelta.t <- (r1+r2)
  res <- list(grad=graddelta.2,utilde=u,delta1=delta1tes,t=t,iter=k,graddelta.t=graddelta.t)
  return(res)
}




loss.pls <- function(t,X,Y,lambda=0){
  Xt <- X%*%diag(t)
  Mt <- t(Xt)%*%Y
  res <- svd(Mt,nu=1,nv=1)
  result <- -(1/dim(X)[1])*res$d[1]+lambda*sum(t)
 # result1 <- (1/dim(X)[1]**2)*power.grad(t,n=10,A,Y,X)$delta1-lambda*sum(t)
  return(result)
}


loss.model <- function(s,X,Y){
  X <- X[,s]
  M <- t(X)%*%Y
  res <- svd(M,nu=1,nv=1)
  result <- -(1/dim(X)[1])*res$d[1]
  return(result)
}


loss.model.penal <- function(s,X,Y,t,lambda){
  n <- dim(X)[1]
  X <- X[,s]
  M <- t(X)%*%Y
  res <- svd(M,nu=1,nv=1)
  result <- (-(1/n)*res$d[1])**2+lambda*sum(t)
  return(result)
}


ADAM.PLS1 <- function(X,Y,lambda,tau=0.5,Niter=200,alpha=0.001,psy=c(0.9,0.999),epoch=10,tol=0.0001,trunc=0.001,t0=NULL){
  n <- dim(X)[1]
  p <- dim(X)[2]
  indice <- 1:p
  Xnew <- X[,indice]
  
  Mnew <- t(Xnew)%*%matrix(Y,ncol=1)
  
  tab.t <- NULL
  if(is.null(t0)) t0 <- rep(0.5,p)
  
  w <- map.t.to.w(t0)
  c <- 10e-8
  #w <- w - alpha*u/sqrt(v+c)
  u <- rep(0,p)
  v <- rep(0,p)
  u.tilde <- rep(0,p)
  v.tilde <- rep(0,p)
  t.new <- t0
  i <- 0
  counter <- 0
  pnew <- p
  wnew <- w
  p.new.iter <- p
  while (i<Niter&counter< epoch){
    
    if(pnew<p.new.iter){
      Xnew <- X[,indice]
      pnew <- length(indice)
      
      Mnew <- t(Xnew)%*%Y
      Anew <- Mnew%*%t(Mnew)
    }
    wnew <- w[indice]
    p.new.iter <- pnew
    i <- i+1
    t.old <- t.new[indice] 
    res.grad <- PLS1.grad(t=c(t.old),Y=Y,X=Xnew,lambda  = lambda )
    gradw <- res.grad$grad
    
    u[indice] <-psy[1]*u[indice]-(1-psy[1])*gradw
    v[indice] <-psy[2]*v[indice]+(1-psy[2])*(gradw*gradw)
    
    u.tilde[indice] <- u[indice]/(1-psy[1]**i)
    v.tilde[indice] <- v[indice]/(1-psy[2]**i)
    w[indice] <- w[indice] + alpha*u.tilde[indice]/sqrt(v.tilde[indice]+c)
    
    t.new <- map.w.to.t(w)
    tab.t <- rbind(tab.t,t.new)
    if(!is.null(trunc)){
      if(any(t.new[indice]<trunc)){
        indice <- which(t.new>trunc)
        t.new[-c(indice)] <- 0
        w[-c(indice)] <- 0
        pnew <- length(indice)
      }
      if(is.null(pnew)) break
    }
    
    maxnorm <- max(abs(t.old-t.new))
    if(maxnorm < tol){counter <- counter + 1} else counter <-0
  }
  
  
  if(i<Niter) convergence <- TRUE else convergence <- FALSE
  t <- map.w.to.t(w)
  s <- as.vector((t>tau))
  
  if(sum(s)==0|sum(s)>n) beta.hat.s <- 0 else{
    beta.hat.s <- 1}
  
  
  result <- list(index.cov=which(s==TRUE),s=s,t=t,w=w,Niter=i,beta.hat.s=beta.hat.s,convergence=convergence,result.t=tab.t,lam=lambda)
  class(result) <- "PLS_BSS"
  return(result)
}


ADAM.PLS.ut <- function(X,Y,delta,lambda,tau=0.5,Niter=200,alpha=0.001,psy=c(0.9,0.999),epoch=10,tol=0.0001,trunc=0.001,t0=NULL){
  n <- dim(X)[1]
  p <- dim(X)[2]
  indice <- 1:p
  Xnew <- X[,indice]
  
  Mnew <- t(Xnew)%*%Y
  Anew <- Mnew%*%t(Mnew)
  tab.t <- NULL
  if(is.null(t0)) t0 <- rep(0.5,p)
  
  w <- map.t.to.w(t0)
  c <- 10e-8
  #w <- w - alpha*u/sqrt(v+c)
  u <- rep(0,p)
  v <- rep(0,p)
  u.tilde <- rep(0,p)
  v.tilde <- rep(0,p)
  t.new <- t0
  i <- 0
  counter <- 0
  pnew <- p
  wnew <- w
  p.new.iter <- p
  
  while (i<Niter&counter< epoch){
   # print(c(i,p))
    if(pnew<p.new.iter){
      Xnew <- X[,indice]
      pnew <- length(indice)
     
      Mnew <- t(Xnew)%*%Y
      Anew <- Mnew%*%t(Mnew)
    }
    wnew <- w[indice]
    p.new.iter <- pnew
    i <- i+1
    t.old <- t.new[indice] 
    res.grad <- power.grad(t=c(t.old),A=Anew,Y=Y,X=Xnew,lambda  = lambda )
    gradw <- res.grad$grad
  
    u[indice] <-psy[1]*u[indice]-(1-psy[1])*gradw
    v[indice] <-psy[2]*v[indice]+(1-psy[2])*(gradw*gradw)
    
    u.tilde[indice] <- u[indice]/(1-psy[1]**i)
    v.tilde[indice] <- v[indice]/(1-psy[2]**i)
    w[indice] <- w[indice] + alpha*u.tilde[indice]/sqrt(v.tilde[indice]+c)
    
    t.new <- map.w.to.t(w)
    tab.t <- rbind(tab.t,t.new)
    if(!is.null(trunc)){
      if(any(t.new[indice]<trunc)){
        indice <- which(t.new>trunc)
        t.new[-c(indice)] <- 0
        w[-c(indice)] <- 0
        pnew <- length(indice)
      }
      if(is.null(pnew)) break
    }
    
    maxnorm <- max(abs(t.old-t.new))
    if(maxnorm < tol){counter <- counter + 1} else counter <-0
  }
  
  
  if(i<Niter) convergence <- TRUE else convergence <- FALSE
  t <- map.w.to.t(w)
  s <- as.vector((t>tau))
  
  if(sum(s)==0|sum(s)>n) beta.hat.s <- 0 else{
    beta.hat.s <- 1}
  
  
  result <- list(index.cov=which(s==TRUE),s=s,t=t,w=w,Niter=i,beta.hat.s=beta.hat.s,convergence=convergence,result.t=tab.t,lam=lambda)
  class(result) <- "PLS_BSS"
  return(result)
}

ADAM.PLS <- function(X,Y,delta,lambda,tau=0.5,Niter=200,alpha=0.001,psy=c(0.9,0.999),epoch=10,tol=0.0001,trunc=0.001,t0=NULL,collect=10){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  indice <- 1:p
  Xnew <- X[,indice]
  
  Mnew <- t(Xnew)%*%Y
  Anew <- Mnew%*%t(Mnew)
  tab.t <- NULL
  if(is.null(t0)) t0 <- rep(0.5,p)
  
  w <- map.t.to.w(t0)
  c <- 10e-8
  #w <- w - alpha*u/sqrt(v+c)
  u <- rep(0,p)
  v <- rep(0,p)
  u.tilde <- rep(0,p)
  v.tilde <- rep(0,p)
  t.new <- t0
  i <- 0
  counter <- 0
  pnew <- p
  wnew <- w
  p.new.iter <- p
  
  while (i<Niter&counter< epoch){
    # print(c(i,p))
    if(pnew<p.new.iter){
      Xnew <- X[,indice]
      pnew <- length(indice)
      
      
    }
    wnew <- w[indice]
    p.new.iter <- pnew
    i <- i+1
    t.old <- t.new[indice] 
    res.grad <- power.grad(t=c(t.old),Y=Y,X=Xnew,lambda  = lambda )
    gradw <- res.grad$grad
    
    u[indice] <-psy[1]*u[indice]-(1-psy[1])*gradw
    v[indice] <-psy[2]*v[indice]+(1-psy[2])*(gradw*gradw)
    
    u.tilde[indice] <- u[indice]/(1-psy[1]**i)
    v.tilde[indice] <- v[indice]/(1-psy[2]**i)
    w[indice] <- w[indice] + alpha*u.tilde[indice]/sqrt(v.tilde[indice]+c)
    
    t.new <- map.w.to.t(w)
    
    if((i%%collect)==0){tab.t <- rbind(tab.t,t.new)}
    
    if(!is.null(trunc)){
      if(any(t.new[indice]<trunc)){
        indice <- which(t.new>trunc)
        t.new[-c(indice)] <- 0
        w[-c(indice)] <- 0
        pnew <- length(indice)
      }
      if(is.null(pnew)) break
    }
    
    maxnorm <- max(abs(t.old-t.new))
    if(maxnorm < tol){counter <- counter + 1} else counter <-0
  }
  tab.t <- rbind(tab.t,t.new)
  
  if(i<Niter) convergence <- TRUE else convergence <- FALSE
  t <- map.w.to.t(w)
  s <- as.vector((t>tau))
  
  if(sum(s)==0|sum(s)>n) beta.hat.s <- 0 else{
    beta.hat.s <- 1}
  
  
  result <- list(index.cov=which(s==TRUE),s=s,t=t,w=w,Niter=i,beta.hat.s=beta.hat.s,convergence=convergence,result.t=tab.t,lam=lambda)
  class(result) <- "PLS_BSS"
  return(result)
}




best.set.combss.PLS <- function(x){
  res <- matrix(FALSE,ncol=length(x),nrow=length(x))
  for(i in 1:length(x)){
    res[(i:length(x)),x[i]] <- TRUE 
  }
  return(res)
}


best.subset.unique.t <- function(result.t){
  #p <- dim(X)[2]
  #############
  ## Add model with t result 
  ############
  res.order <- t(apply(abs(result.t),MARGIN=1,FUN=function(x) order(x,decreasing = TRUE)))
  res.order.unique <- unique(res.order)
  res <- list(result.t=res.order.unique)
  return(res)
  
}


best.subset.combPLS.k <- function(result.t,X,Y,Xtest,Ytest){
  
  p <- dim(X)[2]
  #############
  ## Add model with t result 
  ############
  res.order <- t(apply(abs(result.t),MARGIN=1,FUN=function(x) order(x,decreasing = TRUE)))
  res.order.unique <- unique(res.order)
  print("toto0")
  res.list <- apply(res.order.unique,MARGIN=1,FUN=best.set.combss.PLS,simplify = FALSE)
  res.list <- list.rbind(res.list)
  final.select <- unique(res.list)
  print("toto1")
  #############
  ## Add model with betat  result 
  ############
  
  #res.order <- t(apply(abs(result.betat),MARGIN=1,FUN=function(x) order(x,decreasing = TRUE)))
  #res.order.unique <- unique(res.order)
  #res.list <- apply(res.order.unique,MARGIN=1,FUN=best.set.combss,simplify =FALSE)
  #res.list <- list.rbind(res.list)
  #res.list.unique.betat <- unique(res.list)
  
  #final.select <- rbind(res.list.unique.t,res.list.unique.betat)
  #final.select <- res.list.unique.t
  dim(final.select) -> nb.visited.model
  RSS <- t(apply(final.select,MARGIN=1,FUN=cov.train,Xtrain=X,Y=Y,Xtest=Xtest,Ytest=Ytest))
  RSS.indice <- cbind(RSS,1:dim(RSS)[1])
  RSS.indice.order <- RSS.indice[order(RSS.indice[,1]),]
  print("toto2")
  k.unique <- unique(RSS.indice.order[,1])
  best.subset.res <- NULL
  for (k in k.unique){
    temp <- RSS.indice.order[which(RSS.indice.order[,1]==k),]
    if(class(temp)[1]=="numeric"){
      best.subset.res <- rbind(best.subset.res, temp) 
    }else{
      best.subset.res <- rbind(best.subset.res, temp[order(temp[,2]),][1,])
    }
  }
  print("toto3")
  bestsubset.comb <-cbind(best.subset.res[,1:3],final.select[best.subset.res[,4],])
  colnames(bestsubset.comb) <- c("Dimension","Covariance Train","Covariance Test",paste("X",1:p,sep=""))
  res <- list(bestsubset.comb=bestsubset.comb,nnb.visited.model=nb.visited.model)
  return(res)
  
}



best.subset.BSS.k <- function(result.t,X,Y,K){
  #result.cPLS1$result.Adam -> result.t
  p <- dim(X)[2]
  #############
  ## Add model with t result 
  ############
  res.order <- t(apply(abs(result.t),MARGIN=1,FUN=function(x) order(x,decreasing = TRUE)))
  res.order.unique <- unique(res.order)
  
  
  desired_length <- K#dim(res.order.unique)[2] # or whatever length you want
  empty_list <- vector(mode = "list", length = desired_length)
  empty_list2 <- vector(mode = "list", length = desired_length)
  objective <- rep(0,K)
  mat.var <- matrix(NA,ncol=p,nrow=K)
  for (i in 1:desired_length){
    if(i==1){
      empty_list[[i]] <- unique(res.order.unique[,1:i])}
    else{ empty_list[[i]] <- unique(t(apply(res.order.unique[,1:i],1,sort)))}
    #print(dim(empty_list[[i]]))
    res <- (apply(matrix(empty_list[[i]],ncol=i),MARGIN=1,FUN=cov.train.BSS,Xtrain=X,Y=Y))
    if(i==1) empty_list2[[i]] <- empty_list[[i]][which.min(res)] else{
      empty_list2[[i]] <- empty_list[[i]][which.min(res),]}
    objective[i] <- min(res)
    mat.var[i,empty_list2[[i]]] <- empty_list2[[i]]
    print(c(i,objective[i]))
  }
  
  
  bestsubset.comb <-cbind(1:K,objective,mat.var)
  colnames(bestsubset.comb) <- c("Dimension","Covariance",paste("X",1:p,sep=""))
  res <- list(bestsubset.comb=bestsubset.comb,list.BSS=empty_list2)
  return(res)
  
}






cov.train.BSS <- function(x,Xtrain,Y){
  Xselect <- Xtrain[,x]
  res.svd <- svd(t(Xselect)%*%Y,nu=1,nv=1)
  covtrain <- -(res.svd$d[1])/dim(Xtrain)[1]
  return(covtrain)
}


cov.train <- function(x,Xtrain,Y,Xtest,Ytest){
  Xselect <- Xtrain[,x]
  covtest <- NA
  res.svd <- svd(t(Xselect)%*%Y,nu=1,nv=1)
  covtrain <- -(res.svd$d[1])/dim(Xtrain)[1]
  if(!is.null(Ytest)){
    Xselect.test <- Xtest[,x]
  covtest <- -cov(Xselect.test%*%matrix(res.svd$u,ncol=1),Ytest%*%matrix(res.svd$v,ncol=1))*(dim(Xtest)[1]-1)/dim(Xtest)[1]
  }
  return(c(sum(x),covtrain,covtest))
}

#result.best.subset$best.subset[10,-c(1,2)] -> x
#cov.train(x,X,Y,X.test,Y.test)


PLS.Lambda.Grid <- function(Nlammax,Kmax,X,Y,Xtest,Ytest,delta,tau,Niter,alpha,psy,epoch,tol,trunc,t0){
  ### Set of budget: Nlammax
  ### SIZE max: Kmax  
  ### compute lambda.max
  lambda.max <- (1/(dim(X)[1])**2)*(svd(t(X)%*%Y,nu=1,nv=1)$d[1])**2
  

### Run model at lambda.max
model.Adam <- ADAM.PLS(X=X,Y=Y,delta=delta,lambda=lambda.max,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,trunc=trunc,t0=t0)
result.Adam <- model.Adam$result.t
### check if lambda max is sufficient for selecting the zero model
print("step 1")
while (sum(model.Adam$s)>0 ){
  lambda.max <- 2*lambda.max
  model.Adam <- ADAM.PLS(X=X,Y=Y,delta=delta,lambda=lambda.max,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,trunc=trunc,t0=t0)
  result.Adam <- rbind(result.Adam,model.Adam$result.t)
} 
print("step 2")
vec.modelsize <- sum(model.Adam$s)
vec.lam <- lambda.max  ### define vector of lambda visited 
lam.current <- lambda.max/2 ### first iteration 
### we looking for the lambda with selected a model size greater than the target size Kmax
vec.lam <- c(vec.lam,lam.current)
#print("step 3")
model.Adam <- ADAM.PLS(X=X,Y=Y,delta=delta,lambda=lam.current,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,trunc=trunc,t0=t0)
result.Adam <- rbind(result.Adam,model.Adam$result.t)
vec.modelsize <- c(vec.modelsize,sum(model.Adam$s)) 
#print("step 4")

while (sum(model.Adam$s)<Kmax){
  lam.current <- lam.current/2
  model.Adam <- ADAM.PLS(X=X,Y=Y,delta=delta,lambda=lam.current,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,trunc=trunc,t0=t0)
  result.Adam <- rbind(result.Adam,model.Adam$result.t)
  vec.lam <- c(vec.lam,lam.current)
  vec.modelsize <- c(vec.modelsize,sum(model.Adam$s)) 
  #print(sum(model.Adam$s))
  lam.min <- lam.current
}
#print("step 5")
Nlam <- length(vec.lam)
jj<- 0
while (Nlam < Nlammax ){
  print(paste("step 6= ",jj))
  if(length(unique(vec.modelsize))>=(Kmax+1)) break
  jj <- jj +1
  #print(c("pass=",jj,length(unique(vec.modelsize)),Nlam,unique(vec.modelsize)))   
  k <- 0
  for (j in 1:(length(vec.lam)-1)){
    
    if((vec.modelsize[j+k+1]-vec.modelsize[j+k])>1){
      new.lam <- (vec.lam[j+k]+vec.lam[j+k+1])/2
      vec.lam <- insert(vec.lam,j+k+1,new.lam)
      lam.current <-new.lam  
      model.Adam <- ADAM.PLS(X=X,Y=Y,delta=delta,lambda=lam.current,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,trunc=trunc,t0=t0)
      result.Adam <- rbind(result.Adam,model.Adam$result.t)
      vec.modelsize <- insert(vec.modelsize,j+k+1,sum(model.Adam$s))
      k <- k+1
    }
    Nlam <- length(vec.lam)
    if(Nlam>=Nlammax) break
    if(length(unique(vec.modelsize))>(Kmax+1)) break
  }
} 
#print("end -1")
best.subset <- best.subset.combPLS.k(result.Adam,X,Y,Xtest,Ytest)$bestsubset.comb

#print("end")
result <- list(best.subset=best.subset)
return(result)

}

merge.model <- function(model1,model2,X,Y){
  
  n <- dim(X)[2]
  list1<- model1$best.subset
  list2 <- model2$best.subset
  list1.org <- organize_list(list1,p)
  list2.org <- organize_list(list2,p)
  
  for (j in 1:dim(list1)[1]){
    if((setequal(list1.org[[j]],list2.org[[j]])==FALSE)){
      print(j)
      a <- criterion(x=list1.org[[j]],m=n,data=X,Y=Y)
      b <- criterion(x=list2.org[[j]],m=n,data=X,Y=Y)
      if(b<a){model1$best.subset[j,] <- model2$best.subset[j,]
      list1.org[[j]]<- list2.org[[j]]}
    }
  }
  
  return(list(model=model1,list.model=list1.org))
  
}


PLS.Lambda.new <- function(Grid.lambda,X,Y,Xtest,Ytest,delta,tau,Niter,alpha,psy,epoch,tol,trunc,t0,collect,K){
  result.Adam <- NULL
  j <- 0
  for (lambda in Grid.lambda){
    j <- j+1
    print(paste("lambda",j))
    model.Adam <- ADAM.PLS(X=X,Y=Y,delta=delta,lambda=lambda,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,trunc=trunc,t0=t0,collect=collect)
    res <- best.subset.unique.t(model.Adam$result.t)
    result.Adam <- rbind(result.Adam,res$result.t)
    print(dim(result.Adam))
    
    #  if(j==20){
    #  res.order <- t(apply(abs(result.Adam),MARGIN=1,FUN=function(x) order(x,decreasing = TRUE)))
    #  res.order.unique <- unique(res.order)
    #  res.list <- apply(res.order.unique,MARGIN=1,FUN=best.set.combss.PLS,simplify = FALSE)
    #  res.list <- list.rbind(res.list)
    #  res.list.unique.t <- unique(res.list)
    #final.select <- rbind(res.list.unique.t,res.list.unique.betat)
    #  result.Adam <- unique(res.list.unique.t)
    #  } 
    
    # if(j==35){
    #    res.order <- t(apply(abs(result.Adam),MARGIN=1,FUN=function(x) order(x,decreasing = TRUE)))
    #    res.order.unique <- unique(res.order)
    #   res.list <- apply(res.order.unique,MARGIN=1,FUN=best.set.combss.PLS,simplify = FALSE)
    #  res.list <- list.rbind(res.list)
    #  res.list.unique.t <- unique(res.list)
    #final.select <- rbind(res.list.unique.t,res.list.unique.betat)
    # result.Adam <- unique(res.list.unique.t)
    #} 
    
  }
  print("end")
  print(dim(result.Adam))
  best.subset <- best.subset.BSS.k(result.Adam,X,Y,K=K)
  
  #print("end")
  result <- list(best.subset=best.subset$bestsubset.comb,result.Adam=result.Adam,list.BSS=best.subset$list.BSS)
  return(result)
  
}


PLS.Lambda <- function(Grid.lambda,X,Y,Xtest,Ytest,delta,tau,Niter,alpha,psy,epoch,tol,trunc,t0,collect){
   result.Adam <- NULL
   j <- 0
  for (lambda in Grid.lambda){
    j <- j+1
    print(paste("lambda",j))
  model.Adam <- ADAM.PLS(X=X,Y=Y,delta=delta,lambda=lambda,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,trunc=trunc,t0=t0,collect=collect)
  result.Adam <- rbind(result.Adam,model.Adam$result.t)
  print(dim(result.Adam))
  #best.subset <- best.subset.combPLS.k(result.Adam,X,Y,Xtest,Ytest)$bestsubset.comb
  
  
#  if(j==20){
#  res.order <- t(apply(abs(result.Adam),MARGIN=1,FUN=function(x) order(x,decreasing = TRUE)))
#  res.order.unique <- unique(res.order)
#  res.list <- apply(res.order.unique,MARGIN=1,FUN=best.set.combss.PLS,simplify = FALSE)
#  res.list <- list.rbind(res.list)
#  res.list.unique.t <- unique(res.list)
  #final.select <- rbind(res.list.unique.t,res.list.unique.betat)
#  result.Adam <- unique(res.list.unique.t)
#  } 
  
 # if(j==35){
#    res.order <- t(apply(abs(result.Adam),MARGIN=1,FUN=function(x) order(x,decreasing = TRUE)))
#    res.order.unique <- unique(res.order)
 #   res.list <- apply(res.order.unique,MARGIN=1,FUN=best.set.combss.PLS,simplify = FALSE)
  #  res.list <- list.rbind(res.list)
  #  res.list.unique.t <- unique(res.list)
    #final.select <- rbind(res.list.unique.t,res.list.unique.betat)
   # result.Adam <- unique(res.list.unique.t)
  #} 
  
  }
   print("end")
  best.subset <- best.subset.combPLS.k(result.Adam,X,Y,Xtest,Ytest)$bestsubset.comb
  
  #print("end")
  result <- list(best.subset=best.subset,result.Adam=result.Adam)
  return(result)
  
}





PLS.specific.dimension2 <- function(Nlammax,Kmax,X,Y,Xtest,Ytest,delta,tau,Niter,alpha,psy,epoch,tol,trunc,t0){
  ### Set of budget: Nlammax
  ### SIZE max: Kmax  
  ### compute lambda.max
  lambda.max <- (1/(dim(X)[1])**2)*(svd(t(X)%*%Y,nu=1,nv=1)$d[1])**2
  
  
  ### Run model at lambda.max
  model.Adam <- ADAM.PLS(X=X,Y=Y,delta=delta,lambda=lambda.max,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,trunc=trunc,t0=t0)
  result.Adam <- model.Adam$result.t
  ### check if lambda max is sufficient for selecting the zero model
  #print("step 1")
  while (sum(model.Adam$s)>0 ){
    lambda.max <- 2*lambda.max
    model.Adam <- ADAM.PLS(X=X,Y=Y,delta=delta,lambda=lambda.max,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,trunc=trunc,t0=t0)
    result.Adam <- rbind(result.Adam,model.Adam$result.t)
  } 
  #print("step 2")
  vec.modelsize <- sum(model.Adam$s)
  vec.lam <- lambda.max  ### define vector of lambda visited 
  lam.current <- lambda.max/2 ### first iteration 
  ### we looking for the lambda with selected a model size greater than the target size Kmax
  vec.lam <- c(vec.lam,lam.current)
  #print("step 3")
  model.Adam <- ADAM.PLS(X=X,Y=Y,delta=delta,lambda=lam.current,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,trunc=trunc,t0=t0)
  result.Adam <- rbind(result.Adam,model.Adam$result.t)
  vec.modelsize <- c(vec.modelsize,sum(model.Adam$s)) 
  #print(c("step 4",vec.modelsize))
  
  while (sum(model.Adam$s)!=Kmax){
    lam.current <- lam.current/2
    model.Adam <- ADAM.PLS(X=X,Y=Y,delta=delta,lambda=lam.current,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,trunc=trunc,t0=t0)
    result.Adam <- rbind(result.Adam,model.Adam$result.t)
    vec.lam <- c(vec.lam,lam.current)
    vec.modelsize <- c(vec.modelsize,sum(model.Adam$s)) 
    #print(sum(model.Adam$s))
    lam.min <- lam.current
  }
  #print("step 5")
  Nlam <- length(vec.lam)
  jj<- 0
  while (Nlam < Nlammax ){
    #print(paste("step 6= ",jj))
    if(any(vec.modelsize==Kmax)) break
    jj <- jj +1
    k <- 0
    #print(c("pass=",jj,length(unique(vec.modelsize)),Nlam,unique(vec.modelsize)))   
    for (j in 1:(length(vec.lam)-1)){
      #print(c(j,k,vec.lam))
      
      if(any(vec.modelsize==Kmax)) break
      #print(c(Nlam,vec.modelsize))
      #print(vec.lam)
      
      #print(c(j,k,j+k+1,j+k,vec.modelsize[j+k+1],vec.modelsize[j+k]))
      if((vec.modelsize[j+k+1]-vec.modelsize[j+k])>1){
        #print(c(j,k,j+k+1,j+k,vec.modelsize[j+k+1],vec.modelsize[j+k]))
        new.lam <- (vec.lam[j+k]+vec.lam[j+k+1])/2
        vec.lam <- insert(vec.lam,j+k+1,new.lam)
        lam.current <-new.lam  
        model.Adam <- ADAM.PLS(X=X,Y=Y,delta=delta,lambda=lam.current,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,trunc=trunc,t0=t0)
        result.Adam <- rbind(result.Adam,model.Adam$result.t)
        vec.modelsize <- insert(vec.modelsize,j+k+1,sum(model.Adam$s))
        k <- k+1
      }
      Nlam <- length(vec.lam)
      print(c("Nlam",Nlam))
      print(vec.modelsize)
      if(Nlam>=Nlammax) break
      if(any(vec.modelsize==Kmax)) break
    }
  } 
  #print("end -1")
  best.subset <- best.subset.combPLS.k(result.Adam,X,Y,Xtest,Ytest)$bestsubset.comb
  
  #print("end")
  result <- list(best.subset=best.subset)
  return(result)
  
}


refine.k <- function(k.target,res.lambda.k,X,y,Xtest=NULL,ytest=NULL,tau,Niter,alpha,epoch,tol,CG,trunc,Ntry=100){
  lo <- which.max(res.lambda.k[,1]<k.target)
  up <- lo-1
  ksol <- 0
  lamba.lo <-res.lambda.k[up,2] 
  lamba.up <-res.lambda.k[lo,2] 
  Ntry.t <- 0
  while (ksol!=k.target & Ntry.t<Ntry){
    lambda.test <- (lamba.up+lamba.lo)/2
    
    model.combssR <- ADAM.COMBSS(X,y,delta=dim(X)[1],lambda=lambda.test,tau=tau,Niter=Niter,alpha=alpha,epoch=epoch,tol=tol,CG=CG,trunc=trunc)
    y.pred.train <- as.vector(predict.COMBSS(model.combssR,X))
    RSS.train <- mean((y-y.pred.train)**2)
    if(!is.null(Xtest)){
      y.pred.test <- as.vector(predict.COMBSS(model.combssR,Xtest))
      RSS.test <- mean((ytest-y.pred.test)**2)
    }else RSS.test <- NULL
    lambda.comb <- c(lambda.comb,lam)
    
    ksol <- sum(model.combssR$s)
    #print(c(ksol,lambda.test))
    if(ksol==k.target) break
    if(ksol>k.target) lamba.lo <- lambda.test else lamba.up <- lambda.test
    Ntry.t <- Ntry.t+1
  }
  
  result <- list(lam=lambda.test,X.select=model.combssR$s,k=k.target,RSS.train=RSS.train,RSS.test=RSS.test,Ntry=Ntry.t)
  return(result)
}



PLS.specific.subset.size <- function(Nlammax,Kmax,X,Y,Xtest,Ytest,delta,tau,Niter,alpha,psy,epoch,tol,trunc,t0){
  ### Set of budget: Nlammax
  ### target: Kmax  
  ### compute lambda.max
  k.target <- Kmax
  lambda.max <- (1/(dim(X)[1])**2)*(svd(t(X)%*%Y,nu=1,nv=1)$d[1])**2
  
  ### Run model at lambda.max
  model.Adam <- ADAM.PLS(X=X,Y=Y,delta=delta,lambda=lambda.max,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,trunc=trunc,t0=t0)
  result.Adam <- model.Adam$result.t
  ### check if lambda max is sufficient for selecting the zero model
  print("step 1")
  while (sum(model.Adam$s)>0 ){
    lambda.max <- 2*lambda.max
    model.Adam <- ADAM.PLS(X=X,Y=Y,delta=delta,lambda=lambda.max,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,trunc=trunc,t0=t0)
    result.Adam <- rbind(result.Adam,model.Adam$result.t)
  } 
  print("step 2")
  vec.modelsize <- sum(model.Adam$s)
  vec.lam <- lambda.max  ### define vector of lambda visited 
  lam.current <- lambda.max/2 
  ### first iteration 
  ### we looking for the lambda with selected a model size greater than the target size Kmax
  vec.lam <- c(vec.lam,lam.current)
  print("step 3")
  model.Adam <- ADAM.PLS(X=X,Y=Y,delta=delta,lambda=lam.current,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,trunc=trunc,t0=t0)
  result.Adam <- rbind(result.Adam,model.Adam$result.t)
  vec.modelsize <- c(vec.modelsize,sum(model.Adam$s)) 
  print(c("step 4",vec.modelsize))
  
  while (sum(model.Adam$s)<Kmax){
    lam.current <- lam.current/2
    model.Adam <- ADAM.PLS(X=X,Y=Y,delta=delta,lambda=lam.current,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,trunc=trunc,t0=t0)
    result.Adam <- rbind(result.Adam,model.Adam$result.t)
    vec.lam <- c(vec.lam,lam.current)
    vec.modelsize <- c(vec.modelsize,sum(model.Adam$s)) 
    print(c("step 4bisss",vec.modelsize))
    }
  print(c("step 5",vec.modelsize))
  if(sum(model.Adam$s)!=Kmax){
  
  Nlam <- length(vec.lam)
  
  up <- which.min(vec.modelsize<Kmax)
  lo <- up-1
  ksol <- 0
  lamba.lo <-vec.lam[up] 
  lamba.up <-vec.lam[lo] 
  Ntry.t <- 0
   while (ksol!=Kmax & Ntry.t<Nlammax){
    lambda.test <- (lamba.up+lamba.lo)/2
    
    model.Adam <- ADAM.PLS(X=X,Y=Y,delta=delta,lambda=lambda.test,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,trunc=trunc,t0=t0)
    result.Adam <- rbind(result.Adam,model.Adam$result.t)
    vec.modelsize <- insert(vec.modelsize,lo+1,sum(model.Adam$s))
    ksol <- sum(model.Adam$s)
    vec.lam <- sort(c(vec.lam,lambda.test),decreasing = TRUE)
    
    print(c(ksol,lambda.test))
    if(ksol==k.target) break
    if(ksol>k.target) lamba.lo <- lambda.test else lamba.up <- lambda.test
    Ntry.t <- Ntry.t+1
  }
  }
  #print("end -1")
  best.subset <- best.subset.combPLS.k(result.Adam,X,Y,Xtest,Ytest)$bestsubset.comb
  
  #print("end")
  result <- list(best.subset=best.subset)
  return(result)
  
}




PLS_COMBSS <- function (X, Y, ncomp, mode = "regression", 
                        keepX = rep(ncol(X), ncomp), 
                        scale = TRUE,Nlammax = 30, Kmax=NULL,Xtest=NULL,Ytest=NULL,delta=NULL,tau=0.5,Niter=1000,alpha=0.001,psy=c(0.9,0.999),epoch=10,tol=0.0001,trunc=NULL,t0=NULL) 
  
{
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  q <- ncol(Y)
  p <- ncol(X)
  n <- nrow(X)
  if(is.null(delta)) delta <- n
  if(is.null(Kmax)) Kmax <- round(p/2,digits = 1)
  X.names = dimnames(X)[[2]]
  if (is.null(X.names)) 
    X.names = paste("X", 1:p, sep = "")
  if (dim(Y)[2] == 1) 
    Y.names = "Y"
  else {
    Y.names = dimnames(Y)[[2]]
    if (is.null(Y.names)) 
      Y.names = paste("Y", 1:q, sep = "")
  }
  ind.names = dimnames(X)[[1]]
  if (is.null(ind.names)) {
    ind.names = dimnames(Y)[[1]]
    rownames(X) = ind.names
  }
  if (is.null(ind.names)) {
    ind.names = 1:n
    rownames(X) = rownames(Y) = ind.names
  }
  X.s <- scale(X, scale = scale)
  Y.s <- scale(Y, scale = scale)
  mat.c <- matrix(nrow = p, ncol = ncomp)
  mat.d <- matrix(nrow = q, ncol = ncomp)
  mat.e <- matrix(nrow = q, ncol = ncomp)
  mat.t <- matrix(nrow = n, ncol = ncomp)
  mat.u <- matrix(nrow = n, ncol = ncomp)
  mat.index <- matrix(nrow = ncomp, ncol = p)
  sparsity.x <- dim(X.s)[2] - keepX
  
  result.best.subset <- PLS.Lambda.Grid(Nlammax,Kmax,X.s,Y.s,Xtest,Ytest,delta,tau,Niter,alpha,psy,epoch,tol,trunc,t0)
  
  index.1 <- result.best.subset$best.subset[keepX[1],-c(1:3)]
  index.1 <- as.logical(index.1)
  Xselect <- X.s[,index.1]
  mat.index[1,] <- index.1
  res.svd <- svd(t(Xselect)%*%Y.s,nu=1,nv=1)
  res.deflat <- deflatation.step(X = X.s, Y = Y.s, res.svd$u, 
                                 res.svd$v, mode = mode,index=index.1)
  
  mat.c[, 1] <- res.deflat$c
  if (mode == "regression") 
    mat.d[, 1] <- res.deflat$d
  else mat.e[, 1] <- res.deflat$e
  load.u <- rep(0,p)
  load.u[index.1] <- res.svd$u
  load.v <- res.svd$v
  mat.t[, 1] <- X.s[,index.1] %*% load.u[index.1]
  mat.u[, 1] <- Y.s %*% load.v
  if (ncomp > 1) {
    for (h in 2:ncomp) {
      print(paste("compute the component h = ",h,sep=""))
      res.load <- PLS.Lambda.Grid(Nlammax,Kmax,res.deflat$X.h,res.deflat$Y.h,Xtest,Ytest,delta,tau,Niter,alpha,psy,epoch,tol,trunc,t0)
      
      index.h <- res.load$best.subset[keepX[h],-c(1:3)]
      index.h <- as.logical(index.h)
      mat.index[h,] <- index.h
      Xselect <- res.deflat$X.h[,index.h]
      res.svd <- svd(t(Xselect)%*%res.deflat$Y.h,nu=1,nv=1)
      
      
      load.u.temp <- rep(0,p)
      load.u.temp[index.h] <- res.svd$u
      load.u <- cbind(load.u,load.u.temp)
      load.v <- cbind(load.v, res.svd$v)
      mat.t[, h] <- res.deflat$X.h[,index.h] %*% res.svd$u
      mat.u[, h] <- res.deflat$Y.h %*% res.svd$v
      res.deflat <- deflatation.step(X = res.deflat$X.h, Y = res.deflat$Y.h, res.svd$u, 
                                     res.svd$v, mode = mode,index=index.h)
      mat.c[, h] <- res.deflat$c
      if (mode == "regression") 
        mat.d[, h] <- res.deflat$d
      else mat.e[, h] <- res.deflat$e
      
    }
  }
  else {
    load.u <- matrix(load.u, ncol = 1)
    load.v <- matrix(load.v, ncol = 1)
  }
  rownames(load.u) <- X.names
  rownames(load.v) <- Y.names
  cl = match.call()
  result <- list(call = cl, X = X.s, Y = Y.s, ncomp = ncomp, 
                 mode = mode, keepX = keepX, mat.c = mat.c, 
                 mat.d = mat.d, mat.e = mat.e, loadings = list(X = load.u, 
                                                               Y = load.v), variates = list(X = mat.t, Y = mat.u), 
                 names = list(X = X.names, Y = Y.names, indiv = ind.names), 
                 tol = tol,mat.index=mat.index)
  class(result) = c("sPLS", "spls", "pls")
  return(invisible(result))
}

normv <- function (x) 
  sqrt(sum(x^2))


deflatation.step <- 
  function (X, Y, u.tild.new, v.tild.new, mode,index) 
  {
    xi.h <- X[,index] %*% matrix(c(u.tild.new), ncol = 1)/((normv(u.tild.new))^2)
    w.h <- Y %*% matrix(c(v.tild.new), ncol = 1)/((normv(v.tild.new))^2)
    c.h <- t(X) %*% matrix(xi.h, ncol = 1)/((normv(xi.h))^2)
    d.rh <- t(Y) %*% matrix(xi.h, ncol = 1)/(sum(xi.h * xi.h))
    d.h <- t(Y) %*% matrix(w.h, ncol = 1)/(sum(w.h * w.h))
    X.h <- X - xi.h %*% t(c.h)
    if (mode == "regression") 
      Y.h <- Y - xi.h %*% t(d.rh)
    else Y.h <- Y - w.h %*% t(d.h)
    res <- list(X.h = X.h, Y.h = Y.h, c = c.h, d = d.rh, e = d.h)
    return(res)
  }




predict.PLS.COMBSS <- 
function (object, newdata, ...) 
{
  if (missing(newdata)) 
    stop("No new data available.")
  X = object$X
  Y = object$Y
  q = ncol(Y)
  p = ncol(X)
  if (length(dim(newdata)) == 2) {
    if (ncol(newdata) != p) 
      stop("'newdata' must be a numeric matrix with ncol = ", 
           p, " or a vector of length = ", p, ".")
  }
  if (length(dim(newdata)) == 0) {
    if (length(newdata) != p) 
      stop("'newdata' must be a numeric matrix with ncol = ", 
           p, " or a vector of length = ", p, ".")
    dim(newdata) = c(1, p)
  }
  ncomp = object$ncomp
  a = object$loadings$X
  b = object$loadings$Y
  c = object$mat.c
  means.X = attr(X, "scaled:center")
  means.Y = attr(Y, "scaled:center")
  sigma.X = attr(X, "scaled:scale")
  sigma.Y = attr(Y, "scaled:scale")
  newdata = as.matrix(newdata)
  ones = matrix(rep(1, nrow(newdata)), ncol = 1)
  B.hat = array(0, dim = c(p, q, ncomp))
  Y.hat = array(0, dim = c(nrow(newdata), q, ncomp))
  Y.hat2 = array(0, dim = c(nrow(newdata), q, ncomp))
  t.pred = array(0, dim = c(nrow(newdata), ncomp))
  variates.X = object$variates$X
  betay = list()
  for (h in 1:ncomp) {
    dd = coefficients(lm(Y ~ variates.X[, 1:h, drop = FALSE]))
    if (q == 1) {
      betay[[h]] = (dd[-1])
    }
    if (q >= 2) {
      betay[[h]] = (dd[-1, ])
    }
    W = a[, 1:h, drop = FALSE] %*% solve(t(c[, 1:h, drop = FALSE]) %*% 
                                           a[, 1:h, drop = FALSE])
    B = W %*% drop(betay[[h]])
    Y.temp = scale(newdata, center = means.X, scale = sigma.X) %*% 
      as.matrix(B)
    Y.temp2 = scale(Y.temp, center = FALSE, scale = 1/sigma.Y)
    Y.temp3 = scale(Y.temp2, center = -means.Y, scale = FALSE)
    Y.hat[, , h] = Y.temp3
    t.pred[, h] = scale(newdata, center = means.X, scale = sigma.X) %*% 
      W[, h]
    B.hat[, , h] = B
  }
  rownames(t.pred) = rownames(newdata)
  colnames(t.pred) = paste("dim", c(1:ncomp), sep = " ")
  rownames(Y.hat) = rownames(newdata)
  colnames(Y.hat) = colnames(Y)
  return(invisible(list(predict = Y.hat, variates = t.pred, 
                        B.hat = B.hat, betay = betay)))
}




