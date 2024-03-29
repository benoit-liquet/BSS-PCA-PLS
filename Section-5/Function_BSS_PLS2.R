###########
## Function to upload on GITHUB
###########

library(rlist)
library(R.utils)
library(fields)


compare <- function(a,b){
  a.length <- length(a)
  result <- rep(FALSE, a.length)
  for (i in 1:a.length) {
    result[i] <- all(a[[i]]== b[[i]])
  }
  return(result)}



r1.tilde.t.v <- function(v,t,M){
  Mv <- c(M%*%matrix(v,ncol=1))
  return(2*c(t)*Mv*Mv)
}

map.t.to.w <- function (t) sqrt(-log(1-t))
map.w.to.t <- function (w) 1-exp(-w^2)




cov.train.BSS <- function(x,Xtrain,Y){
  Xselect <- Xtrain[,x]
  res.svd <- svd(t(Xselect)%*%Y,nu=1,nv=1)
  covtrain <- -(res.svd$d[1])/dim(Xtrain)[1]
  return(covtrain)
}



criterion <- function(x,m,data=X,Y=Y){
  Xselect <- data[,as.vector(x)]
  res.svd <- svd(t(Xselect)%*%Y,nu=1,nv=1)
  cov.crit <- -(res.svd$d[1])/m
  return(cov.crit)
}


best.subset.unique.t <- function(result.t){
  res.order <- t(apply(abs(result.t),MARGIN=1,FUN=function(x) order(x,decreasing = TRUE)))
  res.order.unique <- unique(res.order)
  res <- list(result.t=res.order.unique)
  return(res)
  
}



PLS.Lambda.Grid <- function(Nlammax,Kmax,X,Y,tau,Niter,alpha,psy,epoch,tol,t0,collect,Kchoice){
  ### Set of budget: Nlammax
  ### SIZE max: Kmax  
  ### compute lambda.max
  lambda.max <- (1/(dim(X)[1])**2)*(svd(t(X)%*%Y,nu=1,nv=1)$d[1])**2
  
  
  ### Run model at lambda.max
  model.Adam <- ADAM.PLS(X=X,Y=Y,lambda=lambda.max,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,t0=t0,collect=collect)
  result.Adam <- model.Adam$result.t
  ### check if lambda max is sufficient for selecting the zero model
  #print("step 1")
  while (sum(model.Adam$s)>0 ){
    lambda.max <- 2*lambda.max
    model.Adam <- ADAM.PLS(X=X,Y=Y,lambda=lambda.max,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,t0=t0,collect=collect)
    result.Adam <- rbind(result.Adam,model.Adam$result.t)
  } 
  #print("step 2")
  vec.modelsize <- sum(model.Adam$s)
  vec.lam <- lambda.max  ### define vector of lambda visited 
  lam.current <- lambda.max/2 ### first iteration 
  ### we looking for the lambda with selected a model size greater than the target size Kmax
  vec.lam <- c(vec.lam,lam.current)
  #print("step 3")
  model.Adam <- ADAM.PLS(X=X,Y=Y,lambda=lam.current,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,t0=t0,collect=collect)
  result.Adam <- rbind(result.Adam,model.Adam$result.t)
  vec.modelsize <- c(vec.modelsize,sum(model.Adam$s)) 
  #print("step 4")
  
  while (sum(model.Adam$s)<Kmax){
    lam.current <- lam.current/2
    model.Adam <- ADAM.PLS(X=X,Y=Y,lambda=lam.current,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,t0=t0,collect=collect)
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
    if(length(unique(vec.modelsize))>=(Kmax+1)) break
    jj <- jj +1
    #print(c("pass=",jj,length(unique(vec.modelsize)),Nlam,unique(vec.modelsize)))   
    k <- 0
    for (j in 1:(length(vec.lam)-1)){
      
      if((vec.modelsize[j+k+1]-vec.modelsize[j+k])>1){
        new.lam <- (vec.lam[j+k]+vec.lam[j+k+1])/2
        vec.lam <- insert(vec.lam,j+k+1,new.lam)
        lam.current <-new.lam  
        model.Adam <- ADAM.PLS(X=X,Y=Y,lambda=lam.current,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,t0=t0,collect=collect)
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
  best.subset <- best.subset.BSS.k(result.Adam,X,Y,K=Kchoice)
  
  #print("end")
  result <- list(best.subset=best.subset$bestsubset.comb,result.Adam=result.Adam,list.BSS=best.subset$list.BSS)
  return(result)
  
}



PLS.Lambda <- function(Grid.lambda,X,Y,tau,Niter,alpha,psy,epoch,tol,t0,collect,Kchoice){
  result.Adam <- NULL
  j <- 0
  for (lambda in Grid.lambda){
    j <- j+1
    #print(paste("lambda",j))
    model.Adam <- ADAM.PLS(X=X,Y=Y,lambda=lambda,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,t0=t0,collect=collect)
    res <- best.subset.unique.t(model.Adam$result.t)
    result.Adam <- rbind(result.Adam,res$result.t)
    #print(dim(result.Adam))
    
    
  }
  #print("end")
  #print(dim(result.Adam))
  best.subset <- best.subset.BSS.k(result.Adam,X,Y,K=Kchoice)
  
  #print("end")
  result <- list(best.subset=best.subset$bestsubset.comb,result.Adam=result.Adam,list.BSS=best.subset$list.BSS)
  return(result)
  
}


merge.model <- function(model1,model2,X,Y){
  
  n <- dim(X)[2]
  list1.org <- model1$list.BSS
  list2.org <- model2$list.BSS
  
  for (j in 1:length(list1.org)[1]){
    if((setequal(list1.org[[j]],list2.org[[j]])==FALSE)){
      #print(j)
      a <- criterion(x=list1.org[[j]],m=n,data=X,Y=Y)
      b <- criterion(x=list2.org[[j]],m=n,data=X,Y=Y)
      if(b<a){model1$best.subset[j,] <- model2$best.subset[j,]
      list1.org[[j]]<- list2.org[[j]]}
    }
  }
  
  return(list(best.subset=model1$best.subset,list.BSS=list1.org))
  
}

plot.best.subset.PLS <- function(object,Kmax,name.var=NULL){
  result <- object$list.BSS[1:Kmax]
  sort(unique(unlist(result))) -> var.select
  object$best.subset[1:Kmax,-c(1:2)] -> result.BSS
  result.BSS[,var.select]-> BSS
  #pdf(file="BSS-PLS2-APPLI.pdf",width=12)
  col <- rev(gray(seq(0,0.9,length=20),alpha=0.8))
  par(mar = c(5,4,7,2) + 0.8,mfrow=c(1,1))
  BSS[is.na(BSS)] <- 0
  BSS[BSS>0] <- 1
  px <- dim(BSS)[2]
  py <- dim(BSS)[1]
  image(t(BSS),xaxt="n",yaxt="n",x=1:px,y=1:py,col=col,xlab="variable",ylab="Subset Size",main="")
  grid(lwd=1, nx=px, ny=py)
  par(las=2)
  if(!is.null(name.var)) {
    name.var <- name.var[var.select]
    axis(3,at=1:px,labels=name.var)
  }
  axis(2,at=1:py,labels=1:py)
  axis(1,at=1:px,labels=var.select)
}



best.subset.BSS.k <- function(result.t,X,Y,K){
  #result.cPLS1$result.Adam -> result.t
  p <- dim(X)[2]
  #############
  ## Add model with t result 
  ############
  res.order <- t(apply(abs(result.t),MARGIN=1,FUN=function(x) order(x,decreasing = TRUE)))
  res.order.unique <- unique(res.order)
  if(dim(res.order.unique)[1]==1){ res.order.unique <- rbind(res.order.unique,res.order.unique) ## A trick for dealin with a special case !!!
  print("same order")
  }
  desired_length <- K#dim(res.order.unique)[2] # or whatever length you want
  empty_list <- vector(mode = "list", length = desired_length)
  empty_list2 <- vector(mode = "list", length = desired_length)
  objective <- rep(0,K)
  mat.var <- matrix(NA,ncol=p,nrow=K)
  for (i in 1:desired_length){
    #print(i)
    if(i==1){
      empty_list[[i]] <- unique(res.order.unique[,1:i])}
    else{ empty_list[[i]] <- unique(t(apply(res.order.unique[,1:i],1,sort)))}
    #print(dim(empty_list[[i]]))
    res <- (apply(matrix(empty_list[[i]],ncol=i),MARGIN=1,FUN=cov.train.BSS,Xtrain=X,Y=Y))
    if(i==1) empty_list2[[i]] <- empty_list[[i]][which.min(res)] else{
      empty_list2[[i]] <- empty_list[[i]][which.min(res),]}
    objective[i] <- min(res)
    mat.var[i,empty_list2[[i]]] <- empty_list2[[i]]
    #print(c(i,objective[i]))
  }
  
  bestsubset.comb <-cbind(1:K,objective,mat.var)
  colnames(bestsubset.comb) <- c("Dimension","Covariance",paste("X",1:p,sep=""))
  res <- list(bestsubset.comb=bestsubset.comb,list.BSS=empty_list2)
  return(res)
  
}



ADAM.PLS <- function(X,Y,lambda,tau=0.5,Niter=200,alpha=0.001,psy=c(0.9,0.999),epoch=10,tol=0.0001,t0=NULL,collect=10){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  indice <- 1:p
  
  Mnew <- t(X)%*%Y
  Anew <- Mnew%*%t(Mnew)
  tab.t <- NULL
  if(is.null(t0)) t0 <- rep(0.5,p)
  
  w <- map.t.to.w(t0)
  c <- 10e-8
  
  u <- rep(0,p)
  v <- rep(0,p)
  u.tilde <- rep(0,p)
  v.tilde <- rep(0,p)
  t.new <- t0
  i <- 0
  counter <- 0

  while (i<Niter&counter< epoch){
  
    i <- i+1
    t.old <- t.new
    res.grad <- power.grad(t=c(t.old),Y=Y,X=X,lambda  = lambda )
    gradw <- res.grad$grad
    
    u <-psy[1]*u-(1-psy[1])*gradw
    v <-psy[2]*v+(1-psy[2])*(gradw*gradw)
    
    u.tilde <- u/(1-psy[1]**i)
    v.tilde <- v/(1-psy[2]**i)
    w <- w + alpha*u.tilde/sqrt(v.tilde+c)
    
    t.new <- map.w.to.t(w)
    
    if((i%%collect)==0){tab.t <- rbind(tab.t,t.new)}
    
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



power.grad <- function(t,Y,X,tol=0.001,Niter=20,seed=10,lambda=100){
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


pairwise.correlation <- function(Y,label.Y=NULL){
  par(mar=c(5,5,5,7))
  if(is.null(label.Y)) label.Y <- colnames(Y)
  colnames(Y) <- label.Y
  distance <- cor(Y)
  ##  Reverse the columns of the matrix so it will be drawn correctly.
  ncol <- ncol(distance)
  distance2 <- distance[,ncol:1]
  
  #mypalette <- brewer.pal(10,"PiYG")
  mypalette <- c("#8E0152","#C51B7D","#DE77AE","#F1B6DA","#FDE0EF","#E6F5D0","#B8E186","#7FBC41"
                 ,"#4D9221","#276419")
  
  break.seq <- seq(-1,1,by=0.2)
  image(x=1:dim(Y)[2],y=1:1:dim(Y)[2],((distance2)),zlim=c(-1,1),axes=FALSE,xlab="",ylab="",main="Pairwise correlation between the phenotypes",col=mypalette,breaks=break.seq)
  
  
  par(las=2)
  par(cex.axis=1)
  axis(2, at = 1:dim(Y)[2],colnames(distance2))
  par(cex.axis=1)
  axis(1, at = 1:dim(Y)[2],rownames(distance2))
  box()
  image.plot(x=1:dim(Y)[2],y=1:1:dim(Y)[2],((distance2)),diag=FALSE,zlim=c(-1,1),nlevel=10,xlab="",ylab="",main="Pairwise correlation between the phenotypes",col=mypalette,horizontal = FALSE,legend.only=TRUE,lab.breaks=break.seq,breaks=break.seq)
  
  for(i in 1:(ncol-1)) {
    text(i,1:(ncol-i),round(distance2[i,1:(ncol-i)],digits=2))
  }
  
  res <- list(matcor=distance)
}

plotcim.explore <- function(matX,matY) {
  labelY <- colnames(matY)  
  mat.sim <- matrix(NA,ncol=dim(matX)[2],nrow=length(labelY))
  for (i in 1:dim(matX)[2]){
    for (j in 1:length(labelY)){
      mat.sim[j,i] <- cor(matX[,i],matY[,j])
    }}
  cim(mat.sim,row.names = labelY,col.names = 1:dim(matX)[2])
}




PLS.specific.subset.size <- function(Nlammax,Kmax,X,Y,tau,Niter,alpha,psy,epoch,tol,t0,collect=1){
  ### Set of budget: Nlammax
  ### target: Kmax  
  ### compute lambda.max
  k.target <- Kmax
  lambda.max <- (1/(dim(X)[1])**2)*(svd(t(X)%*%Y,nu=1,nv=1)$d[1])**2
  
  ### Run model at lambda.max
  model.Adam <- ADAM.PLS(X=X,Y=Y,lambda=lambda.max,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,t0=t0,collect=collect)
  result.Adam <- model.Adam$result.t
  ### check if lambda max is sufficient for selecting the zero model
  #print("step 1")
  while (sum(model.Adam$s)>0 ){
    lambda.max <- 2*lambda.max
    model.Adam <- ADAM.PLS(X=X,Y=Y,lambda=lambda.max,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,t0=t0,collect=collect)
    result.Adam <- rbind(result.Adam,model.Adam$result.t)
  } 
  #print("step 2")
  vec.modelsize <- sum(model.Adam$s)
  vec.lam <- lambda.max  ### define vector of lambda visited 
  lam.current <- lambda.max/2 
  ### first iteration 
  ### we looking for the lambda with selected a model size greater than the target size Kmax
  vec.lam <- c(vec.lam,lam.current)
  #print("step 3")
  model.Adam <- ADAM.PLS(X=X,Y=Y,lambda=lam.current,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,t0=t0,collect=collect)
  result.Adam <- rbind(result.Adam,model.Adam$result.t)
  vec.modelsize <- c(vec.modelsize,sum(model.Adam$s)) 
  #print(c("step 4",vec.modelsize))
  
  while (sum(model.Adam$s)<Kmax){
    lam.current <- lam.current/2
    model.Adam <- ADAM.PLS(X=X,Y=Y,lambda=lam.current,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,t0=t0,collect=collect)
    result.Adam <- rbind(result.Adam,model.Adam$result.t)
    vec.lam <- c(vec.lam,lam.current)
    vec.modelsize <- c(vec.modelsize,sum(model.Adam$s)) 
    #print(c("step 4bisss",vec.modelsize))
  }
  #print(c("step 5",vec.modelsize))
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
      
      model.Adam <- ADAM.PLS(X=X,Y=Y,lambda=lambda.test,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,t0=t0,collect=collect)
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
  
  best.subset <- best.subset.BSS.k(result.Adam,X,Y,K=k.target)
  
  #print("end")
  result <- list(best.subset=best.subset$bestsubset.comb,result.Adam=result.Adam,list.BSS=best.subset$list.BSS)
  return(result)
  
}

  
BSS.PLS.Object <- function (X, Y, ncomp, mode = "regression", subset=lapply(1:ncomp, function(x) 1:dim(X)[2]), 
                            scale = TRUE) 
{
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  q <- ncol(Y)
  p <- ncol(X)
  n <- nrow(X)
  X.names = dimnames(X)[[2]]
  if (is.null(X.names)) X.names = paste("X", 1:p, sep = "")
  if (dim(Y)[2] == 1) Y.names = "Y" else {
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
  
  Xselect <- X.s[,subset[[1]]]
  res.svd <- svd(t(Xselect)%*%Y.s,nu=1,nv=1)
  
  ###loading
  load.u <- rep(0,p)
  load.u[subset[[1]]] <- res.svd$u
  load.v <- res.svd$v
  
  res.deflat <- deflatation.step(X = X.s, Y = Y.s, load.u, 
                                 load.v, mode = mode)
  mat.c[, 1] <- res.deflat$c
  
  if (mode == "regression") mat.d[, 1] <- res.deflat$d else mat.e[, 1] <- res.deflat$e
  
  mat.t[, 1] <- X.s %*% load.u
  mat.u[, 1] <- Y.s %*% load.v
  if (ncomp > 1) {
    for (h in 2:ncomp) {
      
      Xselect <- res.deflat$X.h[,subset[[h]]]
      res.svd <- svd(t(Xselect)%*%res.deflat$Y.h,nu=1,nv=1)
      
      ###loading
      load.u.temp <- rep(0,p)
      load.u.temp[subset[[h]]] <- res.svd$u
      load.v.temp <- res.svd$v
      
      load.u <- cbind(load.u, load.u.temp)
      load.v <- cbind(load.v, load.v.temp)
      mat.t[, h] <- res.deflat$X.h %*% load.u.temp
      mat.u[, h] <- res.deflat$Y.h %*% load.v.temp
      
      res.deflat <- deflatation.step(X = res.deflat$X.h, Y = res.deflat$Y.h, load.u.temp, 
                                     load.v.temp, mode = mode)
      
      mat.c[, h] <- res.deflat$c
      if (mode == "regression") 
        mat.d[, h] <- res.deflat$d
      else mat.e[, h] <- res.deflat$e
      
    }
  }else {
    load.u <- matrix(load.u, ncol = 1)
    load.v <- matrix(load.v, ncol = 1)
  }
  rownames(load.u) <- X.names
  rownames(load.v) <- Y.names
  cl = match.call()
  result <- list(call = cl, X = X.s, Y = Y.s, ncomp = ncomp, 
                 mode = mode, mat.c = mat.c, 
                 mat.d = mat.d, mat.e = mat.e, loadings = list(X = load.u, 
                                                               Y = load.v), variates = list(X = mat.t, Y = mat.u), 
                 names = list(X = X.names, Y = Y.names, indiv = ind.names) 
  )
  class(result) = c("sPLS", "spls", "pls")
  return(invisible(result))
}
  
deflatation.step <- 
  function (X, Y, u.tild.new, v.tild.new, mode) 
  {
    xi.h <- X %*% matrix(c(u.tild.new), ncol = 1)/((normv(u.tild.new))^2)
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

normv <- function (x) 
  sqrt(sum(x^2))


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
