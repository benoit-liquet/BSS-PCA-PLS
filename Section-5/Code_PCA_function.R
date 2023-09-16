##################3
#### Code for Best Subset Selection PCA
##################


if(F){
  library(rmarkdown)

  render("input.Rmd", md_document(variant = "markdown_github"))
  render("Vignette_PLS2_BSS.Rmd", md_document(variant = "markdown_github"))
  render("Vignette_PCA_COMBSS.Rmd", md_document(variant = "markdown_github"))
}

plot.best.subset.PCA <- function(object,Kmax,name.var=NULL){
  result <- object$list.BSS[2:Kmax]
  sort(unique(unlist(result))) -> var.select
  object$best.subset[2:Kmax,-c(1:2)] -> result.BSS
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
  if(!is.null(name.var)) axis(3,at=1:px,labels=colnames(X[,var.select]))
  axis(2,at=1:py,labels=2:(py+1))
  axis(1,at=1:px,labels=var.select)
}

map.t.to.w <- function (t) sqrt(-log(1-t))
map.w.to.t <- function (w) 1-exp(-w^2)
r1.tilde.t <- function(u,t,A) 2*c(u)*c(A%*%matrix(c(t*u),ncol=1))

library(rlist)
library(R.utils)
library(fields)

if(F){
  X <- rmvnorm(50,mu=c(2,2),sigma=matrix(c(1,-1,-1,2),ncol=2,nrow=2))
  
  X <- rmvnorm(50,mu=rep(2,10),sigma=matrix(c(1,-1,-1,2),ncol=2,nrow=2))
  X <- rmvnorm(500,mu=rep(2,50),sigma=diag(2,50,50),ncol=2,nrow=2)
  
  X <- scale(X)
  p <- dim(X)[2]
  n <- dim(X)[1]
  M <- t(X)%*%X
  A <- M/n
  t0 <- rep(0.5,p)  
  #t0 <- rep(0.99,3)  
  start_time <- Sys.time()
  #power.grad.PCA <- function(t,A,X,tol=0.001,Niter=20,seed=10,lambda){
  Rprof()
  library(profvis)
  profvis({
  res1 <- power.grad.PCA(t=t0,A,X,Niter = 100,lambda = 2)
  
  })
  Rprof(NULL)
  end_time <- Sys.time()
  print(end_time- start_time)
  
  res2 <- numerical.grad.PCA(t=t0,X,h=0.0000000001,pos=2)
  res1$delta1
  res2$delta2
  res2$grad
  res1$graddelta.t[2]
  res1$grad
}


organize_list <- function(obj,p){
  if(dim(obj)[1]!=p) stop("not all model visited")
  list.cPLS <- vector(mode = "list", length = p)
  for (i in 1:p){
    list.cPLS[[i]] <- as.vector(which(obj[i,-c(1:3)]==1))
  }
  return(list.cPLS)
}



numerical.grad.PCA <- function(t,X,h,pos){
  matT <- diag(t)
  Xt <- X%*%matT
  Mt <- t(Xt)%*%Xt
  n <- dim(X)[1]
  At <- Mt/n  
  res <- eigen(At)
  maTh <- matT
  maTh[pos,pos] <- matT[pos,pos]+h
  Xt <- X%*%maTh
  Mt <- t(Xt)%*%Xt
  At <- Mt/n 
  resph <- eigen(At)
  maTh[pos,pos] <- matT[pos,pos]-h
  Xt <- X%*%maTh
  Mt <- t(Xt)%*%Xt
  At <- Mt/n 
  resmh <- eigen(At)
  grad <- (resph$values[1]-resmh$values[1])/(2*h)
  result<- list(delta2=res$values[1],grad=grad)#/sqrt(sum(t**2)))
  return(result)
}



PCA.Lambda.Grid.adapt <- function(Nlammax,Kmax,X,tau,Niter,alpha,psy,epoch,tol,t0,collect=1,Kchoice=NULL){
  ### Set of budget: Nlammax
  ### SIZE max: Kmax  
  ### compute lambda.max
  n <- dim(X)[1]
  lambda.max <- (1/n)*svd((t(X)%*%X)/n,nu=1,nv=1)$d[1]
  #lambda.max <- svd((t(X)%*%X)/n,nu=1,nv=1)$d[1]
  if(is.null(Kchoice)) Kchoice <- dim(X)[2]
  
  ### Run model at lambda.max
  model.Adam <- ADAM.PCA(X=X,lambda=lambda.max,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,t0=t0,collect=collect)
  result.Adam <- model.Adam$result.t
  ### check if lambda max is sufficient for selecting the zero model
  while (sum(model.Adam$s)>0 ){
    lambda.max <- 2*lambda.max
    model.Adam <- ADAM.PCA(X=X,lambda=lambda.max,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,t0=t0,collect=collect)
    result.Adam <- rbind(result.Adam,model.Adam$result.t)
  } 
  
  vec.modelsize <- sum(model.Adam$s)
  vec.lam <- lambda.max  ### define vector of lambda visited 
  lam.current <- lambda.max/2 ### first iteration 
  ### we looking for the lambda with selected a model size greater than the target size Kmax
  vec.lam <- c(vec.lam,lam.current)

  model.Adam <- ADAM.PCA(X=X,lambda=lam.current,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,t0=t0,collect=collect)
  result.Adam <- rbind(result.Adam,model.Adam$result.t)
  vec.modelsize <- c(vec.modelsize,sum(model.Adam$s)) 
 # print(c(Kmax,vec.modelsize))
  while (sum(model.Adam$s)<Kmax){
    lam.current <- lam.current/2
    model.Adam <- ADAM.PCA(X=X,lambda=lam.current,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,t0=t0,collect=collect)
    result.Adam <- rbind(result.Adam,model.Adam$result.t)
    vec.lam <- c(vec.lam,lam.current)
    vec.modelsize <- c(vec.modelsize,sum(model.Adam$s)) 
   # print(sum(model.Adam$s))
    lam.min <- lam.current
  }
  
  Nlam <- length(vec.lam)
  jj<- 0
  while (Nlam < Nlammax ){
    if(length(unique(vec.modelsize))>=(Kmax+1)) break
    jj <- jj +1
    k <- 0
    for (j in 1:(length(vec.lam)-1)){
      
      if((vec.modelsize[j+k+1]-vec.modelsize[j+k])>1){
        new.lam <- (vec.lam[j+k]+vec.lam[j+k+1])/2
        vec.lam <- insert(vec.lam,j+k+1,new.lam)
        lam.current <-new.lam  
        model.Adam <- ADAM.PCA(X=X,lambda=lam.current,tau=tau,Niter=Niter,alpha=alpha,psy=psy,epoch=epoch,tol=tol,t0=t0,collect=collect)
        result.Adam <- rbind(result.Adam,model.Adam$result.t)
        vec.modelsize <- insert(vec.modelsize,j+k+1,sum(model.Adam$s))
     #   print(c(sum(model.Adam$s),vec.modelsize))
        k <- k+1
      }
      Nlam <- length(vec.lam)
      if(Nlam>=Nlammax) break
      if(length(unique(vec.modelsize))>(Kmax+1)) break
    }
  }
  best.subset <- best.subset.BSS.PCA.k(result.Adam,X,K=Kchoice)
  #best.subset <- best.subset.combPLS.k.PCA(result.Adam,X)$bestsubset.comb
  #res <- list(bestsubset.comb=bestsubset.comb,list.BSS=empty_list2)
  result <- list(best.subset=best.subset$bestsubset.comb,list.BSS=best.subset$list.BSS)
  return(result)
  
}

ADAM.PCA <- function(X,lambda,tau=0.5,Niter=200,alpha=0.001,psy=c(0.9,0.999),epoch=10,tol=0.0001,trunc=0.001,t0=NULL,collect=1){
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  Anew <- (t(X)%*%X)/n
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
  pnew <- p
 
  res.time <- NULL
  while (i<Niter&counter< epoch){
    i <- i+1
    t.old <- t.new
    res.grad <- power.grad.PCA(t=c(t.old),A=Anew,X=X,lambda  = lambda )
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
  
  if(i<Niter) convergence <- TRUE else convergence <- FALSE
  t <- map.w.to.t(w)
  s <- as.vector((t>tau))
  
  
  result <- list(index.cov=which(s==TRUE),s=s,t=t,w=w,Niter=i,convergence=convergence,result.t=tab.t,lam=lambda)
  class(result) <- "PCA_BSS"
  return(result)
}


ADAM.PCA.truc <- function(X,delta=n,lambda,tau=0.5,Niter=200,alpha=0.001,psy=c(0.9,0.999),epoch=10,tol=0.0001,trunc=0.001,t0=NULL){
  n <- dim(X)[1]
  p <- dim(X)[2]
  indice <- 1:p
  Xnew <- X[,indice]
  
  # Mnew <- t(Xnew)%*%Y
  Anew <- (t(Xnew)%*%Xnew)/n
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
  res.time <- NULL
  while (i<Niter&counter< epoch){
    
    if(pnew<p.new.iter){
      Xnew <- X[,indice]
      pnew <- length(indice)
      
      #Mnew <- t(Xnew)%*%Y
      #Anew <- Mnew%*%t(Mnew)
      Anew <- (t(Xnew)%*%Xnew)/n
    }
    wnew <- w[indice]
    p.new.iter <- pnew
    i <- i+1
    t.old <- t.new[indice] 
    start_time <- Sys.time()
    res.grad <- power.grad.PCA(t=c(t.old),A=Anew,X=Xnew,lambda  = lambda )
    end_time <- Sys.time()
    res.time <- c(res.time,end_time-start_time)
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
  
  print(sum(res.time))
  result <- list(index.cov=which(s==TRUE),s=s,t=t,w=w,Niter=i,beta.hat.s=beta.hat.s,convergence=convergence,result.t=tab.t,lam=lambda)
  class(result) <- "PCA_BSS"
  return(result)
}


power.grad.PCA <- function(t,A,X,tol=0.001,Niter=20,seed=10,lambda){
  n <- dim(X)[1]
  p <- dim(X)[2]
  # t is the diag of the T matrix
  
  #matT <- diag(t)
  Xt <- t(c(t)%d*%t(X))  # X times T = XT
  
  At <- (t(Xt)%*%Xt)/n
  set.seed(seed)
  u <- rnorm(dim(X)[2])
  u <- u/sqrt(sum(u**2))
  norm.sol <- 1
  Niter <- Niter
  k <-0
  
  
  Bt <- t(t%d*%A) # better than Bt <- A%*%diag(t)
  
  Gt <- matrix(0,ncol=p,nrow=p)
  while((k<Niter) & (norm.sol> tol)) {
    k<- k+1
    temp.At.u <- c(At%*%matrix(u,ncol=1))
    temp <- temp.At.u/sqrt(sum((temp.At.u)**2))
# #    Zt <- t(c(u)%d*%t(Bt))
#     #Zt <- Bt%*%diag(c(u))
#     temp.norm.Atu <- sqrt(sum((At%*%matrix(u,ncol=1))**2))
#     Ft <- (Zt+diag(c(A%*%matrix(t*u,ncol=1))))/temp.norm.Atu
#     dt <- 0.5*At%*%matrix(u,ncol=1)/(temp.norm.Atu^3)
#     gt <- 2*(t(Gt)%*%(At%*%(At%*%matrix(u,ncol=1))))+2*(t(Zt+diag(c(A%*%matrix(t*u,ncol=1)))))%*%(At%*%matrix(u,ncol=1))
#     Ht <- t(outer(as.vector(gt),as.vector(dt)))
#     Gt<- (At%*%Gt)/temp.norm.Atu+Ft-Ht        
    norm.sol <- sqrt(sum((u-temp)**2))
    delta1tes <- sum(u*(At%*%matrix(u,ncol=1)))
    u <- temp
  }
 # Zt <- t(c(u)%d*%t(Bt))
 #  #Zt <- t(Bt)%*%diag(c(u))
 #  temp.norm.Atu <- sqrt(sum((At%*%matrix(u,ncol=1))**2))
 #  Ft <- (Zt+diag(c(Bt%*%u)))/temp.norm.Atu
 #  dt <- 0.5*At%*%matrix(u,ncol=1)/(temp.norm.Atu^3)
 #  gt <- 2*(t(Gt)%*%(At%*%(At%*%matrix(u,ncol=1))))+2*(t(Zt+diag(c(Bt%*%u))))%*%(At%*%matrix(u,ncol=1))
 #  Ht <- t(outer(as.vector(gt),as.vector(dt)))
 #  Gt<- (At%*%Gt)/temp.norm.Atu+Ft-Ht 
 #  
  
  r1 <- r1.tilde.t(u,t,A)
  #r2 <- c(2*t(Gt)%*%(At%*%matrix(u,ncol=1)))
  w <- map.t.to.w(t) 
  #graddelta.2 <- -((1/n^2)*(r1+r2)-lambda)*2*w*exp(-w**2)
  graddelta.2 <- -((1/n)*(r1)-lambda)*2*w*exp(-w**2)
  #graddelta.2 <- -((r1+r2)-lambda)*2*w*exp(-w**2)
  graddelta.t <- (r1)
  res <- list(grad=graddelta.2,utilde=u,delta1=delta1tes,t=t,iter=k,graddelta.t=graddelta.t)
  return(res)
}

PCA.one.component <- function(subset,data,data.d=NULL,center=TRUE,scale=TRUE){
  #X <- scale(data[,subset], center = center, scale = scale)
  if(is.null(data.d)) data.d <- data
  X <- as.matrix(data.d[,subset])
  res <- svd(X,nu=1,nv=1)
  loading <- rep(0,dim(data)[2])
  loading[subset] <- res$v[,1]
  #variates <- drop(X%*%res$v[,1])
  variates <- drop(data%*%matrix(loading,ncol=1))
  result <- list(loading=loading,variates=variates)
  nor2x = sum((data)^2)
  result$var.tot=sum(data^2 / max(1, nrow(X) - 1))
  temp = variates / sum(variates**2)
  result$exp_varX = drop((t(matrix(variates,ncol=1)) %*% data %*% t(data) %*% matrix(temp,ncol=1))/nor2x)
  result$sdev = res$d[1] / sqrt(max(1, nrow(X) - 1))
  result$explained_variance = result$sdev^2 / result$var.tot
  return(result)
}  



plot.CPEV <- function(Object,data,drop=0.1){
  
  p <- dim(data)[2]
  data.CPEV <- data.frame(CPEV=rev(Object[,4]),Sparsity=p-rev(Object[,1]))
  threshold <- max(Object[,4])-max(Object[,4])*drop
  sparsity <- dim(data)[2]-which.min(Object[,4]<threshold)+1
  
  fig <- ggplot(data.CPEV, aes(x=Sparsity, y=CPEV)) +
    geom_linerange(aes(x=Sparsity, ymax=CPEV, ymin=0))+
    geom_hline(aes(yintercept=threshold,linetype="Threshold "),col="red") +
    geom_segment(data = data.frame(x = dim(data)[2]-which.min(Object[,4]<threshold), 
                                   xend = dim(data)[2]-which.min(Object[,4]<threshold), y =0 , yend = Object[which.min(Object[,4]<threshold),4]), mapping = aes(x = x, xend = xend, y = y, 
                                                                                                                                                                yend = yend),col="blue",lwd=2) +
    labs(linetype=paste("Drop of CPEV =",drop*100,"%"))+
    theme(legend.position = c(.5, .5),legend.justification = c("right", "top"), legend.box.just = "right",legend.margin = margin(6, 6, 6, 6))
  
  res <- list(sparsity=sparsity,model.size=p-sparsity,select.var=Object[p-sparsity,-c(1:5)],plot=fig)
   
}



#result.cPCA$best.subset-> object


stat.CPEV <- function(object,deflat.mat=NULL,data,variate.F.PCA,mat.v=NULL,mat.u=NULL,ncomp){  
  res <- NULL 
  res.cor <- NULL
  res.CPEV <- NULL
  for (i in 1:length(object$list.BSS)){
    subset <- object$list.BSS[[i]]
    temp <- PCA.one.component(subset,data,deflat.mat)
    mat.v.new <- cbind(mat.v,temp$loading)
    mat.u.new <- cbind(mat.u,temp$variates)
    temp.cpev <- CPEV(data,mat.v.new,mat.u.new,ncomp=ncomp)
    res.CPEV <- c(res.CPEV,temp.cpev$CPEV[ncomp])
    res <- c(res,temp.cpev$PEV.var[ncomp])
    res.cor <- c(res.cor,cor(temp$variates,variate.F.PCA[,ncomp]))
  }
  
  result <- cbind(object$best.subset[,1:2],signif(100*res,4),signif(100*res.CPEV,4),res.cor,object$best.subset[,-c(1:2)])
  colnames(result)[1:5] <- c("Dimension","Objective fct","PEV(sparse)","CPEV","cor(sPC,PC)")
  return(result)  
  
}  




CPEV <- function(X,mat.v,mat.u,ncomp){
  vect.varX=vector(length=ncomp)
  prop_expl_var=vector(length=ncomp)
  var.tot <- sum(X^2)
  cum.var <- vector(length=ncomp)
for (h in 1:ncomp){  
#h <- ncomp
X.var = X %*% mat.v[,1:h]%*%solve(t(mat.v[,1:h])%*%mat.v[,1:h])%*%t(mat.v[,1:h])
vect.varX[h] = sum(X.var^2)
cum.var[h] <- vect.varX[h]/var.tot
#prop_expl_var <- explained_variance(X, mat.u, ncomp) I do not like this definition of variance explained
prop_expl_var[h] <- sum(mat.u[,h]**2)/var.tot
}

result <- list(CPEV=cum.var,PEV.var=prop_expl_var,var.tot=var.tot)
return(result)
}  
  
  
explained_variance <- function(data, variates, ncomp)
{
  ## pre-allocate output
  expl_var <- vector(mode = 'numeric', length = ncomp)
  names(expl_var) <- paste0('comp', seq_len(ncomp))
  
  norm2.X <- norm(data, type='F')^2 # total variance in the data
  
  for (h in 1:ncomp)
  {
    a <- crossprod(variates[, h, drop=FALSE], data)
    # this is equivalent to calculate the redundancy as detailed in the help file
    expl_var[h] <- tcrossprod(a) / c(crossprod(variates[, h])) / norm2.X
  }
  
  return(expl_var)
}



  
Deflation.step <- function(subset,X,data){
loading <- rep(0,dim(X)[2])
X.select <- X[,subset]  
model <- svd(X.select,nv=1,nu=1)
xih <- X.select%*%model$v[,1] ## definition
ch <- t(X)%*%xih/(sum(xih**2))
X1 <- X - matrix(xih,ncol=1)%*%t(ch)
loading[subset] <- model$v[,1]
variates <- drop(data%*%matrix(loading,ncol=1))
res <- list(Xd=X1,variates=variates,ch=ch,loading=loading)
return(res)
}



#result.cPCA$best.subset-> object







stat.PCA.Comp1 <- function(object,data,var.exp,variate.F.PCA){  
res <- NULL 
res.cor <- NULL
for (i in 1:nrow(object)){
  subset <- as.logical(object[i,-c(1:3)]) 
  temp <- PCA.one.component(subset,data)
  res <- c(res,temp$explained_variance)
  res.cor <- c(res.cor,cor(temp$variates,variate.F.PCA))
}

result <- cbind(object[,1:2],signif(100*res,4),signif(100*res/var.exp,4),res.cor,object[,-c(1:3)])
colnames(result)[1:5] <- c("Dimension","Objective fct","PEV(sparse)","PEV(sparse)/PEV","cor(sPC,PC)")
return(result)  
  
}  


best.subset.BSS.PCA.k <- function(result.t,X,K){
  
  p <- dim(X)[2]
  
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
    res <- (apply(matrix(empty_list[[i]],ncol=i),MARGIN=1,FUN=cov.train.BSS.PCA,Xtrain=X))
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

cov.train.BSS.PCA <- function(x,Xtrain){
  #print(x)
  Xselect <- Xtrain[,x]
  res.svd <- svd((t(Xselect)%*%Xselect)/dim(Xtrain)[1],nu=1,nv=1)
  covtrain <- -(res.svd$d[1])/dim(Xtrain)[1]
  return(covtrain)
}





best.subset.combPLS.k.PCA <- function(result.t,X){
  
  p <- dim(X)[2]
  #############
  ## Add model with t result 
  ############
  res.order <- t(apply(abs(result.t),MARGIN=1,FUN=function(x) order(x,decreasing = TRUE)))
  res.order.unique <- unique(res.order)
  res.list <- apply(res.order.unique,MARGIN=1,FUN=best.set.combss.PCA,simplify = FALSE)
  res.list <- list.rbind(res.list)
  res.list.unique.t <- unique(res.list)
  
  final.select <- unique(res.list.unique.t)
  dim(final.select) -> nb.visited.model
  RSS <- t(apply(final.select,MARGIN=1,FUN=cov.train.PCA,Xtrain=X))
  RSS.indice <- cbind(RSS,1:dim(RSS)[1])
  RSS.indice.order <- RSS.indice[order(RSS.indice[,1]),]
  
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
  
  bestsubset.comb <-cbind(best.subset.res[,1:2],final.select[best.subset.res[,3],])
  colnames(bestsubset.comb) <- c("Dimension","Variance Train",paste("X",1:p,sep=""))
  res <- list(bestsubset.comb=bestsubset.comb,nnb.visited.model=nb.visited.model)
  return(res)
  
}


cov.train.PCA <- function(x,Xtrain){
  Xselect <- Xtrain[,x]
  covtest <- NA
  res.svd <- svd((t(Xselect)%*%Xselect)/dim(Xtrain)[1],nu=1,nv=1)
  covtrain <- -(res.svd$d[1])/dim(Xtrain)[1]
  return(c(sum(x),covtrain))
}

best.set.combss.PCA <- function(x){
  res <- matrix(FALSE,ncol=length(x),nrow=length(x))
  for(i in 1:length(x)){
    res[(i:length(x)),x[i]] <- TRUE 
  }
  return(res)
}
