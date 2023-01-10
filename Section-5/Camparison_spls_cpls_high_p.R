#############
## Simulation section 5.5, effect of dimension d
#############

library(Rfast)
library(mvtnorm)


source("Function_BSS_PLS2.R")
source("generated_true_model_PLS.R")
nsample <- 100
sigma.gamma <- 1
sigma.e.t <- 4#c(2)#,6)#,5,6)
p.t <- c(50,100,200,500)
q <- 10
#sparsity <- 40

Nsimu <- 50
ntruepred <- 10 


tab.cpls <- NULL
tab.spls <- NULL
tabcrit.cpls <- NULL 
tabcrit.spls <- NULL 
tabcrit.true <- NULL 
for (kk in p.t){
  
  spls.true <- 0
  cpls.true <- 0
  res.true <- NULL
  res.spls.crit <- NULL
  res.cpls.crit <- NULL
  sigma.e <- sigma.e.t[1]
  p <- kk
  sparsity <- p - ntruepred
  True_list <- vector(mode = "list", length = 1)
  True_list[[1]] <- (sparsity+1):p
  true.dim <- length(True_list[[1]])
  
  for (jj in 1:Nsimu){
    print(jj)
    set.seed(jj)
    
    data.gen <- simu.true.model(ntrain = nsample,ntest = 500,sigma.gamma,sigma.e,p,q,sparsity)
    #######################################
    X <- data.gen$X.train
    Y <- data.gen$Y.train
    #Xtest <- data.gen$X.test
    
    #Ytest <- data.gen$Y.test
    n <- dim(X)[1]
    #print(X[1,1])
    true.crit <- criterion(x=True_list[[1]],m=n,data=X,Y=Y)
    res.true <- c(res.true,true.crit)
    ############################################################
    #### Compute the true subset by using an exhausitive search
    ##############################################################
    
    
    
    ############
    ## Result with sparse PLS
    ############
    spls_list <- vector(mode = "list", length = 1)
    #for(k in 1:p){
    k <- true.dim
      model.spls <- sgPLS::sPLS(X,Y,ncomp=1,keepX=c(k))
      res.spls <- sgPLS::select.spls(model.spls)$select.X[[1]]
      spls_list[[1]] <- as.vector(res.spls)
    #
      spls.crit <- criterion(x=spls_list[[1]],m=n,data=X,Y=Y)
    
    #######################3
    ## Result with cPLS
    ######################
    if(T){
        result.cPLS <- PLS.specific.subset.size(Nlammax=60,Kmax=true.dim,X,Y,tau=0.5,Niter=1000,alpha=0.001,psy=c(0.9,0.999),epoch=10,tol=0.001,t0=NULL,collect=1) 

      cpls_list <- result.cPLS$list.BSS[[true.dim]]
      cpls.crit <- criterion(x=cpls_list,m=n,data=X,Y=Y)
      res.cpls.crit <- c(res.cpls.crit,cpls.crit)
    }
    ########
    ## Comparison 
    #########
    
    spls.true <- spls.true+ compare(spls_list,True_list)
    cpls.true <- cpls.true + compare(cpls_list,True_list)
    res.spls.crit <- c(res.spls.crit,spls.crit)
    #print(c(spls.crit,cpls.crit,cpls.crit<spls.crit))
  }
  result.int <- list(cpls.true=cpls.true,spls.true=spls.true,tabcrit.cpls=res.cpls.crit,tabcrit.spls=res.spls.crit,tabcrit.true=res.true)
  save(result.int,file=paste("result-high-p-",p,".RData",sep=""))
  
  tab.cpls <- cbind(tab.cpls,cpls.true)  
  tab.spls <- cbind(tab.spls,spls.true)
  tabcrit.cpls <- cbind(tabcrit.cpls,res.cpls.crit)  
  tabcrit.spls <- cbind(tabcrit.spls,res.spls.crit)
  tabcrit.true <- cbind(tabcrit.true,res.true)
  
}  
result <- list(tab.cpls=tab.cpls,tab.spls=tab.spls,tabcrit.cpls=tabcrit.cpls,tabcrit.spls,tabcrit.true=tabcrit.true)
save(result,file="result-high-p-total.RData")


