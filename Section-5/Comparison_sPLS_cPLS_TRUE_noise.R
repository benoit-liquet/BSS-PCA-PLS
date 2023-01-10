#############
## R code for reproducing simulation section 5.2
#############



#############
## Brut Force
#############

library(Rfast)
library(mvtnorm)

###########
## Example Simulation data
##########
source("Function_BSS_PLS2.R")
source("generated_true_model_PLS.R")

sigma.gamma <- 1
sigma.e.t <- c(3,4,5,6,7)
p <- 15
q <- 10
sparsity <- 5

Nsimu <- 100

tab.cpls <- NULL
tab.spls <- NULL
for (kk in sigma.e.t){
  spls.true <- 0
  cpls.true <- 0
  sigma.e <- kk
for (jj in 1:Nsimu){
print(jj)
set.seed(jj)
data.gen <- simu.true.model(ntrain = 50,ntest = 500,sigma.gamma,sigma.e,p,q,sparsity)
#######################################
X <- data.gen$X.train
Y <- data.gen$Y.train
#Xtest <- data.gen$X.test
#Ytest <- data.gen$Y.test
n <- dim(X)[1]


############################################################
#### Compute the true subset by using an exhausitive search
##############################################################
True_list <- vector(mode = "list", length = p)
for(k in 1:p){
  res<-combn(1:p, k)
  sol <- apply(res,MARGIN=2,FUN=criterion,m=n,data=X,Y=Y)
  min.set <- which.min(sol)
  True_list[[k]] <- res[,min.set]
  
}


############
## Result with sparse PLS
############
spls_list <- vector(mode = "list", length = p)
for(k in 1:p){
  model.spls <- sgPLS::sPLS(X,Y,ncomp=1,keepX=c(k))
  res.spls <- sgPLS::select.spls(model.spls)$select.X[[1]]
  spls_list[[k]] <- as.vector(res.spls)
}


#######################3
## Result with cPLS
######################
if(T){
result.cPLS <- PLS.Lambda.Grid(Nlammax=30,Kmax=15,X,Y,tau=0.5,Niter=1000,alpha=0.01,psy=c(0.9,0.999),epoch=10,tol=0.0001,t0=NULL,collect = 1,Kchoice = 15) 
cpls_list <- result.cPLS$list.BSS
}
########
## Comparison 
#########

spls.true <- spls.true+ compare(spls_list,True_list)
cpls.true <- cpls.true + compare(cpls_list,True_list)
}

tab.cpls <- cbind(tab.cpls,cpls.true)  
tab.spls <- cbind(tab.spls,spls.true)
}  

colnames(tab.spls) <- paste("noise = ",sigma.e.t)
colnames(tab.cpls) <- paste("noise = ",sigma.e.t)
result <- list(spls=tab.spls,cpls=tab.cpls)

save(result,file="result_best_ss_noise.Rdata")

