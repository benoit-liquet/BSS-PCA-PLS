###########
## R code for figure 1 (see section 3.1)
##########
library(mvtnorm)
library(rgl)
source("Function_PLS1.R")

#### generate some data for PLS1
set.seed(1)
beta <- c(5,0)
n <- 100
X <- mvtnorm::rmvnorm(n,mean=c(1,1),sigma=matrix(c(3,1,1,2),ncol=2,nrow=2,byrow=T))
Y <- X%*%matrix(beta,ncol=1)+rnorm(n,sd=1)



t0 <- rep(0.5,2)
plot.landscape.path.BSS(lambda=0,t0=t0)
plot.landscape.path.BSS(lambda=150,t0=t0)
plot.landscape.path.BSS(lambda=450,t0=t0)

l




