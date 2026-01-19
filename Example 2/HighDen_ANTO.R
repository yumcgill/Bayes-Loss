library(mvtnorm)
library(glmnet)
library(gtools)
library(devtools)
install_github(repo = "jantonelli111/DoublyRobustHD")
library(DoublyRobustHD)
library(parallel)
options(mc.cores = 23)


high_Anto<-function(seed,n){
  set.seed(seed)
  p = 500
  sig = array(0.3,c(p,p))
  diag(sig) = 1
  x = matrix(rmvnorm(n,mean=rep(0, nrow(sig)), sigma=sig),n,p)
  d = rbinom(n, 1, p=expit(0.45*x[,1] + 0.9*x[,2] - 0.4*x[,5] +1.3*x[,2]*x[,5]+1.8*x[,1]*x[,2]))
  y = rnorm(n, mean=d + 0.5*x[,1] + x[,3] - 0.1*x[,4]-0.2*x[,7] + 1.5*x[,3]*x[,4] +0.6*x[,7]*x[,7] +1.2*x[,1]*x[,3], sd=1)
  
  
  estLinear = DRbayes(y=y, t=d, x=x, nScans=500, nBurn=100, thin=1)
  
  postmean_theta1 = estLinear$TreatEffect
  CIlower = estLinear$TreatEffectCI[1]
  CIupper = estLinear$TreatEffectCI[2] 
  
  
  ci<-0
  
  if(CIlower<1 && CIupper>1) {ci<-1}
  return(list(theta=postmean_theta1,ci=ci))
}

high_Anto_50<-mclapply(sample(c(1:100000000),500),function(x) high_Anto(x,50))
para_high_Anto_50<-unlist(mclapply(high_Anto_50, '[[', "theta"))
coverge_high_Anto_50<-sum(unlist(mclapply(high_Anto_50, '[[', "ci")))/1000
coverge_high_Anto_50
mean(para_high_Anto_50);var(para_high_Anto_50)

high_Anto_100<-mclapply(sample(c(1:100000000),500),function(x) high_Anto(x,100))
para_high_Anto_100<-unlist(mclapply(high_Anto_100, '[[', "theta"))
coverge_high_Anto_100<-sum(unlist(mclapply(high_Anto_100, '[[', "ci")))/1000
coverge_high_Anto_100
mean(para_high_Anto_100);var(para_high_Anto_100)