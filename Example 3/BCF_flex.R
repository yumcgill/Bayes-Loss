library(bcf)
library(parallel)
options(mc.cores = 23) 
BCF_flex<-function(seed,N){
  set.seed(seed)
  #mu_x<-0
  sigma_x<-1
  x1<-rnorm(n=N,mean=1,sd=sigma_x)
  x2<-rnorm(n=N,mean=1,sd=sigma_x)
  x3<-rnorm(n=N,mean=-1,sd=sigma_x)
  x4<-rnorm(n=N,mean=-1,sd=sigma_x)
  
  
  alp0<-0.3
  alp1<-0.9
  alp2<--1.25
  alp3<-1.5
  px<-1/(1+exp(-alp0*x1-alp1*x2-alp2*x3-alp3*x4))
  D<-rbinom(n=N,size=1,prob=px)
  # Now generate the continuous outcome
  
  linp<-D+2*D*x1 + x1+x2+x3+x4 + 0.25*x1^2 + 0.75*x2*x4+0.75*x3*x4
  sigma_y<-1
  Y<-rnorm(N,mean=linp,sd=sigma_y)
  
  
  
  
  ps <- glm(D~x1+x2+x3+x4,family = binomial(link = "logit"))$fitted.values
  
  fit1<-bcf(Y,D,cbind(1,x1,x2,x3,x4),cbind(1,x1,x2,x3,x4),ps,nburn=500, nsim=2000)
  post_theta1<-rowMeans(fit1$tau)
  postmean_theta1<-mean(post_theta1)
  ci<-0 
  if(quantile(post_theta1,prob=c(0.025))<3 && quantile(post_theta1,prob=c(0.975))>3) {ci<-1}
  return(list(theta=postmean_theta1,ci=ci))
  
}


BCF_flex_200<-mclapply(sample(c(1:100000000),1000),function(x) BCF_flex(x,200))
BCF_flex_200_ate<-unlist(mclapply(BCF_flex_200, '[[', "theta"))
coverge_BCF_flex_200<-sum(unlist(mclapply(BCF_flex_200, '[[', "ci")))/1000
coverge_BCF_flex_200
mean(BCF_flex_200_ate);var(BCF_flex_200_ate)
mean(BCF_flex_200_ate) -3
sqrt(mean((BCF_flex_200_ate - 3)^2))

BCF_flex_1000<-mclapply(sample(c(1:100000000),1000),function(x) BCF_flex(x,1000))
BCF_flex_1000_ate<-unlist(mclapply(BCF_flex_1000, '[[', "theta"))
coverge_BCF_flex_1000<-sum(unlist(mclapply(BCF_flex_1000, '[[', "ci")))/1000
coverge_BCF_flex_1000
mean(BCF_flex_1000_ate);var(BCF_flex_1000_ate)
mean(BCF_flex_1000_ate)-3
sqrt(mean((BCF_flex_1000_ate - 3)^2))


BCF_flex_2000<-mclapply(sample(c(1:100000000),1000),function(x) BCF_flex(x,2000))
BCF_flex_2000_ate<-unlist(mclapply(BCF_flex_2000, '[[', "theta"))
coverge_BCF_flex_2000<-sum(unlist(mclapply(BCF_flex_2000, '[[', "ci")))/1000
coverge_BCF_flex_2000
mean(BCF_flex_2000_ate);var(BCF_flex_2000_ate)
mean(BCF_flex_2000_ate)-3
sqrt(mean((BCF_flex_2000_ate - 3)^2))



