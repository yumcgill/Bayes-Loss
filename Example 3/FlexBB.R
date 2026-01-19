library(gtools)
library(parallel)
options(mc.cores = 23)
BB_flex<-function(seed,N){
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

  
  
  post_theta1<-NULL
  for (i in 1:1000){
  w <- as.vector(rdirichlet(1, rep(1,N)))
  ps <- glm(D~x1+x2+x3+x4,family = binomial(link = "logit"),weights=w)$fitted.values
  dataset<-data.frame(cbind(Y,D,x1,x2,x4,ps))
  data.alltrt <- dataset
  data.alltrt$D <- 1
  data.nontrt <- dataset
  data.nontrt$D <- 0
  mod1.lmX <- lm(Y ~ D+I(D*x1)+ps+I(ps*x1),weights=w)
  
  APO.lmX.1 <- mean(predict(mod1.lmX,data.alltrt))
  APO.lmX.0 <- mean(predict(mod1.lmX,data.nontrt))
  post_theta1<-c(post_theta1,APO.lmX.1 - APO.lmX.0)
  }
  
  postmean_theta1<-mean(post_theta1)
  #var_theta1<-var(post_theta1)
  ci<-0 
  if(quantile(post_theta1,prob=c(0.025))<3 && quantile(post_theta1,prob=c(0.975))>3) {ci<-1}
  return(list(theta=postmean_theta1,ci=ci))
  
}


BB_flex_200<-mclapply(sample(c(1:100000000),1000),function(x) BB_flex(x,200))
BB_flex_200_ate<-unlist(mclapply(BB_flex_200, '[[', "theta"))
coverge_BB_flex_200<-sum(unlist(mclapply(BB_flex_200, '[[', "ci")))/1000
coverge_BB_flex_200
mean(BB_flex_200_ate);var(BB_flex_200_ate)
mean(BB_flex_200_ate) -3
sqrt(mean((BB_flex_200_ate - 3)^2))


BB_flex_1000<-mclapply(sample(c(1:100000000),1000),function(x) BB_flex(x,1000))
BB_flex_1000_ate<-unlist(mclapply(BB_flex_1000, '[[', "theta"))
coverge_BB_flex_1000<-sum(unlist(mclapply(BB_flex_1000, '[[', "ci")))/1000
coverge_BB_flex_1000
mean(BB_flex_1000_ate);var(BB_flex_1000_ate)
mean(BB_flex_1000_ate)-3
sqrt(mean((BB_flex_1000_ate - 3)^2))


BB_flex_2000<-mclapply(sample(c(1:100000000),1000),function(x) BB_flex(x,2000))
BB_flex_2000_ate<-unlist(mclapply(BB_flex_2000, '[[', "theta"))
coverge_BB_flex_2000<-sum(unlist(mclapply(BB_flex_2000, '[[', "ci")))/1000
coverge_BB_flex_2000
mean(BB_flex_2000_ate);var(BB_flex_2000_ate)
mean(BB_flex_2000_ate)-3
sqrt(mean((BB_flex_2000_ate - 3)^2))

