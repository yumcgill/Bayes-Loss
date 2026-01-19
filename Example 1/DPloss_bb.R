library(parallel)
options(mc.cores = 23) 
library(gtools)

SI_bb<-function(N){
  #set.seed(seed)
  # First generate covariates X, D
  #N=10
  mu_x<-0
  sigma_x<-1
  x1<-rnorm(n=N,mean=mu_x,sd=sigma_x)
  x2<-rnorm(n=N,mean=mu_x,sd=sigma_x)
  x3<-rnorm(n=N,mean=mu_x,sd=sigma_x)
  x4<-rnorm(n=N,mean=mu_x,sd=sigma_x)
  u1<-abs(x1)/sqrt(1-2/pi)
  
  alp0<-0.4
  alp1<-0.4
  alp2<-0.8
  px<-1/(1+exp(-alp0*u1-alp1*x2-alp2*x3))
  D<-rbinom(n=N,size=1,prob=px)
  # Now generate the continuous outcome
  beta0<-1
  beta1<--1
  beta2<--1
  beta3<--1
  linp<-beta0*D+beta1*u1+beta2*x2+beta3*x4
  sigma_y<-1
  Y<-rnorm(N,mean=linp,sd=sigma_y)
  
  PS_pred <- model.matrix( ~ u1+x2+x3)
  
  
  post_theta1<-NULL
  for (i in 1:10000){
    w<-as.vector(rdirichlet(1,rep(1,N)))
    ps <- glm.fit(PS_pred,D,weights=w,family = binomial(link = "logit"))$fitted.values
    mod<-lm(Y~ D+x1+x2+x4+ps,weights=w)
    post_theta1<-c(post_theta1,mod$coefficients[2])
  }
  post_theta1 = post_theta1[!is.na(post_theta1)]

  postmean_theta1<-mean(post_theta1)
  ci<-0
  if(quantile(post_theta1,prob=c(0.025))<1 && quantile(post_theta1,prob=c(0.975))>1) {ci<-1}
  return(list(theta=postmean_theta1,ci=ci))
}

SenarioI_res_drbb_20<-mclapply(1:1000, function(x) SI_bb(20))
causal_para_drbb_20<-unlist(mclapply(SenarioI_res_drbb_20, '[[', "theta"))
coverge_drbb_20<-sum(unlist(mclapply(SenarioI_res_drbb_20, '[[', "ci")))/1000



SenarioI_res_drbb_50<-mclapply(1:1000, function(x) SI_bb(50))
causal_para_drbb_50<-unlist(mclapply(SenarioI_res_drbb_50, '[[', "theta"))
coverge_drbb_50<-sum(unlist(mclapply(SenarioI_res_drbb_50, '[[', "ci")))/1000


SenarioI_res_drbb_100<-mclapply(1:1000, function(x) SI_bb(100))
causal_para_drbb_100<-unlist(mclapply(SenarioI_res_drbb_100, '[[', "theta"))
coverge_drbb_100<-sum(unlist(mclapply(SenarioI_res_drbb_100, '[[', "ci")))/1000





SenarioI_res_drbb_500<-mclapply(1:1000, function(x) SI_bb(500))
causal_para_drbb_500<-unlist(mclapply(SenarioI_res_drbb_500, '[[', "theta"))
coverge_drbb_500<-sum(unlist(mclapply(SenarioI_res_drbb_500, '[[', "ci")))/1000



### Senario II

SII_bb<-function(N){
  #set.seed(seed)
  # First generate covariates X, D
  #N=10
  mu_x<-0
  sigma_x<-1
  x1<-rnorm(n=N,mean=mu_x,sd=sigma_x)
  x2<-rnorm(n=N,mean=mu_x,sd=sigma_x)
  x3<-rnorm(n=N,mean=mu_x,sd=sigma_x)
  x4<-rnorm(n=N,mean=mu_x,sd=sigma_x)
  u1<-abs(x1)/sqrt(1-2/pi)
  
  alp0<-0.4
  alp1<-0.4
  alp2<-0.8
  px<-1/(1+exp(-alp0*u1-alp1*x2-alp2*x3))
  D<-rbinom(n=N,size=1,prob=px)
  # Now generate the continuous outcome
  beta0<-1
  beta1<--1
  beta2<--1
  beta3<--1
  linp<-beta0*D+beta1*u1+beta2*x2+beta3*x4
  sigma_y<-1
  Y<-rnorm(N,mean=linp,sd=sigma_y)
  
  PS_pred <- model.matrix( ~ x1+x2+x3)
  
  
  
  post_theta1<-NULL
  for (i in 1:1000){
    w<-as.vector(rdirichlet(1,rep(1,N)))
    ps <- glm.fit(PS_pred,D,weights=w,family = binomial(link = "logit"))$fitted.values
    mod<-lm(Y~ D+u1+x2+x4+ps,weights=w)
    post_theta1<-c(post_theta1,mod$coefficients[2])
  }
  post_theta1 = post_theta1[!is.na(post_theta1)]
  ####
  postmean_theta1<-mean(post_theta1)
  ci<-0
  if(quantile(post_theta1,prob=c(0.025))<1 && quantile(post_theta1,prob=c(0.975))>1) {ci<-1}
  return(list(theta=postmean_theta1,ci=ci))
}

SenarioII_res_drbb_20<-mclapply(1:1000, function(x) SII_bb(20))
causal_paraII_drbb_20<-unlist(mclapply(SenarioII_res_drbb_20, '[[', "theta"))
coverge_drbbII_20<-sum(unlist(mclapply(SenarioII_res_drbb_20, '[[', "ci")))/1000


SenarioII_res_drbb_50<-mclapply(1:1000, function(x) SII_bb(50))
causal_paraII_drbb_50<-unlist(mclapply(SenarioII_res_drbb_50, '[[', "theta"))
coverge_drbbII_50<-sum(unlist(mclapply(SenarioII_res_drbb_50, '[[', "ci")))/1000


SenarioII_res_drbb_100<-mclapply(1:1000, function(x) SII_bb(100))
causal_paraII_drbb_100<-unlist(mclapply(SenarioII_res_drbb_100, '[[', "theta"))
coverge_drbbII_100<-sum(unlist(mclapply(SenarioII_res_drbb_100, '[[', "ci")))/1000


SenarioII_res_drbb_500<-mclapply(1:1000, function(x) SII_bb(500))
causal_paraII_drbb_500<-unlist(mclapply(SenarioII_res_drbb_500, '[[', "theta"))
coverge_drbbII_500<-sum(unlist(mclapply(SenarioII_res_drbb_500, '[[', "ci")))/1000



#### Senario III

SIII_bb<-function(N){
  #set.seed(seed)
  # First generate covariates X, D
  #N=10
  mu_x<-0
  sigma_x<-1
  x1<-rnorm(n=N,mean=mu_x,sd=sigma_x)
  x2<-rnorm(n=N,mean=mu_x,sd=sigma_x)
  x3<-rnorm(n=N,mean=mu_x,sd=sigma_x)
  x4<-rnorm(n=N,mean=mu_x,sd=sigma_x)
  u1<-abs(x1)/sqrt(1-2/pi)
  
  alp0<-0.4
  alp1<-0.4
  alp2<-0.8
  px<-1/(1+exp(-alp0*u1-alp1*x2-alp2*x3))
  D<-rbinom(n=N,size=1,prob=px)
  # Now generate the continuous outcome
  beta0<-1
  beta1<--1
  beta2<--1
  beta3<--1
  linp<-beta0*D+beta1*u1+beta2*x2+beta3*x4
  sigma_y<-1
  Y<-rnorm(N,mean=linp,sd=sigma_y)
  
  PS_pred <- model.matrix( ~ x1+x2+x3)
  
  #OR_pred <- model.matrix( ~ D+x1+x2+x4+ps)
  
  post_theta1<-NULL
  for (i in 1:10000){
    w<-as.vector(rdirichlet(1,rep(1,N)))
    ps <- glm.fit(PS_pred,D,weights=w,family = binomial(link = "logit"))$fitted.values
    mod<-lm(Y~ D+x1+x2+x4+ps,weights=w)
    post_theta1<-c(post_theta1,mod$coefficients[2])
  }
  post_theta1 = post_theta1[!is.na(post_theta1)]
  ####
  postmean_theta1<-mean(post_theta1)
  ci<-0
  if(quantile(post_theta1,prob=c(0.025))<1 && quantile(post_theta1,prob=c(0.975))>1) {ci<-1}
  return(list(theta=postmean_theta1,ci=ci))
}

SenarioIII_res_drbb_20<-mclapply(1:1000, function(x) SIII_bb(20))
causal_paraIII_drbb_20<-unlist(mclapply(SenarioIII_res_drbb_20, '[[', "theta"))
coverge_drbbIII_20<-sum(unlist(mclapply(SenarioIII_res_drbb_20, '[[', "ci")))/1000


SenarioIII_res_drbb_50<-mclapply(1:1000, function(x) SIII_bb(50))
causal_paraIII_drbb_50<-unlist(mclapply(SenarioIII_res_drbb_50, '[[', "theta"))
coverge_drbbIII_50<-sum(unlist(mclapply(SenarioIII_res_drbb_50, '[[', "ci")))/1000

SenarioIII_res_drbb_100<-mclapply(1:1000, function(x) SIII_bb(100))
causal_paraIII_drbb_100<-unlist(mclapply(SenarioIII_res_drbb_100, '[[', "theta"))
coverge_drbbIII_100<-sum(unlist(mclapply(SenarioIII_res_drbb_100, '[[', "ci")))/1000



SenarioIII_res_drbb_500<-mclapply(1:1000, function(x) SIII_bb(500))
causal_paraIII_drbb_500<-unlist(mclapply(SenarioIII_res_drbb_500, '[[', "theta"))
coverge_drbbIII_500<-sum(unlist(mclapply(SenarioIII_res_drbb_500, '[[', "ci")))/1000




