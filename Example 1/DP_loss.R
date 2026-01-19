library(parallel)
options(mc.cores = 23) 
library(gtools)

SI_dp<-function(seed,N){
  set.seed(seed)
  # First generate covariates X, D
  #N=20
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
  
  
  ps <- glm.fit(PS_pred,D,family = binomial(link = "logit"))$fitted.values
  
  dataset<-data.frame(cbind(Y,D,x1,x2,x4,ps))
  post_theta1<-NULL
  for (i in 1:1000){
    u<-runif(1000)
    av<-1
    datasetnew<-dataset
    pred<-predict(lm(Y ~ D+x1+x2+x4+ps))
    for(nv in 1:1000){
      if(u[nv] > av/(av+N+nv-1)){
        ind<-sample(1:(N+nv-1),size=1)
        datasetnew<-rbind(datasetnew,datasetnew[ind,])
        pred<-c(pred,pred[ind])
      }else{
        w<-as.vector(rdirichlet(1,rep(1,(N+nv-1))))
        ind<-sample(1:(N+nv-1),prob=w,size=1)
        newdata<-datasetnew[ind,]
        newdata[1]<-rnorm(1,pred[ind],1)
        datasetnew<-rbind(datasetnew,newdata)
        pred<-c(pred,pred[ind])
      }
    }
    
    
    mod<-lm(Y~ D+x1+x2+x4+ps,data=datasetnew[-c(1:N),])
    post_theta1<-c(post_theta1,mod$coefficients[2])
  }
  post_theta1 = post_theta1[!is.na(post_theta1)]
  ####
  postmean_theta1<-mean(post_theta1)
  ci<-0
  if(quantile(post_theta1,prob=c(0.025))<1 && quantile(post_theta1,prob=c(0.975))>1) {ci<-1}
  return(list(theta=postmean_theta1,ci=ci))
}

SenarioI_res_drdp_20<-mclapply(sample(c(1:100000000),1000), function(x) SI_dp(x,20))
causal_para_drdp_20<-unlist(mclapply(SenarioI_res_drdp_20, '[[', "theta"))
coverge_drdp_20<-sum(unlist(mclapply(SenarioI_res_drdp_20, '[[', "ci")))/1000


SenarioI_res_drdp_50<-mclapply(sample(c(1:100000000),1000), function(x) SI_dp(x,50))
causal_para_drdp_50<-unlist(mclapply(SenarioI_res_drdp_50, '[[', "theta"))
coverge_drdp_50<-sum(unlist(mclapply(SenarioI_res_drdp_50, '[[', "ci")))/1000

SenarioI_res_drdp_100<-mclapply(sample(c(1:100000000),1000),function(x) SI_dp(x,100))
causal_para_drdp_100<-unlist(mclapply(SenarioI_res_drdp_100, '[[', "theta"))
coverge_drdp_100<-sum(unlist(mclapply(SenarioI_res_drdp_100, '[[', "ci")))/1000


SenarioI_res_drdp_500<-mclapply(sample(c(1:100000000),1000), function(x) SI_dp(x,500))
causal_para_drdp_500<-unlist(mclapply(SenarioI_res_drdp_500, '[[', "theta"))
coverge_drdp_500<-sum(unlist(mclapply(SenarioI_res_drdp_500, '[[', "ci")))/1000



SII_dp<-function(seed,N){
  set.seed(seed)
  # First generate covariates X, D
  #N=20
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
  
  
  ps <- glm.fit(PS_pred,D,family = binomial(link = "logit"))$fitted.values
  
  dataset<-data.frame(cbind(Y,D,u1,x2,x4,ps))
  post_theta1<-NULL
  for (i in 1:1000){
    u<-runif(1000)
    av<-1
    datasetnew<-dataset
    pred<-predict(lm(Y ~ D+u1+x2+x4+ps))
    for(nv in 1:1000){
      if(u[nv] > av/(av+N+nv-1)){
        ind<-sample(1:(N+nv-1),size=1)
        datasetnew<-rbind(datasetnew,datasetnew[ind,])
        pred<-c(pred,pred[ind])
      }else{
        w<-as.vector(rdirichlet(1,rep(1,(N+nv-1))))
        ind<-sample(1:(N+nv-1),prob=w,size=1)
        newdata<-datasetnew[ind,]
        newdata[1]<-rnorm(1,pred[ind],1)
        datasetnew<-rbind(datasetnew,newdata)
        pred<-c(pred,pred[ind])
      }
    }
    
    
    mod<-lm(Y~ D+u1+x2+x4+ps,data=datasetnew[-c(1:N),])
    post_theta1<-c(post_theta1,mod$coefficients[2])
  }
  post_theta1 = post_theta1[!is.na(post_theta1)]
  ####
  postmean_theta1<-mean(post_theta1)
  ci<-0
  if(quantile(post_theta1,prob=c(0.025))<1 && quantile(post_theta1,prob=c(0.975))>1) {ci<-1}
  return(list(theta=postmean_theta1,ci=ci))
}

SenarioII_res_drdp_20<-mclapply(sample(c(1:100000000),1000), function(x) SII_dp(x,20))
causal_para_drdp_20<-unlist(mclapply(SenarioI_res_drdp_20, '[[', "theta"))
coverge_drdp_20<-sum(unlist(mclapply(SenarioI_res_drdp_20, '[[', "ci")))/1000



SenarioII_res_drdp_50<-mclapply(sample(c(1:100000000),1000), function(x) SII_dp(x,50))
causal_para_drdp_50<-unlist(mclapply(SenarioI_res_drdp_50, '[[', "theta"))
coverge_drdp_50<-sum(unlist(mclapply(SenarioI_res_drdp_50, '[[', "ci")))/1000

SenarioII_res_drdp_100<-mclapply(sample(c(1:100000000),1000),function(x) SII_dp(x,100))
causal_para_drdp_100<-unlist(mclapply(SenarioI_res_drdp_100, '[[', "theta"))
coverge_drdp_100<-sum(unlist(mclapply(SenarioI_res_drdp_100, '[[', "ci")))/1000


SenarioII_res_drdp_500<-mclapply(sample(c(1:100000000),1000), function(x) SII_dp(x,500))
causal_para_drdp_500<-unlist(mclapply(SenarioI_res_drdp_500, '[[', "theta"))
coverge_drdp_500<-sum(unlist(mclapply(SenarioI_res_drdp_500, '[[', "ci")))/1000




SIII_dp<-function(seed,N){
  set.seed(seed)
  # First generate covariates X, D
  #N=20
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
  
  #PS_lin<-glm.fit(PS_pred,D,family = binomial(link = "logit"))$linear.predictors
  ps <- glm.fit(PS_pred,D,family = binomial(link = "logit"))$fitted.values
  #inv_logit<-function(x){exp(x)/(1+exp(x))}
  
  #ps = inv_logit(PS_lin)
  
  #OR_pred <- model.matrix( ~ D+x1+x2+x4+ps)
  dataset<-data.frame(cbind(Y,D,x1,x2,x4,ps))
  post_theta1<-NULL
  for (i in 1:1000){
    u<-runif(1000)
    av<-1
    datasetnew<-dataset
    pred<-predict(lm(Y ~ D+x1+x2+x4+ps))
    for(nv in 1:1000){
      if(u[nv] > av/(av+N+nv-1)){
        ind<-sample(1:(N+nv-1),size=1)
        datasetnew<-rbind(datasetnew,datasetnew[ind,])
        pred<-c(pred,pred[ind])
      }else{
        w<-as.vector(rdirichlet(1,rep(1,(N+nv-1))))
        ind<-sample(1:(N+nv-1),prob=w,size=1)
        newdata<-datasetnew[ind,]
        newdata[1]<-rnorm(1,pred[ind],1)
        datasetnew<-rbind(datasetnew,newdata)
        pred<-c(pred,pred[ind])
      }
    }
    
    
    mod<-lm(Y~ D+x1+x2+x4+ps,data=datasetnew[-c(1:N),])
    post_theta1<-c(post_theta1,mod$coefficients[2])
  }
  post_theta1 = post_theta1[!is.na(post_theta1)]
  ####
  postmean_theta1<-mean(post_theta1)
  ci<-0
  if(quantile(post_theta1,prob=c(0.025))<1 && quantile(post_theta1,prob=c(0.975))>1) {ci<-1}
  return(list(theta=postmean_theta1,ci=ci))
}


SenarioIII_res_drdp_20<-mclapply(sample(c(1:100000000),1000), function(x) SIII_dp(x,20))
causal_para_drdp_20<-unlist(mclapply(SenarioIII_res_drdp_20, '[[', "theta"))
coverge_drdp_20<-sum(unlist(mclapply(SenarioIII_res_drdp_20, '[[', "ci")))/1000


SenarioIII_res_drdp_50<-mclapply(sample(c(1:100000000),1000), function(x) SIII_dp(x,50))
causal_para_drdp_50<-unlist(mclapply(SenarioIII_res_drdp_50, '[[', "theta"))
coverge_drdp_50<-sum(unlist(mclapply(SenarioIII_res_drdp_50, '[[', "ci")))/1000


SenarioIII_res_drdp_100<-mclapply(sample(c(1:100000000),1000),function(x) SIII_dp(x,100))
causal_para_drdp_100<-unlist(mclapply(SenarioIII_res_drdp_100, '[[', "theta"))
coverge_drdp_100<-sum(unlist(mclapply(SenarioIII_res_drdp_100, '[[', "ci")))/1000



SenarioIII_res_drdp_500<-mclapply(sample(c(1:100000000),1000), function(x) SIII_dp(x,500))
causal_para_drdp_500<-unlist(mclapply(SenarioIII_res_drdp_500, '[[', "theta"))
coverge_drdp_500<-sum(unlist(mclapply(SenarioIII_res_drdp_500, '[[', "ci")))/1000
