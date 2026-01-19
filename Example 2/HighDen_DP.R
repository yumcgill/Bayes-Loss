library(mvtnorm)
library(glmnet)
library(gtools)
library(parallel)
options(mc.cores = 23)


expit = function(x) {1/(1+exp(-x))}

stick.breaking<-function(av,Nv){
  u<-rbeta(Nv,1,av)	
  v<-u
  w<-c(1,cumprod(1-u[-Nv]))
  return(v*w)
}

al = 2

high_DP<-function(seed,n){
  set.seed(seed)
  
  p = 20
  sig = array(0.1,c(p,p))
  diag(sig) = 1
  x = matrix(rmvnorm(n,mean=rep(0, nrow(sig)), sigma=sig),n,p)
  d = rbinom(n, 1, p=expit(0.45*x[,1] + 0.9*x[,2] - 0.4*x[,5] +1.3*x[,2]*x[,5]+1.8*x[,1]*x[,2]))
  y = rnorm(n, mean=d + 0.5*x[,1] + x[,3] - 0.1*x[,4]-0.2*x[,7] + 1.5*x[,3]*x[,4] +0.6*x[,7]*x[,7] +1.2*x[,1]*x[,3], sd=1)
  x_des = data.frame(x=x)
  x_int = model.matrix(~.^2-1, data=x_des)
  
  sim.data<-data.frame(x=x_int,d=as.integer(d),y=y)
  Nv<-100
  ps <- glm(d ~ x[,1]+x[,2]+x[,5]+x[,2]*x[,5]+x[,1]*x[,2],family = binomial(link = "logit"))$fitted.values
  pred<-predict(lm(y ~ d+ps))
  sim.data$pred<-pred
  
  start_time<- Sys.time()
  post_theta1<-NULL
  
  
  for (i in 1:1000){
    

    ind<-sample(1:(n),size=Nv,replace = TRUE)
    newdata<-sim.data[ind,]
    newdata$u<-runif(Nv)
    res_ind<-which(newdata$u< al/(al+n))
    newdata[res_ind,]$y<-rnorm(length(res_ind),newdata[res_ind,]$pred,1)
    
    resamp.data<-newdata
    xnew = as.matrix(resamp.data[,-c(212:215)])
    dnew = resamp.data$d
    
    
    w<-stick.breaking(al+n,Nv)
    
    cv_model <- cv.glmnet(xnew, dnew, alpha = 1, family = "binomial")
    resamp.data$ps_est<-predict(cv_model,s="lambda.min", newx = xnew,type="response")
    mod <- lm( y ~ d+ps_est, weights = w, data=resamp.data)
    
    post_theta1<-c(post_theta1,mod$coefficients[2])
    
    
  } 
  ena_time<- Sys.time()
  run_time <- ena_time - start_time
  postmean_theta1<-mean(post_theta1)
  
  ci<-0 
  
  if(quantile(post_theta1,prob=c(0.025))<1 && quantile(post_theta1,prob=c(0.975))>1) {ci<-1}
  return(list(theta=postmean_theta1,ci=ci,run_time=run_time))
}


set.seed(1233)
high_DP_50<-mclapply(sample(c(1:100000000),500),function(x) high_DP(x,50))
para_high_DP_50<-unlist(mclapply(high_DP_50, '[[', "theta"))
coverge_high_DP_50<-sum(unlist(mclapply(high_DP_50, '[[', "ci")))/500
coverge_high_DP_50
mean(para_high_DP_50);var(para_high_DP_50)
mean((para_high_DP_50-1)^2)
mean(unlist(mclapply(high_DP_50, '[[', "run_time")))

high_DP_100<-mclapply(sample(c(1:100000000),500),function(x) high_DP(x,100))
para_high_DP_100<-unlist(mclapply(high_DP_100, '[[', "theta"))
coverge_high_DP_100<-sum(unlist(mclapply(high_DP_100, '[[', "ci")))/500
coverge_high_DP_100
mean((para_high_DP_100-1)^2)
mean(para_high_DP_100);var(para_high_DP_100)
mean(unlist(mclapply(high_DP_100, '[[', "run_time")))

