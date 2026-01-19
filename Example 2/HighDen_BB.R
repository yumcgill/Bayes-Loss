library(mvtnorm)
library(glmnet)
library(gtools)
library(parallel)
options(mc.cores = 23)


expit = function(x) {1/(1+exp(-x))}

high_BB_pois<-function(seed,n){
  set.seed(seed)
  p = 20
  sig = array(0.1,c(p,p))
  diag(sig) = 1
  x = matrix(rmvnorm(n,mean=rep(0, nrow(sig)), sigma=sig),n,p)
  d = rbinom(n, 1, p=expit(0.45*x[,1] + 0.9*x[,2] - 0.4*x[,5] +1.3*x[,2]*x[,5]+1.8*x[,1]*x[,2]))
  y = rnorm(n, mean=d + 0.5*x[,1] + x[,3] - 0.1*x[,4]-0.2*x[,7] + 1.5*x[,3]*x[,4] +0.6*x[,7]*x[,7] +1.2*x[,1]*x[,3], sd=1)
  
  x_des = data.frame(x=x)
  x_int = model.matrix(~.^2-1, data=x_des)
  
  start_time<- Sys.time()
  
  post_theta1<-NULL
  
  for (i in 1:1000){
    w<-as.vector(rdirichlet(1,rep(1,n)))
    cv_model <- cv.glmnet(x_int, d, alpha = 1, family = "binomial",weights=w)
    res_pred <-predict(cv_model,s="lambda.min", newx = x_int,type="response")
    mod = lm(y ~ d+ res_pred,weights=w)
    post_theta1<-cbind(post_theta1,mod$coefficients[2])
  }
  
  ena_time<- Sys.time()
  run_time <- ena_time - start_time
  
  postmean_theta1<-median(post_theta1)
  
  ci<-0
  
  if(quantile(post_theta1,prob=c(0.025))<1 && quantile(post_theta1,prob=c(0.975))>1) {ci<-1}
  return(list(theta=postmean_theta1,ci=ci,run_time=run_time))
}

set.seed(453)

high_BB_pois_50<-mclapply(sample(c(1:100000000),500),function(x) high_BB_pois(x,50))
para_high_BB_pois_50<-unlist(mclapply(high_BB_pois_50, '[[', "theta"))
coverge_high_BB_pois_50<-sum(unlist(mclapply(high_BB_pois_50, '[[', "ci")))/500
coverge_high_BB_pois_50
mean(para_high_BB_pois_50);var(para_high_BB_pois_50)
sqrt(mean((para_high_BB_pois_50-1)^2))
mean(na.omit(as.numeric(unlist(mclapply(high_BB_pois_50, '[[', "run_time")))))





high_BB_pois_100<-mclapply(sample(c(1:100000000),500),function(x) high_BB_pois(x,100))
para_high_BB_pois_100<-unlist(mclapply(high_BB_pois_100, '[[', "theta"))
coverge_high_BB_pois_100<-sum(unlist(mclapply(high_BB_pois_100, '[[', "ci")))/500
coverge_high_BB_pois_100
mean(para_high_BB_pois_100);var(para_high_BB_pois_100)
sqrt(mean((para_high_BB_pois_100-1)^2))
mean(na.omit(as.numeric(unlist(mclapply(high_BB_pois_100, '[[', "run_time")))))


