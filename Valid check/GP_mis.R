library(gtools)
library(parallel)
options(mc.cores = 22) 

### alpha = 0
BB_Valid<-function(seed,N,sigma_out){
  set.seed(seed)
  
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
 
  beta0<-rnorm(1,0,sigma_out)
  beta1<-rnorm(1,-1,sigma_out)
  beta2<-rnorm(1,-1,sigma_out)
  beta3<-rnorm(1,-1,sigma_out)

  linp<-beta0*D+beta1*u1+beta2*x2+beta3*x4
  sigma_y<-1
  Y<-rnorm(N,mean=linp,sd=sigma_y)
  sim.data<-data.frame(x1=x1,x2=x2,x3=x3,x4=x4,u1=u1,D=as.integer(D),Y=Y)
  
  post_theta1<-NULL
  
  for (i in 1:2000){
    w<-as.vector(rdirichlet(1,rep(1,N)))
    ps <-glm(D ~ x1+x2+x3,data=sim.data,family = binomial(link = "logit"),weights =w)$fitted.values
    
    mod <- lm( Y ~ D+u1+x2+x4+ps, weights = w, data=sim.data)
    
    post_theta1<-c(post_theta1,mod$coefficients[2])
  }
 
  ### calculate log likelihood for prior
 
  H_stat<-length(which(post_theta1[201:2000]<beta0))/1800
  
  return(H_stat=H_stat)
}

BB_Valid_100<-mclapply(sample(c(1:100000000),1000),function(x) BB_Valid(x,100,100))
ks.test(unlist(BB_Valid_100),punif,0,1)
BB_Valid_10k<-mclapply(sample(c(1:100000000),1000),function(x) BB_Valid(x,10000,100))
ks.test(unlist(BB_Valid_10k),punif,0,1)


### alpha > 0
stick.breaking<-function(av,Nv){
  u<-rbeta(Nv,1,av)	
  v<-u
  w<-c(1,cumprod(1-u[-Nv]))
  return(v*w)
}

GP_Valid_mis<-function(seed,N,al,sigma_out){
  set.seed(seed)
  
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
  
  beta0<-rnorm(1,0,sigma_out)
  beta1<-rnorm(1,-1,sigma_out)
  beta2<-rnorm(1,-1,sigma_out)
  beta3<-rnorm(1,-1,sigma_out)
  
  linp<-beta0*D+beta1*u1+beta2*x2+beta3*x4
  sigma_y<-1
  Y<-rnorm(N,mean=linp,sd=sigma_y)
  sim.data<-data.frame(x1=x1,x2=x2,x3=x3,x4=x4,u1=u1,D=as.integer(D),Y=Y)
  
  post_theta1<-NULL
  Nv<-40000
  ps <- glm(D ~ u1+x2+x3,data=sim.data,family = binomial(link = "logit"))$fitted.values
  pred<-predict(lm(Y ~ D+ps))
  sim.data$pred<-pred
  
  for (i in 1:1000){ 
    
    #datasetnew<-sim.data
    ind<-sample(1:(N),size=Nv,replace = TRUE)
    newdata<-sim.data[ind,]
    newdata$u<-runif(Nv)
    res_ind<-which(newdata$u<al/(al+N))
    newdata[res_ind,]$Y<-rnorm(length(res_ind),newdata[res_ind,]$pred,1) 
    
    
    resamp.data<-newdata
    w<-stick.breaking(al+N,Nv)
    resamp.data$ps <-glm(D ~ u1+x2+x3,data=resamp.data,family = binomial(link = "logit"),weights =w)$fitted.values
    
    mod <- lm( Y ~ D+ps, weights = w, data=resamp.data)
    
    post_theta1<-c(post_theta1,mod$coefficients[2])
    
  }
  
  
  H_stat<-length(which(post_theta1<beta0))/1000
  
  return(H_stat=H_stat)
}

## n = 100
##alpha = 1
GP_Valid_mis_1<-mclapply(sample(c(1:100000000),1000),function(x) GP_Valid_mis(x,100,1,100))
##alpha = 10
GP_Valid_mis_10<-mclapply(sample(c(1:100000000),1000),function(x) GP_Valid_mis(x,100,10,100))
##alpha = 100
GP_Valid_mis_100<-mclapply(sample(c(1:100000000),1000),function(x) GP_Valid_mis(x,100,100,100))

ks.test(unlist(GP_Valid_mis_1),punif,0,1)

ks.test(unlist(GP_Valid_mis_10),punif,0,1)

ks.test(unlist(GP_Valid_mis_100),punif,0,1)

## n = 10000
##alpha = 1
GP_Valid_mis_1_10kv1<-mclapply(sample(c(1:100000000),1000),function(x) GP_Valid_mis(x,10000,1,100))
##alpha = 10
GP_Valid_mis_10_10kv1<-mclapply(sample(c(1:100000000),1000),function(x) GP_Valid_mis(x,10000,10,100))
##alpha = 100
GP_Valid_mis_100_10kv1<-mclapply(sample(c(1:100000000),1000),function(x) GP_Valid_mis(x,10000,100,100))


ks.test(unlist(GP_Valid_mis_1_10kv1),punif,0,1)

ks.test(unlist(GP_Valid_mis_10_10kv1),punif,0,1)

ks.test(unlist(GP_Valid_mis_100_10kv1),punif,0,1)
