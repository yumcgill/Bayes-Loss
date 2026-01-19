library(gtools)
library(parallel)
options(mc.cores = 23)
stick.breaking<-function(av,Nv=1000){
  u<-rbeta(Nv,1,av)	
  v<-u
  w<-c(1,cumprod(1-u[-Nv]))
  return(v*w)
}

al<-5
Nv <- 5000


DP_flex<-function(seed,N){
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
  sim.data<-data.frame(x1=x1,x2=x2,x3=x3,x4=x4,D=as.integer(D),Y=Y)
  
  
  
  
  
  
  post_theta1<-NULL
  
  for (i in 1:1000){ 
    
    ps <- glm(D~x1+x2+x3+x4,data=sim.data,family = binomial(link = "logit"))$fitted.values
    pred<-predict(lm(Y ~ D+I(D*x1)+ps+I(ps*x1)))
    u<-runif(Nv)
    datasetnew<-sim.data
    
    for(nv in 1:Nv){
      if(u[nv] > al/(al+N)){
        ind<-sample(1:(N),size=1)
        datasetnew<-rbind(datasetnew,sim.data[ind,])
        #pred<-c(pred,pred[ind])
      }else{
        ind<-sample(1:(N),size=1)
        newdata<-sim.data[ind,]
        newdata[1]<-rnorm(1,pred[ind],1)
        datasetnew<-rbind(datasetnew,newdata)
        #pred<-c(pred,pred[ind])
      }
    }  
    
    resamp.data<-datasetnew[-c(1:N),]
    w<-stick.breaking(al+N,Nv)
    resamp.data$ps <-glm(D ~ x1+x2+x3+x4,data=resamp.data,family = binomial(link = "logit"),weights =w)$fitted.values
    
    data.alltrt <- resamp.data
    data.alltrt$D <- 1
    data.nontrt <- resamp.data
    data.nontrt$D <- 0
    
    mod1.lmX <- lm(Y ~ D+I(D*x1)+ps+I(ps*x1),weights=w,data=resamp.data)
    
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


DP_flex_200<-mclapply(sample(c(1:100000000),1000),function(x) DP_flex(x,200))
DP_flex_200_ate<-unlist(mclapply(DP_flex_200, '[[', "theta"))
coverge_DP_flex_200<-sum(unlist(mclapply(DP_flex_200, '[[', "ci")))/1000
coverge_DP_flex_200
mean(DP_flex_200_ate);var(DP_flex_200_ate)
mean(DP_flex_200_ate) -3
sqrt(mean((DP_flex_200_ate - 3)^2))


DP_flex_1000<-mclapply(sample(c(1:100000000),1000),function(x) DP_flex(x,1000))
DP_flex_1000_ate<-unlist(mclapply(DP_flex_1000, '[[', "theta"))
coverge_DP_flex_1000<-sum(unlist(mclapply(DP_flex_1000, '[[', "ci")))/1000
coverge_DP_flex_1000
mean(DP_flex_1000_ate);var(DP_flex_1000_ate)
mean(DP_flex_1000_ate)-3
sqrt(mean((DP_flex_1000_ate - 3)^2))




DP_flex_2000<-mclapply(sample(c(1:100000000),1000),function(x) DP_flex(x,2000))
DP_flex_2000_ate<-unlist(mclapply(DP_flex_2000, '[[', "theta"))
coverge_DP_flex_2000<-sum(unlist(mclapply(DP_flex_2000, '[[', "ci")))/1000
coverge_DP_flex_2000
mean(DP_flex_2000_ate);var(DP_flex_2000_ate)
mean(DP_flex_2000_ate)-3
sqrt(mean((DP_flex_2000_ate - 3)^2))

