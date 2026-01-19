#######################################################################################################################################################

###################### Loading packages ###########################

library(foreign);
library(quantreg);
library(mnormt);
library(gbm);
library(glmnet);
library(MASS);
library(rpart);
library(doParallel)
library(sandwich);
library(hdm);
library(randomForest);
library(nnet)
library(matrixStats)
library(quadprog)
library(xtable)

################ Loading functions and Data ########################


rm(list = ls())  # Clear everything out so we're starting clean
source("ML_Functions.R")  
source("Moment_Functions.R")  
options(warn=-1)
set.seed(1210);
cl <- makeCluster(12, outfile="")

########################### Sample Construction ######################


DML2_flex<-function(seed,N){
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

 data<-data.frame(cbind(GDP=Y,Exprop=D,Latitude=x1,Latitude2=x2,Africa =x3,Asia =x4))

################################ Inputs ##############################

# Outcome Variable
y      <- "GDP";

# Treatment Indicator
d      <- "Exprop";  

# Controls
x      <- "Latitude + Latitude2 + Africa +  Asia "        # use this for tree-based methods like forests and boosted trees
xl     <- "(Latitude + Latitude2 + Africa +  Asia)^2";     # use this for rlasso etc.

# Method names: Boosting, Nnet, RLasso, PostRLasso, Forest, Trees, Ridge, Lasso, Elnet, Ensemble

Boosting     <- list(n.minobsinnode = 1, bag.fraction = .5, train.fraction = 1.0, interaction.depth=2, n.trees=1000, shrinkage=.01, n.cores=1, cv.folds=2, verbose = FALSE, clas_dist= 'adaboost', reg_dist='gaussian')
Forest       <- list(clas_nodesize=3, reg_nodesize=5, ntree=1000, na.action=na.omit, replace=TRUE)
RLasso       <- list(intercept = TRUE)
Nnet         <- list(size=2,  maxit=1000, decay=0.01, MaxNWts=10000,  trace=FALSE)
Trees        <- list(reg_method="anova", clas_method="class")

arguments    <- list(Boosting=Boosting, Forest=Forest,  Nnet=Nnet, Trees=Trees)

ensemble     <- list(methods=c( "Boosting", "Forest"))                       # methods for the ensemble estimation
methods      <- c("Trees", "Boosting", "Forest", "Nnet","Ensemble")          # ML methods that are used in estimation
#methods      <- c("RLasso","RLasso")          # ML methods that are used in estimation

split        <- 50                                                                                                                   # number of splits


################################ Estimation ##################################################

############## Arguments for DoubleML function:

# data:     : data matrix
# y         : outcome variable
# d         : treatment variable
# z         : instrument
# xx        : controls for tree-based methods
# xL        : controls for penalized linear methods
# methods   : machine learning methods
# DML       : DML1 or DML2 estimation (DML1, DML2)
# nfold     : number of folds in cross fitting
# est       : estimation methods (IV, LATE, plinear, interactive)
# arguments : argument list for machine learning methods
# ensemble  : ML methods used ine ensemble method
# silent    : whether to print messages
# trim      : bounds for propensity score trimming


r <- foreach(k = 1:split, .combine='rbind', .inorder=FALSE, .packages=c('MASS','randomForest','neuralnet','gbm', 'sandwich', 'hdm', 'nnet', 'rpart','glmnet')) %dopar% { 
  
  dml <- DoubleML(data=data, y=y, d=d, z=NULL, xx=x, xL=xl, methods=methods, DML="DML2", nfold=2, est="plinear", arguments=arguments, ensemble=ensemble, silent=FALSE, trim=NULL) 
  
  data.frame(t(dml[1,]), t(dml[2,]))
  
}

################################ Compute Output Table ########################################
r<-as.matrix(r)


postmean_theta1<-colQuantiles(r[,1:(length(methods)+1)], probs=0.5)
se1<-colQuantiles(r[,(length(methods)+2):ncol(r)], probs=0.5)
ci<-ifelse((postmean_theta1-1.96*se1)<3 & (postmean_theta1+1.96*se1)>3,1,0)

return(list(theta=postmean_theta1,ci=ci))
}




DML2_flex_200<-mclapply(sample(c(1:100000000),1000),function(x) DML2_flex(x,200))
DML2_flex_200_ate<-do.call(rbind, (mclapply(DML2_flex_200, '[[', "theta")))
coverge_DML2_flex_200<-colSums(do.call(rbind,(mclapply(DML2_flex_200, '[[', "ci"))))/1000
coverge_DML2_flex_200
apply(DML2_flex_200_ate,2,mean)
apply(DML2_flex_200_ate,2,var)
apply(DML2_flex_200_ate,2,mean) -3
apply(DML2_flex_200_ate,2,function(x) sqrt(mean((x-3)^2))) 


DML2_flex_1000<-mclapply(sample(c(1:100000000),1000),function(x) DML2_flex(x,1000))
DML2_flex_1000_ate<-do.call(rbind, (mclapply(DML2_flex_1000, '[[', "theta")))
coverge_DML2_flex_1000<-colSums(do.call(rbind,(mclapply(DML2_flex_1000, '[[', "ci"))))/1000
coverge_DML2_flex_1000
apply(DML2_flex_1000_ate,2,mean)
apply(DML2_flex_1000_ate,2,var)
apply(DML2_flex_1000_ate,2,mean) -3
apply(DML2_flex_1000_ate,2,function(x) sqrt(mean((x-3)^2))) 



DML2_flex_2000<-mclapply(sample(c(1:100000000),1000),function(x) DML2_flex(x,2000))
DML2_flex_2000_ate<-do.call(rbind, (mclapply(DML2_flex_2000, '[[', "theta")))
coverge_DML2_flex_2000<-colSums(do.call(rbind,(mclapply(DML2_flex_2000, '[[', "ci"))))/1000
coverge_DML2_flex_2000
apply(DML2_flex_2000_ate,2,mean)
apply(DML2_flex_2000_ate,2,var)
apply(DML2_flex_2000_ate,2,mean) -3 
apply(DML2_flex_2000_ate,2,function(x) sqrt(mean((x-3)^2))) 

