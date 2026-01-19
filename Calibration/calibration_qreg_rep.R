library(devtools)
install_github("nasyring/GPC", subdir = "GPC")
library(GPC)
library(quantreg)
library(gtools)
library(matrixcalc)
library(expm)
library(parallel)
options(mc.cores = 23)
# Modifying the number of threads used for parallel computation
RcppParallel::setThreadOptions(numThreads = 1)

# data-generating process
rmodel <- function(n, theta) {
  
  X <- rchisq(n, 2) - 2
  Y <- theta[1] + theta[2] * X + rnorm(n, 0, 2)
  return(cbind(X, Y))
  
}


one_iteration_qr <- function(n = 200, 
                             theta.true = c(2,1), 
                             B = 200, 
                             M = 2000, 
                             alpha = 0.05, 
                             eps = 0.05, 
                             w_init = 0.5, seed) {
  set.seed(seed)
  # -----------------------------
  # Data generation
  # -----------------------------
  X <- rchisq(n, 2) - 2
  Y <- theta.true[1] + theta.true[2] * X + rnorm(n, 0, 2)
  data <- cbind(X, Y)
  
  # -----------------------------
  # Initial quantile regression
  # -----------------------------
  theta.hat <- rq(Y ~ X, tau = 0.5)$coef
  
  # -----------------------------
  # Bayesian bootstrap samples
  # -----------------------------
  theta.BB <- matrix(NA, 1000, 2)
  for(bb in 1:1000) {
    w_bb <- as.vector(rdirichlet(1, rep(1, n)))
    theta.BB[bb, ] <- rq(Y ~ X, tau = 0.5, weights = w_bb)$coef
  }
  
  # -----------------------------
  # Bootstrap samples
  # -----------------------------
  data.star <- matrix(0, n, B*2)
  theta.boot <- matrix(0, B, 2)
  for(b in 1:B) {
    id <- sample(n, n, replace = TRUE)
    data.star[, (2*b-1):(2*b)] <- data[id, ]
    theta.boot[b, ] <- as.vector(rq(data.star[, 2*b] ~ data.star[, 2*b-1], tau = 0.5)$coef)
  }
  
  # -----------------------------
  # Initial Gibbs posterior
  # -----------------------------
  w <- w_init
  inital.gibbs <- GibbsMCMC2(n, data, theta.boot,
                             mean(theta.boot[,1]), 0.05,
                             0.1, mean(theta.boot[,2]),
                             alpha, M, w)
  
  diff_mat <- logm(cov(cbind(inital.gibbs$postsamples0,inital.gibbs$postsamples1))) - logm(cov(theta.BB))
  eta.scale <- w * norm(diff_mat, type = "F")
  
  # -----------------------------
  # Learning rate adjustment
  # -----------------------------
  go <- TRUE
  t <- 1
  k <- function(t) (1 + t)^(-0.51)
  
  while(go) {
    cover <- rcpp_parallel_qr(n, data, theta.boot,
                              mean(theta.boot[,1]),
                              mean(theta.boot[,2]),
                              0.1, data.star,
                              alpha, M, B, w)
    diff <- mean(cover) - (1 - alpha)
    if(abs(diff) <= eps || t > 16) {
      go <- FALSE
    } else {
      t <- t + 1
      w <- w + k(t) * diff
    }
  }
  
  # -----------------------------
  # Final Gibbs Posterior
  # -----------------------------
  final.gibbs <- GibbsMCMC2(n, data, theta.boot,mean(theta.boot[,1]),mean(theta.boot[,2]),0.1, 0.1,alpha, M, w)
  
  BB.gibbs <- GibbsMCMC2(n, data, theta.boot,mean(theta.boot[,1]),mean(theta.boot[,2]), 0.1, 0.1, alpha, M, eta.scale)
  
  # -----------------------------
  # Coverage and interval lengths
  # -----------------------------
  cover_SM <- c((final.gibbs$l0 < theta.true[1] & final.gibbs$u0 > theta.true[1]),
                (final.gibbs$l1 < theta.true[2] & final.gibbs$u1 > theta.true[2]))
  
  length_SM <- c(final.gibbs$u0 - final.gibbs$l0,
                 final.gibbs$u1 - final.gibbs$l1)
  
  cover_BBbase <- c((BB.gibbs$l0 < theta.true[1] & BB.gibbs$u0 > theta.true[1]),
                    (BB.gibbs$l1 < theta.true[2] & BB.gibbs$u1 > theta.true[2]))
  
  length_BBbase <- c(BB.gibbs$u0 - BB.gibbs$l0,
                     BB.gibbs$u1 - BB.gibbs$l1)
  
  cover_BB <- c((quantile(theta.BB[,1], 0.025) < theta.true[1] & quantile(theta.BB[,1], 0.975) > theta.true[1]),
                (quantile(theta.BB[,2], 0.025) < theta.true[2] & quantile(theta.BB[,2], 0.975) > theta.true[2]))
  
  length_BB <- c(quantile(theta.BB[,1], 0.975) - quantile(theta.BB[,1], 0.025),
                 quantile(theta.BB[,2], 0.975) - quantile(theta.BB[,2], 0.025))
  
  # -----------------------------
  # Return results
  # -----------------------------
  list(
    w = w,
    eta.scale = eta.scale,
    cover_SM = cover_SM,
    length_SM = length_SM,
    cover_BBbase = cover_BBbase,
    length_BBbase = length_BBbase,
    cover_BB = cover_BB,
    length_BB = length_BB
  )
}

set.seed(121)
qr_res<-mclapply(sample(c(1:100000000),1000), function(x) one_iteration_qr(seed = x))
apply(do.call("rbind", (mclapply(qr_res, '[[', "cover_BBbase"))),2,mean)
apply(do.call("rbind", (mclapply(qr_res, '[[', "length_BBbase"))),2,mean)
mean(apply(do.call("rbind", (mclapply(qr_res, '[[', "cover_BBbase"))),1,sum)==2)

apply(do.call("rbind", (mclapply(qr_res, '[[', "cover_SM"))),2,mean)
apply(do.call("rbind", (mclapply(qr_res, '[[', "length_SM"))),2,mean)
mean(apply(do.call("rbind", (mclapply(qr_res, '[[', "cover_SM"))),1,sum)==2)

apply(do.call("rbind", (mclapply(qr_res, '[[', "cover_BB"))),2,mean)
apply(do.call("rbind", (mclapply(qr_res, '[[', "length_BB"))),2,mean)
mean(apply(do.call("rbind", (mclapply(qr_res, '[[', "cover_BB"))),1,sum)==2)
