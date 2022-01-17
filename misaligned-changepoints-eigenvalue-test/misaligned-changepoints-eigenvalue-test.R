
library(docstring)
library(hdbinseg)
library(InspectChangepoint)
library(changepoint.geo)
library(RMTstat)
 

#----------------------------
# 
# Functions
# 
#----------------------------


sim.cpt.panel <- function(p, n, cpt.region, jump = 5, prob = 0.5) {
  
  #' Simulate a panel with chpt occuring in a region
  #'
  #'@param n, int time points 
  #'@param p, int no. chanels 
  #'@param cpt.region vector, subset of [1,n] where cpt may occur
  #'@param jump float, size of mean change
  #'@param prob float, proportion of chanels which change
  
  X <- matrix(rnorm(p*n),p,n)
  
  for (i in 1:p) {
    
    tau <- ifelse(length(cpt.region) == 1, cpt.region, sample(cpt.region, 1))
    
    i.jump <- abs(rnorm(1)) + jump * rbinom(1,1,prob)
    
    X[i,(tau+1):n] <- X[i,(tau+1):n] + i.jump 
    
  }
  
  return(X)
  
}


cpt.demean <- function(X, cpt.locs) {
  
  #' De-mean panel using changepoint locations
  #' 
  #'@param X matrix, dim p \times n  
  #'@param cpt.locs vector
  
  
  if (length(cpt.locs) == 0) return(t(t(apply(X, 1, function(x) x - mean(x)))))
  
  n <- dim(X)[2]; p <- dim(X)[1]
  
  cpts <- c(0, cpt.locs, n)
  
  for (i in 1:p) for (j in 1:(length(cpts)-1)) X[i,(cpts[j]+1):cpts[j+1]] <-  X[i,(cpts[j]+1):cpts[j+1]] - mean(X[i,(cpts[j]+1):cpts[j+1]])
  
  return(X)
  
}


wishart.eig <- function(X) {
  
  #' Wishart eigenvalue
  #' 
  #' Largest eigenvalue of sample covariance matrix for de-meaned panel. 
  #' Normalized to converge to a Tracy-Widom (\beta = 1) law at the limit. 
  #'
  #'@param X matrix, dim n \times n 
  #'
  #'@references hnps://doi.org/10.1214/aos/1009210544
  
  n <- dim(X)[2]; p <- dim(X)[1]
  
  mu <- (sqrt(n - 1) + sqrt(p)) ** 2
  
  sigma <-  (sqrt(n-1) + sqrt(p)) * ((1/sqrt(n-1) + 1/sqrt(p)) ** (1/3))
  
  lambda.1 <- eigen(X %*% t(X))$values[1]
  
  
  (lambda.1 - mu) / sigma
  
}


#-----------------------
#
# Simulation prelims
#
#-----------------------


set.seed(42)

n <- 500 # no. channel time points 

p <- 500 # no. pannel channels 

K <- 1000 # no. Monte Carlo runs

sparse.jump.size <- 1.5 

sparse.jump.prob <- 2/5 

cpt.loc.alig <- 50 

cpt.region.mis <- 45:50


#---------------------------------------
#
# Max eigenvalues with aligned changes
#
#---------------------------------------

inspect.alig <- numeric(K); sbs.alig <- numeric(K); dcbs.alig <- numeric(K); geom.alig <- numeric(K)

inspect.mis <- numeric(K); sbs.mis <- numeric(K); dcbs.mis <- numeric(K); geom.mis <- numeric(K)

pb <- txtProgressBar(1,K, style = 3)


for (k in 1:K) {
  
  
  X.alig <- sim.cpt.panel(p, n, cpt.loc.alig, sparse.jump.size, sparse.jump.prob)
  
  X.mis <- sim.cpt.panel(p, n, cpt.region.mis, sparse.jump.size, sparse.jump.prob)
  
  
  Y <- cpt.demean(X.alig, inspect(X.alig)$changepoints[,1])
  inspect.alig[k] <- wishart.eig(Y)
  
  Y <- cpt.demean(X.mis, inspect(X.mis)$changepoints[,1])
  inspect.mis[k] <- wishart.eig(Y)


  Y <- cpt.demean(X.alig, sbs.alg(X.alig)$ecp)
  sbs.alig[k] <- wishart.eig(Y)

  Y <- cpt.demean(X.mis, sbs.alg(X.mis)$ecp)
  sbs.mis[k] <- wishart.eig(Y)


  Y <- cpt.demean(X.alig, dcbs.alg(X.alig)$ecp)
  dcbs.alig[k] <- wishart.eig(Y)

  Y <- cpt.demean(X.mis, dcbs.alg(X.mis)$ecp)
  dcbs.mis[k] <- wishart.eig(Y)


  Y <- cpt.demean(X.alig, geomcp(t(X.alig))@dist.cpts)
  geom.alig[k] <- wishart.eig(Y)

  Y <- cpt.demean(X.mis, geomcp(t(X.mis))@dist.cpts)
  geom.mis[k] <- wishart.eig(Y)
  
  
  setTxtProgressBar(pb, k)
  
}


boxplot(cbind(inspect.alig, sbs.alig, dcbs.alig, geom.alig),
        names = c("inspect", "sbs", "dcbs", "geom"),
        outline = F)

abline(h = qtw(0.95), col = "red", lty = 2)



boxplot(cbind(inspect.mis, sbs.mis, dcbs.mis, geom.mis),
        names = c("inspect", "sbs", "dcbs", "geom"),
        outline = F)

abline(h = qtw(0.95), col = "red", lty = 2)


