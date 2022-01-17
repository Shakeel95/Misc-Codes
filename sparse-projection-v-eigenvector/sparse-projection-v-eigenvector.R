library("docstring") # EZ documentation
library("InspectChangepoint") # Wang and Samworth package
library("sparsepca") # sparse PCA


#------------------------------------------------
# Find changepoints by projection onto eigenbasis
#------------------------------------------------

naive.eigen <- function(X, least_varying = FALSE){
  
  #'Locate changepoint with leading eigenvalue projection
  #'
  #'@param X a T x n matrix
  
  C <- cov(X)
  
  if (least_varying) {
    
    n <- ncol(C)
    v <- eigen(C)$vectors[,n]
    
  } else {
    
    v <- eigen(C)$vectors[,1]
    
  }
  
  
  return(list(direction = v, changepoint = which.min(cusum.transform(X %*% v))))
}


fixed.projection <- function(X, v) {
  
  #' Locate changepoint by taking CUSUM of data projected onto a fixed projection direction 
  #' (possibly the oracle direction)
  #' 
  #'@param X a T x n matrix
  #'@param v an n x 1 vector 
  
  which.max(cusum.transform(X %*% v))
  
}


vector.angle <- function(a,b) {
  
  #'Angle between two vectors
  #'
  #'@param a n x 1 vector
  #'@param b n x 1 vector
  
  theta <- abs(acos(sum(a*b) / sqrt(sum(a**2) * sum(b**2))))
  
  if (theta < pi){
    
    return(min(theta, pi - theta))
    
  } else {
    
    theta <- theta - pi
    
    return(min(theta, pi - theta))
    
  }
}

#------------------------------------------------
# Sparse v. dense changes
# 
# visual illustration of the need for sparse PCA
#------------------------------------------------

# set T and n 
set.seed(42)
n = 500

# simulate dense change
X <- matrix(rnorm(n*n),n,n)
X[250:n,] <- X[250:n,] + rep(1,n)
v <- eigen(cov(X))$vectors[,1]

plot.ts(-c(cusum.transform(X %*% v)), 
        ylab = "projection CUSUM",
        main = "dense change")
abline( v = 250, col = "red")


# plot sparse change
Y <- matrix(rnorm(n*n),n,n)
Y[250:n,] <- Y[250:n,] + c(rep(100,.01*n),numeric(n*(1-.01)))
v <- eigen(cov(Y))$vectors[,1]

plot.ts(-c(cusum.transform(Y %*% v)),
        ylab = "projection CUSUM", 
        main = "sparse changepoint")
abline( v = 250, col = "red")


#---------------------------------------
# Simulation:
#
# compare parse projection to eigenvalue
#---------------------------------------

set.seed(42)

n <- 1000 # panel dimension (T = n)
n_runs <- 100 # 
tau <- 500 # changepoint location
norm.const <- c(1,sqrt(0.1),n**(-.25), sqrt(log(n)/n),sqrt(log(log(n))/n))

res.inspect <- matrix(,n_runs,5)
angle.inspect <- matrix(,n_runs,5)

res.eig <- matrix(,n_runs,5)
angle.eig <- matrix(,n_runs,5)

res.oracle <- matrix(,n_runs,5)


for (i in 1:5){
  
  for (j in 1:n_runs){
    
    oracle <- rep(norm.const[i],n) / sqrt(sum(rep(norm.const[i],n)**2))
    X <- matrix(rnorm(n*n),n,n)
    X[tau:n,] <- X[tau:n,] + rep(norm.const[i],n)
    
    inspect.angel <- sparse.svd(cusum.transform(t(X)),.5)
    eigen.sol <- naive.eigen(X)
    
    res.eig[j,i] <- eigen.sol$changepoint
    res.inspect[j,i] <- locate.change(t(X))$changepoint
    res.oracle[j,i] <- fixed.projection(X,oracle)
    
    angle.eig[j,i] <- vector.angle(oracle, eigen.sol$direction)
    angle.inspect[j,i] <- vector.angle(oracle, inspect.angel)
    
    
    cat("\r","finished simulation ", paste0(j), " for sparsity level ", paste0(i))
  }
  
}

#---------------------------
# Inspect simulation results 
#---------------------------

# mean absolute error
apply(res.inspect - tau + 1,2,FUN = function(i) mean(abs(i)))
apply(res.eig - tau +1,2,FUN = function (i) mean(abs(i)))

# median absolute error
apply(res.inspect - tau + 1,2,FUN = function(i) median(abs(i)))
apply(res.eig - tau +1,2,FUN = function (i) median(abs(i)))

# mean angle between estimate and oracle
apply(angle.inspect,2,mean)*(180/pi)
apply(angle.eig,2,mean)*(180/pi)

# median angle between estimate and oracle
apply(angle.inspect,2,median)*(180/pi)
apply(angle.eig,2,median)*(180/pi)


dd <- data.frame(apply(res.inspect - tau + 1,2,FUN = function(i) mean(abs(i))),
                 apply(res.eig - tau + 1,2,FUN = function (i) mean(abs(i))),
                 apply(res.oracle - tau + 1, 2, FUN = function(i) mean(abs(i))),
  
                 apply(res.inspect - tau + 1, 2, FUN = function(i) median(abs(i))),
                 apply(res.eig - tau + 1, 2, FUN = function (i) median(abs(i))),
                 apply(res.oracle - tau + 1, 2, FUN = function(i)median(abs(i))),
                  
                 apply(angle.inspect,2,mean)*(180/pi),
                 apply(angle.eig,2,mean)*(180/pi),
                  
                 apply(angle.inspect,2,median)*(180/pi),
                 apply(angle.eig,2,median)*(180/pi)
)

names(dd) <- c("MAE-inspect","MAE-eig","MAE-oracle",
               "medAE-inspect","medAE-eig","medAE-oracle",
               "MAA-inspect","MAA-eig",
               "medAA-inspect","medAA-eig")
View(round(dd,2))

