library("docstring") # EZ documentation
library("tsDyn") # simulate VAR
library("mvtnorm") # multivariate normal

#----------------------------------------------------------------------
# 
# Function to calculate U-statistic
# 
# Efficient implementation based on reccursion requires O(n^2 * p) time 
# See pseudo-code + explanation in He et al. (2020)
# pp. 18 algorithm 1
# 
#-----------------------------------------------------------------------


U.statistic <- function(X,a,tau) {
  
  #' Two sample U-statistic for changepoint problems
  #'
  #'@param X matrix, data matrix of dimension p x n 
  #'@param a int, order of U-statistics to compute
  #'@param tau int, where to split the sample for two sample test
  
  
  V.k.vec <- function(x,y,a) {
    
    #' Creates a list of intermediate statistics V_k to construct U-stat
    #'
    #'@param x vector, left sample
    #'@param y vector, right sample
    
    
    n.x <- length(x)
    n.y <- length(y)
    V.k <- c()
    
    # build intermediate V-stats
    
    for (i in 1:a) {
      
      v <- 0
      
      for (j in 1:n.x) {
        for (k in 1:n.y) {
          v <- v + ((x[j] - y[k])**i)
        }
      }
      
      V.k <- c(V.k, v)
      
    }
    
    return(V.k)
  }
  
  
  U.l <- function(V.k, a) {
    
    #' Aggregate outputs from V.k.list to build univariate U-statistic
    #'
    #'@param V.k list, output from V.k.list
    #'@param a int, order of U-stats to compute
    
    
    # for a = {1,2} compute values and return
    
    if (a == 1) return(V.k[1])
    
    if (a == 2) return(c(V.k[1], V.k[1]**2 - V.k[2]))
    
    U.vec <- c(V.k[1], (V.k[1]**2) - V.k[2])
    
    # for higher order U-stat build via recursion
    
    for (r in 3:a) {
      
      T.l <- V.k[r]
    
      for (k in rev(1:(r-1))) {
        
        T.l <- V.k[k]*U.vec[r-k] - (r-k)*T.l
        
      }
      
      U.vec <- c(U.vec, T.l)
      
    }
    
    return(U.vec)
    
  }
  
  
  # get data dimensions  
  
  p <- nrow(X)
  n <- ncol(X)
  
  # split sample at estimated changepoint
  
  X.L <- X[,1:tau]
  X.R <- X[,((tau+1):n)]
  
  # sum U-stats over channel coordinates
  
  U.sum <- numeric(a)
  
  for (i in 1:p) {
    
    V.k.vals <- V.k.vec(X.L[i,],X.R[i,] ,a)
    U.sum <- U.sum + U.l(V.k.vals, a)
    
  }
  
  for (i in 1:a) {
    
    U.sum[i] <- U.sum[i] / ((factorial(tau)/factorial(tau - i))*(factorial(n-tau)/factorial(n-tau-i)))
    
  }
  
  return(U.sum)
  
}


#------------------------------------------------
#
# Simulation: Bias under serial correlation
# 
# Investigates three settings: 
#   1. Indipendent Gaussian errors
#   2. Spatially independent AR(1) errors
#   3. Spatially Toeplitz correlated AR(1) errors
#   
#------------------------------------------------


set.seed(42)              # life, the universe and everything
p <- c(100,500,1000)      # no. channels
n <- 100                  # sample size
a <- 10                   # U-stat order 
tau <- 50                 # where to split sample
sim.size <- 100           # number of simulations to run


## Independent standard Gaussian errors

print("Independent standard Gaussian errors")

res.gaussian.errors <- matrix(,3,a)

for (i in 1:3) {
  
  res <- numeric(a)
  
  for (j in 1:sim.size) {
    
    XX <- matrix(rnorm(n*p[i]), nrow = p[i], ncol = n)
    res <- res + U.statistic(XX,a,tau)
    
    cat("\r","finished simulation ", paste0(j), " for sample size ", paste0(p[i]))
    
  }
  
  res.gaussian.errors[i,] <- res / sim.size
  
}



## spatially independent AR(1) errors

print("spatially independent AR(1) errors")

res.AR1.diagonal <- matrix(,3,a)

for (i in 1:3) {
  
  res <- numeric(a)
  
  for (j in 1:sim.size) {
    
    XX <- matrix(,p[i],n)
    for (k in 1:p[i]) {
      XX[k,] <- arima.sim(model = list(ar = 0.5), n)
    }
    
    res <- res + U.statistic(XX,a,tau)
    
    cat("\r","finished simulation ", paste0((j)), "for sample size ", paste0(p[i]))
    
  }
  
  res.AR1.diagonal[i,] <- res / sim.size
  
}


## diagonal VAR(1) serial correlation component with Toeplitz spatial correlation structure

print("diagonal VAR(1) serial correlation component with Toeplitz spatial correlation structure")

res.AR1.toeplitz <- matrix(,3,a)

for (i in 1:3) {
  
  res <- numeric(a) 
  B <- toeplitz(0.5**(0:(p[i]-1)))
  
  for (j in 1:sim.size) {
  
    E <- rmvnorm(n = n, mean = rep(0,p[i]), sigma = B)
    XX <- t(VAR.sim(diag(rep(0.5, p[i])), n = n, lag = 1, include = "none", innov = E))
    res <- res + U.statistic(XX,a,tau)
    
    cat("\r","finished simulation ", paste0((j)), "for sample size ", paste0(p[i]))
    
  }
  
  res.AR1.toeplitz[i,] <- res / sim.size
  
}



## Print results to console

round(res.gaussian.errors,3)
round(res.AR1.diagonal,3)
round(res.AR1.toeplitz,3)



#------------------------------------------------------------
# 
# Simulation: U-statistic dependence under serial correlation
# 
# Investigates two settings: 
#   1. Independent Gaussian errors   
#   2. Spatially uncorrelated AR(1) errors
# 
#-------------------------------------------------------------


## Independent standard Gaussian errors

print("Independent standard Gaussian errors")

res.cov.gaussian.errors <- matrix(,sim.size,a)

for (i in 1:sim.size) {
  
  XX <- matrix(rnorm(n*p[2]), nrow = p[2], ncol = n)
  res.cov.gaussian.errors[i,] <- U.statistic(XX,a,tau)
  
  cat("\r","finished simulation ", paste0(i))
  
}


## Spatially independent AR(1) errors 


print("spatially independent AR(1) errors")

res.cov.AR1 <- matrix(,sim.size,a)

for (i in 1:sim.size) {
  
  XX <- matrix(,p[2],n)
  for (k in 1:p[2]) {
    XX[k,] <- arima.sim(model = list(ar = 0.5), n)
  }
  
  res.cov.AR1[i,] <- U.statistic(XX,a,tau)
  
  cat("\r","finished simulation ", paste0((i)))
  
}
  

## Print result to console 

round(cov(res.cov.gaussian.errors),3)
round(cov(res.cov.AR1),3)
