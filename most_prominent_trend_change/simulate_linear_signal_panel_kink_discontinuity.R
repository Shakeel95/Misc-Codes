library("docstring") # EZ documentation
library("mvtnorm") # multivariate normal
library("tsDyn") # simulate VAR


#------------------------------------------
#
# Functions for various covariance matrices
# 
#------------------------------------------


banded.independence <- function(n, val.1, val.2, val.3){
  
  #'@param n int, matrix size
  #'@param val.1 float, value of diagonals (positive)
  #'@param val.2 float, value of first off diagonal
  #'@param val.3 float, value of second off diagonal 
  
  M <- diag(val.1, n)
  
  M[row(M) == (col(M)-1)] <- val.2
  M[row(M) == (col(M)+1)] <- val.2
  
  M[row(M) == (col(M)-2)] <- val.3
  M[row(M) == (col(M)+2)] <- val.3
  
  return(M)
}


AR.1 <- function(n,val){
  
  #'@param n int, matrix size
  #'@param val float, toeplitz mult. in (0,1)
  
  if (abs(val) >= 1) stop("val should be in (-1,1)")
  
  S <- toeplitz(val**(1:n))
  
  return(S)
}


AR.1.blocks <- function(n, val, block.gen) {
  
  #'@param n int, matrix size
  #'@param val float, toeplitz mult. in (0,1)
  #'@param block.gen function, returns block size from n
  
  m <- round(block.gen(n))
  
  block.list <-seq(from = 1, to = n, by = m)
  block.list <- block.list[-length(block.list)]
  
  M <- diag(n)
  
  for (block in block.list) {
    
    M[block:min(block+m-1,n),block:min(block+m-1,n)] <- AR.1(min(m,n-block), val)
    
  }
  
  return(M)
  
}


#-------------------------------------------------------------
# 
# Function for simulating piecewise linear panel
# containing is single break in s linear signal without jumps
#
#-------------------------------------------------------------



sim.panel <- function(tt,tau,cpt.density,S, noise.type, slope.range) {
  
  #' Simulate pieicewise linear panel with kink discontinuity
  #'
  #'@param tt int, number of time points
  #'@param tau int, location of changepoint 
  #'@param cpt.density float, proportion of time series  
  #'@param S matrix, variance-covariance matrix for error term
  #'@param noise.type string, should be one of 'gaussian', 't.3', 'AR.ind','AR.factor'
  #'@param slope,range vector, range of values slope parameters may take
  
  if (ncol(S) != nrow(S)) stop("S should be a square matrix")
  
  if (slope.range[1] > slope.range[2]) stop("slope.range = c(a,b) should satisfy a < b")
  
  if ((cpt.density > 1) | (cpt.density < 0)) stop("changepoint density should be in (0,1)")
  
  if (!(noise.type %in% c("gaussian","t.3","AR.ind","AR.factor"))) stop("noise.type should be one of 'gaussian', 't.3', 'AR.ind','AR.factor'")
  
  
  n <- ncol(S)
  cpt.series <- sample(1:n,floor(cpt.density*n))
  
  
  X <- matrix(0,tt,n)
  
  for (i in 1:n) {
    
    if (i %in% cpt.series) {
      
      alpha <- runif(1, min = -1, max = 1)
      beta.1 <- runif(1,min = slope.range[1], max = slope.range[2])
      beta.2 <- runif(1,min = slope.range[1], max = slope.range[2])
      
      X[,i] <- c(
        alpha + beta.1*((1:tau) - tau),
        alpha + beta.2*(((tau+1):tt) - tau)
      )
      
    } else {
      
      alpha <- runif(1, min = -1, max = -.5)
      beta.1 <- runif(1, min = -.5, max = .5)
      
      X[,i] <- alpha + beta.1*((1:tt) - tau)
      
    }
  }
  
  
  if (noise.type == "gaussian") E <- rmvnorm(tt, mean = rep(0,n), sigma = S)
  
  if (noise.type == "t.3") E <- rmvt(tt, sigma = S, df = 3)
  
  if (noise.type == "AR.ind") {
    
    E <- matrix(0,tt,n)
    
    for (i in 1:n) E[,i] <- arima.sim(model = list(ar = 0.5), n = tt)
  }
  
  if (noise.type == "AR.factor") E <- VAR.sim(diag(rep(0.5,n)), 
                                              n = tt, 
                                              lag = 1, 
                                              include = "none",
                                              innov = rmvnorm(tt, mean = rep(0,n), sigma = S)
  )
  
  return(X + E)
}




