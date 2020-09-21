library("docstring") # EZ documentation
library("foreach") # nicer loops



#-------------------
# 
# Contrast function
# 
#-------------------


lin.contrast <- function(b, s, e){
  
  #' pcwsLin contrast function for GLR (see NOT paper)
  #' 
  #' 
  #' @param b int, potential changepoint loaction 
  #' @param s int, start
  #' @param e int, end 
  
  l <- e-s
  A <- sqrt(6/(l*(l**2-1)*(1+(e-b+1)*(b-s)+(e-b)*(b-s-1))))
  B <- sqrt(((e-b+1)*(e-b))/((b-s-1)*(b-s)))
  
  contrast <- numeric(l)
  
  contrast[(s+1):b] <- (A*B)*((3*(b-s)+(e-b)-1)*c((s+1):b) - (b*(e-s-1)+2*(s+1)*(b-s)))
  
  contrast[(b+1):e] <- (-1)*(A/B)*((3*(e-b)+(b-s)+1)*c((b+1):e) - (b*(e-s-1)+2*e*(e-b+1)))
  
  return(contrast)
  
}


MAD.pcws.lin <- function(x) {
  
  #' Medain Absolute Deviation estimate of error varaince 
  #'
  #'@param x data vector 
  
  y <- diff(diff(x))/(qnorm(3/4)*sqrt(6))
  
  return(median(abs(y)))
  
}


#------------------------
# 
# GLR aggregation methods
#
#------------------------


get.cpt.GLR <- function(X,S) {
  
  #' Get most likely changepoint location using oracel weighting and naive rescaling
  #'
  #'@param X matrix or data frame
  #'@param S symmetric covaraince matrix for X
  
  tt <- nrow(X)
  n <- ncol(X)
  
  
  CUSUM <- foreach(b = 2:(tt-2), .combine = "cbind") %do% {
    
    apply(X,2, function(x) sum(lin.contrast(b,1,tt)*x))
    
  }
  
  S.diag <- apply(X,2,MAD.pcws.lin)
  
  
  vals.naive <- naive.GLR(CUSUM, S.diag) 
  
  vals.oracle <- oracle.GLR(CUSUM, S)
  
  return(list(
    naive.cpt = which.max(vals.naive) + 2, 
    naive.vals = vals.naive, 
    oracle.cpt = which.max(vals.oracle) + 2, 
    oracle.vals = vals.oracle
  ))
  
}


naive.GLR <- function(CUSUM,S.diag) {
  
  #' aggregate contract functions using naive scaling 
  #'
  #'@param CUSUSM matrix of innter products with linear contrast fu
  #'@param S.diag vector, estimated studentizers
  
  
  if (min(diag(S.diag)) < 0) stop("all series must vary")
  
  n <- nrow(CUSUM)
  
  vals <- sapply(1:ncol(CUSUM), function(i) (sum((CUSUM[,i]/S.diag)**2) - n) / sqrt(2*n) )
  
  return(vals)
  
}


oracle.GLR <- function(CUSUM, S) {
  
  #' Aggregate contrast function using oracle (cov-matrix) weights
  #'
  #'@param CUSUSM matrix of innter products with linear contrast fu
  #'@param S true covaraince matrix for data sequence
  
  
  e <- eigen(S)$values  
  
  if (min(e) < 0) stop("all eigenvalues in S must be positive")
  
  e.1 <- sum(e)
  e.2 <- sum(e**2)
  
  vals <- sapply(1:ncol(CUSUM), 
                 function(i) (sum(CUSUM[,i]**2) - e.1) / sqrt(2*e.2)
                 )
  
  return(vals)
  
}
