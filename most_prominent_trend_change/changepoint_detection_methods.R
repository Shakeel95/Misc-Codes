library("docstring") # EZ documentation
library("robts") # install from R-Forge http://robts.r-forge.r-project.org/


#-----------------------------------------------------------------------------------------------------------------------------------------
#
# Helper functions: 
#
# (1) Contrast function for continuous pcws linear signal from NOT paper
# (2) Generalised Likelihood Ratio (GLR) for kink discontinuity over regions of a signal
# (3) Calculate variance of error component using Mean Absolute Deviation (MAD) estimator for linear signals
# (4) Calculate long run varaince of error by estimating changepoint location chanel by chanel and applying kernel estimate to residuals 
#
# Uses contrast function trick from Baronowski et al. (2019)
#
#-----------------------------------------------------------------------------------------------------------------------------------------


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


vec.GLR <- function(X,b) {
  
  #' Returns a vector of GLR statistics for each chanel
  #'
  #'@param X matrix or data frame, dimension T x n 
  #'@param b int, potential break loaction
  
  
  tt <- nrow(X); n <- ncol(X)
  
  CUSUM <- apply(X,2, function(x) sum(lin.contrast(b,1,tt)*x))
  
  return(CUSUM)
  
}


MAD.pcws.lin <- function(x) {
  
  #' Medain Absolute Deviation estimate of error varaince 
  #'
  #'@param x data vector 
  
  y <- diff(diff(x))/(qnorm(3/4)*sqrt(6))
  
  return(median(abs(y)))
  
}


var.pcws.lin.residuals <- function(x) {
  
  #' Variance of residuals from best fit linear signal with kink discontinuity
  #'  
  #'@param x data vector
  
  tt <- length(x)
  CUSUM <- sapply(2:(length(x)-2), function(b) abs(sum(lin.contrast(b,1,tt)*x)))
  
  k <- which.max(CUSUM) + 2 
  tau <- sapply(1:tt - k, function(t) max(0,t))
  
  res <- lm(x ~ c(1:tt) + tau)$residuals
  
  return(var(res))
  
}


LRV.kernel.lin <- function(x) {
  
  #' Kernel estimate long-run varaince
  #'
  #'@param x data vector 
  
  tt <- length(x)
  CUSUM <- sapply(2:(length(x)-2), function(b) abs(sum(lin.contrast(b,1,tt)*x)))
  
  k <- which.max(CUSUM) + 2 
  tau <- sapply(1:tt - k, function(t) max(0,t))
  
  res <- lm(x ~ c(1:tt) + tau)$residuals
  
  return(max(asymvar.acf(res)$lrv, 
             acf(res, lag.max = 0, type = "covariance", plot = FALSE)$acf / 2 )
         )
  
}


studentize.lin.signal <- function(X, var.est) {
  
  #' Select a procedure for  
  #'
  #'@param var.est string, type of studentizer to calculate
  #'@param X matrix or data frame 
  
  
  if (!(var.est %in% c("MAD","res.var","kernel.LRV"))) stop("var.est must be one of 'MAD','res.var','kernel.LRV'")
  
  
  if (var.est == "res.var") return(apply(X,2,var.pcws.lin.residuals))
  
  if (var.est == "MAD") return(apply(X,2,MAD.pcws.lin))
  
  if (var.est == "kernel.LRV") return(apply(X,2,LRV.kernel.lin))

}


quick.reg <- function(y) {
  
  #' Quickly estimate the slope parameter in a trending regression
  #' 
  #' @param y response vector
  
  
  x <- 1:length(y) 
  
  mu.x <- mean(x)
  mu.y <- mean(y)
  
  slope <- sum((x-mu.x)*(y-mu.y)) / sum((x - mu.x)**2)
  intercept <- mu.y - slope * mu.x
  
  return(c(intercept, slope))
}


#-------------------------------------------------------------------
# 
# Aggreagtion inspired by Enikeeva and Harchaoui (2014)
# "High-dimensional change-point detection with sparse alternatives"
#
# 
#-------------------------------------------------------------------


EH.lin.GLR <- function(CUSUM, X, var.est) {
  
  #' L_2 aggregation of GLR statistics based on Normal likelihood + Chi-square distribution 
  #'
  #'@param CUSUM vector, output from vec.GLR function
  #'@param X matrix or data frame
  #'@param var.est vector, estimated variance / LRV to studentize 
  
  
  CUSUM <- (CUSUM**2) / studentize.lin.signal(X, var.est)
  
  n <- length(CUSUM)
  
  return((sum(CUSUM) - n) / sqrt(2*n))
  
}



EH.scan.GLR <- function(CUSUM, X, var.est) {
  
  #' Scan statistic aggregation of GLR based on marginal likelihood argument 
  #'
  #'@param CUSUM vector, output from vec.GLR function
  #'@param X matrix or data frame
  #'@param var.est vector, estimated variance / LRV to studentize 
  
  
  CUSUM <- (CUSUM**2) / studentize.lin.signal(X, var.est)
  
  CUSUM <- sort(CUSUM, decreasing = TRUE)
  
  CUSUM.scan <- sapply(1:length(CUSUM), function(i) (sum(CUSUM[1:i]) - i) / sqrt(2*i))
  
  k <- which.max(CUSUM.scan)
  
  return(CUSUM.scan[k])
  
}


#--------------------------------------------------------------------------------
# 
# L2 aggegation of GLR statistics + self normalization to avoid tuning parameters
# 
# motivated by: 
#     * Enikeeva and Harchaoui (2014) for likelihood approach
#     * Shao (2015) for time series self normalisation
#
#--------------------------------------------------------------------------------


GLR.self.normalizing <- function(CUSUM, X, b) {
  
  #' Self-normalizing L2 aggregation of univariate GLRs
  #'
  #'@CUSUM vector, output from vec.GLR function
  #'@param X matrix or data frame, dimension T x n 
  #'@param b int, potential break loaction
  
  
  tt <- nrow(X)
  n <- ncol(X)
  U <- (sum(CUSUM**2) - n) / (sqrt(2*n))
  
  W.1 <- 0
  for (k in 3:(b-3)) W.1 <- W.1 + ((sum(vec.GLR(X[1:b,],k)**2) - b) / sqrt(2*b))**2
  
  W.2 <- 0
  for (k in 3:(tt-b-3)) W.2 <- W.2 + ((sum(vec.GLR(X[(b+1):tt,],k)**2) - (tt-b)) / sqrt(2*(tt-b)))**2
  
  W <- (1/tt)*(W.1 + W.2)
  
  
  return((U**2)/W)
  
}


#-------------------------------------------------------
# 
# L infinity (maximum) aggregation of GLR statistics
# 
#   * Motivated by Jirak (2015)
#   * Good for sparse changes
#
#--------------------------------------------------------


GLR.pointwise.max <- function(CUSUM, X, var.est) {
  
  #' Pointwise maximum of squared contrast functions
  #'
  #'@param CUSUM vector, output from vec.GLR function
  #'@param X matrix or data frame
  #'@param var.est vector, estimated variance / LRV to studentize 
  
  CUSUM <- (CUSUM**2) / studentize.lin.signal(X, var.est)
  
  return(max(abs(CUSUM)))
  
}



