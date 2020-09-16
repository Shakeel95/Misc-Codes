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
  
  
  if (!(var.est %in% c("MAD","res.var","kernel.LRV","no.student"))) stop("var.est must be one of 'MAD','res.var','kernel.LRV','no.student'")
  
  
  if (var.est == "res.var") return(apply(X,2,var.pcws.lin.residuals))
  
  if (var.est == "MAD") return(apply(X,2,MAD.pcws.lin))
  
  if (var.est == "kernel.LRV") return(apply(X,2,LRV.kernel.lin))
  
  if (var.est == "no.student") return(rep(1,ncol(X)))  
  
}


is.string <- function(x) {
  
  #' checks if x is string
  
  is.character(x) && length(x) == 1
  
}

#------------------------
# 
# GLR aggregation methods
#
#------------------------


GLR.l2 <- function(CUSUM, X, var.est) {
  
  #' L_2 aggregation of GLR statistics based on Normal likelihood + Chi-square distribution 
  #'
  #'@param CUSUM vector, output from vec.GLR function
  #'@param X matrix or data frame
  #'@param var.est string, type of studentizer to calculate

  
  if (is.string(var.est)) {
    
    CUSUM <- (CUSUM**2) / studentize.lin.signal(X, var.est)
    
  } else {
    
    CUSUM <- (CUSUM**2) / var.est
    
  }

  
  n <- length(CUSUM)
  
  return((sum(CUSUM) - n) / sqrt(2*n))
  
}


GLR.scan <- function(CUSUM, X, var.est) {
  
  #' Scan statistic aggregation of GLR based on marginal likelihood argument 
  #'
  #'@param CUSUM vector, output from vec.GLR function
  #'@param X matrix or data frame
  #'@param var.est string, type of studentizer to calculate
  
  
  if (is.string(var.est)) {
   
    CUSUM <- (CUSUM**2) / studentize.lin.signal(X, var.est) 
    
  } else {
    
    CUSUM <- (CUSUM**2) / var.est
    
  }
  
  
  CUSUM <- sort(CUSUM, decreasing = TRUE)
  
  CUSUM.scan <- sapply(1:length(CUSUM), function(i) (sum(CUSUM[1:i]) - i) / sqrt(2*i))
  
  k <- which.max(CUSUM.scan)
  
  return(CUSUM.scan[k])
  
}


GLR.max <- function(CUSUM, X, var.est, b) {
  
  #' L infinity aggregation of GLR statistics  
  #'
  #'@param CUSUM vector, output from vec.GLR function
  #'@param X matrix or data frame, dimension T x n 
  #'@param b int, potential break loaction
  
  
  if (is.string(var.est)) {
    
    CUSUM <- (CUSUM**2) / studentize.lin.signal(X, var.est)
    
  } else {
    
    CUSUM <- (CUSUM**2) / var.est
    
  }
  
  
  tt <- nrow(X)
  
  return(sqrt((b*(tt-b))/tt)*max(CUSUM))
  
}


#-------------------------------
#
# Contrast function aggregation
# 
#-------------------------------


contrasts.SN <- function(X, b) {
  
  #' Self Normalised contrast function
  #'
  #'@param X matrix or data frame of dimension T x n
  #'@param b int, potential break location 
  
  
  return(0)
  
}


#-------------------------------------------------
#
# Estimate changepoint location using all methods
#
#-------------------------------------------------


estimate.cpt.loc <- function(X, var.est, verbose = T) {
  
  #'Recover changepoint location from each of the three GLR aggreagtion methods
  #'
  #'@param X matrix or data frame
  #'@param var.est 

  
  vals.GLR.l2 <- c()
  vals.GLR.max <- c()
  vals.GLR.scan <- c()
  vals.GLR.SN <- c()
  
  
  for (b in 2:(nrow(X)-2)) {
    
    CUSUM <- vec.GLR(X,b)
    studentizer <- studentize.lin.signal(X,var.est)
    
    vals.GLR.l2 <- c(vals.GLR.l2, GLR.l2(CUSUM, X, studentizer))
    vals.GLR.max <- c(vals.GLR.max, GLR.max(CUSUM, X, studentizer, b))
    vals.GLR.scan <- c(vals.GLR.scan, GLR.scan(CUSUM, X, studentizer))
    vals.GLR.SN <- c(vals.GLR.SN, 0)
    
    if (verbose) cat("\r","Calculated statistic at time point ",paste0(b))
    
  }
  
  if (verbose) cat("\n\n")
  
  return(list(
    res.l2 = list(vals = vals.GLR.l2, cpt = which.max(vals.GLR.l2) + 1), 
    res.max = list(vals = vals.GLR.max, cpt = which.max(vals.GLR.max) + 1), 
    res.scan = list(vals = vals.GLR.scan, cpt = which.max(vals.GLR.scan) + 1), 
    res.SN = list(vals = vals.GLR.SN, cpt = which.max(vals.GLR.SN) + 1)
  ))
  
}
