library("docstring") # EZ documentation
library("mvtnorm") # multivariate normal


#---------------------------
#
# GLR and DIM-SUM functions
#
#---------------------------


GLR <- function(X) {
  
  #' Generalised Likelihood Ratio test
  #' (for mean shift in Guassian data with diagonal covaraince matrix)
  #' 
  #' @param X matrix, dim n x T
  
  
  CUSUM.sqrd <- function(x, tau) {
    
    #' Squared univariate CUSUM statistic
    
    tt <- length(x)
    
    ux.R <- mean(x[1:tau])
    ux.L <- mean(x[(tau+1):tt])
    
    return(((tau*(tt-tau))/tt)*(ux.L - ux.R)**2)
  }
  
  
  n <- dim(X)[1]
  tt <- dim(X)[2]
  
  res <- numeric(tt-2)
  
  for (i in 1:n) {
    
    res <- res + sapply(2:(tt-1), function(tau) CUSUM.sqrd(X[i,], tau))
    
  }
  
  return(list(res = res, cpt = which.max(res)+1))
}


DIM.SUM <- function(X) {
  
  #' Difference-In-Means scaled cuSUMs
  #' 
  #' @param X matrix, dim n x T 
  
  
  dim.sum <- function(x,tau) {
    
    #' Univariate dim-sum statistic
    
    tt <- length(x)
    
    ux.R <- mean(x[1:tau])
    ux.L <- mean(x[(tau+1):tt])
    
    return(sqrt((tau*(tt-tau))/sqrt(tt))*(ux.L - ux.R)**2)
    
  }
  
  n <- dim(X)[1]
  tt <- dim(X)[2]
  
  res <- numeric(tt-2)
  
  for (i in 1:n) {
    
    res <- res + sapply(2:(tt-1), function(tau) dim.sum(X[i,],tau))
    
  }
  
  return(list(res = res, cpt = which.max(res)+1))
  
}


scaled.DIM.SUM <- function(X) {
  
  #' Difference-in-Means scaled cuSUMs
  #' (scaled with inner produce of vector of mean changes)
  #' 
  #' @param X matrix, dim n x T
  
  n <- dim(X)[1]
  tt <- dim(X)[2]
  
  res <- numeric(tt-2)
  
  for (tau in 2:(tt-1)) {
    
    est.proj <- apply(X[,1:tau],1,mean) - apply(X[,(tau):tt-1],1,mean)
    
    res[tau-2] <- sum(est.proj ** 2) * sqrt((tau*(tt-tau))/(tt)) * (norm(est.proj, "2") ** (-1))
    
  }
  
  return(list(res = res, cpt = which.max(res)+2))
   
}


abs.sum.CUSUM <- function(X) {
  
  #' Sum of absolute values of CUSUMs
  #' 
  #' @param X matrix, dim n x T
  
  
  CUSUM.abs <- function(x, tau) {
    
    #' absolute value of univariate CUSUM statistic
    
    tt <- length(x)
    
    ux.R <- mean(x[1:tau])
    ux.L <- mean(x[(tau+1):tt])
    
    return(sqrt((tau*(tt-tau))/tt)*abs(ux.L - ux.R))
  }
  
  
  n <- dim(X)[1]
  tt <- dim(X)[2]
  
  res <- numeric(tt-2)
  
  for (i in 1:n) {
    
    res <- res + sapply(2:(tt-1), function(tau) CUSUM.abs(X[i,], tau))
    
  }
  
  return(list(res = res, cpt = which.max(res)+1))
  
}


#----------------------------------------------------------------
# 
# Simulation study 
# 
# DIM-SUM statistic when the changepoint is close to the bondary. 
# 
#-----------------------------------------------------------------


set.seed(42)

sim.size <- 100
cpt.bal <- 50
cpt.ubal <- 10 


## Balanced chnagepoint simulation

res.DIM.SUM <- c(); res.scaled.DIM.SUM <- c(); res.GLR <- c()

for (k in 1:sim.size) {
  
  Y <- rmvnorm(10, mean = rep(0,100), sigma = diag(100))
  Y[,cpt.bal:100] <- Y[,cpt.bal:100] + 0.5
  
  res.GLR <- c(res.GLR, GLR(Y)$cpt)
  res.DIM.SUM <- c(res.DIM.SUM, DIM.SUM(Y)$cpt)
  res.scaled.DIM.SUM <- c(res.scaled.DIM.SUM, scaled.DIM.SUM(Y)$cpt)
  
  cat("\r","Finished simulation : ", paste0(k))
  
}


res.bal <- data.frame(cbind(res.DIM.SUM, res.scaled.DIM.SUM, res.GLR))

apply(res.bal,2,mean)
apply(res.bal, 2, median)


## Unbalanced changepoint simulation

res.DIM.SUM <- c(); res.scaled.DIM.SUM <- c(); res.GLR <- c()

for (k in 1:sim.size) {
  
  Y <- rmvnorm(10, mean = rep(0,100), sigma = diag(100))
  Y[,cpt.ubal:100] <- Y[,cpt.ubal:100] + 0.5
  
  res.GLR <- c(res.GLR, GLR(Y)$cpt)
  res.DIM.SUM <- c(res.DIM.SUM, DIM.SUM(Y)$cpt)
  res.scaled.DIM.SUM <- c(res.scaled.DIM.SUM, scaled.DIM.SUM(Y)$cpt)
  
  cat("\r","Finished simulation : ", paste0(k))
  
}


res.ubal <- data.frame(cbind(res.DIM.SUM, res.scaled.DIM.SUM, res.GLR))

apply(res.ubal,2,mean)
apply(res.ubal, 2, median)


#------------------------------------------------------------------
#
# Plots: DIM-SUM v. CUSUM with balanced changepoint location
#
#-------------------------------------------------------------------


## DIM-SUM statistic is more peaked around changepoint location than CUSUM. 


set.seed(42)

Y <- rmvnorm(100, mean = rep(0,100), sigma = diag(100))
Y[,50:100] <- Y[,50:100] + 0.5


# DIM-SUM statistic

plot(DIM.SUM(Y)$res,
     type = "l", 
     main = "DIM-SUM", 
     xlab = "time", 
     ylab = "value")

abline(v = 50, col = "red", lty = 2)
abline(v = which.max(DIM.SUM(Y)$res) + 1, col = "blue", lty = 2)


# Summed absolute CUSUMS statistic

plot(abs.sum.CUSUM(Y)$res,
     type = "l", 
     main = "summed CUSUMs", 
     xlab = "time", 
     ylab = "value")

abline(v = 50, col = "red", lty = 2)
abline(v = which.max(DIM.SUM(Y)$res) + 1, col = "blue", lty = 2)


#------------------------------------------------------------------------------------
# 
# Plots: show DIM-SUM v. GLR v. scaled DIM-SUM  with unbalanced changepoint location
# 
#-------------------------------------------------------------------------------------


## DIM-SUM statistic breaks down when the change occurs at the boundary


set.seed(42)

Y <- rmvnorm(100, mean = rep(0,100), sigma = diag(100))
Y[,5:100] <- Y[,5:100] + 0.5


# DIM-SUM statistic

plot(DIM.SUM(Y)$res,
     type = "l", 
     main = "DIM-SUM", 
     xlab = "time", 
     ylab = "value")

abline(v = 10, col = "red", lty = 2)
abline(v = which.max(DIM.SUM(Y)$res) + 1, col = "blue", lty = 2)


# GLR statistic

plot(GLR(Y)$res, 
     type = "l", 
     main = "GLR", 
     xlab = "time", 
     ylab = "value")

abline(v = 10, col = "red", lty = 2)
abline(v = which.max(GLR(Y)$res), col = "blue", lty = 2)


# (scaled) DIM-SUM statistic

plot(scaled.DIM.SUM(Y)$res,
     type = "l", 
     main = "(scaled) DIM-SUM", 
     xlab = "time", 
     ylab = "value")

abline(v = 10, col = "red", lty = 2)
abline(v = which.max(GLR(Y)$res), col = "blue", lty = 2)
