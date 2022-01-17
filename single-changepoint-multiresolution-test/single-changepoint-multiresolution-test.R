

#----------------
#
# Functions 
#
#----------------

require(lpSolveAPI)

minimum_energy_signal <- function(n, rho, cpt_locs) {
  
  #' 
  
  #'@param n int, sample size
  #'@param rho int, numerical constant on signal strength
  #'@param cpt_locs vector, potential cpt locs to sample from
  
  e <- rnorm(n)
  
  tau <- sample(cpt_locs, 1)
  
  jump <- (-1)**rbinom(1,1,0.5) * rho * sqrt(2*log(log(n)) * (n/(tau*(n-tau))))
  
  e + c(rep(0, tau), rep(jump, n-tau))
  
}



multiresolution_endpoints_test <- function(x, alpha = 0.05) {
  
  #' Single chaangepoint test based of multi-resolution fit over all intervals containing one or both endpoints [1,n]
  #'
  #'@param x vector, signal to test 
  
  
  n <- length(x)
  
  cumsum_x <- cumsum(x)
  
  thresh <- sqrt(2*log(log(n))) + (0.5*log(log(log(n))) - log(2*pi) + log(1/log(1/sqrt(1-alpha/2)))) / sqrt(2*log(log(n)))
  
  
  lps.model <- make.lp(nrow = 4*(n-1), ncol = 3)
  
  set.constr.type(lps.model, rep(">=", 4*(n-1)))
  
  set.objfn(lps.model, c(1,0,0))
  
  
  f.rhs <- c(); f.const <- c()
  
  for (i in 1:(n-1)) {
    
    f.rhs <- c(f.rhs, cumsum_x[i]/sqrt(i), (cumsum_x[n]-cumsum_x[i])/sqrt(n-i))
    
    f.const <- c(f.const, sqrt(i), sqrt(n-i))
    
  }
  
  
  f.const <- c(f.const, -f.const)
  
  set.column(lps.model, 1, rep(1, 4*(n-1)))
  
  set.column(lps.model, 2, f.const)
  
  set.column(lps.model, 3, -f.const)
  
  set.rhs(lps.model, c(f.rhs, -f.rhs))
  
  
  solve(lps.model)

  ifelse(get.variables(lps.model)[1] > thresh, 1, 0)
  
}


cusum_test <- function(x, alpha = 0.05) {
  
  #'
  #'
  #'@param x vector, signal to test 
  
  
  cumsum_x <- cumsum(x)
  
  n <- length(x)
  
  thresh <- sqrt(2*log(log(n))) + (0.5*log(log(log(n))) - log(gamma(0.5)) + log(1/log(1/sqrt(1-alpha)))) / sqrt(2*log(log(n)))
  
  res <- c()
  
  for (i in 1:n-1) res <- c(res, sqrt((n-i)/(n*i))*cumsum_x[i] - sqrt(i/(n*(n-i)))*(cumsum_x[n] - cumsum_x[i]))
  
  ifelse(max(abs(res)) > thresh, 1, 0)
  
}


make_power_curves <- function(n, K, rho_seq, trim = 5, alpha = 0.05, prog_bar = TRUE) {
  
  #'
  #'
  #'@param n int, sample size to use
  #'@param K int, no. Monet Carlo runs for 
  #'@param rho_seq vector, sequence of constant to test for minimax energy
  #'@param trim int, don't sample changepoint locations from endpoints
  #'@param prog_bar bool, wheather to show progress bar 
  
  
  nsp_power <- c()
  
  cusum_power <- c()
  
  m <- length(rho_seq)
  
  if (prog_bar) pb <- txtProgressBar(min = 1, max = m, style = 3)
  
  for (i in 1:m) {
    
    nsp_temp <- c()
    cusum_temp <- c()
    
    for (j in 1:K) {
      
      y <- minimum_energy_signal(n, rho_seq[i], trim:(n-trim))
      
      nsp_temp <- c(nsp_temp, multiresolution_endpoints_test(y, alpha))
      cusum_temp <- c(cusum_temp, cusum_test(y, alpha))
      
    }
    
    nsp_power <- c(nsp_power, mean(nsp_temp))
    
    cusum_power <- c(cusum_power, mean(cusum_temp))
    
    if (prog_bar) setTxtProgressBar(pb, i)
    
  }
  
  
  return(list(rho_seq = rho_seq, nsp_power = nsp_power, cusum_power = cusum_power))
  
}


#--------------------------
#
# Power curve simulations
#
#--------------------------

rho_seq <- seq(from = 1, to = 4, length.out = 40)

## simulate

pow.c.100 <- make_power_curves(100, 500, rho_seq)

pow.c.500 <- make_power_curves(500, 500, rho_seq)

pow.c.1000 <- make_power_curves(1000, 500, rho_seq)


#--------------------
#
# Power curve plots
#
#--------------------


plot(pow.c.100$rho_seq, pow.c.100$cusum_power, 
     type = "l", 
     ylim = c(0,1), 
     xlab = "1+rho", 
     ylab = "power", 
     col = "blue")

lines(pow.c.100$rho_seq, pow.c.100$nsp_power, type = "l", col = "red")

lines(pow.c.500$rho_seq, pow.c.500$cusum_power, type = "l", lty = 2, col = "blue")

lines(pow.c.500$rho_seq, pow.c.500$nsp_power, type = "l", lty = 2, col = "red")

lines(pow.c.1000$rho_seq, pow.c.1000$cusum_power, type = "l", lty = 3, col = "blue")

lines(pow.c.1000$rho_seq, pow.c.1000$nsp_power, type = "l", lty = 3 , col = "red")

