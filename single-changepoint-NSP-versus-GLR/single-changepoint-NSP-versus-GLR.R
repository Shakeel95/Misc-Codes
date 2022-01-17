

#-------------------
#
# Functions
#
#-------------------



require(lpSolveAPI)



minimum_energy_signal <- function(n, rho, cpt_loc) {
  
  #' Simulate signal whith changepoint having energy above phase transition 
  #'
  #'@param n int, sample size
  #'@param rho int, numerical constant on signal strength
  #'@param cpt_loc int, changepoint location
  
  e <- rnorm(n)
  
  
  jump <- (-1)**rbinom(1,1,0.5) * (1+rho) * sqrt(2*log(log(n)) * (n/(cpt_loc*(n-cpt_loc))))
  
  e + c(rep(0, cpt_loc), rep(jump, cpt_loc))
  
}


cusum_test <- function(x, alpha = 0.05) {
  
  #' Single changepoint test based on GLR statistic
  #'
  #'@param x vector, signal to test 
  #'@param alpha float, Type I error to control
  
  
  cumsum_x <- cumsum(x)
  
  n <- length(x)
  
  thresh <- sqrt(2*log(log(n))) + (0.5*log(log(log(n))) - log(gamma(0.5)) + log(1/log(1/sqrt(1-alpha)))) / sqrt(2*log(log(n)))
  
  res <- c()
  
  for (i in 1:n-1) res <- c(res, sqrt((n-i)/(n*i))*cumsum_x[i] - sqrt(i/(n*(n-i)))*(cumsum_x[n] - cumsum_x[i]))
  
  ifelse(max(abs(res)) > thresh, 1, 0)
  
}


multiresolution_test <- function(x, alpha = 0.05) {
  
  #' Test a signal for a single pcws mean shift using residuals from multiresolution fit
  #'
  #'@param x vector, signal to test
  #'@param alpha float, Type I error to control

  
  n <- length(x)
  
  thresh <- sqrt(2*log(n)) + (0.5*log(log(n)) + log(0.8/(2*sqrt(pi))) + log(1/log(1/sqrt(1-alpha))))
  
  
  X <- rbind(c(0,0), cbind(rep(1,n), x))
  
  S <- apply(X, 2, cumsum)
  
  n_scans <- floor(log2(n))
  
  
  f_const <- c(); f_rhs <- c()
  
  for (j in 1:(n_scans+1)) {
    
    l <- 2**(j-1)
    
    for (k in 1:(n-l+1)) {
      
      f_const <- c(f_const, (S[k+l,1] - S[k,1]) / sqrt(l))
      
      f_rhs <- c(f_rhs, (S[k+l,2] - S[k,2]) / sqrt(l))
      
    }
  }
  
  
  m <-  2*length(f_rhs)
  
  lps_model <- make.lp(nrow = m, ncol = 3)
  
  set.objfn(lps_model, c(1,0,0))
  
  set.constr.type(lps_model, rep(">=",m))
  
  set.column(lps_model, 1, rep(1,m))
  
  set.column(lps_model, 2, c(f_const, -f_const))
  
  set.column(lps_model, 3, c(-f_const, f_const))
  
  set.rhs(lps_model, c(f_rhs, -f_rhs))
  
  
  solve(lps_model)
  
  ifelse(get.variables(lps_model)[1] > thresh, 1, 0)
  
  
}



simulate_power_curves_rho_seq <- function(n, rho_seq, monte_carlo_reps, cpt_frac = 0.5, alpha = 0.05, seed = 42) {
  
  #' Simulate power curves for different numerical constants on pahse transition boundary for fixed n
  #'
  #'@param n int, signal length
  #'@param rho_seq vector, numerical constants on signal strength
  #'@param monte_carlo_reps int, no. MC reps for calculatin power
  #'@param cpt_frac float, cpt_loc = \floor{n * cpt_frac}
  #'@param alpha float, Tpe I error control
  #'@param seed int
  
  
  set.seed(seed)
  
  K <- length(rho_seq)
  
  pb <- txtProgressBar(min = 1, max = K, style = 3)
  
  multiresolution_power_curve <- numeric(K)
  
  cusum_power_curve <- numeric(K)
  
  
  for (k in 1:K) {
    
    cusum_res <- numeric(monte_carlo_reps)
    
    multiresoluion_res <- numeric(monte_carlo_reps)
    
    for (rep in 1:monte_carlo_reps) {
      
      x <- minimum_energy_signal(n, rho_seq[k], floor(n*cpt_frac))
      
      multiresoluion_res[rep] <- multiresolution_test(x, alpha)
      
      cusum_res[rep] <- cusum_test(x, alpha)
      
    }
    
    setTxtProgressBar(pb, k)
    
    multiresolution_power_curve[k] <- mean(multiresoluion_res)
    
    cusum_power_curve[k] <- mean(cusum_res)
    
  }
  
  return(list(multiresolution_power_curve = multiresolution_power_curve, 
              cusum_power_curve = cusum_power_curve))
  
}



simulate_power_curve_n_seq <- function(rho, n_seq, monte_carlo_reps, cpt_frac = 0.5, alpha = 0.05, seed = 42) {
  
  #' Simulate power curves for different signal lengths given fixed numerical constants on pahse transition boundary 
  #'
  #'@param rho float, numerical constant on energy
  #'@param n_seq vector, signal lengths to test 
  #'@param monte_carlo_reps int, no. MC reps for calculatin power
  #'@param cpt_frac float, cpt_loc = \floor{n * cpt_frac}
  #'@param alpha float, Tpe I error control
  #'@param seed int
  

  set.seed(seed)
  
  K <- length(n_seq)
  
  pb <- txtProgressBar(min = 1, max = K, style = 3)
  
  multiresolution_power_curve <- numeric(K)
  
  cusum_power_curve <- numeric(K)
  
  
  for (k in 1:K) {
    
    cusum_res <- numeric(monte_carlo_reps)
    
    multiresoluion_res <- numeric(monte_carlo_reps)
    
    for (rep in 1:monte_carlo_reps) {
      
      x <- minimum_energy_signal(n_seq[k], rho, floor(n_seq[k]*cpt_frac))
      
      multiresoluion_res[rep] <- multiresolution_test(x, alpha)
      
      cusum_res[rep] <- cusum_test(x, alpha)
      
    }
    
    setTxtProgressBar(pb, k)
    
    multiresolution_power_curve[k] <- mean(multiresoluion_res)
    
    cusum_power_curve[k] <- mean(cusum_res)
    
  }
  
  return(list(multiresolution_power_curve = multiresolution_power_curve, 
              cusum_power_curve = cusum_power_curve))
  
  }



#-------------------------
#
# Simulations and plots
#
#-------------------------


## simulation params

rho_seq <- c(0:4, seq(from = 4.1, to = 5.9, by = 0.1), 6:10)

n_seq <- seq(from = 50, to = 500, by = 12)

monte_carlo_reps <- 100


## simulations (rho seq)

power_n_100 <- get_power_curves_wrapper(100,rho_seq, monte_carlo_reps)

power_n_500 <- get_power_curves_wrapper(500,rho_seq, monte_carlo_reps)



## plots (rho seq)

plot(rho_seq, power_n_100$multiresolution_power_curve,
     type = "l", 
     ylim = c(0,1),
     xlab = "1 + rho",
     ylab = "power",
     col = "red")

lines(rho_seq, power_n_100$cusum_power_curve, type = "l", col = "blue")

lines(rho_seq, power_n_500$multiresolution_power_curve, type = "l", lty = 2, col = "red")

lines(rho_seq, power_n_500$cusum_power_curve, type = "l", lty = 2, col = "blue")


## simulations (n seq)

power_rho_1 <- simulate_power_curve_n_seq(1, n_seq, monte_carlo_reps)

power_rho_2 <- simulate_power_curve_n_seq(2, n_seq, monte_carlo_reps)

power-rho_3 <- simulate_power_curve_n_seq(3, n_seq, monte_carlo_reps) 


## plots (n seq)

plot(n_seq, power_rho_1$multiresolution_power_curve,
     type = "l", 
     ylim = c(0,1),
     xlab = "n",
     ylab = "power",
     col = "red")

lines(n_seq, power_rho_1$cusum_power_curve, type = "l", col = "blue")

lines(n_seq, power_rho_2$multiresolution_power_curve, type = "l", col = "red", lty = 2)

lines(n_seq, power_rho_2$cusum_power_curve, type = "l", col = "blue", lty = 2)

lines(n_seq, power_rho_3$multiresolution_power_curve, type = "l", col = "red", lty = 3)

lines(n_seq, power_rho_3$cusum_power_curve, type = "l", col = "blue", lty = 3)
























