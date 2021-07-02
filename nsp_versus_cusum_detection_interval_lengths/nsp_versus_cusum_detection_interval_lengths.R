
#-------------
#
# Testbed
#
#-------------


make_signal <- function(n, cpt_loc, cpt_energy, pre_change_mean = 0) {
  
  #' Simulate Gaussian signal with changepoint having prescribed energy
  #'
  #'@param n int, sample size
  #'@param cpt_loc int, changepoint location
  #'@param cpt_energy float, changepoint energy (cf. Verzen et al. 2020)
  #'@param pre_change_mean float, pre_change_mean
  
  
  post_change_mean = (-1)**rbinom(1,1,0.5) * (pre_change_mean + cpt_energy*sqrt(n/(cpt_loc*(n-cpt_loc))))
  
  rnorm(n) + c(rep(pre_change_mean, cpt_loc), rep(post_change_mean, n - cpt_loc))
  
}


#------------------------
#
# Changepoint tests
#
#------------------------


require(lpSolveAPI)


make_dyadic_scans <- function(x) {
  
  #' Take averages of signal over all dyadic windows
  #'
  #'@param x vector, signal
  
  
  n <- length(x)
  
  X <- rbind(c(0,0), cbind(rep(1,n), x))
  
  S <- apply(X, 2, cumsum)
  
  n_scans <- floor(log2(n))
  
  scan_sizes <- rep(0,1+n_scans)
  
  dyadic_scans <- array(0, c(n,2,1+n_scans))
  
  for (j in 1:(n_scans+1)) {
    
    l <- 2**(j-1)
    
    for (k in 1:(n-l+1)) dyadic_scans[k,,j] <- (S[k+l,] - S[k,]) / sqrt(l)
    
  }

  return(dyadic_scans)
  
}



multiresolution_test <- function(s, e, dyadic_scans, alpha = 0.05) {
  
  #' Test an interval for pcws mean shift using residuals from multiresolution fit
  #'
  #'@param s int, test window (from)
  #'@param e int, test window (to)
  #'@param dyadic_scans array
  #'@param alpha float, coverage probability control
  
  
  n <- dim(dyadic_scans)[1]
  
  l <- e - s + 1
  
  thresh <- sqrt(2*log(n)) + (0.5*log(log(n)) + log(0.8/(2*sqrt(pi))) + log(1/log(1/sqrt(1-alpha))))
  
  
  f_rhs <- c(); f_const <- c()
  
  J <- floor(log2(e-s+1))
  
  for (j in 0:J) {
    
    f_const <- c(f_const, dyadic_scans[s:(s+l-2**j),1,(j+1)])
    
    f_rhs <- c(f_rhs, dyadic_scans[s:(s+l-2**j),2,(j+1)])
    
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



cusum_test <- function(s, e, cumulative_sum, cusum_thresh) {

  #' Test an interval for pcws mean shift using CUSUM statistic
  #'
  #'@param s int, test window (from)
  #'@param e int, test window (to)
  #'@param cumulative_sum vector, cumulative sum over signal with zero appended to start
  #'@param cusum_thresh float, threshold obtained from `bisect_get_cusum_thresh`
  #'@param alpha float, coverage probability control  
  
  
  n <- length(cumulative_sum) - 1
  
  l <- e - s + 1
  
  
  res <- c()
  
  for (b in s:(e-1)) res <- c(res, sqrt((e-b)/(l*(b-s+1)))*(cumulative_sum[s] - cumulative_sum[b+1]) - sqrt((b-s+1)/(l*(e-b)))*(cumulative_sum[b+1] - cumulative_sum[e+1]))
  
  res <- abs(res)
  
  ifelse(max(res) > cusum_thresh, 1, 0)
  
}




#------------------------
#
# cusum test threshold
#
#------------------------


bisect_get_cusum_thresh <- function(n,thresh_min, thresh_max, desired_prob = 0.05, tollerance = 0.0005) {
  
  #' Get threshold for multiscale cusum test via bisection method
  #'
  #'@param n int, sample size being considered
  #'@param desired_prob float, covergage probability to attain
  #'@param thresh_min float 
  #'@param thresh_max float 
  #'@param tollerance float
  #'
  #'@references see theorem 1.2 in https://arxiv.org/abs/1608.03032
  
  
  thresh_mean <- (thresh_min + thresh_max) / 2
  
  error <- desired_prob - cusum_tail(n, thresh_mean)
  
  
  if (abs(error) <= tollerance) {
    
    return(thresh_mean)
    
  } else if (error > 0) {
    
    bisect_get_cusum_thresh(n, thresh_min, thresh_mean, desired_prob, tollerance)
    
  } else if (error < 0) {
    
    bisect_get_cusum_thresh(n, thresh_mean, thresh_max, desired_prob, tollerance)
    
  }
  
}


cusum_tail <- function(m,b) {
  
  #' Approximate high excursion probability of multiscale cusum test over pure noise
  #'
  #'@references see theorem 1.2 in https://arxiv.org/abs/1608.03032
  
  tail_terms <- combn(1:m, 2, function(u) ifelse(sum(u) <= m, 
                                                 (m - sum(u))/(prod(u)*sum(u))*nu(b*sqrt(u[1]/(u[2]*sum(u))))*nu(b*sqrt(u[2]/(u[1]*sum(u))))*nu(b*sqrt(sum(u)/prod(u))), 
                                                 0))
  b**6 * pnorm(-b) * sum(tail_terms) / 4
  
}


nu <- function(x) {
  
  #'Approximation to Nu function
  #'
  #'@references The Statistics of Gene Mapping (pp. 112)
  
  y <- x/2
  
  (1/y) * (pnorm(y) - 0.5) / (y * pnorm(y) + dnorm(y))
  
}


#--------------------------------------------------
#
# Apply changepoint tests over deterministic grid
#
#--------------------------------------------------


apply_tests_to_signal <- function(x, cpt_loc, cusum_thresh, alpha = 0.05, show_prog = TRUE) {
  
  #' Apply cusum and multiresolution tests to all sub-intervals of signal containig the changepoint
  #'
  #'@param x vector, signal to test
  #'@param cpt_loc float, changepoint location
  #'@param cusum_thresh float, threshold obtained from `bisect_get_cusum_thresh`
  #'@param alpha float, coverage probability control
  #'@param show_prog bool, whetehr to display progress bar
  
  
  n <- length(x)
  
  effective_n <- min(cpt_loc, n - cpt_loc)
  
  
  if (effective_n < 3) return(list(multiresolution_interval = c(0,0), cusum_interval = c(0,0)))
  
  
  pb <- txtProgressBar(min = 1, max = effective_n, style = 3)
  
  dyadic_scans <- make_dyadic_scans(x)
  
  cumulative_sum <- c(0, cumsum(x))
  
  
  multiresolution_flag <- FALSE; cusum_flag <- FALSE
  
  multiresolution_interval <- c(0,0); cusum_interval <- c(0,0)
  
  for (j in 3:effective_n) {
    
    for (k in 1:(j-2)) {
      
      loc <- c(cpt_loc - k, cpt_loc + j - (k+1))
      
      if (!cusum_flag) {
        
        cusum_flag <- cusum_test(loc[1], loc[2], cumulative_sum, cusum_thresh)
        
        if (cusum_flag) cusum_interval <- loc
        
      }
      
      
      if (!multiresolution_flag) {
        
        multiresolution_flag <- multiresolution_test(loc[1], loc[2], dyadic_scans, alpha)
        
        if (multiresolution_flag) multiresolution_interval <- loc
        
      }
      
      if (show_prog) setTxtProgressBar(pb, j)
      
      if (cusum_flag && multiresolution_flag) return(list(multiresolution_interval = multiresolution_interval, cusum_interval = cusum_interval))
      
    }
    
  }
  
  
  return(list(multiresolution_interval = multiresolution_interval, cusum_interval = cusum_interval))
  
  
}



compare_interval_widths <- function(n, cpt_loc, cpt_energy, monte_carlo_reps, cusum_thresh, alpha = 0.05, seed = 42, display_prog_bar = TRUE) {
  
  #' Monte Carlo simulations to compare width of detection intervals from CUSUM test and multiresoolution tests for fixed energy
  #' 
  #' @param n int, length of signal
  #' @param M int, 
  #' @param cpt_loc int, 
  #' @param cpt_energy int
  #' @param monte_carlo_reps int, 
  #' @param cusum_thresh float
  #' @param alpha float
  #' @param seed int
  #' @param display_prog_bar bool, display progress bar
  
  
  set.seed(seed)
  
  if (display_prog_bar) pb <- txtProgressBar(min = 1, max = monte_carlo_reps, style = 3)
  
  
  cusum_res <- matrix(0, nrow = monte_carlo_reps, ncol = 3)
  
  multiresolution_res <- matrix(0, nrow = monte_carlo_reps, ncol = 3)
  
  
  for (j in 1:monte_carlo_reps) {
    
    x <- make_signal(n, cpt_loc, cpt_energy)
    
    test_intervals <- apply_tests_to_signal(x, cpt_loc, cusum_thresh, alpha, show_prog = FALSE)
    
    
    multiresolution_res[j,1:2] <- test_intervals$multiresolution_interval
    
    if (cpt_loc %in% test_intervals[[1]][1]:test_intervals[[1]][2]) multiresolution_res[j,3] <- 1 
    
    
    cusum_res[j,1:2] <- test_intervals$cusum_interval
    
    if (cpt_loc %in% test_intervals[[2]][1]:test_intervals[[2]][2]) cusum_res[j,3] <- 1
    
    
    if (display_prog_bar) setTxtProgressBar(pb, j)
    
  }
  
  
  return(list(multiresolution = multiresolution_res, cusum = cusum_res))
  
}



#---------------------
#
# Wrapper functions
#
#---------------------



run_simulation_wrapper <- function(n, cpt_loc, rho_seq, alpha, cusum_thresh, monte_carlo_reps) {
  
  #' Run simulation to compare average interval widths
  
  
  K <- length(rho_seq)
  
  cusum_widths <- numeric(K)
  
  multiresolution_widths <- numeric(K)
  
  pb <- txtProgressBar(min = 1, max = K, style = 3)
  
  
  for (k in 1:K) {
    
    cpt_energy <- (1 + rho_seq[k]) * sqrt(2*log(n))
    
    RES <- compare_interval_widths(n, cpt_loc, cpt_energy, monte_carlo_reps, cusum_thresh, alpha, display_prog_bar = FALSE)
    

    if (all(RES$multiresolution[,1] == 0)) {
      
      multiresolution_widths[k] <- NaN
      
    } else {
      
      res <- subset(RES$multiresolution, RES$multiresolution[,1] > 0)
      
      multiresolution_widths[k] <- mean(sapply(1:nrow(res), function(i) diff(res[i,1:2]) + 1))
      
    }
    
    
    if (all(RES$cusum[,1] == 0)) {
      
      cusum_widths[k] <- NaN
      
    } else {
      
      res <- subset(RES$cusum, RES$cusum[,1] > 0)
     
      cusum_widths[k] <- mean(sapply(1:nrow(res), function(i) diff(res[i,1:2]) + 1))
       
    }
    
    setTxtProgressBar(pb, k)
    
  }
  
  
  return(list(cusum_widths = cusum_widths, multiresolution_widths = multiresolution_widths))
  
}



plot_simulation_wrapper <- function(rho_seq, cusum_widths, multiresolution_widths) {
  
  #' Plot outputs from `run_simulation_wrapper`
  
  
  plot(rho_seq, multiresolution_widths,
       ylim = c(0, max(max(multiresolution_widths, na.rm = TRUE), max(cusum_widths, na.rm = TRUE))),
       type = "l",
       lty = 2,
       xlab = "rho",
       ylab = "mean interval width"
  )
  
  lines(rho_seq, cusum_widths, type = "l", lty = 2, col = "red")
  
}



#--------------------------
#
# Simulations and plots
#
#--------------------------


## simulation params


rho_seq <- seq(from = 0, to = 6, by = 0.1)
K <- length(rho_seq)
M <- 10**3
monte_carlo_reps <- 100
alpha  <- 0.05


## n = 100

n <- 100; cpt_loc <- n/2

cusum_thresh_100 <- bisect_get_cusum_thresh(n, 4, 10, alpha) # 4.070312


res_n_100 <- run_simulation_wrapper(n, cpt_loc, rho_seq, alpha, cusum_thresh_100, monte_carlo_reps)

plot_simulation_wrapper(rho_seq, res_n_100$cusum_widths, res_n_100$multiresolution_widths)


## n = 500 

n <- 500; cpt_loc <- n/2
cusum_thresh_500 <- bisect_get_cusum_thresh(n, 4, 10, alpha) # 4.644531

res_n_100 <- run_simulation_wrapper(n, cpt_loc, rho_seq, alpha, cusum_thresh_500, monte_carlo_reps, M)

plot_simulation_wrapper(rho_seq, res_n_100$cusum_widths, res_n_100$multiresolution_widths)


# n = 1000 


n <- 1000; cpt_loc <- n/2
cusum_thresh_1000 <- bisect_get_cusum_thresh(n, 4, 10, alpha) # 4.84375

res_n_100 <- run_simulation_wrapper(n, cpt_loc, rho_seq, alpha, cusum_thresh_1000, monte_carlo_reps, M)

plot_simulation_wrapper(rho_seq, res_n_100$cusum_widths, res_n_100$multiresolution_widths)

