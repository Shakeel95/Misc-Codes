#--------------------------------------------
#
# Inverting the runs test (Davies and Kovac)
#
#--------------------------------------------

runs.lower.bound <- function(yy, max.runs)
{
  #' Lower bound on isotonic reg function by inverting runs test
  #'
  #'@param yy data vector
  #'@param max.runs 
  #'
  #'@references DOI: 10.1214/aos/996986501
  
  nn <- length(yy)
  max.runs <- floor(max.runs)
  
  ll <- numeric(nn)
  ll[1:max.runs] <- -Inf
  
  for (ii in (max.runs+1):nn) ll[ii] <- max(c(ll[ii-1], min(yy[(ii-max.runs):ii])))

  ll
}


runs.thresh <- function(nn, alpha = .1)
{
  #' Threshold for runs test (own derivation)
  #'
  #'@param nn int, sample size 
  #'@param alpha float, test size
  
  log2(nn-1) + (1/log(2)) * log(1/log(1/(1-alpha)))
}



#-------------------------------------------------------------
#
# Inverting multi-resolution sign tests (Dumbgen and Johns)
#
#--------------------------------------------------------------

centered.multi.res.lower.bound <- function(yy, kappa) 
{
  #' Lower bound on isotonic reg function by inverting multi-resolution sign tests with scale correction
  #'
  #' Corresponds to test T^{0} in the paper (pp. 522)
  #'
  #'@param yy data vector
  #'@param kappa critical value (pp. 523, table 1)
  #'
  #'@references https://doi.org/10.1198/1061860043506

  nn <- length(yy)
  ll <- numeric(nn)
  
  for (kk in 1:nn)
  {
    if (kk == 1)
    {
      ll[kk] <- -Inf
    } else {
      ll[kk] <- ll[kk-1]
    }
    
    jj <- kk 
    
    while(jj > 0)
    {
      if (jj == kk)
      {
        SS <- 0 
        r.new <- Inf
      }
      
      if (yy[min(jj+1,nn)] > ll[kk])
      {
        SS <- SS + 1
        r.new <- min(r.new, yy[min(jj+1,nn)])
      } else {
        SS <- SS - 1
      }
      
      dd <- kk-jj+1
      
      if (SS <= sqrt(dd) * (kappa + sqrt(2*log(exp(1)*nn/dd))))
      {
        jj <- jj - 1
      } else {
        ll[kk] <- r.new
        jj <- kk 
      }
    }
  }
  
  ll
}


weighted.multi.res.lower.bound <- function(yy, kappa) 
{
  #' Lower bound on isotonic reg function by inverting weighted multi-resolution sign tests
  #'
  #' Corresponds to test T^{1} in the paper (pp. 522)
  #'
  #'@param yy data vector
  #'@param kappa critical value (pp. 523, table 1)
  #'
  #'@references https://doi.org/10.1198/1061860043506
  
  nn <- length(yy)
  ll <- numeric(nn)
  
  for (kk in 1:nn)
  {
    if (kk == 1)
    {
      ll[kk] <- -Inf
    } else {
      ll[kk] <- ll[kk-1]
    }
    
    jj <- kk
    
    while(jj > 0)
    {
      if (jj == kk)
      {
        S.0 <- 0
        S.1 <- 0
        r.new <- Inf
      }
      
      if (yy[min(jj+1,nn)] > ll[kk])
      {
        S.0 <- S.0 + 1
        r.new <- min(r.new, yy[min(jj+1,nn)])
      } else {
        S.0 <- S.0 - 1
      }
      
      S.1 <- S.0 + S.1
      dd <- kk-jj+1
      beta.d <- sqrt(dd*(dd+1)*(2*dd+1)/6)
      
      if (S.1 <= beta.d*(kappa + sqrt(2*log(exp(1*nn/dd)))))
      {
        jj <- jj - 1
      } else {
        ll[kk] <- r.new
        jj <- kk 
      }
    }
  }
  
  ll
}


#---------------------------------------------------------------
#
# Inverting 'vanilla' multi-resolution sign tests (Fryzlewicz)
#
#---------------------------------------------------------------


thresh_kab_bern <- function(nn, alpha = 0.1) 
{
  #' Threshold for vanilla multi-resolution test
  #'
  #'@param nn int, sample size 
  #'@param alpha float, test size
  #'
  #'@references https://doi.org/10.1016/j.spa.2014.03.015
  #'@references https://arxiv.org/abs/2109.02487
  
  alpha <- alpha / 2
  
  an <- sqrt(2 * log(nn*log(nn)^(-1/2)))
  tau <- -log(-1/(0.2740311) * log(1 - alpha))
  
  an + tau / an
}

multi.res.lower.bound <- function(yy, kappa) 
{
  #' Lower bound on isotonic reg function by inverting multi-resolution sign tests
  #'
  #' Specifically, inverts test in robust NSP paper
  #'
  #'@param yy data vector
  #'@param kappa critical value (pp. 523, table 1)
  #'
  #'@references https://arxiv.org/abs/2109.02487
  
  nn <- length(yy)
  ll <- numeric(nn)
  
  for (kk in 1:nn)
  {
    if (kk == 1)
    {
      ll[kk] <- -Inf
    } else {
      ll[kk] <- ll[kk-1]
    }
    
    jj <- kk 
    
    while(jj > 0)
    {
      if (jj == kk)
      {
        SS <- 0 
        r.new <- Inf
      }
      
      if (yy[min(jj+1,nn)] > ll[kk])
      {
        SS <- SS + 1
        r.new <- min(r.new, yy[min(jj+1,nn)])
      } else {
        SS <- SS - 1
      }
      
      dd <- kk-jj+1
      
      if (SS <= sqrt(dd)*kappa)
      {
        jj <- jj - 1
      } else {
        ll[kk] <- r.new
        jj <- kk 
      }
    }
  }
  
  ll
}
