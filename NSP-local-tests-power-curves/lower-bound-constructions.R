#=====================================================================================
#
# Robust lower bounds on isotonic regression functions by inverting different tests
#
#======================================================================================

##
## Scale corrected multi-resolution sign test: Dumbgen and Johns

thresh.dumbgen.T0 <- function(nn)
{
  #' Critical values for T^{0} in Dumbgen and Johns at \alpha = 0.1
  #'
  #'@param nn int, sample size
  
  if (nn <= 100) return(1.120)
  if (nn <= 200) return(1.198)
  if (nn <= 300) return(1.246)
  if (nn <= 700) return(1.297)
  
  if (nn <= 1000) return(1.325)
  if (nn <= 2000) return(1.387)
  if (nn <= 5000) return(1.424)
  
  return(1.467)
}

corrected.multi.res.lower.bound <- function(yy, kappa) 
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
  NN <- kappa[1]
  lambda <- kappa[2]
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
      
      if (SS <= sqrt(dd) * (lambda + sqrt(2*log(exp(1)*NN/dd))))
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


##
## Multi-resolution sign tests: Fryzlewicz


thresh.kab.bern <- function(nn, alpha = 0.1) 
{
  #' Threshold for vanilla multi-resolution test
  #'
  #'@param nn int, sample size 
  #'@param alpha float, test size
  #'
  #'@references https://doi.org/10.1016/j.spa.2014.03.015
  #'@references https://arxiv.org/abs/2109.02487
    
  an <- sqrt(2 * log(nn*log(nn)^(-1/2)))
  tau <- -log(-1/(2*0.2740311) * log(1 - alpha))
  
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
