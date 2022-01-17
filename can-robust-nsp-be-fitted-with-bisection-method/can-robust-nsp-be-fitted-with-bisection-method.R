
#===============
#
# Functions
#
#================


multi.sup.norm <- function(x)
{
  #' Multi-resolution sup norm of a vector
  #'
  #'@param x vector 
  
  cumsum.x <- c(0, cumsum(x))
  n <- length(x)
  
  loc.max <- 0 
  run.max <- 0 
  
  for (s in 1:n) for (e in (s+1):(n+1)) 
  {
    loc.max <- abs(cumsum.x[e] - cumsum.x[s]) / sqrt(e-s)
    
    if (loc.max > run.max) run.max <- loc.max
  }
  
  return(run.max)
  
}


add.midpoints.sort <- function(x)
{
  #' Adds midpoints to vector and returns sorted values
  #'
  #'@param x vector
  
  n <- length(x)
  
  x.mid <- numeric(n-1)
  
  for (ii  in 1:(n-1)) x.mid[ii] <- (x[ii] + x[ii+1]) / 2
  
  x <- c(x, x.mid)
  
  return(sort(x))
  
}


plot.wrapper <- function(x)
{
  #' Plots sign multi-resolution norm of residuals against candidate fit
  #' for all 2*n - 1 candidate fits leading to different residual sign sequence
  #' combinations
  #'
  #'@param x
  #'@references 

  n <- length(x)
    
  xx <- add.midpoints.sort(x)

  xx.sup.norm <- c()
  
  for (yy in xx) 
  {
    xx.sup.norm <- c(
      xx.sup.norm, 
      multi.sup.norm(sign(x - rep(yy,n)))
    )
  }
  
  plot(xx.sup.norm ~ xx, 
       type = "l",
       xlab = "candidate fit",
       ylab = "sup-norm of residuals"
       )
  
}

#========================
#
# Simultion study
#
#========================


# global params 

set.seed(100)

n <- 100 


# gaussian noise 

e.gauss <- rnorm(n)

plot.wrapper(e.gauss)


# cuachy noise 

e.cauchy <- rcauchy(n)

plot.wrapper(e.cauchy)


# poison noise

e.pois <- rpois(n, 1)

plot.wrapper(e.pois)


# hypergeometric noise

e.hyper <- rhyper(n, 10, 5, 4)

plot.wrapper(e.hyper)
