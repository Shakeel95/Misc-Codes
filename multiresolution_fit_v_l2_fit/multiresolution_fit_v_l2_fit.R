
#-----------------------------------
#
# Multi-resolution fit functions
#
#-----------------------------------


library(lpSolve)


fit.with.multiresolution <- function(y,x, include.intercept = TRUE)
{
  #' Fit linear model params using multiresolution norm
  #'
  #'@param y vector, dim (n x 1)
  #'@param x matrix, dim (n x k)
  #'@param include.intercept bool 

  n <- length(y)
  
  if (include.intercept) x <- cbind(rep(1,n),x)
  if (is.vector(x)) x <- matrix(x)
  k <- ncol(x)-1
  
  y.cumsum <- c(0, cumsum(y))
  x.cumsum <- rbind(rep(0,k+1), apply(x, 2, cumsum))

  lhs.constraints <- c()
  rhs.constraints <- c()
  
  for (ii in 1:n) for (jj in (ii+1):(n+1))
  {
    rhs.constraints <- c(rhs.constraints, (y.cumsum[jj] - y.cumsum[ii]) / sqrt(jj-ii))
    lhs.constraints <- rbind(lhs.constraints, (x.cumsum[jj,] - x.cumsum[ii,]) / sqrt(jj-ii))
  }
  
  rhs.constraints <- c(rhs.constraints, -rhs.constraints)
  
  q <- length(rhs.constraints)
  
  lhs.constraints <- cbind(lhs.constraints, -lhs.constraints)
  lhs.constraints <- rbind(lhs.constraints, -lhs.constraints)
  lhs.constraints <- cbind(rep(1,q), lhs.constraints)
  
  f.obj <- c(1, rep(0, 2*(k+1)))
  f.dir <- rep(">=",q)
  
  lp.solution <- lp(
    direction = "min",
    objective.in = f.obj, 
    const.mat = lhs.constraints, 
    const.dir = f.dir,
    const.rhs = rhs.constraints
  )
  
  return(lp.solution$solution[-1])
  
}


get.multiresolution.fitted.values <- function(mr.fit, x)
{
  #' Get fitted values from multi-resolution sup-norm regression
  #'
  #'@param mr.fit vector, output from `fit.with.multiresolution`
  #'@param x matrix, dim (n x k)

  
   if (is.vector(x))
   {
     if (any(mr.fit > 0))
     {
       ind <- which(mr.fit > 0)
       param <- mr.fit[ind] * ifelse(ind == 1, 1, -1)
       return(param * x)
     }
     
     return(x*0)
   }
   
   q <- length(mr.fit) / 2
   res <- numeric(nrow(x))
   
   for (ii in 1:q) {
     
     if (any(mr.fit[c(ii,q+ii)] > 0))
     {
      ind <- which(mr.fit[c(ii,q+ii)] > 0)
      param <- mr.fit[c(ii,q+ii)][ind] * ifelse(ind == 1, 1, -1)
      res <- res + x[,ii] * param
     }
   }
   
   return(res)
}


#-----------------------------------------------------------------------------
#
# Multiresolution fit and least squares fit in pieicewise polynomial settings
#
#------------------------------------------------------------------------------

n <- 100
x <- 1:n / n 
e <- rnorm(n)

par(mfrow = c(1,2))

b.loc <- floor(n/2)
x.lhs <- x[1:b.loc]
x.rhs <- x[(b.loc+1):n]
b <- x[b.loc]


## pieicewise constant

s <- c(
  rep(1, b.loc),
  rep(2, n - b.loc)
)

lm.fit.noiseless <- lm(s ~ 1)
mr.fit.noiseless <- fit.with.multiresolution(s,x**0, include.intercept = FALSE)

plot(s ~ x, type = "l", col = "grey", main = "noiseless fit")
lines(lm.fit.noiseless$fitted.values ~ x , col = "red")
lines(get.multiresolution.fitted.values(mr.fit.noiseless, x**0) ~ x, col = "blue")



y <- s + e

lm.fit.noisy <- lm(y ~ 1)
mr.fit.noisy <- fit.with.multiresolution(y,x**0, include.intercept = FALSE)

plot(y ~ x, type = "l", col = "grey", main = "Noisy fit")
lines(s ~ x, col = "grey")
lines(lm.fit.noisy$fitted.values ~ x, col = "red")
lines(get.multiresolution.fitted.values(mr.fit.noisy, x**0) ~ x, col = "blue")


## piecewise linear (kink)

s <- c(
  1 + 0.5*(x.lhs - b),
  1 + 5*(x.rhs - b)
)

lm.fit.noiseless <- lm(s ~ x)
mr.fit.noiseless <- fit.with.multiresolution(s, x)

plot(s ~ x, type = "l", col = "grey", main = "noiseless fit")
lines(lm.fit.noiseless$fitted.values ~ x, col = "red", lty = 2)
lines(get.multiresolution.fitted.values(mr.fit.noiseless, cbind(x**0, x)) ~ x, col = "blue")


y <- s + e

lm.fit.noisy <- lm(y ~ x)
mr.fit.noisy <- fit.with.multiresolution(y,x)

plot(y~x, type = "l", col = "grey", main = "noisy fit")
lines(s ~ x, col = "grey")
lines(lm.fit.noisy$fitted.values ~ x, col = "red")
lines(get.multiresolution.fitted.values(mr.fit, cbind(x**0, x)) ~ x, col = "blue")


## piecewise linear (jump)

s <- c(
  1 + 0.5*(x.lhs - b),
  1.5 + 1.5*(x.rhs - b)
)

lm.fit.noiseless <- lm(s ~ x)
mr.fit.noiseless <- fit.with.multiresolution(s,x)

plot(s ~ x, type = "l", col = "grey", main = "noiseless fit")
lines(lm.fit.noiseless$fitted.values ~ x, col = "red")
lines(get.multiresolution.fitted.values(mr.fit.noiseless, cbind(x**0, x)) ~ x, col = "blue")

y <- s + e

lm.fit.noisy <- lm(y ~ x)
mr.fit.noisy <- fit.with.multiresolution(y,x)

plot(y~x, type = "l", col = "grey", main = "noisy fit")
lines(s ~ x, col = "grey")
lines(lm.fit.noisy$fitted.values ~ x, col = "red")
lines(get.multiresolution.fitted.values(mr.fit, cbind(x**0, x)) ~ x, col = "blue")


## pieicewise quadratic

s <- c(
  1 + 0.5*(x.lhs - b) + 20*(x.lhs - b)**2,
  1 + 5*(x.rhs - b) + 2*(x.rhs - b)**2
)

lm.fit.noiseless <- lm(s ~ x + I(x**2))
mr.fit.noiseless <- fit.with.multiresolution(s, cbind(x, x**2))

plot(s ~ x, type = "l", col = "grey", main = "noiseless fit")
lines(lm.fit.noiseless$fitted.values ~ x, col = "red")
lines(get.multiresolution.fitted.values(mr.fit.noiseless, cbind(x**0, x, x**2)) ~ x, col = "blue")

y <- s + e 

lm.fit.noisy <- lm(y ~ x + I(x**2))
mr.fit.noisy <- fit.with.multiresolution(y, cbind(x, x**2))


plot(y ~ x, type = "l", col = "grey", main = "noisy fit")
lines(s ~ x, col = "grey")
lines(lm.fit.noiseless$fitted.values ~ x, col = "red")
lines(get.multiresolution.fitted.values(mr.fit.noiseless, cbind(x**0, x, x**2)) ~ x, col = "blue")

