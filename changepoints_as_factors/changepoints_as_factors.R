#---------------------------------------------------------------------#
#                                                                     #
# Explore whether changepoints in a panel can be mistake for factors  #
#                                                                     #    
#   (1) simulate panel with changepoints + uncorrelated noise         #
#   (2) identify number of factoors with eigenvalue ratio test        #
#   (3) check whether number of factors corresponsd to number of      #
#                                                                     #
#---------------------------------------------------------------------#


library("mvtnorm") # multivariate normal
library("docstring") # EZ documentation


#----------------------------------------#
# Functions for simulation and detection #
#----------------------------------------#


eigenvalue.ratio <- function(S) {
  
  #' Ratio of successive eigenvalues
  #' 
  #' @param S symmetric matrix
  
  if (!isSymmetric.matrix(S)) stop("S should be a symmetric matrix.")
  
  eig.vals <- eigen(S)$values
  
  k <- ncol(S)
  
  sapply(1:(k-1), function(i) eig.vals[i] / eig.vals[i+1])
  
}


pcws.lin.panel <- function(n, tt, n.cpt, d.cpt, min.change = -5, max.change = 5) {
  
  #' SImulate a panel with piecewise linear mean
  #'
  #'@param n int, number of channels
  #'@param tt int, number of time points
  #'@param n.cpt int, number of changepoints
  #'@param d.cpt float, proportion of chanels with change 
  #'@param min.change float, min jump in mean at changepoint
  #'@param max.change float, max jump in mean at changepoint
  
  E <- rmvnorm(tt, rep(0,n), diag(rep(1,n)))
  A <- matrix(0,tt,n)
  
  m <- ceiling(n*d.cpt)
  
  cpt.locs <- sort(sample(2:(tt-1), n.cpt))
  
  for (tau in cpt.locs) {
    
    chanels <- sample(1:n,m)
    delta <- runif(m, min.change, max.change)
    
    A[(tau+1):tt,chanels] <- A[(tau+1):tt,chanels] + t(replicate(tt-tau,delta)) 
    
  }
  
  return(A + E)
  
}


#--------------------------# 
# Simulation               # 
#--------------------------#

set.seed(42)
sim.size <- 100

n <- 50
tt <- 100 
n.cpt <- 20
d.cpt <- 0.8

res.eig <- matrix(,n,sim.size)
res.eig.rartio <- matrix(,n-1,sim.size)


for (i in 1:sim.size) {
  
  X <- pcws.lin.panel(n,tt,n.cpt, d.cpt)
  S <- cov(X)
  
  res.eig[,i] <- eigen(S)$values
  res.eig.rartio[,i] <- eigenvalue.ratio(S)
  
  cat("\r","finished simulation ", paste0(i))
  
}

matplot(res.eig, type = "l")

plot(apply(res.eig,1,mean), type = "b")
abline(v = n.cpt, col = "red", lty = 2)

matplot(res.eig.rartio, type = "l")

plot(apply(res.eig.rartio,1,mean), type = "b")
abline(v = n.cpt, col = "red", lty = 2)


