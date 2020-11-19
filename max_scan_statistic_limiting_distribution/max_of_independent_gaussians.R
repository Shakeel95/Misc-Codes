library(foreach)


#-------------------------------------------------------------------------
#
# Exact plots to verify limit distribution of max of independent gaussians
#
#-------------------------------------------------------------------------


# normalizing constants

a <- function(n) sqrt(2*log(n)) - (0.5*log(log(n)) + log(2*sqrt(pi))) / sqrt(2*log(n))

b <- function(n) 1 / sqrt(2*log(n))
  

# exact & limit CDFs

exact.distr <- function(x,n) pnorm(a(n) + b(n)*x) ** n 

limit.distr <- function(x) exp(-exp(-x))


# grid 

x <- seq(from = -2, to = 3, length.out = 1000)

n.seq <- seq(from = 10, to = 100, by = 30)


# CDF plots

plot(x, limit.distr(x), type = "l", ylab = "CDF")

for (n in n.seq) lines(x, sapply(x, function(j) exact.distr(j,n)), col = n, type = "l")



#-------------------------------------------------------------------------
#
# Simulations to verify limit distribution of max of independent gaussians
#
#-------------------------------------------------------------------------


# empirical CDFs

CDFs <- list()

for (e in 1:length(n.seq)) {
  
  n <- n.seq[e]
  
  max.sim <- foreach(i=1:1000, .combine = c) %do% (max(rnorm(n)) - a(n)) / b(n)
  
  CDFs[[e]] <- ecdf(max.sim)
  
}



# plot limit distribution and ecdfs

plot(x, limit.distr(x), type = "l")

for (e in 1:length(n.seq)) lines(x, CDFs[[e]](x), col = e)
