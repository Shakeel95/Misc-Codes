
#============================================================================
#
# Exact plots to verify limit distribution of max of independent gaussians
#
#============================================================================

library(foreach)


#--------
#
# Setup 
#
#--------


# normalizing constants

a <- function(n) sqrt(2*log(n)) - (0.5*log(log(n)) + log(2*sqrt(pi))) / sqrt(2*log(n))

b <- function(n) 1 / sqrt(2*log(n))
  

# exact & limit CDFs

exact.cdf <- function(x,n) pnorm(a(n) + b(n)*x) ** n 

limit.cdf <- function(x) exp(-exp(-x))


# exact & limit PDFs

limit.pdf <- function(x) exp(-x)*exp(-exp(-x))

exact.pdf <- function(x,n) n * (pnorm(a(n) + b(n)*x)**(n-1)) * dnorm(a(n) + b(n)*x) * b(n)


# grid

steps <- 3

K <- 1000 

x <- seq(from = -3, to = 6, length.out = 1000)

n.seq <- c(seq(from = 5, to = 10, length.out = steps),1000)

cols <- colorRampPalette(c("blue", "red"))(steps)



#----------------------------
#
# Exact PDF & CDF plots
#
#----------------------------


# exact CDF v. limit CDF

plot(x, limit.cdf(x), type = "l", lty = 2, ylab = "CDF")

for (n in n.seq) lines(x, exact.cdf(x,n), col = n.cols[which(n.seq == n)], type = "l")

legend(-3, 1, 
       legend = sapply(n.seq, function(i) paste0("n = ", round(i))), 
       col = n.cols, 
       lty = rep(1,length(n.seq)))


# exact PDF v. limit PDF

plot(x, limit.pdf(x), type = "l", lty = 2, ylab = "PDF")

for (n in n.seq) lines(x, exact.pdf(x,n), col = n.cols[which(n.seq == n)])

legend(-3, 0.35, 
       legend = sapply(n.seq, function(i) paste0("n = ", round(i))), 
       col = n.cols, 
       lty = rep(1,length(n.seq)))



#-------------------------------------------------------------------------
#
# Simulations to verify limit distribution of max of independent gaussians
#
#-------------------------------------------------------------------------


# empirical CDFs

CDFs <- list()

for (e in 1:length(n.seq)) {
  
  n <- n.seq[e]
  
  max.sim <- foreach(i=1:K, .combine = c) %do% (max(rnorm(n)) - a(n)) / b(n)
  
  CDFs[[e]] <- ecdf(max.sim)
  
}



# plot limit distribution and ecdfs

plot(x, limit.cdf(x), type = "l", lty = 2)

for (e in 1:length(n.seq)) lines(x, CDFs[[e]](x), col = n.cols[e])

legend(-3, 1, 
       legend = sapply(n.seq, function(i) paste0("n = ", round(i))), 
       col = n.cols, 
       lty = rep(1,length(n.seq)))
