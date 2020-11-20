library(foreach)
library(parallel)
library(doParallel)


#-----------------------------------------
#
# Simulations to check scan stat distr...
#
#-----------------------------------------


# register backend + set seed

set.seed(42)

cl <- makeCluster(detectCores() - 1)

registerDoParallel(cl)


# normalizing constants

H <- 0.858547 

a <- function(n) sqrt(2*log(n)) - (0.5*log(log(n)) + log(2*sqrt(pi))) / sqrt(2*log(n))

b <- function(n) 1 / sqrt(2*log(n))


# limiting CDF

limit.distr <- function(x) exp(-exp(-x))


# grids 

x <- seq(from = -5, to = 5, length.out = 1000)

steps <- 5

tt.seq <- seq(from = 10, to = 50, length.out = steps)

tt.cols <- colorRampPalette(c("blue","red"))(steps)

n.seq <- seq(from = 10, to = 50, length.out = steps)

n.cols <- colorRampPalette(c("blue","red"))(steps)


#-------------------------------------------------------------
#
# Simulations for increasing number of channels (n -> \infty)
#
#-------------------------------------------------------------


tt <- 100


# get empirical CDFs

CDF.n.increasing <- list()


for (q in 1:steps) {
  
  
  n <- n.seq[q]
  
  L.infty.scans <- foreach(i=1:100, .combine = c) %do% {
    
    
    channel.scans <- foreach(j=1:n, .combine = c, .inorder = F) %dopar% {
      
      e <- rnorm(tt)
      
      std.increments <- combn(1:tt,2, function(k) (sum(e[k[1]:k[2]])) / sqrt(k[2] + 1 - k[1]))
      
      max(std.increments)
      
    }
    
    (max(channel.scans) - a(H*n*log(n)*tt)) / b(H*n*log(n)*tt)
    
  }
  
  
  cat(paste0("\r", "finished calculating ecdf for n = ", n.seq[q]))
  
  CDF.n.increasing[[q]] <- ecdf(L.infty.scans)
  
}



# plot limiting cdf and ecdf

plot(x, limit.distr(x), type = "l", lty = 2)

for (q in 1:steps) lines(x, CDF.n.increasing[[q]](x), col = n.cols[q])

legend(-5, 1, 
       legend = sapply(n.seq, function(i) paste0("n = ", i)), 
       col = n.cols, 
       lty = rep(1,steps))


#----------------------------------------------------------------
#
# Simulations for increasing number of time points (T -> \infty)
#
#----------------------------------------------------------------


n <- 100


# get empirical CDFs

CDF.tt.increasing <- list()


for (q in 1:steps) {
  
  
  tt <- tt.seq[q]
  
  L.infty.scans <- foreach(i=1:100, .combine = c) %do% {
    
    
    channel.scans <- foreach(j=1:n, .combine = c, .inorder = F) %dopar% {
      
      e <- rnorm(tt)
      
      std.increments <- combn(1:tt,2, function(k) (sum(e[k[1]:k[2]])) / sqrt(k[2] + 1 - k[1]))
      
      max(std.increments)
      
    }
    
    (max(channel.scans) - a(H*n*log(n)*tt)) / b(H*n*log(n)*tt)
    
  }
  
  
  cat(paste0("\r", "finished calculating ecdf for n = ", tt.seq[q]))
  
  CDF.tt.increasing[[q]] <- ecdf(L.infty.scans)
  
}



# plot limiting cdf and ecdf

plot(x, limit.distr(x), type = "l", lty = 2)

for (q in 1:steps) lines(x, CDF.tt.increasing[[q]](x), col = n.cols[q])

legend(-5, 1, 
       legend = sapply(n.seq, function(i) paste0("T = ", i)), 
       col = n.cols, 
       lty = rep(1,steps))



# stop paralell backend

stopCluster(cl)
