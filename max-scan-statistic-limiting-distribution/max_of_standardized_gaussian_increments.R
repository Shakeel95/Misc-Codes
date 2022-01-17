library(docstring)
library(foreach)
library(parallel)
library(doParallel)


#---------------------------------------------------------------------------------
#
# Simulations to verify limit distribution maximum standardized gaussian increment
#
# Paper: Kabluchko (2018) https://arxiv.org/abs/0706.1849
#
#---------------------------------------------------------------------------------


# normalizing constants

# H value from Kabluchko's email to Piotr

H <- 0.858547 

a <- function(n) sqrt(2*log(n)) + (0.5*log(log(n)) + log(H) - log(2*sqrt(pi))) / sqrt(2*log(n))

b <- function(n) 1 / sqrt(2*log(n))



# limiting CDF

limit.distr <- function(x) exp(-exp(-x))


# grid 

x <- seq(from = -5, to = 5, length.out = 1000)

n.seq <- seq(from = 100, to = 500, by = 100)

n.cols <- colorRampPalette(c("blue", "red"))(5)



# get empirical CDFs

CDFs <- list()

cl <- makeCluster(8)

registerDoParallel(cl)



for (q in 1:length(n.seq)) {
  
  
  n <- n.seq[q]

  std.increments <- foreach(i=1:1000, .combine = c) %dopar% {

    e <- rnorm(n)

    std.increments <- combn(1:n,2, function(k) (sum(e[k[1]:k[2]])) / sqrt(k[2] + 1 - k[1]))

    (max(std.increments) - a(n)) / b(n)

  }
  
  cat(paste0("\r", "finished calculating ecdf for n = ", n.seq[q]))
  
  CDFs[[q]] <- ecdf(std.increments)
  
}


# plot limiting cdf and ecdf

plot(x, limit.distr(x), type = "l")

for (q in 1:length(n.seq)) lines(x, CDFs[[q]](x), col = n.cols[q])

legend(-5, 1, 
       legend = sapply(n.seq, function(i) paste0("n = ", i)), 
       col = n.cols, 
       lty = rep(1,length(n.seq)))


# stop paralell backend

stopCluster(cl)
