
#======================================
#
# Numerical studies: execution time
#
#======================================


library(microbenchmark)


rnsp.timer <- function(n.seq, ff, reps = 1, digits = 2)
{
  #' Get execution times for RNSP on specified distribution
  #'
  #'@param n.seq int, length of noise vector
  #'@param ff function, generates noise 
  #'@param reps int, number of function calls to time
  #'@param digits int, numer of digits to round to

  res <- matrix(, ncol = 2, nrow = length(n.seq))
  colnames(res) <- c("multiresolution", "runs")
  rownames(res) <- n.seq

  for (ii in 1:length(n.seq))
  {
    ee <- ff(n.seq[ii])
    
    res[ii,1] <- microbenchmark(nsp_robust(ee,method = "multiresolution"), times = reps)$time / 10**9
    res[ii,2] <- microbenchmark(nsp_robust(ee,method = "runs"), times = reps)$time / 10**9
  }
  
  return(round(res, digits))
  
}


##
## Time function calls
##

set.seed(42)

n.seq <- c(100,500,1000)

ff.pois <- function(nn) rpois(nn, lambda = 1)

ff.geom <- function(nn) rgeom(nn, prob = .1)


call.times <- list(
  gaussian.noise = rnsp.timer(n.seq, ff = rnorm), 
  cauchy.noise = rnsp.timer(n.seq, ff = rcauchy), 
  poisson.noise = rnsp.timer(n.seq, ff = ff.pois), 
  geometric.noise = rnsp.timer(n.seq, ff = ff.geom)
)



for (mm in call.times)
{
  barplot(t(mm), 
          xlab = "signal length",
          ylab = "execuition time (seconds)",
          beside = TRUE)
}
