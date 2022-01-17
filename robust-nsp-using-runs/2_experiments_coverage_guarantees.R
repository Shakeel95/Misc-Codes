
#========================================
#
# Numerical studies: coverage guarantees
#
#========================================


##
## All models from RNSP paper
## 


model.list <- list(
  plain.gauss = function() rnorm(100), 
  plain.poisson = function() as.numeric(rpois(200, 1)), 
  heterogeneous.gauss = function() c(rep(1, 100), rep(8, 50), rep(1, 100)) * rnorm(250), 
  symmetric.bernoulli = function() as.numeric(rbinom(200, 1, .5)), 
  plain.cauchy = function() rcauchy(100), 
  mix.1 = function()
  {
    sample(3, size=300, replace = TRUE, prob=c(.35, .3, .35)) -> xx
    xx[xx != 2] <- rnorm(sum(xx !=2 ))
    xx
  }, 
  mix.2 = function() rpois(200, 5) + rnorm(200)/30
)


##
## Numerical study
##


set.seed(42)

n.reps <- 100

pb <- txtProgressBar(min = 1, max = n.reps, style = 3)

coverage.probabilities <- matrix(0, ncol = 2, nrow = length(model.list))
rownames(coverage.probabilities) <- names(model.list)
colnames(coverage.probabilities) <- c("multiresolution", "runs")


for (ii in 1:length(model.list))
{
  
  multiresolution.detections <- numeric(n.reps)
  runs.detections <- numeric(n.reps)
  
  cat("\r", names(model.list)[ii])
  
  for (jj in 1:n.reps)
  {
    ee <- model.list[[ii]]()
    multiresolution.detections[jj] <- nrow(nsp_robust(ee, method = "multiresolution")$intervals) == 0
    runs.detections[jj] <- nrow(nsp_robust(ee, method = "runs")$intervals) == 0
    
    setTxtProgressBar(pb, jj)
  }
  
  coverage.probabilities[ii,1] <- 100 * mean(multiresolution.detections)
  coverage.probabilities[ii,2] <- 100 * mean(runs.detections)
}



