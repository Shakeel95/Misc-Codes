
#-------------------------------------------------
#
# Check coverage guarantees for new NSP method 
# based on inverting multi-resolution tests
#
#-------------------------------------------------

library(magrittr)

##
## Models from RNSP paper

model.list <- list(
  plain.gauss = function() rnorm(100),
  plain.poisson = function() as.numeric(rpois(200, 1)),
  heterogeneous.gauss = function() c(rep(1, 100), rep(8, 50), rep(1, 100)) * rnorm(250),
  symmetric.bernoulli = function() as.numeric(rbinom(200, 1, .5)),
  plain.cauchy = function() rcauchy(100),
  mix.2 = function() rpois(200, 5) + rnorm(200)/30
)


##
## Simulation study

set.seed(42)

n.reps <- 100

pb <- txtProgressBar(min = 1, max = n.reps, style = 3)
coverage.probabilities <- matrix(0, ncol = 4, nrow = length(model.list))
rownames(coverage.probabilities) <- names(model.list)

colnames(coverage.probabilities) <- c("multi.res", "corrected.multi.res")


for (ii in 3:length(model.list))
{
  multi.res.detections <- numeric(n.reps)
  corrected.multi.res.detections <- numeric(n.reps)
  
  cat("\r", names(model.list)[ii])
  
  for (jj in 1:n.reps)
  {
    ee <- model.list[[ii]]()
    
    multi.res.detections[jj] <- try(
      nrow(nsp_robust(ee, method = "multi.res", thresh = thresh.kab.bern(length(ee)))$intervals) == 0
    )
    if ("try-error" %in% class(multi.res.detections[jj])) multi.res.detections[jj] <- NaN
      
    corrected.multi.res.detections[jj] <- try(
      nrow(nsp_robust(ee, method = "corrected.multi.res", thresh = thresh.dumbgen.T0(length(ee)))$intervals) == 0
    )    
    if ("try-error" %in% class(corrected.multi.res.detections[jj])) corrected.multi.res.detections[jj] <- NaN

    setTxtProgressBar(pb, jj)
  }
  
  coverage.probabilities[ii,1] <- 100 * mean(multi.res.detections, na.rm = TRUE)
  coverage.probabilities[ii,2] <- 100 * mean(corrected.multi.res.detections, na.rm = TRUE)
}

