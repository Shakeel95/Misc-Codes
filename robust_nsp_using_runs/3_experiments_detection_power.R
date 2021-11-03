
#======================================
#
# Numerical studies: detection power
#
#======================================


##
## Changepoint models 
##


model.blocks <- list(
  name = "blocks",
  simulate = function() c(rep(0,205),rep(14.64,62),rep(-3.66,41),rep(7.32,164),rep(-7.32,40),rep(10.98,308),rep(-4.39,82),rep(3.29,430),rep(19.03,225),rep(7.68,41),rep(15.37,61),rep(0,389)) + rnorm(2048, sd = 10), 
  cpt.locs = c(205, 267, 308, 472, 512, 820, 902, 1332, 1557, 1598, 1659)
)

model.cauchy <- list(
  name = "Cauchy",
  simulate = function() c(rcauchy(100, 1), rcauchy(100, 2), rcauchy(100, 1)), 
  cpt.locs = c(100, 200, 300)
)

model.bursts <- list(
  name = "Bursts",
  simulate = function() (c(rep(1, 200), rep(3, 80), rep(1, 200), rep(3, 80), rep(1, 200),rep(4, 40)) * rnorm(800))**2, 
  cpt.locs = c(200, 280, 480, 560, 760)
)

model.poisson <- list(
  name = "Poisson",
  simulate = function() as.numeric(rpois(350, c(rep(1, 50), rep(4, 50), rep(10, 50), rep(2, 200)))), 
  cpt.locs = c(50, 100, 150)
)


all.changepoint.models <- list(blocks = model.blocks, 
                               cauchy = model.cauchy, 
                               bursts = model.bursts, 
                               poisson = model.poisson
                               )

##
## Evaluation metrics
##


check.for.spurious.intervals <- function(nsp.obj, model.obj)
{
  #' Checks whether any spurious intervals were detected
  #'
  #'@param nsp.obj
  #'@param model.obj
  
  cpt.locs <- model.obj$cpt.locs
  detection.intervals <- nsp.obj$intervals
  
  if (nrow(detection.intervals) == 0) return(FALSE)
  
  for (ii in 1:nrow(detection.intervals))
  {
    if(!any(cpt.locs %in% detection.intervals[ii,1]:detection.intervals[ii,2])) return(TRUE)
  }
  
  return(FALSE)
}


get.prop.genuine.intervals <- function(nsp.obj, model.obj)
{
  #' Returns proportion of intervals which contain a changepoint
  #'
  #'@param nsp.obj
  #'@param model.obj

  cpt.locs <- model.obj$cpt.locs
  detection.intervals <- nsp.obj$intervals
  
  if (nrow(detection.intervals) == 0) return(NaN)
  
  correct.detections <- 0
  
  for (ii in 1:nrow(detection.intervals))
  {
    if(any(cpt.locs %in% detection.intervals[ii,1]:detection.intervals[ii,2])) correct.detections <- correct.detections + 1
  }
  
  return(correct.detections / nrow(detection.intervals))
}


get.num.genuine.intervals <- function(nsp.obj, model.obj)
{
  #' Returns the number of intervals which contain a changepoint
  #'
  #'@param nsp.obj
  #'@param model.obj
  
  cpt.locs <- model.obj$cpt.locs
  detection.intervals <- nsp.obj$intervals
  
  if (nrow(detection.intervals) == 0) return(0)
  
  correct.detections <- 0
  
  for (ii in 1:nrow(detection.intervals))
  {
    if(any(cpt.locs %in% detection.intervals[ii,1]:detection.intervals[ii,2])) correct.detections <- correct.detections + 1
  }
  
  return(correct.detections)
}


get.average.interval.length <- function(nsp.obj, model.obj)
{
  #' Returns the number of intervals which contain a changepoint
  #'
  #'@param nsp.obj
  #'@param model.obj

  detection.intervals <- nsp.obj$intervals
  
  if (nrow(detection.intervals) == 0) return(NaN)
  
  mean(apply(detection.intervals[,1:2],1,diff))
}


all.metrics <- list(spurious = check.for.spurious.intervals, 
                    prop.genuine = get.prop.genuine.intervals,
                    num.genuine = get.num.genuine.intervals,
                    average.length = get.average.interval.length)


update.metrics <- function(cur.metrics, all.metrics, nsp.obj, model.obj)
{
  #'
  #'
  #'@param cur.metrics vector or matrix
  #'@param all.metrics list of functions
  #'@param nsp.obj 
  #'@param model.obj 

  new.metrics <- c()
  
  for (metric in all.metrics) new.metrics <- c(new.metrics, metric(nsp.obj, model.obj))
  
  cbind(cur.metrics, new.metrics)
}


##
## Numerical study
##

set.seed(42)

n.reps <- 100

pb <- txtProgressBar(min = 1, max = n.reps, style = 3)

cpt.models.metrics <- list()

for (changepoint.model in all.changepoint.models)
{

  multires.metrics <- c()
  runs.metrics <- c()
  
  for (jj in 1:n.reps)
  {
    cat("\r", changepoint.model$name)
    
    ee <- changepoint.model$simulate()
    
    rnsp.multires.obj <- nsp_robust(ee, method = "multiresolution")
    multires.metrics <- update.metrics(multires.metrics, all.metrics, rnsp.multires.obj, changepoint.model)
    
    rnsp.runs.obj <- nsp_robust(ee, method = "runs")
    runs.metrics <- update.metrics(runs.metrics, all.metrics, rnsp.runs.obj, changepoint.model)
    
    setTxtProgressBar(pb, jj)
  }
  
  metrics.averages <- data.frame(
    multiresolution = apply(multires.metrics, 1, function(ii) mean(ii, na.rm = TRUE)), 
    runs = apply(runs.metrics, 1, function(ii) mean(ii, na.rm = TRUE))
  )
  
  rownames(metrics.averages) <- names(all.metrics)
  
  cpt.models.metrics[[changepoint.model$name]] <- metrics.averages 
}

