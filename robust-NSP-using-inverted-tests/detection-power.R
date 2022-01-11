
#======================================
#
# Detection Power Simulation
#
#======================================


#----------------------
#
# Changepoint Models
#
#----------------------

changepoint.models <- list( 

  model.blocks = list(
    name = "blocks",
    signal = c(rep(0,205),rep(14.64,62),rep(-3.66,41),rep(7.32,164),rep(-7.32,40)),
    simulate = function() c(rep(0,205),rep(14.64,62),rep(-3.66,41),rep(7.32,164),rep(-7.32,40)) + rnorm(512, sd = 10),
    cpt.locs = c(205, 267, 308, 472)
  ),

  model.cauchy = list(
    name = "Cauchy",
    signal = c(rep(2.5,100), rep(5,100), rep(2.5,100)),
    simulate = function() c(rcauchy(100, 1), rcauchy(100, 2), rcauchy(100, 1)),
    cpt.locs = c(100, 200, 300)
  ),

  model.bursts = list(
    name = "Bursts",
    signal = c(rep(1, 200), rep(3, 80), rep(1, 200), rep(3, 80), rep(1, 200),rep(4, 40)) ** 2,
    simulate = function() (c(rep(1, 200), rep(3, 80), rep(1, 200), rep(3, 80), rep(1, 200),rep(4, 40)) * rnorm(800))**2,
    cpt.locs = c(200, 280, 480, 560, 760)
  ),

  model.poisson = list(
    name = "Poisson",
    signal = c(rep(1, 50), rep(4, 50), rep(10, 50), rep(2, 200)),
    simulate = function() as.numeric(rpois(350, c(rep(1, 50), rep(4, 50), rep(10, 50), rep(2, 200)))),
    cpt.locs = c(50, 100, 150)
  )
)



#------------------------
#
# Evaluation metrics
#
#------------------------


metrics <- list(
  
  spurious = function(detection.intervals, model.obj)
  {
    #' Checks whether any spurious intervals were detected
    #'
    #'@param detection.intervals
    #'@param model.obj
    
    cpt.locs <- model.obj$cpt.locs
    
    if (nrow(detection.intervals) == 0) return(FALSE)
    
    for (ii in 1:nrow(detection.intervals))
    {
      if(!any(cpt.locs %in% detection.intervals[ii,1]:detection.intervals[ii,2])) return(TRUE)
    }
    
    return(FALSE)
  },
  
  
  prop.genuine = function(detection.intervals, model.obj)
  {
    #' Returns proportion of intervals which contain a changepoint
    #'
    #'@param detection.intervals
    #'@param model.obj
    
    cpt.locs <- model.obj$cpt.locs
    
    if (nrow(detection.intervals) == 0) return(NaN)
    
    correct.detections <- 0
    
    for (ii in 1:nrow(detection.intervals))
    {
      if(any(cpt.locs %in% detection.intervals[ii,1]:detection.intervals[ii,2])) correct.detections <- correct.detections + 1
    }
    
    return(correct.detections / nrow(detection.intervals))
  }, 
  
    
  num.genuine = function(detection.intervals, model.obj)
  {
    #' Returns the number of intervals which contain a changepoint
    #'
    #'@param nsp.obj
    #'@param model.obj
    
    cpt.locs <- model.obj$cpt.locs
    
    if (nrow(detection.intervals) == 0) return(0)
    
    correct.detections <- 0
    
    for (ii in 1:nrow(detection.intervals))
    {
      if(any(cpt.locs %in% detection.intervals[ii,1]:detection.intervals[ii,2])) correct.detections <- correct.detections + 1
    }
    
    return(correct.detections)
  },
  
  
  average.length = function(detection.intervals, model.obj)
  {
    #' Returns the number of intervals which contain a changepoint
    #'
    #'@param detection.intervals
    #'@param model.obj
    
    if (nrow(detection.intervals) == 0) return(NaN)
    
    mean(apply(detection.intervals,1,diff))
  }
  )


update.metrics <- function(cur.metrics, metrics, detection.intervals, model.obj)
{
  #'
  #'
  #'@param cur.metrics vector or matrix
  #'@param metrics list of functions
  #'@param nsp.obj 
  #'@param model.obj 
  
  new.metrics <- c()
  for (metric in metrics) new.metrics <- c(new.metrics, metric(detection.intervals, model.obj))
  cbind(cur.metrics, new.metrics)
}



#--------------------------
#
# Numerical study
#
#--------------------------

## compare with MQS: https://doi.org/10.1080/01621459.2020.1859380

if (!require(mqs)) devtools::install_github("ljvanegas/mqs")


set.seed(42)
n.reps <- 100
pb <- txtProgressBar(min = 0, max = n.reps, style = 3)
cpt.models.metrics <- list()


for (changepoint.model in changepoint.models)
{
  
  multiresolution.metrics <- c()
  mqs.metrics <- c()
  runs.metrics <- c()
  multi.res.metrics <- c()
  corrected.multi.res.metrics <- c()

  for (jj in 1:n.reps)
  {
    cat("\r", changepoint.model$name)
    
    ee <- changepoint.model$simulate()
    nn <- length(ee)
    
    rnsp.obj <- nsp_robust(ee, method = "multiresolution", thresh = thresh.kab.bern(nn))$intervals[,1:2]
    multiresolution.metrics <- update.metrics(multiresolution.metrics, metrics, rnsp.obj, changepoint.model)
    
    rnsp.obj <- nsp_robust(ee, method = "runs", thresh = runs.thresh(nn))$intervals[,1:2]
    runs.metrics <- update.metrics(runs.metrics, metrics, rnsp.obj, changepoint.model)

    mqs.obj <- mqs::mqse(ee, alpha = .1)$confInt
    mqs.metrics <- update.metrics(mqs.metrics, metrics, mqs.obj, changepoint.model)

    rnsp.obj <- try(
      nsp_robust(ee, method = "multi.res", thresh = thresh.kab.bern(nn))$intervals[,1:2]
    )
    if ("try-error" %in% class(rnsp.obj)) rnsp.obj <- matrix(0,0,2)
    multi.res.metrics <- update.metrics(multi.res.metrics, metrics, rnsp.obj, changepoint.model)

    rnsp.obj <- try(
      nsp_robust(ee, method = "corrected.multi.res", thresh = thresh.dumbgen.T0(nn))$intervals[,1:2]
    )
    if ("try-error" %in% class(rnsp.obj)) rnsp.obj <- matrix(0,0,2)
    corrected.multi.res.metrics <- update.metrics(corrected.multi.res.metrics, metrics, rnsp.obj, changepoint.model)
    
    setTxtProgressBar(pb, jj)
  }
  
  metrics.averages <- data.frame(
    multiresolution = apply(multiresolution.metrics, 1, function(ii) mean(ii, na.rm = TRUE)),
    runs = apply(runs.metrics, 1, function(ii) mean(ii, na.rm = TRUE)),
    multi.res = apply(multi.res.metrics, 1, function(ii) mean(ii, na.rm = TRUE)),
    corrected.multi.res = apply(corrected.multi.res.metrics, 1, function(ii) mean(ii, na.rm = TRUE))
  )
  
  rownames(metrics.averages) <- names(metrics)
  
  cpt.models.metrics[[changepoint.model$name]] <- metrics.averages 
}

cpt.models.metrics
