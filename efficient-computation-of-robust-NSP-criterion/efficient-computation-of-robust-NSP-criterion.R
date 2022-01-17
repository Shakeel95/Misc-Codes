
#==================
#
# Functions
#
#==================


sign.multi.sup.norm <- function(x)
{
  #' Sign multi-resolution sup norm of a vector
  #'
  #'@param x vector 
  
  x <- sign(x)
  cumsum.x <- c(0, cumsum(x))
  n <- length(x)
  
  loc.max <- 0 
  run.max <- 0 
  
  for (s in 1:n) for (e in (s+1):(n+1)) 
  {
    loc.max <- abs(cumsum.x[e] - cumsum.x[s]) / sqrt(e-s)
    
    if (loc.max > run.max) run.max <- loc.max
  }
  
  return(run.max)
  
}


get.candidate.fits <- function(y)
{
  #' Get candidate fits as decribed in RNSP paper
  #'
  #'@param y vector 
  #'@references https://stats.lse.ac.uk/fryzlewicz/nsp/rnsp.pdf
  
  sort(c(y, zoo::rollmean(y,2), min(y)-1, max(y)+1))
  
}


brute.force.search <- function(yy, candidate.fits)
{
  #' Find best candidtae fit using brute force search
  #'
  #'@param yy vector
  #'@param candidate.fits vector
  
  
  nn <- length(yy)
  
  res <- c()
  
  for (loc.fit in candidate.fits) res <- c(res, sign.multi.sup.norm(yy - rep(loc.fit,nn)))
  
  return(list(all.fits = res, norm.min = min(res)))
}


golden.section.search <- function(yy, candidate.fits, min.width = 10)
{
  #' Find best candidtae fit using golden section search
  #'
  #'@param yy vector
  #'@param candidate.fits vector
  #'@param min.width int
  
  mm <- length(candidate.fits)
  nn <- length(yy)
  
  RR <- 0.61803399
  CC <- 1 - RR 
  
  x1 <- 1 
  x2 <- floor(mm/2) 
  x3 <- floor(mm/2 + CC * mm/2)
  x4 <- mm 
  
  all.sections <- list()
  all.sections[[1]] <- c(x1,x2,x3,x4)
  
  f2 <- sign.multi.sup.norm(yy - rep(candidate.fits[x2],nn))
  f3 <- sign.multi.sup.norm(yy - rep(candidate.fits[x3],nn))
  
  while(abs(x4 - x1) > min.width)
  {
    if (f2 > f3)
    {
      
      x1 <- x2
      x2 <- x3
      x3 <- floor(RR * x3 + CC * x4) 
      x4 <- x4
      
      f2 <- f3  
      f3 <- sign.multi.sup.norm(yy - rep(candidate.fits[x3],nn))
      
    } else {
      
      x4 <- x3 
      x3 <- x2  
      x2 <- floor(RR * x2 + CC * x1)
      x1 <- x1
      
      f3 <- f2
      f2 <- sign.multi.sup.norm(yy - rep(candidate.fits[x2],nn))
      
    }
    
    all.sections[[length(all.sections)+1]] <- c(x1,x2,x3,x4)
    
  }
  
  
  res <- c()
  for (ii in x1:x4) res <- c(res, sign.multi.sup.norm(y - rep(candidate.fits[ii],nn)))
  
  return(list(all.sections = all.sections, norm.min = min(res)))
  
}



corrected.golden.section.search <- function(x1, x4, yy, candidate.fits, min.width = 10)
{
  #' Find best candidtae fit using golden section search with recursion on flat regions 
  #'
  #'@param yy vector
  #'@param candidate.fits vector
  #'@param min.width int
  
  nn <- length(yy)
  RR <- 0.61803399
  CC <- 1 - RR 
  
  min.width.condition <- (x4-x1) > min.width
  flat.interval.conditioin <- FALSE
  
  if (min.width.condition)
  {
    res <- c()
    for (ii in x1:x4) res <- c(res, sign.multi.sup.norm(yy - rep(candidate.fits[ii],nn)))
    return(min(res))
  }
  
  x2 <- floor(x4/2) 
  x3 <- floor(x4/2 + CC * x4/2)
  
  f2 <- sign.multi.sup.norm(yy - rep(candidate.fits[x2],nn))
  f3 <- sign.multi.sup.norm(yy - rep(candidate.fits[x3],nn))
  
  while(min.width.condition | flat.interval.conditioin)
  {
    if (f2 > f3)
    {
      
      x1 <- x2
      x2 <- x3
      x3 <- floor(RR * x3 + CC * x4) 
      x4 <- x4
      
      f2 <- f3  
      f3 <- sign.multi.sup.norm(yy - rep(candidate.fits[x3],nn))
      
      min.width.condition <- (x4 - x1) < min.width
      
    } else if (f3 > f2) {
      
      x4 <- x3 
      x3 <- x2  
      x2 <- floor(RR * x2 + CC * x1)
      x1 <- x1
      
      f3 <- f2
      f2 <- sign.multi.sup.norm(yy - rep(candidate.fits[x2],nn))
      
      min.width.condition <- (x4 - x1) < min.width
      
    } else {
      
      flat.interval.conditioin <- TRUE
      
    }
  }
  
  
  if (flat.interval.conditioin) return(min(
    corrected.golden.section.search(x1, floor((x2+x3)/2), yy, candidate.fits, min.width),
    corrected.golden.section.search(ceiling((x2+x3)/2), x4, yy, candidate.fits, min.width)
  ))
  
  if (min.width.condition)
  {
    res <- c()
    for (ii in x1:x4) res <- c(res, sign.multi.sup.norm(yy - rep(candidate.fits[ii],nn)))
    return(min(res))
  }
  
}


#========================
#
# Numerical experiments
#
#========================

set.seed(42)
nn <- 300
y <- rep(10, nn) + rcauchy(nn)
candidates <- get.candidate.fits(y)


## Accuracy

"""
The naive implementation of golden section search does not always find the best candidate fit. 
"""


res.brute.force <- brute.force.search(y, candidates)

res.golden.section <- golden.section.search(y, candidates)


for (ii in res.golden.section$all.sections)
{
  plot(res.brute.force$all.fits, type = "l")
  abline(v = ii, col = "red", lty = 2)
}


corrected.golden.section.search(
  x1 = 1, 
  x4 = length(candidates), 
  yy = y, 
  candidate.fits = candidates, 
  min.width = 30)


## Speed

"""
Recursive function calls in R are very slow; the adjusted golden section method can be slower than brute force. 
The naive implementation is (as expected) significantly faster than brute force search
"""

microbenchmark::microbenchmark(
  brute.force.search(y, candidates), 
  times = 1 
)

microbenchmark::microbenchmark(
  golden.section.search(y, candidates), 
  times = 1 
)

microbenchmark::microbenchmark(
  corrected.golden.section.search(
    x1 = 1, 
    x4 = length(candidates), 
    yy = y, 
    candidate.fits = candidates, 
    min.width = 10), 
  times = 10 
)

