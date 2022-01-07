## Implements robust NSP variants to be compared 
## 
##
## IMPORTANT! 
##
## Before running the code load the file NSP_for+Github***.R from: https://github.com/pfryz/nsp/
##


#-------------------------------------------------------------------
#
# Generic functions to test for change points using inverted tests
#
#-------------------------------------------------------------------


get.bounds.non.decreasing <- function(yy, bound.construction, thresh)
{
  #' Construct confidence set assuming regression function is non-decreasing
  #'
  #'@param yy data vector
  #'@param bound.construiction function for constructing lower bound, see `lower-bound-constructions`
  #'@param thresh threshold for corresponding method
  
  ll <- bound.construction(yy, thresh)
  uu <- bound.construction(-rev(yy), thresh)
  
  list(ll = ll, uu = -rev(uu))
}


get.bounds.non.increasing <- function(yy, bound.construction, thresh)
{
  #' Construct confidence set assuming regression function is non-decreasing
  #'
  #'@param yy data vector
  #'@param bound.construiction function for constructing lower bound, see `lower-bound-constructions`
  #'@param thresh threshold for corresponding method
  
  ll <- bound.construction(rev(yy), thresh)
  uu <- bound.construction(-yy, thresh)
  
  list(ll = rev(ll), uu = -uu)
}



confidence.set.test <- function(yy, method, thresh)
{
  #' Test for change in piecewise constant mean using confidence set
  #'
  #'@param yy data vector
  #'@param method which test to invert, choices are: "runs", "multi.res", and "corrected.multi.res"
  #'@param thresh corresponding threshold
  
  if (method == "runs") bound.construction <- runs.lower.bound
  if (method == "corrected.multi.res") bound.construction <- corrected.multi.res.lower.bound
  if (method == "corrected.multi.res.dyadic") bound.construction <- corrected.multi.res.dyadic.lower.bound
  if (method == "multi.res") bound.construction <- multi.res.lower.bound
  if (method == "multi.res.dyadic") bound.construction <- multi.res.dyadic.lower.bound
  
  bb.non.decreasing <- get.bounds.non.decreasing(yy, bound.construction, thresh)
  non.decreasing <- (min(bb.non.decreasing$uu) >= max(bb.non.decreasing$ll))
  
  bb.non.increasing <- get.bounds.non.increasing(yy, bound.construction, thresh)
  non.increasing <- (min(bb.non.increasing$uu) >= max(bb.non.decreasing$ll))
  
  as.numeric(!(non.increasing & non.decreasing))
}


#-------------------------------------------------------------------
#
# Extending the RNSP function to accept our function
#
#-------------------------------------------------------------------



rnsp <- nsp_robust <- function(x, method = "multiresolution", M = 1000, thresh = NULL, alpha = 0.1, thresh.type = "bern", thresh.sim.times = 100, deg = 0, max.length = 3000, overlap = FALSE, zeros = TRUE)
{
  #'
  #'
  #'@param x vector, the data 
  #'@param method string, one of "multiresolution", "runs", "multi.res", or "corrected.multi.res"
  #'@param M int, 
  #'@param thresh float,
  #'@param alpha float, 
  #'@param thresh.type string
  #'@param thresh.sim.times 
  #'@param deg 
  #'@param max.length 
  #'@param overlap
  #'@param zeros
  
  n <- length(x)
  
  if (is.null(thresh)) stop("must set threshold manually")
  
  
  res <- iter_random_checks_scan_signed(c(1, n), x, M, thresh, max.length, overlap, 0, deg, zeros, method)
  
  intervals <- data.frame(t(order_chron(res)))
  colnames(intervals) <- c("starts", "ends", "values")
  
  list(intervals=intervals, threshold.used=thresh)
  
}


thresh_kab_bern <- function(n, alpha = 0.1) {
  
  an <- sqrt(2 * log(n*log(n)^(-1/2)))
  
  tau <- -log(-1/(2*0.2740311) * log(1 - alpha))
  
  an + tau / an
  
}


iter_random_checks_scan_signed <- function(ind, x, M, thresh, max.length, overlap = FALSE, buffer = 0, deg = 0, zeros = TRUE, method) {
  
  s <- ind[1]
  e <- ind[2]
  n <- e - s + 1
  x <- x[s:e]
  
  if (n > 1) {
    
    next.int <- random_checks_scan_signed_2stage(c(1,n), x, M, thresh, max.length, deg, zeros, method)
    
    if (!is.na(next.int$selected.ind))  {
      
      if (!overlap) {
        
        if (next.int$selected.val[1,1]-buffer >= 2) left <- iter_random_checks_scan_signed(c(1, next.int$selected.val[1,1]-buffer), x, M, thresh, max.length, overlap, buffer, deg, zeros, method) else left <- matrix(NA, 3, 0)
        if (n - next.int$selected.val[1,2]-buffer >= 1) {
          
          right <- iter_random_checks_scan_signed(c(next.int$selected.val[1,2]+buffer, n), x, M, thresh, max.length, overlap, buffer, deg, zeros, method)
          if (dim(right)[2]) right <- right + c(rep(next.int$selected.val[1,2]-1+buffer, 2), 0)
          
        }
        else right <- matrix(NA, 3, 0)
      }
      
      else {
        
        
        if (floor(mean(next.int$selected.val[1,1:2]))-buffer >= 2) left <- iter_random_checks_scan_signed(c(1, floor(mean(next.int$selected.val[1,1:2]))-buffer), x, M, thresh, max.length, overlap, buffer, deg, zeros, method) else left <- matrix(NA, 3, 0)
        if (n - floor(mean(next.int$selected.val[1,1:2]))-buffer >= 2) {
          right <- iter_random_checks_scan_signed(c(floor(mean(next.int$selected.val[1,1:2]))+1+buffer, n), x, M, thresh, max.length, overlap, buffer, deg, zeros, method)
          if (dim(right)[2]) right <- right + c(rep(floor(mean(next.int$selected.val[1,1:2]))+buffer, 2), 0)
          
        }
        else right <- matrix(NA, 3, 0)
        
        
      }
      
      
      return(cbind(t(next.int$selected.val), left, right))
      
      
    }
    
    else(return(matrix(NA, 3, 0)))
    
    
  }
  
  else(return(matrix(NA, 3, 0)))
  
}


random_checks_scan_signed_2stage <- function(ind, x, M, thresh, max.length, deg = 0, zeros = TRUE, method) {
  
  s1 <- random_checks_scan_signed(ind, x, M, thresh, max.length, deg, zeros, method)
  
  if (!is.na(s1$selected.ind)) {
    
    s <- s1$selected.val[1,1] + ind[1] - 1
    e <- s1$selected.val[1,2] + ind[1] - 1
    
    s2 <- random_checks_scan_signed(c(s,e), x, M, thresh, max.length, deg, zeros, method)
    
    if (!is.na(s2$selected.ind)) {
      
      replacement <- s2$selected.val
      replacement[1,1:2] <- replacement[1,1:2] + s - ind[1]
      s1$selected.val <- replacement
      
    }
    
  }
  
  s1	
  
}


random_checks_scan_signed <- function(ind, x, M, thresh, max.length, deg = 0, zeros = TRUE, method) {
  
  s <- ind[1]
  e <- ind[2]
  n <- e - s + 1
  
  if (n > 1) {
    
    x <- x[s:e]
    
    M <- min(M, (n-1)*n/2)
    
    ind <- grid_intervals_sorted(n, M)      # record where this function can be found
    
    M <- dim(ind)[2]
    
    res <- matrix(0, M, 3)
    
    res[,1:2] <- t(ind)
    
    zero.check <- TRUE
    j <- 1
    
    while (zero.check && (j <= M)) {
      
      res[j,3] <- check_interval_robust_sup_norm(res[j,1:2], x, thresh, max.length, deg, zeros, method)
      zero.check <- (res[j,3] == 0)
      j <- j + 1
      
    }
    
    if (zero.check) {
      
      selected.ind <- NA
      selected.val <- matrix(0, 0, 3)
      
    }
    
    else {
      
      selected.ind <- j-1
      selected.val <- res[selected.ind,,drop=FALSE]
      
    }
    
    
  }
  
  else {
    
    selected.val <- matrix(0, 0, 3)
    selected.ind <- NA
    M <- 0
    
  }
  
  list(selected.ind = selected.ind, selected.val = selected.val, M.eff=max(1, M))
  
}


check_interval_robust_sup_norm <- function(ind, x, thresh, max.length, deg, zeros = TRUE, method) {
  
  if (ind[2] - ind[1] + 1 > max.length) return(0) 
  
  if (method == "multiresolution") 
  {
    sup.nm <- robust_sup_norm_rv_wwozr(x[ind[1]:ind[2]], zeros)$sup.norm
    sup.nm <- sup.nm * (sup.nm > thresh) 
    
  } else {
    
    sup.nm <- confidence.set.test(x[ind[1]:ind[2]], method, thresh)
  }
  
  return(sup.nm)
  
}



robust_sup_norm_rv_wwozr <- function(data, zeros = TRUE) {
  
  # rv stands for "rank version"
  # wwozr stands for "with or without zeros"
  # this version uses rank rather than order
  # the user can set zeros=FALSE if they are sure there are no masses at any local medians (e.g. if the distribution of the data is continuous)
  # requires plyr::mapvalues !!
  
  n <- length(data)
  
  if (n <= 1) {
    
    sup.norm <- 0
    d <- (zeros + 1) * n + 1
    left.max <- right.max <- comb.max <- rep(0, d)
    centre <- data
    sig <- data		
    
  } else {	
    
    data.r <- rank(data, ties.method = "min")
    data.r.s <- sort(unique(data.r))
    how.many <- length(data.r.s)
    data.sr <- plyr::mapvalues(data.r, data.r.s, 1:how.many)
    
    if (zeros) {
      
      pseudo.data <- matrix(1, n, 2*how.many+1)
      
      for (i in 1:n) {
        pseudo.data[i, (2*data.sr[i]+1):(2*how.many+1)] <- -1
        pseudo.data[i, 2*data.sr[i]] <- 0
      }
      
    }
    else {
      
      pseudo.data <- matrix(1, n, how.many+1)
      
      for (i in 1:n) pseudo.data[i, (data.sr[i]+1):(how.many+1)] <- -1
      
    }
    
    left.cumsum	<- apply(pseudo.data, 2, cumsum)
    left.max <- apply(abs(left.cumsum) / sqrt(1:n), 2, max)
    
    right.cumsum <- apply(pseudo.data[n:1,], 2, cumsum)
    right.max <- apply(abs(right.cumsum) / sqrt(1:n), 2, max) 
    
    comb.max <- pmax(left.max, right.max)
    
    which.comb.min <- which.min(comb.max)
    
    sup.norm <- comb.max[which.comb.min]
    
    if (zeros) {
      
      which.comb.trim.min <- which.min(comb.max[2:(2*how.many)])
      centre.ind.up <- data.r.s[floor((which.comb.trim.min + 1)/2)]
      centre.ind.down <- data.r.s[ceiling((which.comb.trim.min + 1)/2)]
      centre <- 1/2 * (data[which(data.r == centre.ind.up)][1] + data[which(data.r == centre.ind.down)][1])
      
    }
    else {
      
      if (how.many > 1) {
        which.comb.trim.min <- which.min(comb.max[2:how.many])
        centre.ind.up <- data.r.s[floor(which.comb.trim.min + 1/2)]
        centre.ind.down <- data.r.s[ceiling(which.comb.trim.min + 1/2)]
        centre <- 1/2 * (data[which(data.r == centre.ind.up)][1] + data[which(data.r == centre.ind.down)][1])
        
      }
      else centre <- data[1]
      
      
    }
    
    sig <- rep(centre, n)
    
  }
  
  return(list(sup.norm = sup.norm, left.max = left.max, right.max = right.max, comb.max = comb.max, centre = centre, sig = sig))
  
}
