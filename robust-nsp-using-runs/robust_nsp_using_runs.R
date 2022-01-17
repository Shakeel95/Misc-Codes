
## IMPORTANT! 
##
## Before runningthe code load the file NSP_for+Github***.R from: 
##
## https://github.com/pfryz/nsp/
##


#=================
#
# Runs test
#
#=================

get.bounds.non.increasing <- function(yy, max.runs)
{
  #' Get upper and lower bounds for admissible class of non-increasing function
  #'
  #'@param yy vector
  #'@param max.runs int
  #'
  #'@references DOI: 10.1214/aos/996986501
  
  nn <- length(yy)
  
  uu <- numeric(nn)
  uu[1:max.runs] <- Inf
  
  for (ii in (max.runs+1):nn)
  {
    if (all(yy[(ii-max.runs-1):(ii-1)] > uu[(ii-max.runs-1):(ii-1)])) return(list(shape.check = FALSE))
    uu[ii] <- min(c(uu[ii-1], max(yy[(ii-max.runs):ii])))
  }
  
  
  ll <- numeric(nn)
  ll[(nn + 1 - max.runs):nn] <- -Inf
  
  for (ii in (nn - max.runs):1)
  {
    ll[ii] <- max(ll[ii+1], min(yy[ii:(ii+max.runs)]))
  }
  
  return(list(shape.check = TRUE, uu = uu, ll = ll))
  
}


get.bounds.non.decreasing <- function(yy, max.runs)
{
  #' Get upper and lower bounds for admissible class of non-increasing function
  #'
  #'@param yy vector
  #'@param max.runs int
  #'
  #'@references DOI: 10.1214/aos/996986501
  
  nn <- length(yy)
  
  ll <- numeric(nn)
  ll[1:max.runs] <- -Inf
  
  for (ii in (max.runs+1):nn)
  {
    if (all(yy[(ii-max.runs-1):(ii-1)] < ll[(ii-max.runs-1):(ii-1)])) return(list(shape.check = FALSE)) # checking for a change in monotonicity 
    ll[ii] <- max(c(ll[ii-1], min(yy[(ii-max.runs):ii])))
  }
  
  uu <- numeric(nn)
  uu[(nn + 1 - max.runs):nn] <- Inf
  
  for (ii in (nn - max.runs):1)
  {
    uu[ii] <- min(uu[ii+1], max(yy[ii:(ii+max.runs)]))
  }
  
  return(list(shape.check = TRUE, uu = uu, ll = ll))
  
}


check.if.p0.seperable <- function(ll, uu)
{
  #' Check if bounds are seperable by degree 0 polynomial  
  #'
  #'@param ll vector
  #'@param uu vector
  
  return(
    ifelse(min(uu) >= max(ll), FALSE, TRUE)
  )
}


check.if.p1.seperable <- function(ll, uu)
{
  #'Check if bounds are seperable by degree 1 polynomial  
  #'
  #'@param ll vector
  #'@param uu vector
  
  nn <- length(ll)
  
  ll.inf.ind <- which(abs(ll) == Inf)
  ll.hull <- cbind(setdiff(1:nn, ll.inf.ind), ll[-ll.inf.ind])
  ll.hull.ind <- chull(ll.hull)
  
  uu.inf.ind <- which(abs(uu) == Inf)
  uu.hull <- cbind(setdiff(1:nn, uu.inf.ind), uu[-uu.inf.ind])
  uu.hull.ind <- chull(uu.hull)
  
  ll.hull <- ll.hull[ll.hull.ind,]
  uu.hull <- uu.hull[uu.hull.ind,]
  
  if (is.null(dim(ll.hull))|is.null(dim(uu.hull)))
  {
    
    return(FALSE) 
    
  } else if ((nrow(ll.hull)>2)&(nrow(uu.hull)>2)) {
    
    return(convex.hulls.intersect(uu.hull, ll.hull))
    
  } else {
    
    return(FALSE)
    
  } 
}


check.if.pq.seperable <- function(ll, uu)
{
  #' Check if bounds are seperable by degree q > 1 polynomial  
  #'
  #'@param ll vector
  #'@param uu vector
  
  return(NULL)
}

runs.thresh <- function(n, a = .1)
{
  #' Threshold for coverage control of runs test 
  #'
  #'@param n int 
  #'@param a float
  
  log2(n-1) + (1/log(2)) * log(1/log(1/(1-a)))
  
}


runs.test <- function(yy, thresh, deg)
{
  #' Test for changepoint using the runs test
  #'
  #'@param yy vector 
  #'@param thresh float 
  
  max.runs <- floor(thresh)
  
  if (max.runs >= length(yy)) return(0)
  
  
  no.fits <- TRUE
  
  test.non.increasing <- get.bounds.non.increasing(yy, max.runs)
  
  if (test.non.increasing$shape.check) 
  {
    
    no.fits <- FALSE
    
    if (deg == 0)
    {
      
      linearly.seperable <- check.if.p0.seperable(test.non.increasing$ll, test.non.increasing$uu)
      if (linearly.seperable) return(1) 
      
    } else if (deg == 1) {
      
      test.linearly.seperable <- check.if.p1.seperable(test.non.increasing$ll, test.non.increasing$uu)
      if (test.linearly.seperable) return(1) 
      
    } else {
      
      return(runs.test.experimental(yy,thresh))
      
    }
    
  }
  
  
  test.non.decreasing <- get.bounds.non.decreasing(yy, max.runs)
  
  if (test.non.decreasing$shape.check)
  {
    
    no.fits <- FALSE
    
    if (deg == 0)
    {
      
      test.linearly.seperable <- check.if.p0.seperable(test.non.decreasing$ll, test.non.decreasing$uu)
      if (test.linearly.seperable) return(1)
      
    } else if (deg == 1) {
      
      test.linearly.seperable <- check.if.p1.seperable(test.non.decreasing$ll, test.non.decreasing$uu)
      if (test.linearly.seperable) return(1) 
      
    } else {
      
      return(runs.test.experimental(yy,thresh))
      
    }
  }
  
  
  if (no.fits) return(1)
  
  return(0)
  
}


runs.test.experimental <- function(yy, thresh) 
{
  #' Test for changepoint using the runs test
  #'
  #'@param yy vector 
  #'@param thresh float 
  
  max.runs <- floor(thresh)
  
  if (max.runs >= length(yy)) return(0)
  
  test.non.increasing <- get.bounds.non.increasing(yy, max.runs)
  
  test.non.decreasing <- get.bounds.non.decreasing(yy, max.runs)
  
  ifelse(!any(test.non.decreasing$shape.check, test.non.increasing$shape.check), 1, 0)
}


#===============================
#
# Extending the RNSP function
#
#===============================



rnsp <- nsp_robust <- function(x, method = "multiresolution", M = 1000, thresh = NULL, alpha = 0.1, thresh.type = "bern", thresh.sim.times = 100, deg = 0, max.length = 3000, overlap = FALSE, zeros = TRUE)
{
  #'
  #'
  #'@param x vector, the data 
  #'@param method string, one of "multiresolution" or "runs"
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
  
  if (is.null(thresh)) 
  {
    if (method == "runs") {
      thresh <- runs.thresh(n, alpha)
    } else {
      if (thresh.type == "bern") thresh <- thresh_kab_bern(n, alpha)
      if (thresh.type == "gauss") thresh <- 0.9 * thresh_kab(n, alpha)
      if (thresh.type == "sim")
      {
        msn <- sim_multi_sup_norm(n, thresh.sim.times)
        thresh <- as.numeric(quantile(msn, 1-alpha))
      } 
    }
  }
  
  
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
    
  } else if (method == "runs") {
    
    sup.nm <- runs.test(x[ind[1]:ind[2]], thresh, deg)
    
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


multi_sup_norm <- function(eps) {
  
  
  n <- length(eps)
  
  eps.cum <- c(0, cumsum(eps))
  
  max.use <- max.abs.use <- 0
  
  for (i in 0:(n-1)) for (j in (i+1):n) {
    
    scan.stat <- (eps.cum[j+1] - eps.cum[i+1]) / sqrt(j-i)
    
    if (abs(scan.stat) > max.abs.use) max.abs.use <- abs(scan.stat)
    
  }
  
  max.abs.use
  
}
