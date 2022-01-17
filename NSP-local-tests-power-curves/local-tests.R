
#--------------------------------
#
# Local test from NSP paper
#
#--------------------------------

rnsp.fryz <- function(xx, method, thresh = NULL)
{
  #' Local test for change points from RNSP paper
  #'
  #'@param xx vector, data
  #'@param thresh float, threshold
  #'@param method string,
  #'
  #'@references https://arxiv.org/abs/2109.02487
  #'@references https://github.com/pfryz/nsp

  nn <- length(xx)
  if (is.null(thresh)) thresh <- thresh.kab.bern(nn)
  
  xx.r <- rank(xx, ties.method = "min")
  xx.r.s <- xx.r |> unique() |> sort()
  nn.unique <- length(xx.r.s)
  
  xx.sr <- plyr::mapvalues(xx.r, xx.r.s, 1:nn.unique)
  
  xx.pseudo <- matrix(1, nn, 2*nn.unique+1)
  
  for (ii in 1:nn) 
  {
    xx.pseudo[ii, (2*xx.sr[ii]+1):(2*nn.unique+1)] <- -1
    xx.pseudo[ii, 2*xx.sr[ii]] <- 0
  }
  
  if (method == "anchored")
  {
    
    left.cumsum	<- apply(xx.pseudo, 2, cumsum)
    left.max <- apply(abs(left.cumsum) / sqrt(1:nn), 2, max)
    
    right.cumsum <- apply(xx.pseudo[nn:1,], 2, cumsum)
    right.max <- apply(abs(right.cumsum) / sqrt(1:nn), 2, max) 
    
    comb.max <- pmax(left.max, right.max)
    which.comb.min <- which.min(comb.max)
    
    sup.norm <- comb.max[which.comb.min]
    
  } else if (method == "complete") {
    
    sup.norm <- apply(xx.pseudo, 2, multi.sup.norm) |> min()
    
  }
  
  (sup.norm > thresh) |> as.numeric()
  
}


multi.sup.norm <- function(xx)
{
  #' Compute multi-resolution norm of a vector
  #'
  #'@param xx vector
  
  nn <- length(xx)
  xx.cumsum <- c(0,cumsum(xx))
  
  loc.sum <- 0
  max.sum <- 0
  
  for (ss in 1:nn) for (ee in (ss+1):(nn+1))
  {
    loc.sum <- abs(xx.cumsum[ee]-xx.cumsum[ss])/sqrt(ee-ss)
    if (loc.sum > max.sum) max.sum <- loc.sum
  }
  
  max.sum
}


#------------------------------------------------
#
# Local test by inverting multi-resolution test
#
#------------------------------------------------


rnsp.inv <- function(xx, method, thresh)
{
  #' Local test for change points based on inverting multi-resolution sign test
  #'
  #'@param xx vector, data
  #'@param thresh float, threshold
  #'@param method string,
  #'
  #'@references https://doi.org/10.1198/1061860043506
  
  nn <- length(xx)
  
  if (method == "mr.inv")
  {
    bb.non.decreasing <- get.bounds.non.decreasing(xx, multi.res.lower.bound, thresh)
    non.decreasing <- (min(bb.non.decreasing$uu) >= max(bb.non.decreasing$ll))
    
    bb.non.increasing <- get.bounds.non.increasing(xx, multi.res.lower.bound, thresh)
    non.increasing <- (min(bb.non.increasing$uu) >= max(bb.non.decreasing$ll))
    
    return(as.numeric(!(non.increasing & non.decreasing)))
    
  } else if (method == "mr.inv.corrected") {
        
    bb.non.decreasing <- get.bounds.non.decreasing(xx, corrected.multi.res.lower.bound, thresh)
    non.decreasing <- (min(bb.non.decreasing$uu) >= max(bb.non.decreasing$ll))
    
    bb.non.increasing <- get.bounds.non.increasing(xx, corrected.multi.res.lower.bound, thresh)
    non.increasing <- (min(bb.non.increasing$uu) >= max(bb.non.decreasing$ll))
    
    return(as.numeric(!(non.increasing & non.decreasing)))
  }
}


get.bounds.non.decreasing <- function(xx, bound.construction, thresh)
{
  #' Construct confidence set assuming regression function is non-decreasing
  #'
  #'@param xx data vector
  #'@param bound.construiction function for constructing lower bound, see `lower-bound-constructions`
  #'@param thresh threshold for corresponding method
  
  ll <- bound.construction(xx, thresh)
  uu <- bound.construction(-rev(xx), thresh)
  
  list(ll = ll, uu = -rev(uu))
}


get.bounds.non.increasing <- function(xx, bound.construction, thresh)
{
  #' Construct confidence set assuming regression function is non-decreasing
  #'
  #'@param xx data vector
  #'@param bound.construiction function for constructing lower bound, see `lower-bound-constructions`
  #'@param thresh threshold for corresponding method
  
  ll <- bound.construction(rev(xx), thresh)
  uu <- bound.construction(-xx, thresh)
  
  list(ll = rev(ll), uu = -uu)
}

