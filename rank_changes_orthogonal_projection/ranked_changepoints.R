library("docstring")

#---------------------
#
# Contrast functions
#
#---------------------


lin.contrast <- function(s, b, e){
  
  #' pcwsLin contrast function
  #' 
  #' 
  #' @param b int, potential changepoint loaction 
  #' @param s int, start
  #' @param e int, end 
  
  l <- e-s
  A <- sqrt(6/(l*(l**2-1)*(1+(e-b+1)*(b-s)+(e-b)*(b-s-1))))
  B <- sqrt(((e-b+1)*(e-b))/((b-s-1)*(b-s)))
  
  contrast <- numeric(l)
  
  contrast[(s+1):b] <- (A*B)*((3*(b-s)+(e-b)-1)*c((s+1):b) - (b*(e-s-1)+2*(s+1)*(b-s)))
  
  contrast[(b+1):e] <- (-1)*(A/B)*((3*(e-b)+(b-s)+1)*c((b+1):e) - (b*(e-s-1)+2*e*(e-b+1)))
  
  return(contrast)
  
}


const.contrast <- function(s, b, e){
  
  #' pcwsConst contrast function
  #' 
  #' 
  #' @param b int, potential changepoint loaction 
  #' @param s int, start
  #' @param e int, end 
  
  l <- e-s

  contrast <- numeric(l)
  
  contrast[(s+1):b] <- sqrt((e-b)/(l*(b-s)))
  
  contrast[(b+1):e] <- -sqrt((b-s)/(l*(e-b)))
  
  return(contrast)
  
}

#-----------------------
#
# Projection functions
#
#-----------------------


orth.space.proj <- function(x,v) {
  
  #' Project x onto the orthogonal compliment of the space spanned by v 
  #' 
  #' @param x vector
  #' @param v vector 
  
  (diag(length(v)) - outer(v,v)) %*% x
  
}



#------------------------------------------------------
#
# Functions for discovering changepoints by importance
#
#------------------------------------------------------


ranked.pcws.const <- function(x,N) {
  
  #' Ranked changepoint in pcsw constant mean model
  #' 
  #' @param x vector, observations from piecewise constant model
  #' @param N int, max changepoints to detect
  
  
  n <- length(x)
  
  cpt.locs <- list()
  cpt.locs[[1]] <- 2 
  
  x.proj <- list()
  x.proj[[1]] <- x
  
  
  while(length(cpt.locs) < N+1) {
    
    contras.vals <- sapply(2:(n-1), function(b) abs(sum(const.contrast(1,b,n)*x)))
    
    b <- which.max(contras.vals[-c(unlist(cpt.locs)-1)]) + 1
    
    x <- orth.space.proj(x, const.contrast(1,b,n))
    
    
    x.proj[[length(x.proj)+1]] <- x 
    
    cpt.locs[[length(cpt.locs)+1]] <- b
    
  }
  
  
  return(list(cpt.locs = cpt.locs, x.proj = x.proj))
  
}


ranked.pcws.lin <- function(x,N) {
  
  #' ....
  #' 
  #' @param x vector, observations from piecewise constant model
  #' @param N int, max changepoints to detect
  
  
  n <- length(x)
  
  cpt.locs <- list()
  cpt.locs[[1]] <- 1 
  
  x.proj <- list()
  x.proj[[1]] <- x
  
  
  while(length(cpt.locs) < N+1) {
    
    contras.vals <- sapply(2:(n-1), function(b) abs(sum(const.contrast(1,b,n)*x)))
    
    b <- which.max(contras.vals) + 1
    
    x <- orth.space.proj(x, const.contrast(1,b,n))
    
    
    x.proj[[length(x.proj)+1]] <- x 
    
    cpt.locs[[length(cpt.locs)+1]] <- b
    
  }
  
  
  return(list(cpt.locs = cpt.locs, x.proj = x.proj))
  
}


s <- c(rep(0,20),rep(1,40),rep(-2,40))

e <- rnorm(100, sd = 0.6)

x <- s + e



teeth <- rep(c(rep(0,10),rep(1,20)),10)

saw <- rep(c(1:11,10:2),)

x <- teeth + e 
