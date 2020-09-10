library("docstring") # EZ documentation

#-------------------------------------------------
# U-statistic for panel slope changepoint problem
#
# (use small samples, runs SLOW)
#--------------------------------------------------



U.slopes <- function(X,a,tau){
  
  #' U-statistic for detecing kink discontinuities in linear panel
  #'
  #'@param X matrix, dimension p x n
  #'@param a int, order of the U-statistic
  #'@param tau int, location of changepoint to investigate 
  
  
  V.k.list <- function(x,y,a) {
    
    #' Compute V statistic for two sample difference in slopes problem (component-wise)
    #'
    #'@param x vector, first sample 
    #'@param y vector, second sample
    #'@param a int, order of the U-statistic
    
    n.x <- length(x)
    n.y <- length(y)

    x.slopes <- combn(1:n.x,2,function(i) return((x[i[2]]-x[i[1]]) / (i[2] - i[1])))
    y.slopes <- combn(1:n.y,2,function(j) return((y[j[2]]-y[j[1]]) / (j[2] - j[1])))
    
    V.k <- as.list(numeric(a))
    
    
    for (j in 1:choose(n.x,2)) {
      for (k in 1:choose(n.y,2)) {
        
        v <- (x.slopes[j] - y.slopes[k])
        
        for (i in 1:a) V.k[[i]] <- V.k[[i]] + v**i
        
      }
    }
    
    return(V.k)
    
  }
  
  
  
  perm <- function(n,k){
    
    #' Permutations with large numbers
    #'
    #'@param n int, set size 
    #'@param k int, subset size 
    
    choose(n,k) * factorial(k)
    
  }
  
  
  
  U.l <- function(V.k, a) {
    
    #' Aggregate outputs from V.k.list to build univariate U-statistic (component-wise)
    #'
    #'@param V.k list, output from V.k.list
    #'@param a int, order of U-stats to compute
    
    
    # for a = {1,2} compute values and return
    
    if (a == 1) return(V.k[[1]])
    
    if (a == 2) return(c(V.k[[1]], V.k[[1]]**2 - V.k[[2]]))
    
    U.vec <- c(V.k[[1]], (V.k[[1]]**2) - V.k[[2]])
    
    # for higher order U-stat build via recursion
    
    for (r in 3:a) {
      
      T.l <- V.k[[r]]
      
      for (k in rev(1:(r-1))) {
        
        T.l <- V.k[[k]]*U.vec[[r-k]] - (r-k)*T.l
        
      }
      
      U.vec <- c(U.vec, T.l)
      
    }
    
    return(U.vec)
    
  }
  

  x.L <- x[1:tau]
  x.R <- x[(tau+1):length(x)]
  
  n.L <- choose(length(x.L),2)
  n.R <- choose(length(x.R),2)
  
  
  V.k.vals <- V.k.list(x.L,x.R,a)
  U.l.vals <- U.l(V.k.vals,a)
  

  return(U.l.vals / sapply(1:a,function(b) perm(n.L,b)*perm(n.R,b)))
  
}

