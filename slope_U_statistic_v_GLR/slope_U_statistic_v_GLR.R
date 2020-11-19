library("docstring") # EZ documentation 


#--------------------------------------------
# 
# GLR function for detecting change in slope 
# 
#--------------------------------------------


GLR.slope <- function(x) {
  
  #' GLR approach for detecting location of change in slope
  #'
  #'@param x vector, univariate signal 
  
  
  lin.contrast <- function(s, b, e){
    
    #' pcwsLin contrast function for GLR (see NOT paper)
    #' 
    #' @param b int, change
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
  
  
  n <- length(x)
  
  res <- sapply(2:(n-1), function(b) sum(x * lin.contrast(1,b,n))**2)
  
  return(list(res = res, cpt = which.max(res) + 1))
  
}



#-------------------------------------------
# 
# Two-sample t-test for slope parameter 
# 
#-------------------------------------------


t.test.slope <- function(x) {
  
  #' Rolling two sample t-test for changepoint location
  #' 
  #' @param x vector, univariate signal
  
  
  b <- function(y,s,e) {
    
    #' Quick! estimate slope parameter
    #' 
    #' @param y vector, signal
    #' @param s int, start here 
    #' @param e int, end here
    
    y <- y[s:e]
    x <- 1:length(y)
    
    u.x <- mean(x)
    u.y <- mean(y)
    
    
    return(sum((y - u.y)*(x - u.x)) / sum((x - u.x)**2))
    
  }
  
  
  v.correct <- function(s,b,e) {
  
    #' Proportional to variance of diff-in-slopes stat
    #' 
    #'@param s int, start here 
    #'@param b int, break here
    #'@param e int, end here  
    
    return(
      
      sum(((s:b) - mean(s:b))**2)**(-1) + sum((((b+1):e) - mean((b+1):e))**2)**(-1)
      
    )
    
  }
  
  
  n <- length(x)
  
  res <- sapply(3:n-2, function(j) abs(b(x,1,j) - b(x,j+1,n)) / v.correct(1,j,n))
  
  return(list(res = res, cpt = which.max(res)))
  
}



#-----------------------------------------------------
# 
# Two-sample U statistic for detecting change in slope
# 
#-----------------------------------------------------


U.slope <- function(x) {
  
  #' Two-sample U-statistic approach for detecting location of a change in slope
  #'   
  #' N.B. naive implementation == super slow!    
  #'   
  #' @param x vector, univaraite signal
  
  
  n 
  
  for (b in 2:)
  
  L.combn <- 
  
}

j <- 0 

start <- Sys.time()

for (i in 1:(100**4)) j <- j + 1 

end <- Sys.time()

end - start



#-------------
# 
# Simulations
#
#-------------

sim.size <- 500

set.seed(42)

cpt.loc <- 50 

GLR.slope(y)

t.test.slope(y)



## balanced changepoint location

set.seed(42)
cpt.loc <- 50

t.test.res.bal <- c()
GLR.res.bal <- c()

for (i in 1:sim.size) {
  
  e <- 5*rnorm(100)
  s <- 0.5*(1:100) + 0.505*ifelse(1:100 - cpt.loc > 0 , 1:100 - cpt.loc, 0)
  y <- s + e 
  
  t.test.res.bal <- c(t.test.res.bal, abs(t.test.slope(y)$cpt - cpt.loc))
  GLR.res.bal <- c(GLR.res.bal, abs(GLR.slope(y)$cpt - cpt.loc))
  
  cat("\r","Finished simulation: ", paste0(i))
  
}


## unbalanced changepoint location

set.seed(42)
cpt.loc <- 20

t.test.res.ubal <- c()
GLR.res.ubal <- c()

for (i in 1:sim.size) {
  
  e <- 5*rnorm(100)
  s <- (1:100) + 0.505*ifelse(1:100 - cpt.loc > 0 , 1:100 - cpt.loc, 0)
  y <- s + e 
  
  t.test.res.ubal <- c(t.test.res.ubal, abs(t.test.slope(y)$cpt - cpt.loc))
  GLR.res.ubal <- c(GLR.res.ubal, abs(GLR.slope(y)$cpt - cpt.loc))
  
  cat("\r","Finished simulation: ", paste0(i))
  
}


y <- rnorm(100) + .01*(1:100)


sum(combn(1:100, 2, function(i) diff(y[i]) / diff(i))) / choose(100,2)


