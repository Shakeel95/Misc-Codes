max.cpt.alignment <- function(tt,n,N,cardinality,runs = 100) {
  
  #' @param tt int, time points to simulate
  #' @param n int, panel size
  #' @param N int, number of changepoints
  #' @param cardinality int, cardinality of common changepoint set
  #' @param runs int, number of runs over which to  average
  
  max.cpt.random <- c()
  max.cpt.factor <- c()
  
  for (r in 1:runs) {
    
    X <- matrix(0,n,tt)
    for (i in 1:n) X[i,sample(1:tt,min(N,tt))] <- 1 
    max.cpt.random <- c(max.cpt.random, max(apply(X,2,sum)))
    
    X <- matrix(0,n,tt)
    cpt.common <- sample(1:tt,min(cardinality,tt))
    for (i in 1:n) X[i,sample(cpt.common,min(N,length(cpt.common)))] <- 1
    max.cpt.factor <- c(max.cpt.factor, max(apply(X,2,sum)))
    
  }
  
  return(list(max.random = median(max.cpt.random), 
              max.factor = median(max.cpt.factor)))
}



#--------------------------------------
# N. changepoints = diverging linearly
# set cardinatity = diverging linearly
#--------------------------------------

params <- as.list(data.frame(rbind(
  seq(from = 100, to = 1050, by = 50), # panel dimensions
  round(0.01*seq(from = 100, to = 1050, by = 50)), # numer of changepoints
  round(0.5*seq(from = 100, to = 1050, by = 50)) # cardinality of 
  )))

max.random <- c()
max.factor <- c()

for (param in params){
  
  res <- max.cpt.alignment(tt = param[1], n = param[1], N = param[2], cardinality = param[3])
  
  max.random <- c(max.random, res$max.random)
  max.factor <- c(max.factor, res$max.factor)
  
  cat("\r","calculated T = ", paste0(param[1]))
  
}

matplot(cbind(max.factor,max.random), 
        type = "b",
        xlab = "T,n",
        ylab = "maximum aligned changepoitns", 
        xaxt = "n"
)
axis(1, at = 1:20, seq(from = 100, to = 1050, by = 50))



#--------------------------------------
# N. changepoints = fixed
# set cardinatity = diverging linearly
#--------------------------------------

params <- as.list(data.frame(rbind(
  seq(from = 100, to = 1050, by = 50), # panel dimensions
  rep(10,20), # numer of changepoints
  round(0.5*seq(from = 100, to = 1050, by = 50)) # cardinality of common cpt set
)))

max.random <- c()
max.factor <- c()

for (param in params){
  
  res <- max.cpt.alignment(tt = param[1], n = param[1], N = param[2], cardinality = param[3])
  
  max.random <- c(max.random, res$max.random)
  max.factor <- c(max.factor, res$max.factor)
  
  cat("\r","calculated T = ", paste0(param[1]))
  
}

matplot(cbind(max.factor,max.random), 
        type = "b",
        xlab = "T,n",
        ylab = "maximum aligned changepoitns", 
        xaxt = "n"
)
axis(1, at = 1:20, seq(from = 100, to = 1050, by = 50))



#---------------------------------------------
# N. changepoints = fixed
# set cardinatity = diverging lOgaraithmically
#---------------------------------------------

params <- as.list(data.frame(rbind(
  seq(from = 100, to = 1050, by = 50), # panel dimensions
  rep(10,20), # numer of changepoints
  round(10*log(seq(from = 100, to = 1050, by = 50))) # cardinality of common cpt set
)))

max.random <- c()
max.factor <- c()

for (param in params){
  
  res <- max.cpt.alignment(tt = param[1], n = param[1], N = param[2], cardinality = param[3])
  
  max.random <- c(max.random, res$max.random)
  max.factor <- c(max.factor, res$max.factor)
  
  cat("\r","calculated T = ", paste0(param[1]))
  
}

matplot(cbind(max.factor,max.random), 
        type = "b",
        xlab = "T,n",
        ylab = "maximum aligned changepoitns", 
        xaxt = "n"
)
axis(1, at = 1:20, seq(from = 100, to = 1050, by = 50))



#-----------------------------------
# Finding the number of factors
# by ranking changepoint alignments
#-----------------------------------

N <- 10
n <- 500
tt <- 500
error <- 50
runs <- 100
cardinality <- 50


# no factor structure

res <- matrix(0,runs,tt)

for (r in 1:runs) {
  
  X <- matrix(0,n,tt)
  for (i in 1:n) X[i,sample(1:tt,N+rpois(1,error))] <- 1
  
  res[r,] <- sort(apply(X,2,sum))
}

matplot(t(res),
        xlab = "column", 
        ylab = "No. aligned changepoints",
        type = "l")

# factor structure

error <- 0 
res <- matrix(0,runs,tt)

for (r in 1:runs) {
  
  common <- sample(1:tt,cardinality)
  X <- matrix(0,n,tt)
  
  for (i in 1:n) {
    cpt <- c(sample(common,N),
             sample(setdiff(1:tt,common),rpois(1,error)))
    X[i,cpt] <- 1
  }
  
  res[r,] <- sort(apply(X,2,sum))
}

matplot(t(res),
        main = "Factor structure", 
        xlab = "column", 
        ylab = "No. aligned changepoints",
        type = "l")

abline( v = tt - cardinality, col = "red")



# factor structure + idiosyncratic changepoints

error <- 50 
res <- matrix(0,runs,tt)

for (r in 1:runs) {
  
  common <- sample(1:tt,cardinality)
  X <- matrix(0,n,tt)
  
  for (i in 1:n) {
    cpt <- c(sample(common,N),
             sample(setdiff(1:tt,common),rpois(1,error)))
    X[i,cpt] <- 1
  }
  
  res[r,] <- sort(apply(X,2,sum))
}

matplot(t(res),
        xlab = "column", 
        ylab = "No. aligned changepoints",
        type = "l")

abline( v = tt - cardinality, col = "red")
