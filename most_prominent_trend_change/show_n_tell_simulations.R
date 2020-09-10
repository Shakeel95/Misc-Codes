#-------------------------
# Setup: 
#-------------------------

w.dir <- "" # set wd to project folder! 
setwd(w.dir)

source("changepoint_detection_methods.R") # changepoint detection code
source("simulate_linear_signal_panel_kink_discontinuity.R") # panel simulation code

cpt.locs <- 3:97


#-------------------------------------------
# * Serially: uncorrelated
# * cross-sectionally: uncorrelated
# * small dense change 
#-------------------------------------------


set.seed(123)
cpt.loc <- 50 

res.lin <- c(); res.max <- c(); res.scan <- c(); res.SN <- c()
X <- sim.panel(tt = 100,tau = cpt.loc, cpt.density = 0.8, S = diag(rep(1,100)), noise.type = "gaussian", slope.range = c(-0.025,0.025))

for (i in cpt.locs) {
  
  CUSUM <- vec.GLR(X,i)
  res.lin <- c(res.lin, EH.lin.GLR(CUSUM,X,"MAD"))
  res.scan <- c(res.scan, EH.scan.GLR(CUSUM,X,"MAD"))
  res.SN <- c(res.SN, GLR.self.normalizing(CUSUM, X, i))
  res.max <- c(res.max, J.max.GLR(CUSUM, X, "MAD",b))
  
  cat("\r",paste0(i))
  
}

par(mfrow = c(2,2)) 

plot(cpt.locs, res.lin,
     main = paste("L2 aggregation:",which.max(res.lin)+3), 
     type = "l", 
     xlab = "time", 
     ylab = "value")

abline(v = cpt.loc, col = "red", lty = 2)
abline(v = which.max(res.lin)+3, col = "blue", lty = 2)

plot(cpt.locs, res.scan, 
     main = paste("scan aggregation:", which.max(res.scan) + 3),
     type = "l",
     xlab = "time", 
     ylab = "value")

abline(v = cpt.loc, col = "red", lty = 2)
abline(v = which.max(res.scan)+3, col = "blue", lty = 2)

plot(cpt.locs, res.max, 
     main = paste("L infinity aggregation:", which.max(res.max)+3),
     type = "l",
     xlab = "time", 
     ylab = "value")

abline(v = cpt.loc, col = "red", lty = 2)
abline(v = which.max(res.max)+3, col = "blue", lty = 2)

plot(cpt.locs, res.SN, 
     main = paste("L2 + self-normalizing:", which.max(res.SN) + 3), 
     type = "l",
     xlab = "time", 
     ylab = "value")

abline(v = cpt.loc, col = "red", lty = 2)
abline(v = which.max(res.SN)+3, col = "blue", lty = 2)




#-------------------------------------------
# * Serially: uncorrelated
# * cross-sectionally: uncorrelated
# * large sparse change
#-------------------------------------------

set.seed(123)
cpt.loc <- 20 

res.lin <- c(); res.max <- c(); res.scan <- c(); res.SN <- c()
X <- sim.panel(tt = 100,tau = cpt.loc, cpt.density = 0.05, S = diag(rep(1,100)), noise.type = "ar.factor", slope.range = c(-0.1,0.1))

for (i in cpt.locs) {
  
  CUSUM <- vec.GLR(X,i)
  res.lin <- c(res.lin, EH.lin.GLR(CUSUM,X,"MAD"))
  res.scan <- c(res.scan, EH.scan.GLR(CUSUM,X,"MAD"))
  res.SN <- c(res.SN, GLR.self.normalizing(CUSUM, X, i))
  res.max <- c(res.max, J.max.GLR(CUSUM, X, "MAD",b))
  
  cat("\r",paste0(i))
  
}

par(mfrow = c(2,2)) 


plot(cpt.locs, res.lin,
     main = paste("L2 aggregation:",which.max(res.lin)+3), 
     type = "l", 
     xlab = "time", 
     ylab = "value")

abline(v = cpt.loc, col = "red", lty = 2)
abline(v = which.max(res.lin)+3, col = "blue", lty = 2)

plot(cpt.locs, res.scan, 
     main = paste("scan aggregation:", which.max(res.scan) + 3),
     type = "l",
     xlab = "time", 
     ylab = "value")

abline(v = cpt.loc, col = "red", lty = 2)
abline(v = which.max(res.scan)+3, col = "blue", lty = 2)

plot(cpt.locs, res.max, 
     main = paste("L infinity aggregation:", which.max(res.max)+3),
     type = "l",
     xlab = "time", 
     ylab = "value")

abline(v = cpt.loc, col = "red", lty = 2)
abline(v = which.max(res.max)+3, col = "blue", lty = 2)

plot(cpt.locs, res.SN, 
     main = paste("L2 + self-normalizing:", which.max(res.SN) + 3), 
     type = "l",
     xlab = "time", 
     ylab = "value")

abline(v = cpt.loc, col = "red", lty = 2)
abline(v = which.max(res.SN)+3, col = "blue", lty = 2)



#-------------------------------------------
# * Serially: correlated
# * cross-sectionally: correlated
# * small dense change
# * balanced changepoint (50)
#-------------------------------------------


set.seed(123)
cpt.loc <- 50 

res.lin <- c(); res.max <- c(); res.scan <- c(); res.SN <- c()
X <- sim.panel(tt = 100,tau = cpt.loc, cpt.density = 0.8, S = AR.1(100, 0.9), noise.type = "AR.factor", slope.range = c(-0.05,0.05))

for (i in cpt.locs) {
  
  CUSUM <- vec.GLR(X,i)
  
  # res.lin <- c(res.lin, EH.lin.GLR(CUSUM,X,"kernel.LRV"))
  # res.scan <- c(res.scan, EH.scan.GLR(CUSUM,X,"kernel.LRV"))
  # res.SN <- c(res.SN, GLR.self.normalizing(CUSUM, X, i))
  # res.max <- c(res.max, J.max.GLR(CUSUM, X, "kernel.LRV",b))
  
  res.lin <- c(res.lin, EH.lin.GLR(CUSUM,X,"MAD"))
  res.scan <- c(res.scan, EH.scan.GLR(CUSUM,X,"MAD"))
  res.SN <- c(res.SN, SN.new(CUSUM, X, i))
  res.max <- c(res.max, J.max.GLR(CUSUM, X, "MAD",b))
  
  cat("\r",paste0(i))
  
}

par(mfrow = c(2,2)) 

plot(cpt.locs, res.lin,
     main = paste("L2 aggregation:",which.max(res.lin)+3), 
     type = "l", 
     xlab = "time", 
     ylab = "value")

abline(v = cpt.loc, col = "red", lty = 2)
abline(v = which.max(res.lin)+3, col = "blue", lty = 2)

plot(cpt.locs, res.scan, 
     main = paste("scan aggregation:", which.max(res.scan) + 3),
     type = "l",
     xlab = "time", 
     ylab = "value")

abline(v = cpt.loc, col = "red", lty = 2)
abline(v = which.max(res.scan)+3, col = "blue", lty = 2)

plot(cpt.locs, res.max, 
     main = paste("L infinity aggregation:", which.max(res.max)+3),
     type = "l",
     xlab = "time", 
     ylab = "value")

abline(v = cpt.loc, col = "red", lty = 2)
abline(v = which.max(res.max)+3, col = "blue", lty = 2)

plot(cpt.locs, res.SN, 
     main = paste("L2 + self-normalizing:", which.max(res.SN) + 3), 
     type = "l",
     xlab = "time", 
     ylab = "value")

abline(v = cpt.loc, col = "red", lty = 2)
abline(v = which.max(res.SN)+3, col = "blue", lty = 2)


