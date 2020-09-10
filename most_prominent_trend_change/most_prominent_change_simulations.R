#-------------------------
# Setup: 
#-------------------------

w.dir <- "" # set wd to project folder! 
setwd(w.dir)

source("changepoint_detection_methods.R") # changepoint detection code
source("simulate_linear_signal_panel_kink_discontinuity.R") # panel simulation code

  
#-------------------------------------
# 
# Simulation (1)
# 
#   * Spatially uncorrelated errors
#   * Serially unocrrelated errors
# 
#--------------------------------------


X <- sim.panel(tt = 100,
               tau = 50,
               cpt.density = 1, 
               diag(rep(1,100)), 
               noise.type = "AR.ind", 
               slope.range = c(-0.05,0.05))

to.check <- 3:97

res.lin.MAD <- c()
res.scan.MAD <- c()

res.lin.kern <- c()
res.scan.kern <- c()

for (i in to.check) {
  
  CUSUM <- vec.GLR(X,i)
  
  res.scan.MAD <- c(res.scan.MAD, EH.scan.GLR(CUSUM,X, var.est = "MAD"))
  res.lin.MAD <- c(res.lin.MAD, EH.lin.GLR(CUSUM,X,var.est = "MAD"))
  
  res.lin.kern <- c(res.lin.kern, EH.lin.GLR(CUSUM, X, var.est = "kernel.LRV"))
  res.scan.kern <- c(res.scan.kern, EH.scan.GLR(CUSUM, X, var.est = "kernel.LRV"))
  
  cat("\r","finished calculating stat for b = ", paste0(i))
  
}

par(mfrow = c(2,2))

plot(to.check,res.scan.MAD, type = "l", main = "scan CUSUM (MAD)")
abline(v = 50, col = "red")

plot(to.check, res.lin.MAD, type = "l", main = "linear CUSUM (MAD)")
abline(v = 50, col = "red")

plot(to.check,res.scan.kern, type = "l", main = "scan CUSUM (kernel)")
abline(v = 50, col = "red")

plot(to.check, res.lin.kern, type = "l", main = "linear CUSUM (kernel)")
abline(v = 50, col = "red")



S <- diag(runif(200, min = 0, max = 3))
X <- sim.panel(tt = 100, tau = 20, cpt.density = 0.1, S, noise.type = "AR.factor", slope.range = c(-0.5,0.5))

res.SN <- c()
res.lin <- c()
to.check <- 3:97

for (i in to.check) {
  
  CUSUM <- vec.GLR(X,i)
  res.SN <- c(res.SN,GLR.self.normalizing(CUSUM,X,i))
  cat("\r","just finished calculating ",paste0(i),"-th statistic")
  
}

par(mfrow = c(1,2))

plot(to.check, res.lin, type = "l", main = "linear statistic")
abline(v = 20, col = "red")

plot(to.check, res.SN, type = "l", main = "SN statistic")
abline(v = 20, col = "red")
