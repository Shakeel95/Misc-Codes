#---------------------
#
# Paralellize 
# 
#----------------------


library("doParallel")
library("foreach")

# register backend
cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)



#------------------------------
#
# Global simulation parameters
# 
#------------------------------ 

set.seed(42)

sim.size <- 100 # number of simulations to run
n <- 100 # number of channels 
tt <- 100 # number of time points
cpt.density <- .75 # density at which all methods work well
cpt.loc <- floor(.5*tt) # balanced changepoint location always gives most power
slope.vals <- seq(.005, .075, length.out =  10)



#-------------------------------------------
#
# Simulation (1): identity dispersion matrix
#
#-------------------------------------------


res.simulation.1 <-  matrix(,4,length(slope.vals))


for (i in 1:length(slope.vals)) {
  
  
  res <- foreach (j = 1:sim.size, .combine = "rbind") %dopar% {
    
    source("changepoint_detection_methods.R") 
    source("simulate_linear_signal_panel_kink_discontinuity.R")
    
    X <- sim.panel(tt = tt,
                   tau = cpt.loc, 
                   cpt.density = cpt.density, 
                   S = diag(runif(n,1,3)), 
                   noise.type = "gaussian", 
                   slope.range = c(-slope.vals[i], slope.vals[i]))
    
    
    est.cpt <- estimate.cpt.loc(X,"MAD", verbose = F)
    
    c(
      abs(est.cpt$res.l2$cpt - cpt.loc),
      abs(est.cpt$res.max$cpt - cpt.loc),
      abs(est.cpt$res.scan$cpt - cpt.loc),
      abs(est.cpt$res.SN$cpt - cpt.loc)
    )
    
  }
  
  
  cat("\r", "Finished simulating slope change magnitude", paste0(round(slope.vals[i],3)))
  
  res.simulation.1[,i] <- apply(res, 2, mean)
  
}


rownames(res.simulation.1) <- c("L2", "max", "scan", "SN")
colnames(res.simulation.1) <- round(slope.vals, 3)
View(res.simulation.1)

write.csv(res.simulation.1, file = "data/slope_change_mag_1.csv")



#--------------------------------------------------------------
#
# Simulation (2): toeplitz (autoregressive) spatial dependence
#
#---------------------------------------------------------------


res.simulation.2 <-  matrix(,4,length(slope.vals))


for (i in 1:length(slope.vals)) {
  
  
  res <- foreach (j = 1:sim.size, .combine = "rbind") %dopar% {
    
    source("changepoint_detection_methods.R") 
    source("simulate_linear_signal_panel_kink_discontinuity.R")
    
    X <- sim.panel(tt = tt,
                   tau = cpt.loc, 
                   cpt.density = cpt.density, 
                   S = 3*AR.1(n, (-1)**rbinom(1,1,0.5) *runif(1,(1/3),1)), 
                   noise.type = "gaussian", 
                   slope.range = c(-slope.vals[i], slope.vals[i]))
    
    
    est.cpt <- estimate.cpt.loc(X,"MAD", verbose = F)
    
    c(
      abs(est.cpt$res.l2$cpt - cpt.loc),
      abs(est.cpt$res.max$cpt - cpt.loc),
      abs(est.cpt$res.scan$cpt - cpt.loc),
      abs(est.cpt$res.SN$cpt - cpt.loc)
    )
    
  }
  
  
  cat("\r", "Finished simulating slope change magnitude", paste0(round(slope.vals[i],3)))
  
  res.simulation.2[,i] <- apply(res, 2, mean)
  
}


rownames(res.simulation.2) <- c("L2", "max", "scan", "SN")
colnames(res.simulation.2) <- round(slope.vals, 2)
View(res.simulation.2)

write.csv(res.simulation.2, file = "data/slope_change_mag_2.csv")



#----------------------------------------------------------
#
# Simulation (3): toeplitz spatial dependence + heavy tails
#
#----------------------------------------------------------


res.simulation.3 <-  matrix(,4,length(slope.vals))


for (i in 1:length(slope.vals)) {
  
  
  res <- foreach (j = 1:sim.size, .combine = "rbind") %dopar% {
    
    source("changepoint_detection_methods.R") 
    source("simulate_linear_signal_panel_kink_discontinuity.R")
    
    X <- sim.panel(tt = tt,
                   tau = cpt.loc, 
                   cpt.density = cpt.density, 
                   S = 3*AR.1(n, (-1)**rbinom(1,1,0.5) *runif(1,(1/3),1)), 
                   noise.type = "t.3", 
                   slope.range = c(-slope.vals[i], slope.vals[i]))
    
    
    est.cpt <- estimate.cpt.loc(X,"MAD", verbose = F)
    
    c(
      abs(est.cpt$res.l2$cpt - cpt.loc),
      abs(est.cpt$res.max$cpt - cpt.loc),
      abs(est.cpt$res.scan$cpt - cpt.loc),
      abs(est.cpt$res.SN$cpt - cpt.loc)
    )
    
  }
  
  
  cat("\r", "Finished simulating slope change magnitude", paste0(round(slope.vals[i],3)))
  
  res.simulation.3[,i] <- apply(res, 2, mean)
  
}


rownames(res.simulation.3) <- c("L2", "max", "scan", "SN")
colnames(res.simulation.3) <- round(slope.vals, 2)
View(res.simulation.3)

write.csv(res.simulation.3, file = "data/slope_change_mag_3.csv")



#-------------------------------------------------------------------------------------
#
# Simulation (4): toeplitz spatial dependence + AR(1) (coef = 0.5) serial correlation
#
#-------------------------------------------------------------------------------------


res.simulation.4 <-  matrix(,4,length(slope.vals))


for (i in 1:length(slope.vals)) {
  
  
  res <- foreach (j = 1:sim.size, .combine = "rbind") %dopar% {
    
    source("changepoint_detection_methods.R") 
    source("simulate_linear_signal_panel_kink_discontinuity.R")
    
    X <- sim.panel(tt = tt,
                   tau = cpt.loc, 
                   cpt.density = cpt.density, 
                   S = 3*AR.1(n, (-1)**rbinom(1,1,0.5) *runif(1,(1/3),1)), 
                   noise.type = "AR.factor", 
                   slope.range = c(-slope.vals[i], slope.vals[i]))
    
    
    est.cpt <- estimate.cpt.loc(X,"MAD", verbose = F)
    
    c(
      abs(est.cpt$res.l2$cpt - cpt.loc),
      abs(est.cpt$res.max$cpt - cpt.loc),
      abs(est.cpt$res.scan$cpt - cpt.loc),
      abs(est.cpt$res.SN$cpt - cpt.loc)
    )
    
  }
  
  
  cat("\r", "Finished simulating slope change magnitude", paste0(round(slope.vals[i],3)))
  
  res.simulation.4[,i] <- apply(res, 2, mean)
  
}


rownames(res.simulation.4) <- c("L2", "max", "scan", "SN")
colnames(res.simulation.4) <- round(slope.vals, 2)
View(res.simulation.4)

write.csv(res.simulation.4, file = "data/slope_change_mag_4.csv")



#-----------------------
#
# End paralellize 
#
#-----------------------


stopCluster(cl)