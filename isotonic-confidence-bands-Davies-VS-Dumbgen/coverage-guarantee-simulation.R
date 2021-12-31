set.seed(42)

coverage.guarantee <- function(ff, KK, thresh.vals)
{
  #'
  #'
  #'@param ff 
  #'@param KK
  #'@param thresh.vals
  
  res <- data.frame(matrix(0,KK,4))
  names(res) <- c("runs", "multires", "multires (centered)", "multires (weighted)")
  
  pb <- txtProgressBar(min = 0, max = KK, style = 3)
  
  for (ii in 1:KK)
  {
    yy <- ff + rcauchy(length(ff))
    
    f.fit <- runs.lower.bound(yy,thresh.vals[1])
    res[ii,1] <- all(ff > f.fit)
    
    f.fit <- multi.res.lower.bound(yy, thresh.vals[2])
    res[ii,2] <- all(ff > f.fit)
    
    f.fit <- centered.multi.res.lower.bound(yy, thresh.vals[3])
    res[ii,3] <- all(ff > f.fit)
    
    f.fit <- weighted.multi.res.lower.bound(yy, thresh.vals[4])
    res[ii,4] <- all(ff > f.fit)
    
    setTxtProgressBar(pb,ii)
  }
  
  res
}

#--------------------------------
#
# Constant median function
#
#---------------------------------

## sample size = 100
## 

const.ff.coverage.100 <- coverage.guarantee(
  ff = rep(1,100), 
  KK = 100, 
  thresh.vals = c(runs.thresh(100),
                  thresh_kab_bern(100),
                  1.335,
                  1.209)
)

apply(const.ff.coverage.100,2,mean)

## sample size = 500
## 

const.ff.coverage.500 <- coverage.guarantee(
  ff = rep(1,500), 
  KK = 100, 
  thresh.vals = c(runs.thresh(500),
                  thresh_kab_bern(500),
                  1.502,
                  1.359)
)

apply(const.ff.coverage.500,2,mean)


#---------------------------------
#
# Linear median function
#
#------------------------------------

## sample size = 100
## 

lin.ff.coverage.100 <- coverage.guarantee(
  ff = (1:100 / sqrt(100)), 
  KK = 100, 
  thresh.vals = c(runs.thresh(100),
                  thresh_kab_bern(100),
                  1.335,
                  1.209)
)

apply(lin.ff.coverage.100,2,mean)

## sample size = 500
## 

lin.ff.coverage.500 <- coverage.guarantee(
  ff = (1:500 / sqrt(500)), 
  KK = 100, 
  thresh.vals = c(runs.thresh(500),
                  thresh_kab_bern(500),
                  1.502,
                  1.359)
)

apply(lin.ff.coverage.500,2,mean)

