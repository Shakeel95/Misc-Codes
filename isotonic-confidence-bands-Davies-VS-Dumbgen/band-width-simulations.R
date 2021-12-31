set.seed(42)

median.widths <- function(ff, KK, thresh.vals)
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
    res[ii,1] <- median(abs(ff - f.fit))
    
    f.fit <- multi.res.lower.bound(yy, thresh.vals[2])
    res[ii,2] <- median(abs(ff - f.fit))
    
    f.fit <- centered.multi.res.lower.bound(yy, thresh.vals[3])
    res[ii,3] <- median(abs(ff - f.fit))
    
    f.fit <- weighted.multi.res.lower.bound(yy, thresh.vals[4])
    res[ii,4] <- median(abs(ff - f.fit))
    
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

const.ff.widths.100 <- median.widths(
  ff = rep(1,100), 
  KK = 100, 
  thresh.vals = c(runs.thresh(100),
                  thresh_kab_bern(100),
                  1.335,
                  1.209)
)

apply(const.ff.widths.100,2,mean)

## sample size = 500
## 

const.ff.widths.500 <- median.widths(
  ff = rep(1,500), 
  KK = 100, 
  thresh.vals = c(runs.thresh(500),
                  thresh_kab_bern(500),
                  1.502,
                  1.359)
)

apply(const.ff.widths.500,2,mean)


#---------------------------------
#
# Linear median function
#
#------------------------------------

## sample size = 100
## 

lin.ff.widths.100 <- median.widths(
  ff = (1:100 / sqrt(100)), 
  KK = 100, 
  thresh.vals = c(runs.thresh(100),
                  thresh_kab_bern(100),
                  1.335,
                  1.209)
)

apply(const.ff.widths.100,2,mean)

## sample size = 500
## 

const.ff.widths.500 <- median.widths(
  ff = (1:500 / sqrt(500)), 
  KK = 100, 
  thresh.vals = c(runs.thresh(500),
                  thresh_kab_bern(500),
                  1.502,
                  1.359)
)

apply(const.ff.widths.500,2,mean)


#---------------------------------
#
# Quadratic median function
#
#------------------------------------

## sample size = 100
## 

lin.ff.widths.100 <- median.widths(
  ff = (1:100 / sqrt(100))**2, 
  KK = 100, 
  thresh.vals = c(runs.thresh(100),
                  thresh_kab_bern(100),
                  1.335,
                  1.209)
)

apply(const.ff.widths.100,2,mean)

## sample size = 500
## 

const.ff.widths.500 <- median.widths(
  ff = (1:500 / sqrt(500))**2, 
  KK = 100, 
  thresh.vals = c(runs.thresh(500),
                  thresh_kab_bern(500),
                  1.502,
                  1.359)
)

apply(const.ff.widths.500,2,mean)
