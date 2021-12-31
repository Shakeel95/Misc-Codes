
set.seed(42)

#--------------------------------
#
# Constant median function
#
#---------------------------------

nn <- 100
ff <- rep(5,nn)
yy <- ff + rcauchy(nn)

plot(yy, ylim = c(-5,15))
lines(ff, type = "l", lty = 2)

lines(
  runs.lower.bound(yy,runs.thresh(n)), 
)

lines(
  multi.res.lower.bound(yy, thresh_kab_bern(n)),
  col = "red"
)

lines(
  centered.multi.res.lower.bound(yy, 1.335), 
  col = "blue"
)

lines(
  weighted.multi.res.lower.bound(yy, 1.209), 
  col = "green"
)

legend(1,15,
       c("truth", "runs", "multires","centered","weighted"),
       col = c(rep("black",2),"red","blue","green"),
       lty = c(2,1,1,1,1))


#--------------------------------
#
# Linear median function
#
#---------------------------------


nn <- 100
ff <- 1:nn / sqrt(nn)
yy <- ff + rcauchy(nn)

plot(yy, ylim = c(-5,15))
lines(ff, type = "l", lty = 2)

lines(
  runs.lower.bound(yy,runs.thresh(n)), 
)

lines(
  multi.res.lower.bound(yy, thresh_kab_bern(n)),
  col = "red"
)

lines(
  centered.multi.res.lower.bound(yy, 1.335), 
  col = "blue"
)

lines(
  weighted.multi.res.lower.bound(yy, 1.209), 
  col = "green"
)

legend(1,15,
       c("truth", "runs", "multires","centered","weighted"),
       col = c(rep("black",2),"red","blue","green"),
       lty = c(2,1,1,1,1))
