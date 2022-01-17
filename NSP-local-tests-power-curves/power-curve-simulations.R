
#----------------------------
#
# Power curve simulations
#
#----------------------------

set.seed(42)

monte.carlo.reps <- 100
jump.seq.len <- 10

pb <- txtProgressBar(min = 0, max = jump.seq.len, style = 3)

source("lower-bound-constructions.R")
source("local-tests.R")

local.tests = list(
  function(xx) rnsp.fryz(xx, "anchored"),
  function(xx) rnsp.fryz(xx, "complete"),
  function(xx) rnsp.inv(xx, "mr.inv", length(xx) |> thresh.kab.bern()),
  function(xx) rnsp.inv(xx, "mr.inv.corrected", c(length(xx), length(xx) |> thresh.dumbgen.T0()))
)

local.test.10mult = list(
  function(xx) rnsp.fryz(xx, "anchored", (10*length(xx)) |> thresh.kab.bern()),
  function(xx) rnsp.fryz(xx, "complete", (10*length(xx)) |> thresh.kab.bern()),
  function(xx) rnsp.inv(xx, "mr.inv", (10*length(xx)) |> thresh.kab.bern()),
  function(xx) rnsp.inv(xx, "mr.inv.corrected", c(10*length(xx),(10*length(xx)) |> thresh.dumbgen.T0()))
)


plot.power.curves <- function(res, jumps)
{
  #' Plot the power curves
  #' 
  #'@param res array
  #'@param jumps vector

  res.mean <- c()
  for (ii in 1:length(jumps)) res.mean <- cbind(res.mean, apply(res[,,ii],1,mean))
  
  plot(res.mean[1,]~jumps,
       type = "b",
       ylim = c(0,1),
       xlab = "jump size",
       ylab = "power")
  
  abline(h = 1, lty = 2)
  
  for (ii in 2:4) lines(
    res.mean[ii,]~jumps, 
    type = "b",
    pch = ii, 
    col = ii
    )
}


##
## Gaussian noise

res.gauss <- array(1, c(4,monte.carlo.reps,jump.seq.len))
res.gauss.10mult <- array(1, c(4,monte.carlo.reps,jump.seq.len))
jumps <- seq(from = 0.5, to = 1.5, length.out = jump.seq.len)

for (ii in 1:jump.seq.len)
{
  ff <- c(rep(0,50), rep(jumps[ii],50))
  for (jj in 1:monte.carlo.reps)
  {
    ee <- ff + rnorm(100)
    for (kk in 1:4) res.gauss[kk,jj,ii] <- local.tests[[kk]](ee)
  }
  setTxtProgressBar(pb, ii)
}

png("gaussian-noise")
plot.power.curves(res.gauss, jumps)
dev.off()

for (ii in 1:jump.seq.len)
{
  ff <- c(rep(0,50), rep(jumps[ii],50))
  for (jj in 1:monte.carlo.reps)
  {
    ee <- ff + rnorm(100)
    for (kk in 1:4) res.gauss.10mult[kk,jj,ii] <- local.test.10mult[[kk]](ee)
  }
  setTxtProgressBar(pb, ii)
}

png("gaussian-noise-10mult")
plot.power.curves(res.gauss.10mult, jumps)
dev.off()


##
## Cauchy noise

res.cuachy <- array(1, c(4,monte.carlo.reps,jump.seq.len))
res.cuachy.10mult <- array(1, c(4,monte.carlo.reps,jump.seq.len))
jumps <- seq(from = 0.5, to = 2.5, length.out = jump.seq.len)

for (ii in 1:jump.seq.len)
{
  ff <- c(rep(0,50), rep(jumps[ii],50))
  for (jj in 1:monte.carlo.reps)
  {
    ee <- ff + rcauchy(100)
    for (kk in 1:4) res.cuachy[kk,jj,ii] <- local.tests[[kk]](ee)
  }
  setTxtProgressBar(pb, ii)
}


png("cauchy-noise")
plot.power.curves(res.cuachy, jumps)
dev.off()


for (ii in 1:jump.seq.len)
{
  ff <- c(rep(0,50), rep(jumps[ii],50))
  for (jj in 1:monte.carlo.reps)
  {
    ee <- ff + rcauchy(100)
    for (kk in 1:4) res.cuachy.10mult[kk,jj,ii] <- local.test.10mult[[kk]](ee)
  }
  setTxtProgressBar(pb, ii)
}

png("cauchy-noise-10mult")
plot.power.curves(res.cuachy.10mult, jumps)
dev.off()


##
## Chi-square noise

res.chi.sqr <- array(1, c(4,monte.carlo.reps,jump.seq.len))
res.chi.sqr.10mult <- array(1, c(4,monte.carlo.reps,jump.seq.len))
jumps <- seq(from = 0, to = 20, length.out = jump.seq.len)

for (ii in 1:jump.seq.len)
{
  ff <- c(rep(1,50), rep(1+jumps[ii],50))
  for (jj in 1:monte.carlo.reps)
  {
    ee <- (rnorm(100)**2)*ff
    for (kk in 1:4) res.chi.sqr[kk,jj,ii] <- local.tests[[kk]](ee)
  }
  setTxtProgressBar(pb, ii)
}

png("chi-square-noise")
plot.power.curves(res.chi.sqr, jumps)
dev.off()

for (ii in 1:jump.seq.len)
{
  ff <- c(rep(1,50), rep(1+jumps[ii],50))
  for (jj in 1:monte.carlo.reps)
  {
    ee <- (rnorm(100)**2)*ff
    for (kk in 1:4) res.chi.sqr.10mult[kk,jj,ii] <- local.test.10mult[[kk]](ee)
  }
  setTxtProgressBar(pb, ii)
}

png("chi-square-noise-10mult")
plot.power.curves(res.chi.sqr.10mult, jumps)
dev.off()


#
# Poisson noise

res.poisson <- array(1, c(4,monte.carlo.reps,jump.seq.len))
res.poisson.10mult <- array(1, c(4,monte.carlo.reps,jump.seq.len))
jumps <- seq(from = 0, to = 3, length.out = jump.seq.len)

for (ii in 1:jump.seq.len)
{
  for (jj in 1:monte.carlo.reps)
  {
    ee <- c(rpois(50,1),rpois(50,1+jumps[ii]))
    for (kk in 1:4) res.poisson[kk,jj,ii] <- local.tests[[kk]](ee)
  }
  setTxtProgressBar(pb, ii)
}

png("poisson-noise")
plot.power.curves(res.poisson, jumps)
dev.off()

for (ii in 1:jump.seq.len)
{
  for (jj in 1:monte.carlo.reps)
  {
    ee <- c(rpois(50,1),rpois(50,1+jumps[ii]))
    for (kk in 1:4) res.poisson.10mult[kk,jj,ii] <- local.test.10mult[[kk]](ee)
  }
  setTxtProgressBar(pb, ii)
}

png("poisson-noise-10mult")
plot.power.curves(res.poisson.10mult, jumps)
dev.off()

##
## Bernoulli noise

res.bern <- array(1, c(4,monte.carlo.reps,jump.seq.len))
res.bern.10mult <- array(1, c(4,monte.carlo.reps,jump.seq.len))
jumps <- seq(from = 0, to = .9, length.out = jump.seq.len)

for (ii in 1:jump.seq.len)
{
  for (jj in 1:monte.carlo.reps)
  {
    ee <- c(rbinom(50,1,.1),rbinom(50,1,.1 + jumps[ii]))
    for (kk in 1:4) res.bern[kk,jj,ii] <- local.tests[[kk]](ee)
  }
  setTxtProgressBar(pb, ii)
}

png("bernoulli-noise")
plot.power.curves(res.bern, jumps)
dev.off()

for (ii in 1:jump.seq.len)
{
  for (jj in 1:monte.carlo.reps)
  {
    ee <-  c(rbinom(50,1,.1),rbinom(50,1,.1 + jumps[ii]))
    for (kk in 1:4) res.bern.10mult[kk,jj,ii] <- local.test.10mult[[kk]](ee)
  }
  setTxtProgressBar(pb, ii)
}

png("bernoulli-noise-10mult")
plot.power.curves(res.bern.10mult, jumps)
dev.off()

##
## Geometric noise 

res.geom <- array(1, c(4,monte.carlo.reps,jump.seq.len))
res.geom.10mult <- array(1, c(4,monte.carlo.reps,jump.seq.len))
jumps <- seq(from = 0, to = .5, length.out = jump.seq.len)

for (ii in 1:jump.seq.len)
{
  for (jj in 1:monte.carlo.reps)
  {
    ee <- c(rgeom(50,.1),rgeom(50,.1+jumps[ii]))
    for (kk in 1:4) res.geom[kk,jj,ii] <- local.tests[[kk]](ee)
  }
  setTxtProgressBar(pb, ii)
}

png("geometric-noise")
plot.power.curves(res.geom, jumps)
dev.off()

for (ii in 1:jump.seq.len)
{
  for (jj in 1:monte.carlo.reps)
  {
    ee <-  c(rbinom(50,1,.1),rbinom(50,1,.1 + jumps[ii]))
    for (kk in 1:4) res.geom.10mult[kk,jj,ii] <- local.test.10mult[[kk]](ee)
  }
  setTxtProgressBar(pb, ii)
}

png("geometric-noise-10mult")
plot.power.curves(res.geom.10mult, jumps)
dev.off()
