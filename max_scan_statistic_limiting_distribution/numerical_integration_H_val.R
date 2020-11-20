
#---------------------------------------------------------------
#
# Numerical integration for estimating the constant H in L2 scan
#
# (Piotr's code)
#
#---------------------------------------------------------------


# Running the chunk gives: H = 1.57609

  
H.integrand <- function(y, up.sum = 10000) {

    k <- 1:up.sum

    exp(-4 * sum(1/k * pnorm(-sqrt(k/(y)))))

}


H.value <- function(up.sum = 10000, grid.size = 0.01, grid.endpoint = 10) {

  dy <- seq(from = grid.size, by = grid.size, to = grid.endpoint)

  n <- length(dy)

  H.int <- rep(0, n)
  
  for (i in 1:n) H.int[i] <- H.integrand(dy[i], up.sum)
      
  sum(H.int * grid.size)
    
}


H.value(10**5, .001, 20)

