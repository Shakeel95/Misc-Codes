
#--------------------------------------
#
# Piecewise constant signal functions
#
#--------------------------------------

teeth <- function(n.teeth,block.size, jump) {
  
  #' Teeth function from WBS paper
  #' 
  #' @param n.teeth int, number of blocks
  #' @param block.size int, size of teeth
  #' @param jump float, jump size
  
  
  block <- c(rep(0,block.size),rep(jump,block.size))
  
  rep(block, n.teeth)
  
}

rand.teeth <- function(n.teeth, spacing, min.jump, max.jupm) {
  
  #'Teeth function 
  #'
  #'@param n.teeth 
  #'
  
}



#-----------------------------------
#
# Piecewise linear signal functions
#
#-----------------------------------



#-----------------------------------
#
# Mixed polynomial signal functions
#
#-----------------------------------


