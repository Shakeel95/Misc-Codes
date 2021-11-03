
#===========================================
#
# Test whether two convex hulls intersect
#
#===========================================


convex.hulls.intersect <- function(PP,QQ)
{
   #' Determin whether two convex hulls intersect
   #'
   #'@param PP matrix, convex hull 
   #'@param QQ matrix, convex hull

  nn <- nrow(PP)
  mm <- nrow(QQ)
  
  for (ii in 2:(nn-1)) for (jj in 2:(mm-1))
  {
    if (triangle.intersection(PP[c(1,ii,ii+1),],QQ[c(1,jj,jj+1),])) return(TRUE)
  }
  
  return(FALSE)
  
}


#===============================================================
#
# Helper functions
#
# *mostly from section 1.6 of "Computational Geometry in C"
#
#===============================================================


triangle.intersection <- function(T1, T2)
{
  #' Check whethere T1 intersects T2
  #'
  #'@param T1 matrix, edge coordinates 
  #'@param T2 matrix, edge coordinates
  
  edge.combinations <- utils::combn(3,2)
  for (ii in 1:3)
  {
    aa <- T1[edge.combinations[1,ii],]
    bb <- T1[edge.combinations[2,ii],]
    for (jj in 1:3)
    {
      cc <- T2[edge.combinations[1,jj],]
      dd <- T2[edge.combinations[2,jj],]
      
      if (line_segment_intersection(aa,bb,cc,dd)) return(TRUE)
    }
  }
  
  
  for (ii in 1:3) if (point_contained_in_triangle(T2[ii,], T1)) return(TRUE) 
  
  return(FALSE)
  
}

point_contained_in_triangle <- function(pp, TT)
{
  #' Checks if point pp is contained in the triangle TT
  #'
  #'@param pp vector, cooridnates of point
  #'@param TT matrix, coordinates of trinagle edges
  #'
  #'@references https://blackpawn.com/texts/pointinpoly/default.html
  
  #  Compute vectors        
  v0 <- TT[3,] - TT[1,]
  v1 <- TT[2,] - TT[1,]
  v2 <- pp - TT[1,]
  
  ## Compute dot products
  dot00 <- sum(v0*v0)
  dot01 <- sum(v0*v1)
  dot02 <- sum(v0*v2)
  dot11 <- sum(v1*v1)
  dot12 <- sum(v1*v2)
  
  ## Compute barycentric coordinates
  
  invDenom <- 1 / (dot00 * dot11 - dot01 * dot01)
  u <- (dot11 * dot02 - dot01 * dot12) * invDenom
  v <- (dot00 * dot12 - dot01 * dot02) * invDenom
  
  ## Check if point is in triangle
  
  return((u >= 0) & (v >= 0) & (u + v < 1))
  
}


line_segment_intersection <- function(aa, bb, cc, dd)
{
  #' Checks for intersection between line segments ab and cd
  #' 
  #' @param aa vector
  #' @param bb vector 
  #' @param cc vector 
  #' @param dd vector
  
  if (propper_intersection(aa,bb,cc,dd))
  {
    return(TRUE)
    
  } else if (impropper_intersection_between(aa,bb,cc)|
             impropper_intersection_between(aa,bb,dd)|
             impropper_intersection_between(cc,dd,aa)|
             impropper_intersection_between(cc,dd,bb)) {
    return(TRUE)
  } else {
    return(FALSE) 
  }
}


propper_intersection <- function(aa, bb, cc, dd)
{
  #' Propper intersection between line segments 
  #' 
  #' Checks for a proper intersection between ab and cd; see fig 1.23 in "Computational Geometry in C"
  
  if (Colinear(aa,bb,cc)|
      Colinear(aa,bb,dd)| 
      Colinear(cc,dd,aa)|
      Colinear(cc,dd,bb)
  ) return(FALSE)
  
  return(
    xor(Left(aa,bb,cc), Left(aa,bb,dd)) & xor(Left(cc,dd,aa), Left(cc,dd,bb))
  )
}


impropper_intersection_between <- function(aa,bb,cc)
{
  #' Imropper intersection between line segments 
  #' 
  #' Checks for an improper intersection between ab and cd; see fig 1.24 in "Computational Geometry in C"
  
  if (!Colinear(aa,bb,cc)) return(FALSE)
  
  if (aa[1] != bb[1])
  {
    return(
      ((aa[1] <= cc[1]) & (cc[1] <= bb[1])) | ((aa[1] >= cc[1]) & (cc[1] >= bb[1]))
    )
  } else {
    return(
      ((aa[2] <= cc[2]) & (cc[2] <= bb[2])) | ((aa[2] >= cc[2]) & (cc[2] >= bb[2]))
    )
  }
}


##
## Undocumented helpers
##

Area2 <- function(aa,bb,cc) (bb[1] - aa[1]) * (cc[2] - aa[2]) - (cc[1] - aa[1]) * (bb[2] - aa[2])

Left <- function(aa,bb,cc) Area2(aa, bb, cc) > 0

Right <- function(aa,bb,cc) !Left(aa,bb,cc)

LeftOn <- function(aa,bb,cc) Area2(aa, bb, cc) >= 0 

Colinear <- function(aa,bb,cc) Area2(aa, bb, cc) == 0 

LinearInterpolant <- function(p1,p2) list(slope = (p2[2]-p1[2])/(p2[1]-p1[1]), intercept = p2[2] - p2[1]*(p2[2]-p1[2])/(p2[1]-p1[1]))

PointPosition <- function(pp, halfspace) sign(pp[2] - (halfspace$intercept + halfspace$slope * pp[1]))

HalfspaceIntersection <- function(halfspace1, halfspace2) NULL 