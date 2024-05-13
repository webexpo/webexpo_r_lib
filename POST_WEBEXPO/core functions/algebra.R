# Author: Patrick Bélisle
#
# Version 1.0 (Jun 2023)


# -----------------------------------


Quadratic.solution <- function(A, B, C)
{
  # Note: there is another version with a similar name -- quadratic.solution -- below
  delta <- B**2 - 4*A*C
  negative.delta <- delta < 0
  symmetry.axis <- -B/(2*A)

  if (negative.delta)
  {
    left <- NA
    right <- NA
  }
  else
  {
    m <- sqrt(delta)/(2*abs(A))
    left  <- symmetry.axis - m
    right <- symmetry.axis + m
  }

  list(left=left, right=right, negative.delta=negative.delta, symmetry.axis=symmetry.axis)
} # end of Quadratic.solution


# modif_0.11 (added target)
quadratic.solution <- function(theta, target=0, l=-Inf, u=Inf)
{
  # Note: there is another version with a similar name -- Quadratic.solution -- above

  C <- theta[1] - target
  B <- theta[2]
  A <- theta[3]

  delta <- B^2 - 4*A*C

  if (delta < 0)
  {
    soln <- NA
  }
  else
  {
    soln <- (-B + c(-1,1)*sqrt(delta))/2/A
    soln <- soln[soln > l & soln < u]
    if (length(soln) == 0) soln <- NA
  }

  soln
} # end of quadratic.solution


real.cubic.roots <- function(theta, target=0, l=numeric(0), u=numeric(0), epsilon=1e-8) 
{ 
  a <- theta[4] 
  b <- theta[3] 
  c <- theta[2] 
  d <- theta[1] 
  
  if (a == 0) 
  { 
    if (length(l) == 0) l <- -Inf 
    if (length(u) == 0) u <- +Inf 
    roots <- quadratic.solution(b, c, d, l=l, u=u) 
  } 
  else 
  { 
    roots <- numeric(0) 
    f <- function(x, k=theta){sum(k*c(1,x^seq(3)))} 
    f.prime <- function(x, k=theta){sum(k[-1]*seq(3)*c(1,x^seq(2)))} 

    any.soln <- list(L=FALSE, M=FALSE, R=FALSE)
    start <- list(L=numeric(0), M=numeric(0), R=numeric(0)) 
    inflexion.point <- list(x=NA, y=NA) 
    local.minmax <- list(x=rep(NA, 2), y=rep(NA, 2)) 
    
    inflexion.point$x <- -b/3/a 
    inflexion.point$y <- f(inflexion.point$x) 
    
    A <- 3*a 
    B <- 2*b 
    C <- c 
    delta <- B^2 - 4*A*C 
    
    any.middle.section <- delta > 0 
    
    if (any.middle.section) 
    { 
      local.minmax$x <- inflexion.point$x + c(-1, 1)*sign(A)*sqrt(delta)/2/A 
      for (i in 1:2) local.minmax$y[i] <- f(local.minmax$x[i]) 
    } 
    else 
    { 
      local.minmax$x <- rep(inflexion.point$x, 2) 
      local.minmax$y <- rep(inflexion.point$y, 2) 
    } 
    
    
    any.soln$L <- (a > 0 & local.minmax$y[1] > target) | (a < 0 & local.minmax$y[1] < target) 
    any.soln$R <- (a > 0 & local.minmax$y[2] < target) | (a < 0 & local.minmax$y[2] > target) 
    
    if (any.middle.section) 
    { 
      tmp <- sort(local.minmax$y) # vector of length 2 
      any.soln$M <- tmp[1] <= target & target <= tmp[2] 
      if (any.soln$M) start$M <- inflexion.point$x 
    } 
    
    
    # find points with +/- 1 slope (depending on sign of a)
    # (if necessary)
    
    pm1slope <- list(x=numeric(0), helpfull=any.middle.section) 
    
    if (!pm1slope$helpfull) 
    { 
      s <- f.prime(inflexion.point$x) 
      pm1slope$helpfull <- abs(s) < 1 
    } 
    
    if (pm1slope$helpfull) 
    { 
      s <- sign(a) 
      pm1slope$x <- (-B + c(-1,1)*s*sqrt(B^2-4*A*(C-1)))/2/A 
    }  
    
    
    # Look for root in the left(L) section
    
    start$L <- local.minmax$x[1] 
    
    if (any.soln$L & length(u) > 0) 
    { 
      if (u < local.minmax$x[1]) 
      { 
        y <- f(u) 
        any.soln$L <- (a > 0 & y >= target) | (a < 0 & y <= target) 
        start$L <- u 
      } 
    } 
    
    
    if (any.soln$L) 
    { 
      slope <- f.prime(start$L) 
      if (abs(slope) < 1) start$L <- pm1slope$x[1] 
      
      if (length(l) > 0) 
      { 
        if (l > local.minmax$x[1]) 
        { 
          any.soln$L <- FALSE
        } 
        else 
        { 
          y <- f(l) 
          any.soln$L <- (a > 0 & y <= target) | (a < 0 & y >= target) 
          if (l > start$L) start$L <- l 
        } 
      } 
    } 
    
    
    if (any.soln$L) 
    { 
      # run Newton-Raphson algorithm to find local root
      
      continue <- TRUE
      x <- start$L 
      y <- f(x) 
      
      while (continue) 
      { 
        change <- (y-target)/f.prime(x) 
        x <- x - change 
        y <- f(x) 
        continue <- abs(y-target) > epsilon 
      } 
      
      roots <- x 
    } 
    
    
    # Look for root in the middle(M) section
    
    
    if (any.soln$M & length(l) > 0) 
    { 
      any.soln$M <- l < local.minmax$x[2] 
      
      if (any.soln$M) 
      { 
        y <- f(l) 
        any.soln$M <- (a > 0 & y >= target) | (a < 0 & y <= target) 
      } 
    } 
    

    if (any.soln$M & length(u) > 0) 
    { 
      any.soln$M <- u > local.minmax$x[1] 
      
      if (any.soln$M) 
      { 
        y <- f(u) 
        any.soln$M <- (a > 0 & y <= target) | (a < 0 & y >= target) 
      } 
    }     
    
    
    if (any.soln$M) 
    { 
      continue <- TRUE
      x <- start$M 
      y <- f(x) 
      
      while (continue) 
      { 
        change <- (y-target)/f.prime(x) 
        x <- x - change 
        y <- f(x) 
        continue <- abs(y-target) > epsilon 
      } 
      
      roots <- c(roots, x) 
    } 
    
    
    # Look for root in the right(R) section
    
    start$R <- local.minmax$x[2] 
    
    if (any.soln$R & length(l) > 0) 
    { 
      if (l > local.minmax$x[2]) 
      { 
        y <- f(l) 
        any.soln$R <- (a > 0 & y <= target) | (a < 0 & y >= target) 
        start$R <- l 
      } 
    } 
    
    
    if (any.soln$R) 
    { 
      slope <- f.prime(start$R) 
      if (abs(slope) < 1) start$R <- pm1slope$x[2] 
      
      if (length(u) > 0) 
      { 
        if (u < local.minmax$x[2]) 
        { 
          any.soln$R <- FALSE
        } 
        else 
        { 
          y <- f(u) 
          any.soln$R <- (a > 0 & y >= target) | (a < 0 & y <= target) 
          if (u < start$R) start$R <- u 
        } 
      } 
    } 
    
    
    if (any.soln$R) 
    { 
      # run Newton-Raphson algorithm to find local root
      
      continue <- TRUE
      x <- start$R 
      y <- f(x) 
      
      while (continue) 
      { 
        change <- (y-target)/f.prime(x) 
        x <- x - change 
        y <- f(x) 
        continue <- abs(y-target) > epsilon 
      } 
      
      roots <- c(roots, x) 
    }  
  } 
  
  
  roots 
} # end of real.cubic.roots 


rotate <- function(x, y, x0, y0, theta) 
{ 
  # find coordinates of a set of points (x,y) after rotation of an angle theta around the pivot (x0, y0)

  # Args:
  # -----
  # x, y: two vectors of same length
  # (x0, y0): point around which the rotation is performed
  # theta: rotation angle

  rho <- sqrt((x-x0)^2+(y-y0)^2) 
  phi <- atan2(y-y0, x-x0) 
  new.phi <- phi + theta 
  new.x <- x0 + rho*cos(new.phi) 
  new.y <- y0 + rho*sin(new.phi) 

  list(x=new.x, y=new.y) 
} # end of rotate 
