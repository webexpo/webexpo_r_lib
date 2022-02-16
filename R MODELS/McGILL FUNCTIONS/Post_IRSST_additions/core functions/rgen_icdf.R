# Version 0.2 (Apr 2021)
# [distributed]

                                                                               
# Change Log *****************************************************************
#                                                      
#   When updated, look for comments with new_* and modif_*  
#   to rapidly identify new/modified code.                                      
#
#
# Version 0.2 (Apr 2021)
# ----------------------- 
#   The fct uses a completely rewritten version of ref.points (much easier to understand & maintain).
#     That function was relegated to file ref_points.R
#
#
# Version 0.1 (Mar 2021)
# -----------------------  
#   Initial values in the search for the distrn remote endpoints is possible through
#     the argument remote.start
#     - Alternatively, an approximation to the density can be done to find initial values in the above search
#     through the argument dens.approx (1: normal, 2: log-normal, 3: beta, 4: square root inverted gamma)
#       [this argument is ignored if the above-mentioned remote.start is provided]
#     - In the absence of the above two, the function log.f.inv.remote will be used if defined,
#       otherwise a quadratic approximation will be applied locally) around the mode to find good starting points (hence the fct log.f.inv.remote, sometimes difficult to derive, is not necessary anymore)
#
#   The function now returns the remote distrn endpoints,
#     which may prove very useful when used as starting values in the 
#     search for remote endpoints in the next call to rgen.icdf (in an MCMC perspective)
#
#   The code searching for mode & remote limits of the distrn was copied for dens.gen.icdf,
#   but the ICDF section was completely rewritten.
#     - $approx was added to the output of ref.points


rgen.icdf <- function(o, u=runif(1), 
                      start=numeric(0), remote.start=numeric(0), 
                      precision=1e-6, f.ratio.remote=1e-8, 
                      save.objects=character(0), unlink=T)
{             
  # o is a list of:
  #   f, log.f, log.f.prime, log.f.second: functions
  #   A: includes all arguments to the above-mentioned functions
  
    
  # In the search for the distrn remote endpoints, an approximation to the distrn can be done
  # to find *initial values* in that search. The following approximating distrns are available
  # dens.approx = 1 -> normal approx
  #             = 2 -> log-normal approx
  #             = 3 -> beta approx
  #             = 4 -> square root inverted gamma
  #
  # That approximation density will also be used to find initial values in the ICDF section.
  
  
  if (length(save.objects) > 0)
  {
    # save.objects was useful while debugging the development version of this software:
    # it may not be necessary to have it in other languages' versions.
    # It was left here in case we need it again.
    
    save(o, u, precision, start, remote.start, f.ratio.remote, file=save.objects)
  }
      
  
  # !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!**!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
  


  ref <- ref.points(o, start, remote.start=remote.start, epsilon=precision, f.ratio.remote=f.ratio.remote)
    o$A$M <- ref$M


  area <- function(a, b, obj=o){integrate(obj$f, a, b, obj$A)$value}
  
  area.tot <- area(ref$remote[1], ref$remote[2])
  
  
  if (area.tot > diff(ref$remote))
  {
    # There is clearly a higher mode than that found earlier.
    # It would be better to find it before going on.
    
    # We look for a higher mode than the local mode found before by starting the search from
    # the right-hand side remote endpoint (followed by a search starting from the left-hand side 
    # remote endpoint if necessary)
        
    
    # Start from right-side remote
    tmp <- ref.points(o, ref$remote[2], epsilon=precision, find.remote=F, second.guessing.allowed=F)
    
    # Look again for the higher mode, starting from left-side remote, if the above did not succeed
    if (tmp$M < ref$M) tmp <- ref.points(o, ref$remote[1], epsilon=precision, find.remote=F, second.guessing.allowed=F)
    
    if (tmp$M > ref$M)
    {
      # Recalculate remote endpoints
      remote <- ref.points(o, tmp$mode, remote.start=ref$remote, 
                           known=list(mode=tmp$mode, M=tmp$M, lower.local.max=ref$mode), 
                           f.ratio.remote=f.ratio.remote)
                           
      ref <- tmp
        ref$bimodal <- T
        ref$remote <- remote$remote
        ref$approx <- remote$approx
        o$A$M <- ref$M
      
      area.tot <- area(ref$remote[1], ref$remote[2])
    }
  }

  
  
  target <- u*area.tot
  
  bounds <- ref$remote
  
  area.mode <- area(ref$remote[1], ref$mode)
  u.mode <- area.mode / area.tot
  
  
  # Use the mode as a starting point *pour couper la poire en deux*
  
  x <- ref$mode
  
  if (target <= area.mode)
  { 
    bounds[2] <- ref$mode
    a <- area.mode
  }
  else
  {
    bounds[1] <- ref$mode
    target <- target - area.mode
    a <- 0
  }
  
  
  # If the function was called with a suggested distrn approx (through o$dens.approx)
  # then use that distrn quantile as a starting point (if on the right side of the mode)
  
  
  if (o$dens.approx > 0)
  {
    if (any(is.na(ref$approx))) ref$approx <- Dens.approx(ref$mode, o)$approx
  
    if      (o$dens.approx == 1) approx.x <-  qnorm(u, ref$approx$mu, ref$approx$sigma)
    else if (o$dens.approx == 2) approx.x <- qlnorm(u, ref$approx$mu, ref$approx$sigma)
    else if (o$dens.approx == 3) approx.x <-  qbeta(u, ref$approx$alpha, ref$approx$beta)
    else                         approx.x <-  sqrt(1/qgamma(1-u, ref$approx$alpha, ref$approx$beta))
    
    # If the starting point suggested above lies on the right side of the mode (and within bounds), 
    # we accept it; otherwise we'll simply start from the mode
    # (x is currently at the mode)
      
      # If the beta approx is bad, then ref$approx returns NAs and approx.x will also be NA, hence the condition below.
    
    if (!is.na(approx.x))
    {
      if ((u < u.mode & approx.x < x & approx.x > bounds[1]) | (u > u.mode & approx.x > x & approx.x < bounds[2]))
      { 
        x <- approx.x
        a <- area(bounds[1], x)
        
        if (a < target)
        {
          bounds[1] <- x
          target <- target - a
          a <- 0
        }
        else bounds[2] <- x
      }
    }
  }


  continue <- T
  

  while (continue)
  {
    previous.x <- x
    f <- o$f(x, o$A)
    change <- (target - a) / f
    x <- x + change
    
    if (x > bounds[2])
    {
      fp <- f * o$log.f.prime(previous.x, o$A) / area.tot
      t <- Quadratic.solution(fp/2, f, -target)
      
      if (t$negative.delta) x <- previous.x + t$symmetry.axis
      else x <- previous.x + t$right
      
      if (x > bounds[2]) x <- (previous.x + bounds[2]) / 2 
    }
    else if (x < bounds[1])
    {
      fp <- f * o$log.f.prime(previous.x, o$A)  / area.tot      
      t <- Quadratic.solution(fp/2, f, a - target)
      x <- previous.x + t$left
      
      
      if (x < bounds[1]) x <- (previous.x + bounds[1]) / 2
    }
    
    
    converged <- abs(x - previous.x) < precision
    continue <- !converged
    
    
    if (continue)
    {
      a <- area(bounds[1], x)
      
      if (a < target)
      {
        bounds[1] <- x
        target <- target - a
        a <- 0
      }
      else bounds[2] <- x
    }  
  }
  
   
  if (length(save.objects) > 0 & unlink) unlink(save.objects)
  

  list(x=x, ref=list(mode=ref$mode, remote=ref$remote))
} # end of rgen.icdf