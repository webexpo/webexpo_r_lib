# Version 1.2 (Dec 2023)
#         Shared/distributed: no
#		      Last shared version: 1.1
# -----------------------------------


                                                                               
  # Change Log *****************************************************************
  #
  #
  # Version 1.2 (Dec 2023)
  # ----------------------
  #
  #
  # Version 1.1 (Dec 2023)
  # ----------------------
  #   Added an alternate starting point to ref.points search (from left side)
  #   to avoid starting from a point with null slope.
  #
  #
  # Version 1.0 (Jun 2023)
  # -----------------------
  #   T/F tokens were changed for TRUE/FALSE.
  #
  #
  # Version 0.26 (Mar 2021)
  # -----------------------  
  #   Initial values in the search for the distrn remote endpoints is now possible through
  #     the argument remote.start
  #     - Alternatively, an approx to the density can be done to find initial values in the above search
  #     through the argument dens.approx (1: normal, 2: log-normal, 3: beta)
  #       [this argument is ignored if the above-mentioned remote.start is provided]
  #     - In the absence of the above two, the function log.f.inv.remote will be used if defined,
  #       otherwise a quadratic approximation will be applied locally) around the mode to find good starting points (hence the fct log.f.inv.remote, sometimes difficult to derive, is not necessary anymore)
  #
  #   The argument return.remote was added, to include the remote endpoints in the returned object:
  #     that may prove very useful as using the actual remote endpoint values as starting values in the 
  #     search for remote endpoints in the next call to dens.gen.icdf (in an MCMC perspective)
  #
  #   (Also added a convergence criterion for when the remote limit is at the edge of the variable's domain)
  #
  #
  # Version 0.25 (Feb 2021)
  # -----------------------
  #   Added inestimable.lower.limit to a call to secured.NR.search
  #   Added na.rm to a call to min()
  #   Added a search for real mode on both sides of (potentially local) first mode found when potientally.bimodal=T
  #   The fct higher.local.max was revisited, removing the need for a difficult-to-derive o$start() function
  #   The dimension o$potentially.bimodal is now obsolete: the potential bimodality of any posterior distrn will always be assessed by dens.gen.icdf
  #
  # Version 0.23
  # ------------
  #   Relaxed the precautions added in version 0.22
  #
  # Version 0.22
  # ------------
  #   Added some protections following changes in version 0.21
  #
  # Version 0.21
  # ------------
  #   Brought slight modifications to dens.gen.icdf
  #
  # Version 0.20
  # ------------
  #   Added a convergence criterion to secured.NR.search
  #
  # Version 0.18
  # ------------
  #   Added argument RData to each function calling dens.gen.icdf
  #     -> useful at development stage only, for monitoring parameters values when dens.gen.icdf is crashing
  #
  #
  # Version 0.15 (March 2019)
  # -------------------------
  #
  # - Using a quadratic approximation (around the mode) to log.f to find starting points in the search
  #   for remote limits (in fct ref.points) when log.f.inv.remote is not defined
  # - Added $M to A after finding the mode (= the value of log.f at its mode)
  #  (it may be used by log.f.inv.remote to find a starting point in the search for 'remote' limits)
  # - Also moved the following functions from file fcts.R to this file (that is, they are now embedded in dens.gen.icdf)
  #   [these functions remain unchanged, though]
  #   -> Potentially.bimodal
  #   -> quadratic.solution
  #
  # Version 0.14 (Jan 2018)
  # -----------------------
  #    [not documented]
  #
  # Version 0.13 (Aug 2017)
  # -----------------------
  #
  # -  modified fct split.starting.values.leftRight
  # -  added function unexplored.range
  # -  added a convergence criterion to secured.NR.search
  #
  #
  # -  inestimable.logf.at.limits &
  #    inestimable.logfprime.at.limits
  #    were dropped as arguments to be provided by user (too much error-prone)
  #    -> they are now computed within function
  #
  # - function Potentially.bimodal was added, to define default value for potentially.bimodal
  #   -> potentially.bimodal could still be passed through o$potentially.bimodal dimension
  #      as it was in earlier versions, or directly through potentially.bimodal function's argument
  #
  # Version 0.12 (Jun 2017)
  # -----------------------
  #
  # - inestimable.lower.limit (logical scalar) was changed for
  #    i) inestimable.logf.at.limits      (logical vector of length 2)
  #   ii) inestimable.logfprime.at.limits (logical vector of length 2)
  #   > look for inestimable.logf in code
  #
  # - math.lower.limit (scalar) was changed for math.limits (vector of length 2)
  #   > look for math.limits in code
  #
  # - a few more modifications / additions
  #
  # changed NR.max.iter to 200


dens.gen.icdf <- function(o, A=o$A, range=numeric(0), u=runif(1), start=mean(range), precision=1e-6, inestimable.lower.limit=FALSE, NR.max.iter=200, save.objects=character(0), unlink=TRUE, remote.start=numeric(0), dens.approx=0, return.remote=FALSE)
{ 
  # ICI remettre save.objects=temp(my.save.objects, time.stamp=T) pour simulation
             
  # o is a list of:
  #   f, log.f, log.f.prime, log.f.second: functions
  
  #   start:                               function
  #   logNormal.distrn:                    logical
  #   A:                                   the list of functions arguments
  
  # In earlier versions, o also included  the function
  #   log.f.inv.remote
  # which is not necessary anymore (but will be used if older versions calling dens.gen.icdf provide it)
  
  # In the search for the distrn remote endpoints, an approximation to the distrn can be done
  # to find *initial values* in that search. The following approximating distrns are available
  # dens.approx = 1 -> normal approx
  #             = 2 -> log-normal approx
  #             = 3 -> beta approx
  
  
  if (length(save.objects) > 0) 
  { 
    # save.objects was useful while debugging the development version of this software:
    # it may not be necessary to have it in other languages' versions.
    # It was left here in case we need it again.
    
    save(o, A, range, u, precision, inestimable.lower.limit, start, remote.start, dens.approx, 
         file=save.objects) 
  } 
      
  
  ping.pong.tie.breaker <- function(area, A, start, area.x, target,  
                                    math.lower.limit, ref.range, inestimable.lower.limit, epsilon,  
                                    distrn.leftSide=FALSE, hp.mult=ifelse(distrn.leftSide, -1, 1),
                                    max.count=100) 
  { 
    continue <- TRUE

    count <- 0 
    x <- start 
    visited.x <- x 

    while (continue) 
    { 
      hp <- A$f(x, A) * hp.mult 
      change <- (target - area.x)/hp 
      x <- x + change 
      
      if (x < ref.range[1])  
      { 
        x <- ifelse(ref.range[1] == math.lower.limit & inestimable.lower.limit, (x - change + math.lower.limit) / 2, ref.range[1]) 
      } 
      
      
      area.x <- smoothed.area.estimate(area, x, A, math.lower.limit, ref.range, inestimable.lower.limit, hp.mult=hp.mult) 

      count <- count + 1 
      converged <- abs(area.x - target) < epsilon 
      visited.x <- c(visited.x, x) 
      caught.in.loop <- any(duplicated(visited.x) & is.finite(visited.x)) 
      continue <- !converged & count <= max.count & !caught.in.loop 
    } 
    
    
    if (!converged & caught.in.loop) 
    { 
      # We have found the series of points that are repeatedly visited:
      # recalibrate and try Newton-Raphson again
    
      w <- which(visited.x == x) 
      w <- min(w) # repeated value (loop) start on visited.x[w] 
      loop.x <- visited.x[seq(from=w+1, to=count+1)] 

      if (distrn.leftSide) x <- max(loop.x) 
      else x <- min(loop.x) 

      area.x <- smoothed.area.estimate(area, x, A, math.lower.limit, ref.range, inestimable.lower.limit, hp.mult=hp.mult) 

      if (distrn.leftSide) A$upper.limit <- x 
      else A$lower.limit <- x 

      target <- target - area.x 
      area.x <- 0 

      count <- 0 
      continue <- TRUE

      while (continue) 
      { 
        hp <- A$f(x, A) * hp.mult 
        change <- (target - area.x)/hp 
        x <- x + change 
        area.x <- area(x, A) 
        count <- count + 1 
        converged <- abs(area.x - target) < epsilon 
        continue <- !converged & count <= max.count 
      } 
    } 


    if (!converged) stop("Newton-Raphson algorithm did not converge. Sorry.\n") 

    x 
  } # end of ping.pong.tie.breaker 
 
  
  ref.points <- function(o, A, start, range, inestimable.lower.limit=FALSE, f.ratio.remote=1e-8, epsilon=1e-6, max.niter=NR.max.iter, find.remote=TRUE)
  { 
    # Arguments:
    # ----------------------------------------------
    # o: same nature as o in dens.gen.icdf arguments
    # start: EITHER a vector of length 2 [two potential starting points]
    #        or a scalar
    # f.ratio.remote: scalar
    # epsilon: scalar
    # range: vector of length 2

    higher.local.max <- function(o, A, dip.x, range, hs, epsilon)  
    { 
      # Find a starting point (for N-R to work smoothly) on each side of dip.x

      #i) on right side

      x <- dip.x 
      fp <- 0 
      fs <- hs 
      continue <- TRUE

      while (continue) 
      { 
        change <- 1/fs 
        x <- x + change 
        if (x > range[2]) x <- (x - change + range[2]) / 2 

        fp <- o$log.f.prime(x, A) 
        fs <- o$log.f.second(x, A) 
        continue <- fp > 0 & fs > 0 
      } 

      start.right <- x 


      #ii) on left side

      x <- dip.x 
      fp <- 0 
      fs <- hs 
      continue <- TRUE

      while (continue) 
      { 
        change <- - 1/fs 
        x <- x + change 
        if (x < range[1]) x <- (x - change + range[1]) / 2 

        fp <- o$log.f.prime(x, A) 
        fs <- o$log.f.second(x, A) 
        continue <- fp < 0 & fs > 0 
      } 
       
      start.left <- x 

      
      # Run pure Newton-Raphson to find local max on each side
      
      # i) left side

      x <- start.left 
      continue <- TRUE
      hp <- o$log.f.prime(x, A) 

      while (continue) 
      { 
        change <- hp/o$log.f.second(x, A) 
        x <- x - change 
        hp <- o$log.f.prime(x, A)  
        continue <- abs(hp) > epsilon 
      } 

      local.mode.left <- x 
      

      # ii) right-side

      x <- start.right 
      continue <- TRUE
      hp <- o$log.f.prime(x, A) 

      while (continue) 
      { 
        change <- hp/o$log.f.second(x, A) 
        x <- x - change 
        hp <- o$log.f.prime(x, A)  
        continue <- abs(hp) > epsilon 
      } 

      local.mode.right <- x 
      
      h.left  <- o$log.f(local.mode.left, A) 
      h.right <- o$log.f(local.mode.right, A) 
      
      x <- ifelse(h.left > h.right, local.mode.left, local.mode.right) 
      x 
    } # end of higher.local.max 

    
    max.fastTrack <- function(o, A, lower.limit, rightside.x, rightside.hp, epsilon=1e-4, epsilon.x=1e-20) 
    {     
      x <- rightside.x 
      h <- o$log.f(x, A) 
      hp <- rightside.hp 
        
      b <- list(x=c(NA, x), h=c(NA, h), hp=c(NA, hp)) 
      
      x <- (x + lower.limit) / 2 
      continue <- TRUE

 
      # Find a point with positive slope (that is, to the left of the mode)      
      while (continue) 
      { 
        hp <- o$log.f.prime(x, A) 
        
        if (hp > 0) 
        { 
          continue <- is.infinite(hp) 
          
          if (continue) 
          { 
            next.x <- (3*x - lower.limit) / 2 
            lower.limit <- x 
            
            continue <- abs(next.x - x) > epsilon.x 
            if (continue) x <- next.x 
          } 
        } 
        else 
        { 
          x <- (x + lower.limit) / 2 
        }  
      }  
      
      
      hp <- o$log.f.prime(x, A) 
      converged <- is.finite(hp) & hp > 0 
      continue <- converged 
      
      
      while (continue) 
      { 
        side <- ifelse(hp > 0, 1, 2) 
        b$x[side] <- x 
        b$h[side] <- o$log.f(x, A) 
        b$hp[side] <- hp 
      
        b.diff <- diff(b$hp) 
          
        if (b.diff == 0) x <- mean(b$x) 
        else  
        { 
          x <- (diff(b$hp*b$x) - diff(b$h)) / b.diff 
          if (x <= b$x[1] | x >= b$x[2]) x <- mean(b$x) 
        }  
          
        hp <- o$log.f.prime(x, A) 
        continue <- abs(hp) > epsilon 
      } 
      
      h <- o$log.f(x, A) 
      
      list(converged=converged, x=x, h=h) 
    } # end of max.fastTrack 
    

    split.starting.values.leftRight <- function(o, target, A, range, start.values) 
    { 
      point.with.slope.closest.to.secant <- function(start, o, A, target, range, right.side=TRUE)
      { 
        # y <- o$log.f(start, A) 
        y <- rep(NA, length(start)) 
        for (i in seq(along=start)) y[i] <- o$log.f(start[i], A) 
        
        sup <- start[y>target] 
        inf <- start[y<target] 

        if (right.side) 
        { 
          if (length(sup) > 0) 
          { 
            sup <- max(sup) 
            range[1] <- sup 
          } 


          if (length(inf) > 0) 
          { 
            inf <- min(inf) 
            range[2] <- inf 
          } 
        } 
        else 
        { 
          if (length(sup) > 0) 
          { 
            sup <- min(sup) 
            range[2] <- sup 
          } 

          if (length(inf) > 0) 
          { 
            inf <- max(inf) 
            range[1] <- inf 
          } 
        } 


        x <- c(inf, sup) # one or two values 

        if (length(x) > 1) 
        { 
          # pick the one with the slope closest to straight line cutting area in two equal parts
          
          #slope <- o$log.f.prime(x, A)
          slope <- rep(NA, length(x)) 
          for (i in seq(along=x)) slope[i] <- o$log.f.prime(x[i], A) 
          
          angle <- atan(slope) 

          #y <- o$log.f(x, A)
          y <- rep(NA, length(x)) 
          for (i in seq(along=x)) y[i] <- o$log.f(x[i], A) 
          
          secant.slope <- diff(y)/diff(x) 
          secant.angle <- atan(secant.slope) 

          angle.diff <- abs(angle-secant.angle) 
          w <- which.min(angle.diff)[1] 
          x <- x[w] 
        } 


        list(x=x, range=range) 
      } # end of point.with.slope.closest.to.secant 


      start.values <- start.values[!is.na(start.values)]

      ls.start <- start.values[start.values <= A$mode & start.values >= range[1]]
      rs.start <- start.values[start.values >= A$mode & start.values <= range[2]]

      start <- list(left=NA, right=NA) 
      range <- list(left=c(range[1], A$mode), right=c(A$mode, range[2])) 

      tmp <- point.with.slope.closest.to.secant(ls.start, o, A, target, range$left, right.side=FALSE)
      start$left <- tmp$x 
      range$left <- tmp$range 

      tmp <- point.with.slope.closest.to.secant(rs.start, o, A, target, range$right) 
      start$right <- tmp$x 
      range$right <- tmp$range 


      if (length(start$left) == 0)  start$left  <- NA 
      if (length(start$right) == 0) start$right <- NA 

      list(start=start, range=range)
    } # end of split.starting.values.leftRight 

    
    # --------------------------------------------------------------------------
    # If range is finite and log.f not inestimable at lower end, 
    # then we first look for values at both ends before searching for a mode
    
    mode <- list(x=numeric(0), h=numeric(0), found=FALSE)

    bimodal <- FALSE

    
    if (all(is.finite(range)) & !inestimable.lower.limit) 
    { 
      log.f.prime <- rep(NA, 2) 
      for (i in 1:2) log.f.prime[i] <- o$log.f.prime(range[i], A) 
      
      # If slope (f') has the same sign at both ends, then there is no mode between the 2 endpoints (assuming a unimodal distrn)
      # and we then assume that the highest end value for log.f is hence the local mode
      
      log.f.prime.sign <- sign(log.f.prime) 
      
      if (prod(log.f.prime.sign) > 0) 
      { 
        log.f <- rep(NA, 2) 
        for (i in 1:2) log.f[i] <- o$log.f(range[i], A) 
        w <- which.max(log.f) 
        mode <- list(x=range[w], h=log.f[w], found=TRUE)
      } 
    } 
 
    
    if (!mode$found) 
    { 
      # Choose a starting point among suggested starting points, if more than one point was suggested
    
      if (length(start) > 1) 
      { 
        start <- c(start, mean(start)) # Add the middle point to the two potential starting points 
        start.len <- length(start) 
        log.f <- rep(NA, start.len) 
        for (i in 1:start.len) log.f[i] <- o$log.f(start[i], A=A) 
        w <- which.max(log.f) 
        start <- start[w] 
      } 
    

      # Find mode by Newton-Raphson
        
      o.mode <- list(h=o$log.f.prime, hp=o$log.f.second, range=range) 
      nr <- secured.NR.search(o.mode, A, start, epsilon, range=range, max.niter=max.niter, inestimable.lower.limit=inestimable.lower.limit)
       
      if (nr$converged) 
      { 
        hs <- o$log.f.second(nr$x, A) 
        if (hs > 0) 
        { 
          # We have found a mimimum point (between two local modes):
          # we redo the search on both sides of this dip
          mode$x <- higher.local.max(o, A, nr$x, range, hs, epsilon)
          bimodal <- TRUE
        } 
        else 
        { 
          mode$x <- nr$x 
        } 
        
        mode$found <- TRUE
      } 
      else 
      { 
        x <- nr$bounds$x 
        h <- nr$bounds$h 
        h.on.both.sides <- prod(sign(h)) < 0 
        diff.x <- abs(diff(x)) 
        if (diff.x < epsilon & h.on.both.sides) 
        {  
          mode$x <- mean(x) 
          mode$found <- TRUE
        } 
        else  
        { 
          converged <- FALSE
          
          if (diff(nr$bounds$x) < epsilon & all(nr$bounds$h < 0)) 
          { 
            if (is.finite(range[1])) 
            { 
              # We try a little bit further to find the actual mode, 
              # since it looks like we are almost there...
            
              converged <- TRUE
              tmp <- max.fastTrack(o, A, range[1], nr$bounds$x[1], nr$bounds$h[1]) # yes: 4th argument is bounds$h (and not hp!) 
            
              if (tmp$converged) 
              {  
                mode$x <- tmp$x 
                mode$found <- TRUE
              } 
              else 
              { 
                mode$x <- range[1] 
                mode$found <- !inestimable.lower.limit  
              } 
            } 
            else 
            { 
              x <- nr$bounds$x[1] 
              h <- nr$bounds$h[1] 
              tmp <- abs(h/o$log.f.second(x, A)) 
              if (tmp < 1e-12) 
              {  
                mode <- list(x=x, h=h, found=TRUE)
                converged <- TRUE
              } 
            } 
          } 
        
          if (!converged) stop("Algorithm did not converge.\n", save.objects, '\n') 
        } 
      } 
      
      
      if (mode$found) mode$h <- o$log.f(mode$x, A=A) 
    } 
    

    if (find.remote) 
    { 
      log.f.ratio.remote <- log(f.ratio.remote)
      if (mode$found) target <- mode$h + log.f.ratio.remote 
      A$mode <- mode$x # some 'remote.*' functions need that information 


      # Find remote points on both sides of mode
    
      if (mode$found) 
      { 
        x <- rep(NA, 2) 
        f <- rep(NA, 2) 
        fp <- rep(NA, 2) 

        fs <- o$log.f.second(mode$x, A=A) # scalar 
        d <- 10/abs(fs) 


        # remote.start <- o$log.f.inv.remote(target, A=A) # vector
        
        
        if (length(remote.start) == 0) 
        { 
          if (dens.approx == 1) 
          { 
            # Normal distrn approx
            approx.sd <- sqrt(-1/fs) 
            remote.start <- mode$x + 3*c(-1,1)*approx.sd 
            
            if      (remote.start[1] < range[1]) remote.start[1] <- (mode$x + range[1]) / 2 
            else if (remote.start[2] > range[2]) remote.start[2] <- (mode$x + range[2]) / 2 
          } 
          else if (dens.approx == 2) 
          { 
            # log-Normal distrn approx
            approx.sd2 <- -1/fs/mode$x 
            approx.mu <- log(mode$x) - 2*approx.sd2 
            remote.start <- qlnorm(c(0.001, 0.999), approx.mu, sqrt(approx.sd2)) 
          } 
          else if (dens.approx == 3) 
          { 
            # Beta distrn approx
            approx.u <- (1-mode$x)**2 
            approx.v <- mode$x**2 
            approx.A <- matrix(c(approx.u, approx.v, mode$x-1, mode$x), byrow=TRUE, ncol=2)
            approx.y <- c(-fs*approx.u*approx.v+approx.u+approx.v, 2*mode$x-1) 
            approx.parms <- as.vector(solve(approx.A) %*% approx.y) 
            remote.start <- qbeta(c(0.001, 0.999), approx.parms[1], approx.parms[2]) 
          } 
          else if (!is.na(match("log.f.inv.remote", names(o)))) 
          { 
            remote.start <- o$log.f.inv.remote(target, A=A) 
          } 
          else 
          { 
            # Quadratic approx around the mode (locally)
            approx.m <- sqrt(-3/fs) 
            remote.start <- mode$x + c(-1,1)*approx.m 
            
            if      (remote.start[1] < range[1]) remote.start[1] <- (mode$x + range[1]) / 2 
            else if (remote.start[2] > range[2]) remote.start[2] <- (mode$x + range[2]) / 2 
          } 
        } 
                
        
        # start <- split.starting.values.leftRight(o, target, A, range, remote.start) # vector of length 2
        search <- split.starting.values.leftRight(o, target, A, range, remote.start) # a list with dimensions $start and $range
        start <- c(search$start$left, search$start$right)

        mode.on.border <- range == mode$x
        j.seq <- which(!mode.on.border & is.na(start))


        for (j in j.seq) 
        {
          if (is.na(start[j])) 
          { 
            dir <- ifelse(j == 1, -1, 1) 

            tmp <- mode$x + dir*d 

            accepted <- tmp > range[1] & tmp < range[2] 
            if (!accepted) 
            { 
              my.d <- d 
              while (!accepted) 
              { 
                my.d <- my.d/2 
                tmp <- mode$x + dir*my.d 
                accepted <- tmp > range[1] & tmp < range[2] 
              } 
            } 


            fp[1] <- o$log.f.prime(tmp, A=A) 
            accepted <- is.finite(fp[1]) 

            while (!accepted) 
            { 
              tmp <- (tmp + mode$x) / 2 
              fp[1] <- o$log.f.prime(tmp, A=A) 
              accepted <- is.finite(fp[1]) 
            } 

            a <- fp[1]/d # estimated rate of slope acceleration 
            d2 <- sqrt(2*abs(log(f.ratio.remote)/a)) # estimated distance we need to get from mode,  
                                         # at the rate above, to reach the f.ratio.remote distance

            x2 <- mode$x + dir*d2 
            accepted <- x2 > range[1] & x2 < range[2] 

            while (!accepted) 
            { 
              d2 <- d2/2 
              x2 <- mode$x + dir*d2 
              accepted <- x2 > range[1] & x2 < range[2] 
            } 

            fp2 <- o$log.f.prime(x2, A=A)  
            accepted <- is.finite(fp2) 

            while (!accepted) 
            { 
              d2 <- d2/2 
              x2 <- mode$x + dir*d2 
              fp2 <- o$log.f.prime(x2, A=A)  
              accepted <- is.finite(fp2)           
            } 

            x <- c(tmp, x2) 
            f <- rep(NA, 2) 
            for (i in 1:2) f[i] <- o$log.f(x[i], A=A) 
            fp[2] <- fp2  

            # Approximating the shape of log.f.prime by a straight line, 
            # we can estimate the point where the distance
            # prescribed by f.ratio.remote is reached 

            X <- matrix(c(x, 1, 1), ncol=2) # 2 x 2 matrix 
            theta <- as.vector(solve(X) %*% fp) 

            qA <- theta[1]  
            qB <- theta[2]  

            inner.side <- ifelse(xor(j==1, x[1]>x[2]), 2, 1) 
            inner.x <- x[inner.side] 

            qC <- -qA*inner.x^2 - qB*inner.x - (target - f[inner.side]) 

            if (j == 1) 
            { 
              l <- range[1] 
              u <- mode$x 
            } 
            else 
            { 
              l <- mode$x 
              u <- range[2] 
            } 


            roots <- quadratic.solution(c(qC, qB, qA), l=l, u=u)

            if (length(roots) == 1) 
            { 
              if (!is.na(roots))  
              { 
                start[j] <- roots 
              } 
              else 
              { 
                start[j] <- x[3 - inner.side] 
              } 
            } 
            else if (j == 1) 
            { 
              start[j] <- max(roots) 
            } 
            else 
            { 
              start[j] <- min(roots) 
            } 
          } 
        } 


        # Suggest another starting point (on both sides)
        # based on the 2nd degree polynomial approximation (Taylor series devpmt) of log.f

        if (!any(mode.on.border)) 
        { 
          delta <- sqrt(2*log.f.ratio.remote/fs) 

          # Left side

          x <- mode$x - delta 

          if (x > range[1]) 
          { 
            alt.start <- c(start[1], x) 
            f <- rep(NA, 2) 
            for (i in 1:2) f[i] <- abs(o$log.f(alt.start[i], A) - target) 
            d <- abs(f - target) 
            w <- which.min(d)[1] 
            start[1] <- alt.start[w] 
          } 

          # Right side

          x <- mode$x + delta 

          if (x < range[2]) 
          { 
            alt.start <- c(start[2], x) 
            f <- rep(NA, 2) 
            for (i in 1:2) f[i] <- abs(o$log.f(alt.start[i], A) - target) 
            d <- abs(f - target) 
            w <- which.min(d)[1] 
            start[2] <- alt.start[w] 
          } 
        } 


        oh <- list(h=o$log.f, hp=o$log.f.prime, hs=o$log.f.second)

        if (is.na(start[1])) 
        { 
          remote.left <- mode$x 
        } 
        else 
        { 
          nr <- secured.NR.search(oh, A, start[1], epsilon, target=target,  
                   range=search$range$left, max.point=mode, 
                   inestimable.lower.limit=inestimable.lower.limit, max.niter=max.niter)

          if (nr$converged) 
          { 
            remote.left <- nr$x 
          } 
          else if (is.finite(range[1])) 
          { 
            remote.left <- range[1] 
          } 
          else 
          { 
            stop('fin tempo ICI\n') 
          } 
        } 


        if (is.na(start[2])) 
        { 
          remote.right <- mode$x 
        } 
        else 
        { 
          nr <- secured.NR.search(oh, A, start[2], epsilon, target=target,  
                   range=search$range$right, max.point=mode, max.niter=max.niter) 

          if (nr$converged) 
          { 
            remote.right <- nr$x 
          } 
          else 
          { 
            stop("Algorithm did not converge.\n", save.objects, '\n') # should not happen 
          } 
        } 
      }  
      else 
      {
        # mode was not found
        remote.left <- NA 
        remote.right <- NA 
      }   
    } 
    else 
    {
      remote.left <- NA 
      remote.right <- NA 
      log.f.ratio.remote <- NA 
    }  
    

    list(mode=mode, remote=c(remote.left, remote.right), log.f.ratio.remote=log.f.ratio.remote, bimodal=bimodal)   
  } # end of ref.points 
      
    
  secured.NR.search <- function(o, A, start, epsilon, target=0, range=numeric(0),  
                                max.point=list(x=numeric(0), h=numeric(0)), 
                                inestimable.lower.limit=FALSE, max.niter=200, expected.hpsign=-1, cat.x=FALSE)
  {             
    # Throughout this function, 'bounds' is a list with the following elements/dimensions:
    # x, h, hp, range:   numeric (vectors of length 2)
    # higher.side:       numeric (1)
    # include.soln:      logical (1)

    
    cubic.extrapolation <- function(target, bounds, monotonic, range) 
    { 
      # target: scalar
      # bounds: list
     
      higher.side <- bounds$higher.side  

      # Approximate h by a cubic curve
      # (with appropriate values and slopes at both ends)
      
      X <- matrix(c(1, 1, 0, 0, 
                    bounds$x, 1, 1, 
                    bounds$x^2, 2*bounds$x,  
                    bounds$x^3, 3*bounds$x^2), nrow=4) 
      
      y <- c(bounds$h, bounds$hp) 
      
      theta <- Xinverse.y(X, y)  
      
      
      if (any(is.nan(theta)) | any(is.infinite(theta))) 
      { 
        use.linear.extrapolation <- TRUE
      }      
      else 
      { 
        tmp <- real.cubic.roots(theta, target) 
        use.linear.extrapolation <- length(tmp) == 0
      }                   
      
      
      if (use.linear.extrapolation) 
      { 
        # use a linear extrapolation instead since X is probably not invertible
        X <- matrix(c(bounds$x, 1, 1), nrow=2) # 2 x 2 matrix 
        theta <- Xinverse.y(X, bounds$h) 
        tmp <- (target - theta[2]) / theta[1] 
      } 


      inbound.solns <- tmp[tmp > bounds$x[1] & tmp <= bounds$x[2]] 
      inrange.solns <- tmp[tmp > range[1] & tmp <= range[2]] 
      
      
      if (bounds$include.soln) 
      { 
        m <- length(inbound.solns) 
        
        if (m == 0) 
        { 
          # Use a linear extrapolation instead
          X <- matrix(c(bounds$x, 1, 1), nrow=2) # 2 x 2 matrix 
          theta <- Xinverse.y(X, bounds$h) 
          soln <- (target - theta[2]) / theta[1] 
        } 
        else if (m == 1) 
        { 
          soln <- inbound.solns 
        } 
        else  
        { 
          soln <- mean(inbound.solns) 
        } 
      } 
      else 
      { 
        tmp <- inrange.solns 
        m <- length(tmp) 
        
        if (m == 0) 
        { 
          soln <- NA 
        } 
        else if (m == 1) 
        { 
          soln <- tmp 
        } 
        else if (monotonic) 
        { 
          if (higher.side == 1) 
          { 
            tmp <- tmp[tmp < bounds$x[1]] 
            if (length(tmp) > 0) 
            { 
              soln <- max(tmp) 
            } 
            else 
            { 
              soln <- NA 
            } 
          } 
          else 
          { 
            tmp <- tmp[tmp > bounds$x[2]] 
            if (length(tmp) > 0) 
            { 
              soln <- min(tmp) 
            } 
            else 
            { 
              soln <- NA 
            } 
          } 
        } 
        else 
        { 
          # keep the closest solution to any bound (or center of bounds, equivalently)
          d <- abs(tmp - mean(bounds$x)) 
          w <- which.min(d) 
          soln <- tmp[w] 
        } 
      } 
      
      soln  
    } # end of cubic.extrapolation 
    

    leftSide.fastScan <- function(bounds, o, A, target, lowerside.limit, h2modeless, ls.start=list(x=numeric(0), h=numeric(0), hp=numeric(0)), epsilon=1e-8, lower.side=1) 
    {     
      # in input bounds, the slopes (hp values) must be negative & positive, respectively
      # (thus pointing towards a 'dip', or local minimum)
    
      dip.characteristics <- function(b, o, A, epsilon=1e-8) 
      { 
        continue <- TRUE

        while (continue) 
        { 
          # intersection between the tangents at both ends of bounds (b)
          
          b.diff <- diff(b$hp) 
          
          if (b.diff == 0)  
          { 
            x <- mean(b$x) 
          } 
          else 
          { 
            x <- (diff(b$hp*b$x) - diff(b$h)) / b.diff  
            if (x <= b$x[1] | x >= b$x[2]) x <- mean(b$x) 
          } 
          
          h <- o$h(x, A) 
          hp <- o$hp(x, A) 

          side <- ifelse(hp < 0, 1, 2) 
          b$x[side] <- x 
          b$h[side] <- h 
          b$hp[side] <- hp 
          continue <- diff(b$x) > epsilon 
        } 

        list(x=x, h=h, hp=hp) 
      } # end of dip.characteristics 


      leftbump.characteristics <- function(o, start, hp.start, target, A, lowerside.limit, h2modeless, ls.start=list(x=numeric(0), h=numeric(0), hp=numeric(0)), epsilon=1e-8) 
      { 
        out <- list(x=numeric(0), h.max=numeric(0)) 

        # Find maximum point

        lim <- c(lowerside.limit, start) 
        if (length(ls.start$x) == 1) lim[1] <- ls.start$x 
        
        x <- start 
        hp <- hp.start 
        continue <- TRUE

        while (continue) 
        { 
          if (hp < 0) lim[2] <- x 
          else lim[1] <- x 

          hs <- -abs(o$hs(x, A)) 
          change <- hp/hs 
          x <- x - change 

          if (x > lim[2])      x <- (x + change + lim[2]) / 2 
          else if (x < lim[1]) x <- (x + change + lim[1]) / 2 

          hp <- o$hp(x, A) 
          converged <- abs(hp) < epsilon 
          continue <- !converged 
          
          if (continue) 
          { 
            h <- o$h(x, A) 
            continue <- h < h2modeless 
          } 
        } 


        if (converged) 
        { 
          out$h.max <- o$h(x, A) 
          continue <- out$h.max > target 
        } 


        if (continue) 
        { 
          # Find bump-left-side solution

          lim <- c(lowerside.limit, x) 
          
          if (length(ls.start$x) == 1) 
          { 
            x <- ls.start$x 
            h <- ls.start$h 
          } 
          else 
          { 
            hs <- o$hs(x, A) 
            change <- (hp-1)/hs 
            x <- x - change 
            if (x < lim[1]) x <- (x + change + lim[1]) / 2 

            h <- o$h(x, A) 
          } 

          while (continue) 
          { 
            hp <- o$hp(x, A) 
            change <- (h-target)/hp 
            x <- x - change 

            if (x > lim[2])      x <- (x + change + lim[2]) / 2 
            else if (x < lim[1]) x <- (x + change + lim[1]) / 2 

            h <- o$h(x, A) 
            continue <- abs(h-target) > epsilon 
          } 

          out$x <- x 
          out$h <- h 
        } 

        out$converged <- converged 
        out  
      } # end of leftbump.characteristics 
      

      higher.side <- 3 - lower.side 
      dip <- dip.characteristics(bounds, o, A, epsilon=epsilon) 
      leftbump <- leftbump.characteristics(o, bounds$x[lower.side], bounds$hp[lower.side], target, A, lowerside.limit, h2modeless, ls.start=ls.start)  
      
      
      if (leftbump$converged) 
      { 
        found.soln <- leftbump$h.max > target  
        continue <- !found.soln 
              
        if (found.soln) 
        { 
          x <- leftbump$x 
          h <- leftbump$h 
          hp <- NA 

        } 
        else 
        { 
          x <- dip$x 
          h <- dip$h 
          hp <- (bounds$x[higher.side] - h) / (bounds$x[higher.side] - x)  
        } 
      } 
      else 
      { 
        x <- lowerside.limit 
        h <- target 
        hp <- NA 
        continue <- FALSE
      } 
      

      list(x=x, h=h, hp=hp, continue=continue) 
    } # end of leftSide.fastScan 
  
  
    new.bound <- function(start, bounds, monotonic, o, A, target, range, above=wb$above) 
    { 
      x <- start 
      h <- o$h(x, A) 
      hp <- o$hp(x, A) 
      accepted <- is.finite(h) & is.finite(hp) 

      
      if (!accepted) 
      { 
        # Try a cubic extrapolation if we have two points in bounds
        
        n.bounds <- sum(!is.na(bounds$x)) 
        
        if (n.bounds == 2) 
        { 
          tmp <- cubic.extrapolation(target, bounds, monotonic, range) 
          
          if (!is.na(tmp)) 
          { 
            x <- tmp 
            h <- o$h(x, A) 
            hp <- o$hp(x, A) 
            accepted <- is.finite(h) & is.finite(hp) 
          } 
        } 
      } 
      
      
      if (!accepted) 
      { 
        if (above) 
        { 
          bound <- max(bounds$x, na.rm=TRUE)
        } 
        else 
        { 
          bound <- min(bounds$x, na.rm=TRUE)
        } 
      } 
      
      
      while (!accepted) 
      { 
        x <- (x + bound) / 2 
        h <- o$h(x, A) 
        hp <- o$hp(x, A) 
        accepted <- is.finite(h) & is.finite(hp) 
      } 

      list(x=x, h=h, hp=hp) 
    } # end of new.bound 


    unexplored.range <- function(visited.x, epsilon) 
    { 
      L <- length(visited.x) 
      j.odd  <- seq(from=1, to=L, by=2) 
      j.even <- seq(from=2, to=L, by=2) 

      x.odd  <- visited.x[j.odd] 
      x.even <- visited.x[j.even] 

      outer.range <- max(visited.x) - min(visited.x) 

      mean.x.odd  <- mean(x.odd) 
      mean.x.even <- mean(x.even) 

      if (mean.x.odd < mean.x.even) 
      { 
        inner.lim <- c(max(x.odd), min(x.even)) 
      } 
      else 
      { 
        inner.lim <- c(max(x.even), min(x.odd)) 
      } 

      inner.range <- diff(inner.lim) 

      continue <- inner.range > epsilon & inner.range/outer.range > 0.5 
      x <- ifelse(continue, mean(inner.lim), NA) 

      list(continue=continue, x=x) 
    } # end of unexplored.range 
    
    
    unreached.target <- function(target, h, higher.side, monotonic, mode.search, left.side=FALSE)
    {       
      if (!monotonic | mode.search) 
      { 
        unreached <- FALSE
      } 
      else 
      { 
        pos.slope <- higher.side == 2 
        unreached <- xor(pos.slope, h > target) 
        if (left.side) unreached <- !unreached 
      } 
      
      unreached 
    } # end of unreached.target 


    within.bounds <- function(x, bounds) 
    {       
      nbounds <- sum(!is.na(bounds$x)) 
      
      if (nbounds == 2) 
      {  
        above <- x > bounds$x[2] 
        below <- x < bounds$x[1] 
      } 
      else 
      { 
        w <- which(!is.na(bounds$x)) 
        bound <- bounds$x[w] 
        above <- x > bound 
        below <- !above 
      } 
      
      within <- !above & !below 
      out <- list(within=within, below=below, above=above) 

      out # list with 3 dimensions: within, above & below 
    } # end of within.bounds 
        
    
    # --- Fct starts (secured.NR.search) -----------------------------------------------------------
    
    mode.search <- length(max.point$x) == 0 
    monotonic <- !mode.search  
    
    limit <- list(x=range, checked=rep(FALSE,2), usable=rep(FALSE,2), h=rep(NA,2), hp=rep(NA,2))
    limit$checked[1] <- inestimable.lower.limit
    
    bounds <- list(x=rep(NA, 2), h=rep(NA, 2), hp=rep(NA, 2), include.soln=FALSE)

    
    if (length(max.point$x) > 0) 
    { 
      higher.side <- ifelse(start < max.point$x, 2, 1) 
      
      bounds$x[higher.side]  <- max.point$x 
      bounds$h[higher.side]  <- max.point$h 
      bounds$hp[higher.side] <- 0 
      found2bounds <- TRUE
      
      expected.hpsign <- 2*higher.side - 3 
      correct4hp.sign <- start > max.point$x 
      
      h2modeless <- 2*max.point$h - target 
    } 
    else 
    { 
      found2bounds <- FALSE
      higher.side <- 1 
      correct4hp.sign <- TRUE
    } 
    
    
    lower.side <- 3 - higher.side 
    x <- start 

    # make sure that x is not out of range
    if (x < limit$x[1]) x <- (max.point$x + limit$x[1])/2  

        
    bounds$range <- range 
    bounds$higher.side <- higher.side 


    # Register initial x into bounds

    h <- o$h(x, A) 
    hp <- o$hp(x, A) 

    
    if (length(max.point$x) > 0) 
    { 
      bounds$include.soln <- h < target 
      side <- lower.side 
    } 
    else 
    { 
      side <- 1 
    } 
    
    
    bounds$x[side] <- x 
    bounds$h[side] <- h 
    bounds$hp[side] <- hp 
    visited.x <- x 
    

    # Do a first step

    if (correct4hp.sign) hp <- expected.hpsign * abs(hp) 
    change <- (h-target)/hp                               
    x <- x - change                                     


    # continue until convergence

    unreachable.mode <- FALSE
    continue <- TRUE
    count <- 0 

    while (continue) 
    { 
      computed.h <- FALSE
      accepted <- FALSE
      count <- count + 1 
      
      while (!accepted) 
      { 
        wb <- within.bounds(x, bounds) 
        
        if (bounds$include.soln) 
        { 
          if (wb$within) 
          { 
            accepted <- TRUE
          } 
          else 
          { 
            x <- cubic.extrapolation(target, bounds, monotonic, range) 
            accepted <- !is.na(x) 
          } 
        } 
        else if (wb$above) 
        { 
          # bounds do not include solution and we are looking above bounds
          
          if (x > limit$x[2]) 
          { 
            if (limit$checked[2]) 
            { 
              accepted <- TRUE
              h <- limit$h[2] 
              
              
              if (unreached.target(target, h, higher.side, monotonic, mode.search)) 
              { 
                # force end of while-loop, as target is not reached even at this end
                x <- limit$x[2] 
                h <- target 
                computed.h <- TRUE
              } 
              else if (limit$usable[2]) 
              { 
                x <- bounds$range[2] 
                hp <- limit$hp[2] 
                computed.h <- TRUE
              } 
              else 
              { 
                x <- (x + change + bounds$range[2])/2 
              } 
            } 
            else 
            { 
              limit$checked[2] <- TRUE
              h <- o$h(limit$x[2], A) 

              if (unreached.target(target, h, higher.side, monotonic, mode.search)) 
              { 
                # force end of while-loop, as target is not reached even at this end
                x <- limit$x[2] 
                h <- target 
                accepted <- TRUE
                computed.h <- TRUE
              } 
              else 
              { 
                limit$h[2] <- h 
              
                if (is.finite(h)) 
                { 
                  hp <- o$hp(limit$x[2], A) 
                  limit$hp[2] <- hp 
                  limit$usable[2] <- is.finite(hp) 
                  accepted <- limit$usable[2] 
                } 
                else 
                { 
                  accepted <- FALSE
                } 
              
                if (accepted) 
                { 
                  x <- bounds$range[2] 
                  computed.h <- TRUE
                  unreachable.mode <- mode.search & h > 0
                } 
                else 
                { 
                  x <- (x + change + bounds$range[2])/2 
                } 
              } 
            } 
          } 
          else 
          { 
            tmp <- new.bound(x, bounds, monotonic, o, A, target, range) 
          
            x <- tmp$x 
            h <- tmp$h 
            hp <- tmp$hp 
            computed.h <- TRUE
            accepted <- TRUE
          } 
        } 
        else 
        { 
          # bounds do not include solution and we are looking below bounds
          
          if (x <= limit$x[1]) 
          { 
            # we are looking out of the variable domain (bad!)
            
            if (!limit$checked[1]) 
            { 
              limit$checked[1] <- TRUE
              
              h <- o$h(limit$x[1], A) 
              limit$h[1] <- h 
              
              if (is.finite(h)) 
              { 
                if (unreached.target(target, h, higher.side, monotonic, mode.search, left.side=TRUE))
                { 
                  # force end of while-loop, as target is not reached even at this end
                  x <- limit$x[1] 
                  h <- target 
                  accepted <- TRUE
                } 
                else 
                { 
                  hp <- o$hp(limit$x[1], A)  
                  limit$hp[1] <- hp 
                  limit$usable[1] <- is.finite(hp) 
                  accepted <- limit$usable[1] 
                } 
                
                computed.h <- accepted 
                if (accepted) x <- limit$x[1] 
              } 
              else 
              { 
                limit$usable[1] <- FALSE
                x <- (min(bounds$x, na.rm=TRUE) + limit$x[1]) / 2
              } 
            } 
            else if (limit$usable[1]) 
            { 
              x <- limit$x[1] 
              h <- limit$h[1] 
              hp <- limit$hp[1] 
              computed.h <- TRUE
              accepted <- TRUE
            } 
            else 
            { 
              # not limit$usable[1]
              x <- (min(bounds$x, na.rm=TRUE) + limit$x[1]) / 2
            } 
          } 
          else 
          { 
            tmp <- new.bound(x, bounds, monotonic, o, A, target, range) 
          
            x <- tmp$x 
            h <- tmp$h 
            hp <- tmp$hp 
            computed.h <- TRUE
            accepted <- TRUE
            
            unreachable.mode <- mode.search & ((x == limit$x[2] & h > 0) | (x == limit$x[1] & h < 0)) 
          } 
        } 
      } 

      
      if (!computed.h) 
      { 
        h <- o$h(x, A) 
        hp <- o$hp(x, A) 
      } 
        
      #converged <- abs(h - target) < epsilon | unreachable.mode
      converged <- abs(h - target) < epsilon | unreachable.mode | (!mode.search & x == limit$x[side] & h > target)

      if (!converged & bounds$include.soln) 
      { 
        converged <- diff(bounds$x) < epsilon 
        if (converged) x <- mean(bounds$x) 
      } 

      continue <- !converged 


      if (continue) 
      {         
        # register results in bounds
        
        h.lt.target <- h < target 

        if (bounds$include.soln) 
        { 
          side <- ifelse(h.lt.target, lower.side, higher.side) 
        } 
        else 
        { 
          bounds$include.soln <- xor(h.lt.target, bounds$h[1] < target)  
        
          if (found2bounds) 
          { 
            if (mode.search) 
            { 
              side <- ifelse(h.lt.target, lower.side, higher.side) 
              change.opposite.side.bounds <- TRUE
            } 
            else if (bounds$include.soln) 
            { 
              # bounds newly include solution (happening for the first time with current x)
              wb <- within.bounds(x, bounds) 
              
              if (wb$within) 
              { 
                if (higher.side == 2) 
                { 
                  side <- lower.side 
                  change.opposite.side.bounds <- FALSE

                  
                  if (hp < 0) 
                  { 
                    bounds$x[lower.side] <- x 
                    bounds$h[lower.side] <- h 
                    bounds$hp[lower.side] <- hp 
                    
                    tmp <- leftSide.fastScan(bounds, o, A, target, range[lower.side], h2modeless) 
                    
                    x <- tmp$x 
                    h <- tmp$h 
                    hp <- tmp$hp 
                    continue <- tmp$continue    
                  } 
                } 
                else 
                { 
                  stop("Scenario imprevu /1a --- ICI\n", save.objects, '\n') 
                } 
              } 
              else if (wb$below) 
              { 
                side <- 1 
                change.opposite.side.bounds <- TRUE
              } 
              else 
              { 
                # wb$above
                side <- 2 
                change.opposite.side.bounds <- TRUE
              } 
            } 
            else 
            { 
              # bounds still do not include solution
              
              wb <- within.bounds(x, bounds) 
              
              if (wb$within) 
              { 
                change.opposite.side.bounds <- FALSE
                
                if (h > max(bounds$h)) 
                { 
                  side <- higher.side 
                } 
                else if (h < min(bounds$h)) 
                { 
                  side <- ifelse(hp < 0, lower.side, higher.side) 
                } 
                else 
                { 
                  side <- higher.side 
                } 
              } 
              else if (wb$below) 
              { 
                if (higher.side == 2) 
                { 
                  if (hp < 0) 
                  { 
                    bounds$x[lower.side] <- x 
                    bounds$h[lower.side] <- h 
                    bounds$hp[lower.side] <- hp         
                    
                    tmp <- leftSide.fastScan(bounds, o, A, target, range[lower.side], h2modeless) 
                    
                    x <- tmp$x 
                    h <- tmp$h 
                    hp <- tmp$hp 
                    continue <- tmp$continue 
                    change.opposite.side.bounds <- TRUE
                  } 
                  else if (h < min(bounds$h)) 
                  { 
                    side <- lower.side 
                    change.opposite.side.bounds <- TRUE
                  } 
                  else 
                  { 
                    leftside.potential.start.point <- list(x=x, h=h, hp=hp) 
                  
                    lim <- c(x, bounds$x[lower.side]) 
                    x <- mean(lim) 
                    hp <- o$hp(x, A) 
                    
                    while (hp > 0) 
                    { 
                      h <- o$h(x, A) 
                      side <- ifelse(h > bounds$h[lower.side], 1, 2) 
                      lim[side] <- x 
                      x <- mean(lim) 
                      hp <- o$hp(x, A) 
                    } 
                    
                    h <- o$h(x, A) 
                    
                    bounds$x[higher.side] <- bounds$x[lower.side] 
                    bounds$h[higher.side] <- bounds$h[lower.side] 
                    bounds$hp[higher.side] <- bounds$hp[lower.side] 
                    
                    bounds$x[lower.side] <- x 
                    bounds$h[lower.side] <- h 
                    bounds$hp[lower.side] <- hp 
                  
                    tmp <- leftSide.fastScan(bounds, o, A, target, range[lower.side], h2modeless, ls.start=leftside.potential.start.point) 
                    
                    x <- tmp$x 
                    h <- tmp$h 
                    hp <- tmp$hp 
                    continue <- tmp$continue    
                  
                    side <- higher.side 
                    change.opposite.side.bounds <- FALSE
                  } 
                } 
                else 
                { 
                  stop("Scenario imprevu /1b --- ICI\n", save.objects, '\n') 
                } 
              } 
              else 
              { 
                # wb$above
                side <- lower.side 
                change.opposite.side.bounds <- TRUE
              }               
            } 
            
            
            if (change.opposite.side.bounds) 
            { 
              opposite.side <- 3 - side 

              bounds$x[opposite.side]  <- bounds$x[side] 
              bounds$h[opposite.side]  <- bounds$h[side] 
              bounds$hp[opposite.side] <- bounds$hp[side] 
            } 
          } 
          else 
          { 
            # !found2bounds
            found2bounds <- TRUE
            
            side <- ifelse(h < bounds$h[1], lower.side, higher.side) 
            
            if (side == 1) 
            { 
              bounds$x[2]  <- bounds$x[1] 
              bounds$h[2]  <- bounds$h[1] 
              bounds$hp[2] <- bounds$hp[1] 
            } 
          } 
        } 
        

        bounds$x[side]  <- x 
        bounds$h[side]  <- h 
        bounds$hp[side] <- hp 
        visited.x <- c(visited.x, x) 
        
        # converged <- abs(h - target) < epsilon | unreachable.mode
        converged <- abs(h - target) < epsilon | unreachable.mode | (!mode.search & x == limit$x[side] & h > target)

        
        if (!converged) 
        { 
          continue <- count < max.niter
          
          if (!continue) 
          { 
            # Before we give up, we check one last thing:
            # if the last few steps were all leaning in the same direction,
            # then convergence may only be a question of time & patience!
            # Give it (yet) another chance!
            
            tmp <- sign(diff(visited.x)) 
            last.change.direction <- tmp[length(tmp)] 
            iter.in.opposite.direction <- which(tmp==-last.change.direction) 
            
            if (length(iter.in.opposite.direction) == 0) 
            { 
              continue <- TRUE
            } 
            else 
            { 
              last.iter.in.opposite.direction <- max(iter.in.opposite.direction) 
              if ((count - last.iter.in.opposite.direction) > 20) continue <- TRUE
            } 


            # Alternatively, if range was unsatisfactorily explored, take a point
            # in the middle of unexplored range and redo

            if (!continue) 
            { 
              tmp <- unexplored.range(visited.x, epsilon) 
              continue <- tmp$continue 
              if (continue) x <- tmp$x 
            } 

            
            if (continue) 
            { 
              count <- 0 
              visited.x <- numeric(0) 
            } 
          } 
        } 
               
        
        if (continue) 
        { 
          found2bounds <- TRUE
          abs.hp <- abs(hp) 
          if (correct4hp.sign) hp <- expected.hpsign * abs.hp 
          
          change <- (h - target)/hp 
          pct.change <- abs(change)/diff(bounds$x) 
          
          if (pct.change < 1e-4 & abs.hp > 1000) 
          { 
            if (change < 0) p <- c(0.9, 0.1) 
            else p <- c(0.1, 0.9)  
            x <- sum(p * bounds$x) 
          } 
          else x <- x - change             
        } 
        else if (!converged & diff(bounds$x) < epsilon) 
        {
          if (bounds$include.soln) 
          { 
            converged <- TRUE
            continue <- FALSE
            x <- mean(bounds$x) 
          } 
          else 
          { 
            change <- - (max(bounds$h)-target)/mean(bounds$hp) 
            x <- mean(bounds$x) + change + epsilon*sign(change) 
            continue <- TRUE
            count <- 0 
          } 
        } 


        if (cat.x) 
        { 
          cat('count=', count, '\n') 
          cat('x = ', x, '\n') 
          cat('h = ', h, '\n') 
          cat('target = ', target, '\n') 
          cat('abs(h-target)', abs(h-target), '\n') 
          cat('bounds$x = ', bounds$x, '\n') 
          cat('diff(bounds$x) = ', diff(bounds$x), '\n') 
          if (diff(bounds$x) < 0) cat('***************\n') 
        } 
      } 
    } 
      
    list(x=x, converged=converged, bounds=bounds) 
  } # end of secured.NR.search 


  smoothed.area.estimate <- function(area, x, A, math.lower.limit, ref.range, inestimable.lower.limit, d=min(c(1e-4, diff(ref.range)/1000)), max.count=10000, hp.mult=1, remote.right=FALSE)
  {  
    f.x <- A$f(x, A) 
     
    a <- rep(NA, 3) 
    # values around x
    x.star <- x + c(-1, 0, 1)*d 
    
    
    if (x.star[1] < math.lower.limit)  
    { 
      if (inestimable.lower.limit) 
      { 
        x.star[1] <- (x + math.lower.limit)/2 
      } 
      else 
      { 
        x.star[1] <- math.lower.limit 
      } 
    } 
    
    
    for (i in 1:3) a[i] <- area(x.star[i], A) 
    a.x <- a[2] 
    
    observed.slopes <- hp.mult*diff(a)/diff(x.star) 
    expected.slope <- f.x 
    
    observed.slope.angles <- atan(observed.slopes) 
    expected.slope.angle  <- atan(expected.slope) 
    angle.diff <- max(abs(expected.slope.angle-observed.slope.angles)) 

    estimate.ok <- angle.diff < 0.035 & all(observed.slopes > 0) 
    # 0.035 radians is approximately 2 degrees
    
    continue <- !estimate.ok 
    
    
    if (continue) 
    { 
      # We look for a point on the left side of x and one on its right side
      # where the integral [area] seems to be correctly estimated
      # (with positive slopes of value close to expected value [density])
      
      # If we are estimating the area on the remote right end of the distrn (remote.right=T)
      # then we will be happy if an estimate to the left of x is obtained 
      
      found <- list(ok=rep(FALSE,2), x=rep(NA,2), a=rep(NA,2), f=rep(NA,2))
      side <- 0 
      
        # Store starting points
        tmp.x <- x.star 
        tmp.a <- a 
        
      
      while (continue) 
      { 
        side <- side + 1 # side = 1 -> Left  of x 
                         # side = 2 -> Right of x
        
        direction <- ifelse(side==1, -1, 1) 
        
        
        if (side == 2) 
        { 
          x.star <- tmp.x 
          a <- tmp.a 
        } 
        
        j <- ifelse(side==1, 1, 3) 
        z <- x.star[j] 
        
        count <- 0 
        on.border <- FALSE
        
        while (continue) 
        { 
          z <- z + direction * d 
        
          # make sure that we don't step to the left of mathematical lower limit
          if (z <= math.lower.limit) 
          { 
            if (inestimable.lower.limit) 
            { 
              z <- (math.lower.limit + z + d) / 2 
            } 
            else 
            { 
              z <- math.lower.limit 
              on.border <- TRUE
            } 
          }  
          
          
          if (side == 1 & z == math.lower.limit) 
          { 
            # This is still possible even with the above 'protection' against
            # this scenario, due to numerical imprecision
            a.z <- 0 
            continue <- FALSE
            estimate.ok <- TRUE
          } 
          else 
          { 
            a.z <- area(z, A) 
          
          
            if (side == 1) 
            { 
              x.star <- c(z, x.star[-3]) 
              a <- c(a.z, a[-3]) 
            } 
            else 
            { 
              x.star <- c(x.star[-1], z) 
              a <- c(a[-1], a.z) 
            }  
        
    
            observed.slopes <- hp.mult*diff(a)/diff(x.star) 
            expected.slope <- A$f(z, A) 
        
            observed.slope.angles <- atan(observed.slopes) 
            expected.slope.angle  <- atan(expected.slope) 
            angle.diff <- max(abs(expected.slope.angle-observed.slope.angles)) 

            estimate.ok <- (angle.diff < 0.035 & all(observed.slopes > 0)) | angle.diff < 1e-6 
            continue <- !estimate.ok & !on.border  
          } 

          
          if (continue) 
          { 
            if (remote.right & side == 1) 
            { 
              continue <- expected.slope/f.x < 10000 
            } 
            else 
            { 
              count <- count + 1 
              continue <- count <= max.count 
            } 
          } 
        } 
        
        
        if (side == 1 & remote.right) 
        { 
          continue <- !estimate.ok 
          if (estimate.ok) a.x <- a.z 
        } 
        else continue <- side == 1 
        

        # Store results
        found$ok[side] <- estimate.ok 
        found$x[side]  <- z 
        found$a[side]  <- a.z 
        found$f[side]  <- expected.slope         

      } # end of while-continue 
      
      
      if (!remote.right | !found$ok[1]) 
      { 
        if (found$ok[1] & found$ok[2]) 
        { 
          f.diff <- diff(found$f) 
          extrapolate <- f.diff == 0 
        
          if (!extrapolate) 
          { 
            x.intersect <- - (diff(found$a) - diff(found$f*found$x)) / f.diff 
            extrapolate <- x.intersect < found$x[1] | x.intersect > found$x[2] 
          
            if (!extrapolate) 
            { 
              j <- ifelse(x <= x.intersect, 1, 2) 
              a.x <- found$a[j] + found$f[j]*(x-found$x[j]) 
            } 
          } 
        
        
          if (extrapolate) a.x <- found$a[1] + diff(found$a)/diff(found$x) * (x - found$x[1]) 
        } 
        else if (found$ok[1] | found$ok[2]) 
        { 
          j <- which(found$ok) 
          a.x <- found$a[j] + found$f[j]*(x-found$x[j]) 
        } 
        else 
        { 
          stop('Could not get a smoothed cumulative density estimate.\n', save.objects, '\n') 
        } 
      } 
    } 


    a.x 
  } # end of smoothed.area.estimate 
       

  ref <- ref.points(o, A, start, range, inestimable.lower.limit=inestimable.lower.limit, epsilon=precision) 
    bimodal <- ref$bimodal


  if (ref$mode$found & !bimodal) 
  { 
    ref.alt.prec <- 1e-6 

    ref.alt.range.l <- c(ref$remote[1], ref$mode$x) 
    ref.alt.range.r <- c(ref$mode$x, ref$remote[2]) 

    # right-side 

    ref.alt <- ref.points(o, A, start=ref.alt.range.r[2], range=ref.alt.range.r, epsilon=precision, find.remote=FALSE)
   
    if (abs(ref.alt$mode$x - ref$mode$x) > ref.alt.prec)  
    {  
      # If it is a new global mode, run ref.points again to update the value of the remote limits

      if (ref.alt$mode$h > ref$mode$h) ref <- ref.points(o, A, start=ref.alt$mode$x, range, inestimable.lower.limit=inestimable.lower.limit, epsilon=precision) 
      bimodal <- TRUE
    } 
    

    # ICI js ai ajoute diff(cond ref.alt.range.l) > 0
    if (!bimodal && diff(ref.alt.range.l) > 0)
    { 
      # left-side

      # ICI js
      ref.alt.start <- ref.alt.range.l[1]
      ref.alt.start.slope <- o$log.f.prime(ref.alt.start, A)

      if (ref.alt.start.slope == 0)
      {
        midpoint <- mean(ref.alt.range.l) # 1st possibility

        my.h <- o$log.f(ref.alt.start, A)
        bottom2top.straightline.slope <- (ref$mode$h - my.h) / diff(ref.alt.range.l)
        dx <- bottom2top.straightline.slope / o$log.f.second(ref.alt.start, A)
        tmp <- ref.alt.start + dx
        if (tmp < ref.alt.start) ref.alt.start <- midpoint
        else                     ref.alt.start <- tmp
      }


      ref.alt <- ref.points(o, A, start=ref.alt.start, range=ref.alt.range.l, epsilon=precision, find.remote=FALSE)
   
      if (abs(ref.alt$mode$x - ref$mode$x) > ref.alt.prec)  
      {    
        # If it is a new global mode, run ref.points again to update the value of the remote limits 
        
        if (ref.alt$mode$h > ref$mode$h) ref <- ref.points(o, A, start=ref.alt$mode$x, range, inestimable.lower.limit=inestimable.lower.limit, epsilon=precision) 
        bimodal <- TRUE
      } 
    } 
  } 


  if (ref$mode$found) 
  { 
    math.lower.limit <- o$math.lower.limit 

    A$f <- o$f 
    A$lower.limit <- ref$remote[1] 
    A$M <- ref$mode$h # (kind of a) standardizing constant 
    area <- function(x, A){integrate(A$f, A$lower.limit, x, A=A)$value} 
    
    area.tot <- smoothed.area.estimate(area, ref$remote[2], A, math.lower.limit, ref$remote, inestimable.lower.limit, remote.right=TRUE)
    target <- u*area.tot 
    precision <- precision*area.tot 

    x <- ref$mode$x 
    area.x <- 0

    # ICI js ai ajoute condition et valeur alternative a area.mode
    if (x != ref$remote[1]) area.mode <- smoothed.area.estimate(area, x, A, math.lower.limit, ref$remote, inestimable.lower.limit)
    else area.mode <- 0

    u.mode <- area.mode/area.tot  
  
    converged <- abs(area.mode - target) < precision 
    continue <- !converged 
  } 
  else 
  { 
    stop('mode not found\n') 
  } 
  
  
  if (continue) 
  { 
    distrn.leftSide <- u < u.mode 
  
    if (distrn.leftSide) 
    { 
      area <- function(x, A){integrate(A$f, x, A$upper.limit, A=A)$value} 
      A$upper.limit <- ref$mode$x 
      range <- c(ref$remote[1], ref$mode$x) 
    } 
    else 
    { 
      A$lower.limit <- ref$mode$x 
      range <- c(ref$mode$x, ref$remote[2]) 
    } 
  } 
  
  

  if (continue & bimodal) 
  { 
    # Slower but safer algorithm for posterior distrns that are potentially bimodal
        
    f.mode <- A$f(ref$mode$x, A) 
    start <- ref$mode$x  + (target-area.mode)/f.mode 
    
    continue <- FALSE
    converged <- TRUE
    
    if (distrn.leftSide) 
    { 
      highest.point <- ref$remote[1] 
      h.max <- area.mode 
      target <- (u.mode - u)/u.mode * h.max 
    } 
    else 
    { 
      highest.point <- ref$remote[2] 
      h.max <- smoothed.area.estimate(area, ref$remote[2], A, math.lower.limit, range, inestimable.lower.limit, remote.right=TRUE)
      target <- (u-u.mode)/(1-u.mode) * h.max 
    } 
    
    
    o.cum <- list(h=area, hp=o$f) 
    
    nr <- secured.NR.search(o.cum, A, start, precision, target=target,  
                             range=range, max.point=list(x=highest.point, h=h.max), 
                             max.niter=NR.max.iter) 
                             
    if (nr$converged)  
    { 
      x <- nr$x 
    } 
    else 
    { 
      # it seems we were caught in an infinite loop, had we not limited the number of iterations;
      # see if we can rule that problem out
      
      x <- nr$bounds$x[1] 
      h <- nr$bounds$h[1] 
      
      x <- ping.pong.tie.breaker(area, A, x, h, target, math.lower.limit, range, inestimable.lower.limit, precision, distrn.leftSide=distrn.leftSide) 
    } 
  } 
    
    

  if (continue) 
  { 
    count <- 0 
    direction <- ifelse(distrn.leftSide, -1, 1)     

    if (distrn.leftSide) 
    { 
      target <- (u.mode - u)/u.mode * area.mode   
    } 
    else 
    { 
      area.remote.right <- smoothed.area.estimate(area, ref$remote[2], A, math.lower.limit, range, inestimable.lower.limit, remote.right=TRUE)
      target <- (u-u.mode)/(1-u.mode) * area.remote.right 
    } 
    
    
    # New starting point (based on a normal approximation)
    z.max <- sqrt(-2*ref$log.f.ratio.remote) 
    z <- qnorm(u) 
    if (abs(z) > z.max) z <- sign(z) * z.max * 0.99 
    if (distrn.leftSide) 
    { 
      z <- min(z, 0) 
      side <- 1 
    } 
    else 
    { 
      z <- max(z, 0) 
      side <- 2 
    } 
    side.len <- abs(ref$mode$x - ref$remote[side]) 
    side.sigma <- side.len / z.max 
    x <- ref$mode$x + z * side.sigma 
    area.x <- area(x, A)     
  } 
  
  

  
  if (continue) 
  { 
    # The second while-block below existed before, but is now embedded in this if-block, and
    # preceded by an almost identical while-block, which is more cautious since we do not start
    # at mode anymore, which ensured convergence to solution without further precautions
    # (on possible stepping out of bounds)

    continue <- area.x > target
    
    while (continue) 
    { 
      hp <- o$f(x, A) 
      change <- (target - area.x)/hp*direction 
      x <- x + change 
        
      if (direction > 0) 
      { 
        if (x < A$lower.limit) x <- (x - change + A$lower.limit) / 2 
      } 
      else if (x > A$upper.limit) x <- (x - change + A$upper.limit) / 2 

    
      area.x <- area(x, A) 
      count <- count + 1 
      converged <- abs(area.x - target) < precision 
      continue <- !converged  & count <= NR.max.iter & area.x > target 
    } 
    
    continue <- !converged & count <= NR.max.iter
    
    # same as above, but without the direction-checks, which slow down the process
    while (continue) 
    { 
      hp <- o$f(x, A) 
      change <- (target - area.x)/hp*direction 
      x <- x + change 
            
      area.x <- area(x, A) 
      count <- count + 1 
      converged <- abs(area.x - target) < precision 
      continue <- !converged  & count <= NR.max.iter 
    } 
  } 
  
  
  if (!converged) 
  { 
    if (x > ref$remote[2] | x < ref$remote[1]) 
    { 
      stop("stepped out of bounds -- due to numerical imprecision?", save.objects, '\n') 
    } 
    else 
    { 
      # it seems we were caught in an infinite loop, had we not limited the number of iterations;
      # see if we can rule that problem out

      x <- ping.pong.tie.breaker(area, A, x, area.x, target, math.lower.limit, range, inestimable.lower.limit, precision, distrn.leftSide=distrn.leftSide) 
    } 
  } 

   
  if (length(save.objects) > 0 & unlink) unlink(save.objects) 


  if (return.remote) list(x=x, mode=ref$mode$x, remote=ref$remote) 
  else x 
} # end of dens.gen.icdf 
