                                                                                 
                                                                                         
#  Version 0.25 (Feb 2022)
                                                                                                                              
                                                                    
# Change Log *******************************************************************          
#
#
# Version 0.25 (Feb 2022)
# -----------------------
#   We have included the area of each riskband in the output, under the name geo$area
#
#
# Version 0.24 (Feb 2022)
# -----------------------
#   Added the function SEG.riskband.posterior.plot to show a few informative
#   graphics - scatter plot of (mu, sigma), and histograms of P95, mu & sigma
# 
#
# Version 0.23 (Feb 2022)                      
# -----------------------   
#   The function smooth() was revisited --- the function go.east.young.man was dropped
#   as its code was simplified and incorporated in smooth(): the calculations of CPL
#   now start at h02 (the height at entry in riskband #2) rather than at h01 (the height at entry in riskband #1).
#
#   The fct last.riskband was modified to also treat first riskband
#   -> hence it was renamed closing.riskband 
#   For the same reason, revisit.last2riskbands was renamed revisit.closing2riskbands
#
#
# Version 0.22 (Feb 2022)                      
# -----------------------
#   We found that the option (c) in version 0.21b was the best and dropped the other two.
#   We have also added the fct revisit.last2riskbands at the very end of CPL computation,
#   making the next-to-last riskband convex when possible.
#   We have also moved the functions rb.init, rb.update & better.dens to a separate file location
#   as they are generic functions which we could easily be used in other projects, and renamed
#   them bs.init, bs.update & better.result
#
#                  
# Version 0.21b (Jan 2022)
# ------------------------
#   Back to version 0.21, with a slight addition for model comparison,
#   giving preseance to endregions with either
#     a) min difference between the slopes of their two segments.
#     b) min angle between these two segments
#     c) the lower angle form by these two segments and the last segment of the preceding region 
#        and the 1st segment of the next riskband.
#   We have also added a last 'brush' to each inner riskband (independently), relaxing the breakpoint location of the two lines in it
#   to get a distribution as smooth as possible at a low runtime cost (revisit.inner.riskbands);
#   finally, once we get out of the loops looking for an optimal distribution curve, we try again to fond better breakpoint locations
#   for the inner three regions, this time not independently, but taking into account the new curve defined by each combination of
#   the revisited breakpoints: this is done only with the final distrn, as this is (slightly) more computer intensive 
#   and would hence be too long to do at each step in the search for an optimal distribution curve.  
#
#
# Version 0.21a (Jan 2022)
# -----------------------
#   Modified the way two candidates for CPL distrns are compared
#    (the way they are 'measured' through () function)
#    -> the results are bad, many slopes are very sharp, which we try to avoid 
#    => We step back to version 0.21
#
#                         
# Version 0.21 (Jan 2022)
# -----------------------
#   Revisited the way we close with last riskband, where the line.break's
#     are free to occur anywhere in the region (and not only a x95=b[5])
#
#                                                      
# Version 0.20 (Jan 2022)
# -----------------------  
#   Revisited the way optimal CPL is assessed along the way
#     We have simplified the way s4 (the slope in riskband #4) is estimated
#     from h4 (the height of the distrn fct at the beginning of the same riskband)
#     by calculating the regression coefficients at each each iteration from cumulative moments
#     [rather than re-calculating the coefficient estimates from scratch at each iteration from 'history']
#
#
# Version 0.19 (Jan 2022)
# -----------------------
#   Added the function last.riskband, revisiting the way we 'close' CPL
#   The output of () is made multi-dimensional for a more flexible way of comparing CPL distrn fcts
#
#
# Version 0.18 (Jan 2022)
# -----------------------
#   Dropped the parm1.is.intercept approach 
#     -> parm1 is ALWAYS the intercept now
#   Added the use of 'history' to make convergence to CPL faster
#
#
# Version 0.17 (Jan 2022)
# -----------------------
#   Simplified how the objects were manipulated within smooth() function 
#     (in the computation of CPL prior distrn)
#
#   Two ideas were explored in the search for the path to CPL prior density estimation
#     - in v. 0.17a, a simple (programmatically) approach was chosen, but it obviously 
#       involved the repetition of the same paths along the way
#     - in v. 0.17b, the approach is a little bit more sophisticated but yet easy to understand
#       We think v. 0.17b, however it necessitated more work, will run faster (less repetition, but more "ifs")
#       -> And it did, indeed! We adopt that code and move on to version 0.18 
#           
#    
# Version 0.16 (Jan 2022)
# -----------------------
#   Added a function (SEG.riskband.plot.CPL.prior)
#   to draw the prior distrn used when the option use.continuous.piecewise.linear.prior=T
#   is used.
#   - The code within SEG.riskband was also slighlty modified in order for the output to include
#     all the information necessary to do that plot.
#
#
# Version 0.15 (Jan 2022)
# -----------------------
#   Now allows the possibility to use a Continuous Piecewise-Linear prior distrn
#     (through argument use.continuous.piecewise.linear.prior, defaulting to FALSE)
#
#
# Version 0.14 (Nov 2021)
# -----------------------
#   Now returns x95 sample
#     as well as its posterior median and 95%-UCL
#
#
# Version 0.13 (Mar 2021)
# -----------------------
#   When region.prior.prob is left undefined, two options are possible:
#     either a) regions are equally probable if equally.probable.regions=T (the default)
#         or b) a uniform prior is used on mu & sigma (if equally.probable.regions=F).
#               In earlier versions, the argument equally.probable.regions did not exist,
#               but the algorithm assumed a uniform prior (b) on mu & sigma 
#               when region.prior.prob was not defined.
# 
#   Important note: equally.probable.regions is ignored/irrelevant when region.prior.prob is specified.
#
#   Also:
#     - The default cut-off values were changed for the most-used values 0.01 / 0.1 / 0.5 / 1
#     - A correction was brought to the calculation of area with shape = 45 
#      (which will [slightly] affect prior the density calculated for each riskband)
#     - Correction: now use the new definition of region.prior.logdens rather than region.prior.logprob
#       in the calculation of segment posterior probabilities (that was incorrect!)
#
#
# Version 0.12
# ------------
#   Added argument save.RData
#     -> useful at development stage only, for monitoring parameters values when dens.gen.icdf is crashing
#
# Version 0.11
# ------------
#  - using out.logout.moments to compute 'moments'
#  - using default.inits
#  - also using argument logNormal.distrn in call to data.summary
#
# 
# Version 0.10 [replaces model.Banerjee]:
# ------------
#  - dropped gen.sum and gen.sum2
#  - now using:
#      empty.matrices
#      me.gen
#      me.gen.object
#      mu.truncatedData.gen
#      mu.truncatedData.gen.object
#      sigma.gen.object
#      sigma.truncatedData.gen
#      sigma.truncatedData.gen.object
#      sqrt.invertedGamma.gen
#      y.gen
#      y.gen.inits
#
#  - functions logp.from.logpcum and pnorm.logpcum were relegated to file fcts.R
#  - NOT using sigma.gen.icdf anymore
#  - cv.range and me.sd.range were added 
#     only one of me.sd.range or cv.range must be specified [of length 2, if any]
# - objects mu.sample & sd.sample were regrouped in list 'sample'
#     while burnin.mu & burnin.sd were regrouped in list 'burnin'
#     > both are wrapped through function out.sample (defined in fcts.R) at the end of function
               
                                       
SEG.riskband <- function(y=numeric(0), lt=numeric(0), gt=numeric(0),
  interval.lower=numeric(0), interval.upper=numeric(0),
  n.iter=15000, n.burnin=500, n.thin=1, n.chains=1, monitor.burnin=F,
  A=c(0.01, 0.1, 0.5, 1),
  region.prior.prob=rep(NA, R),
  gm.min=exp(mu.lower), gm.max=exp(mu.upper),
    mu.lower=-3, mu.upper=6.2,
  gsd.min=1.05, gsd.max=4, 
    sigma.lower=log(gsd.min), sigma.upper=log(gsd.max), 
  init.mu=rep(default.inits$mu, n.chains),
  init.sigma=rep(default.inits$sigma, n.chains),
  outcome.is.logNormally.distributed=TRUE,
  equally.probable.regions=TRUE,  # ignored when region.prior.prob is provided
  use.continuous.piecewise.linear.prior=F,
  quantile=0.95,
  me.sd.range=numeric(0), cv.range=numeric(0),
  save.RData=FALSE, RData.dir='c:/users/jerome')
{
  # IMPORTANT:
  # Selon la valeur de outcome.is.logNormally.distributed (TRUE/FALSE)
  # nous devrions demander à l'usager 
  #  i) seulement (gm.min, gm.max, gsd.min, gsd.max) si outcome.is.logNormally.distributed = TRUE
  # ii) seulement (mu.min, mu.max, sd.min, sd.max)   sinon


  # ___ Useful function(s) ____________________________________________________________________________
  
  
  # new_0.15
  cpld.fmu.given.sigma <- function(sigma, o, lim, z)
  {
    # Function to calculate the conditional distrn f(mu|sigma) in a Continuous Piecewise Linear Distrn
  
    # o is a riskband.cpld() output
    
    h0 <- o$h[-length(o$h)]
    s <- o$s
    r <- o$r

    mu.lim <- lim$mu


    # Find boundaries for each interval
    
    breakpoints <- o$x95 - z * sigma

      mu.lower <- breakpoints[-length(breakpoints)]
      mu.upper <- breakpoints[-1]

    
    # Remove regions that are out of bounds for current sigma

    w <- which(mu.upper <= mu.lim[1] | mu.lower >= mu.lim[2])

    if (length(w) > 0)
    {
      h0 <- h0[-w]
      s <- s[-w]
      r <- r[-w]
      mu.lower <- mu.lower[-w]
      mu.upper <- mu.upper[-w]
    }


    # Distance between a point (mu, sigma) and the lower bound (given by line mu+z*sigma = x95) to its left
    # is given by t = (mu - mu_lower) * cos(theta)
    #
    # The density f(mu|sigma) is given by
    # f = h0 + s * t 
    #   = h0 + s * cos(theta) * (mu - mu_lower)
    #   = h0 - mu_lower * s * cos(theta) + s * cos(theta) * mu
    #   = b + m * mu
  
    b <- h0 - mu.lower * s * o$geometry$cos.theta
    m <- s * o$geometry$cos.theta


    mu.lower <- pmax(mu.lower, mu.lim[1])
    mu.upper <- pmin(mu.upper, mu.lim[2])
  
    list(m=m, b=b, r=r, lower=mu.lower, upper=mu.upper)
  } # end of cpld.fmu.given.sigma
  

  # new_0.15
  cpld.fsigma.given.mu <- function(mu, o, lim, z)
  {
    # Function to calculate the conditional distrn f(sigma|mu) in a Continuous Piecewise Linear Distrn
    # o: a riskband.cpld() output

    h0 <- o$h[-length(o$h)]
    s <- o$s
    
    sigma.lim <- lim$sigma

  
    # Find boundaries for each interval

    breakpoints <- (o$x95 - mu) / z

      sigma.lower <- breakpoints[-length(breakpoints)]
      sigma.upper <- breakpoints[-1]
      

    # Remove regions that are out of bounds for current sigma

    w <- which(sigma.upper <= sigma.lim[1] | sigma.lower >= sigma.lim[2])

    if (length(w) > 0)
    {
      h0 <- h0[-w]
      s <- s[-w]

      sigma.lower <- sigma.lower[-w]
      sigma.upper <- sigma.upper[-w]
    }
    
    
    # Distance between a point (mu, sigma) and the lower bound (given by line mu+z*sigma = x95)
    # is given by t = (sigma - sigma_lower) * sin(theta)
    #
    # The density f(sigma|mu) is given by
    # f = h0 + s * t = h0 + s * (sigma - sigma_lower) * sin(theta)
    #   = h0 - s * sigma_lower * sin(theta) + s * sin(theta) * sigma
    #   = b + m * sigma
  
  
    b <- h0 - s * sigma.lower * o$geometry$sin.theta
    m <- s * o$geometry$sin.theta
    
    
    sigma.lower <- pmax(sigma.lower, sigma.lim[1])
    sigma.upper <- pmin(sigma.upper, sigma.lim[2])
  
    list(m=m, b=b, lower=sigma.lower, upper=sigma.upper)
  } # end of cpld.fsigma.given.mu
  
  
  CPL.prior <- function(region.prior.prob, lim, z, A, outcome.is.logNormally.distributed=T, precision=1e-4)
  {     
    # --- Useful functions -----------------------------------------------------  


    better.result <- function(result1, result2)
    {
      b1 <- c(result1 < result2, T)
      b2 <- c(result1 > result2, T)
      
      wb1 <- which(b1)[1]
      wb2 <- which(b2)[1]
      
      wb1 <= wb2 # include equality as being better -- even if not 'better' strictly speaking
    } # end of better.result
 
    
    geometry <- function(region.prior.prob, lim, z, A, outcome.is.logNormally.distributed=T, max.deepness=1/3)
    {
      if (outcome.is.logNormally.distributed) A <- log(A)
      
      theta <- atan(z)
      cos.theta <- cos(theta)
      R <- length(A) + 1

      mu.lim    <- lim$mu
      sigma.lim <- lim$sigma

      H <- diff(sigma.lim)
      L <- H / cos.theta
      Bt <- H*z               # base of corner triangle
      d <- Bt * cos.theta     # length of chord going from (mu, sigma) domain outer corner to triangle diagonal

      m <- L / d

        mu0 <- mu.lim[1]
        mu1 <- mu.lim[2]
        sigma1 <- sigma.lim[2]
        mu.top <- c(mu0, A - z*sigma1, mu1 - Bt)

      B <- diff(mu.top)
      b <- B*cos.theta

      area <- B*H
        w <- c(1,5)
        area[w] <- area[w] + Bt * H / 2

      u <- region.prior.prob / area
      p <- region.prior.prob / L
      hbar <- p / b  # meaningless for the left-end and right-end regions
                
          
      # regression coefficients that permit to calculate s1 from s2 in paralleliped regions
      
      s1.intercept <- 8/b * hbar
      s1.h0 <- -8/b
         
         
      # Calculate the floor level to respect max deepness
      # and the highest height (on either side of the riskband)
      # to be able to respect max deepness
      
      lo <- u * max.deepness    # max deepness = 1/3 of the piecewise-uniform level
      hi <- 4 * hbar - 3 * lo
        hi[c(1,5)] <- Inf # meaningless for the left-end and right-end regions
        
        # (Note that even in the absence of floor/deepness consideration, the height
        #  of the distrn fct at the border of any riskband is still constrained
        #  by the fact that it cannot go in the negative values [at the riskband mid-section 
        #  or on the other side of the riskband] to compensate.)
 
      
      # At the border between two riskbands, the lower value allowed for the distrn height
      # is the lower of the two adjascent riskand 'lo' levels
      
      ground <- pmin(lo[-R], lo[-1])
        
      lo[c(1,5)] <- 0    


      # See whether or not it is possible to design a distrn curve that does not go 'deep' in at least one riskband
      
      impossible <- any(hi[-1] < lo[-R])
      impossible <- impossible | any(hi[-R] < lo[-1])   
          
      # inner-region lower points 
      #   (for the calculation of density curves appropriateness, 
      #   adding a cost to curves going too low when compared to piecewise-uniform prior inner-regions)
      
      u.low <- pmin(u[-R], u[-1])
      u.low <- c(u.low[1], u[2], u.low[2], 
                           u[3], u.low[3],
                           u[4], u.low[4])

      # Region bounds

      x95b <- c(mu.lim[1] + z*sigma.lim[1], A, mu.lim[2] + z*sigma.lim[2])
      
   
      R1b1 <- c(x95b[1], x95b[2])
      R1b2 <- c(x95b[1] + c(0, Bt), x95b[2])
      
      
      inner.b <- rep(NA, 2*R-3)
        tmp.inner.b <- x95b[-c(1,R+1)]
        inner.b[2*seq(R-1)-1] <- tmp.inner.b
        # inner riskbands midpoints
        inner.b[2*seq(R-2)] <- c(mean(tmp.inner.b[c(1,2)]),
                                 mean(tmp.inner.b[c(2,3)]),
                                 mean(tmp.inner.b[c(3,4)]))
      
      
      R5b1 <- c(x95b[R], x95b[R+1])
      R5b2 <- c(x95b[R], x95b[R+1] - c(Bt, 0))
      
      
      endregion.is.local.min <- rep(NA, 2)
      uniq.u <- u[!duplicated(u)]
      endregion.is.local.min[1] <- uniq.u[1] < uniq.u[2]
      
      uniq.u <- rev(u)
      uniq.u <- uniq.u[!duplicated(uniq.u)]
      endregion.is.local.min[2] <- uniq.u[1] < uniq.u[2]
      
      j <- c(1,5)
      s.fade.out <- u[j] * max.deepness / (d + b[j]) * c(1, -1)
           

      list(z=z, H=H, L=L, Bt=Bt, d=d, m=m, cos.theta=cos.theta,
           mu.lim=mu.lim, sigma.lim=sigma.lim, R=R,
           B=B, b=b, area=area, u=u, u.low=u.low,
           R1b1=R1b1, R1b2=R1b2, inner.b=inner.b,
           R5b1=R5b1, R5b2=R5b2, 
           region.prior.prob=region.prior.prob, hbar=hbar, p=p,
           endregion.is.local.min=endregion.is.local.min, s.fade.out=s.fade.out,
           hi=hi, lo=lo, ground=ground, impossible=impossible,        
           outcome.is.logNormally.distributed=outcome.is.logNormally.distributed,
           x95b=x95b, # useful for call to *plot fct only
           s1.intercept=s1.intercept, s1.h0=s1.h0, reverse=F, 
           max.deepness=max.deepness)
    } # end of geometry
    
    
    geometry.completeDefn <- function(g)
    {
      if (g$endregion.is.local.min[2]) 
      {
        # Correction to higher entry level in region R-1 when region R is a local min
        # (to ensure a density curve going down at entry of region R, at a higher point than u[R])
        
        r <- g$R - 1
        u <- g$u[g$R]
        hbar <- g$hbar[r]
        
        new.hi <- 4*hbar - 3*u
        if (g$hi[r] > new.hi) g$hi[r] <- new.hi
      } 
      
      
      if (!g$impossible & g$endregion.is.local.min[1])
      {
        d <- g$d
        b <- g$b[1]
        p <- g$p[1]
        
        k1 <- d^2/3 + d*b
        
        #           __
        # plateau  /
        
        s <- p / k1
        plateau.h <- s * d
        
        if (plateau.h > g$hi[2]) g$impossible <- T
      }
      
      g
    } # end of geometry.completeDefn
  
    
    geometry.reverse <- function(g)
    {
      g$reverse <- !g$reverse
      
      g$B  <- rev(g$B)
      g$b  <- rev(g$b)
        
      g$hbar    <- rev(g$hbar)
      g$hi      <- rev(g$hi)
      g$lo      <- rev(g$lo)
      g$ground  <- rev(g$ground)
      g$p       <- rev(g$p)
        
      g$u                       <- rev(g$u)
      g$u.low                   <- rev(g$u.low)
      g$region.prior.prob       <- rev(g$region.prior.prob)
      g$endregion.is.local.min  <- rev(g$endregion.is.local.min)
        
      g$s1.intercept <- rev(g$s1.intercept)
      g$s1.h0        <- rev(g$s1.h0)
      
      g$inner.b <- -rev(g$inner.b)
      g$x95b    <- -rev(g$x95b)
      
  
            
      tmp.R5b1 <- g$R5b1
      tmp.R5b2 <- g$R5b2
        
      g$R5b1 <- -rev(g$R1b1)
      g$R5b2 <- -rev(g$R1b2)
        
      g$R1b1 <- -rev(tmp.R5b1)
      g$R1b2 <- -rev(tmp.R5b2)
      
        
      g
    } # end of geometry.reverse 
    
  
    smooth <- function(g, epsilon=precision)
    {
      # -- Useful fcts ---------------------------------------------------------
      
      angle.between.slope.changes <- function(s)
      {
        s.module <- sqrt(1+s^2)
          
        S <- length(s)
        scalar.product <- 1 + s[-1]*s[-S]
        tmp <- pmin(scalar.product / s.module[-1] / s.module[-S], 1) # protection against numerical imprecision
    
        acos(tmp)
      } # end of angle.between.slope.changes 
      
      
      augmented.body <- function(o, appendix, right.side=T) 
      {
        trim.appendix <- length(appendix$h) > length(appendix$s)

        if (right.side)
        {
          if (trim.appendix) appendix$h <- appendix$h[-1]
          
          o$h <- c(o$h, appendix$h)
          o$s <- c(o$s, appendix$s)
        }
        else
        {
          if (trim.appendix) o$h <- o$h[-1] # we trim o$h [same result as trimming new.h]

          o$h <- c(appendix$h, o$h)
          o$s <- c(appendix$s, o$s)
        }
        
        o
      } # end of augmented.body
    
    
      best.guess <- function(h0, reg4)
      {
        n <- reg4$n
        
        if (n > 2)
        {
          # Do a linear extrapolation
          
          b.denom <- n*reg4$x2.sum - reg4$x.sum**2
          
          if (b.denom == 0)
          {
            theta <- Inf
            margin <- 0 # irrelevant
          }
          else
          {
            b <- (n*reg4$xy.sum - reg4$x.sum*reg4$y.sum) / b.denom
            a <- (reg4$y.sum - b*reg4$x.sum) / n
            s <- a + b*h0 
            theta <- atan(s)
          
            sd2 <- (reg4$y2.sum - 2*a*reg4$y.sum - 2*b*reg4$xy.sum + n*a*a + 2*a*b*reg4$x.sum + reg4$x2.sum*b*b) / (n-1)
            if (sd2 < 0) sd2 <- 0 # protection against numerical imprecision
            
          
            if (sd2 > 0)
            {
              sd <- sqrt(sd2)
              s.upper <- s + sign(s)*sd
              margin <- abs(atan(s.upper) - theta)
            }
            else margin <- 0.01
          }
        }
        else
        {
          theta <- Inf
          margin <- 0 # irrelevant
        }
        
        
        list(init=theta, margin=margin)
      } # end of best.guess 
                
      
      closing.riskband <- function(g, G, trimmed.body, foldable=F)
      {        
        # G: closing riskband geometrics
        
        h0 <- ifelse(G$is.last.riskband, trimmed.body$h[length(trimmed.body$h)], trimmed.body$h[1]) # the height of the density on the inner side of
        # the closing riskband

  
        if (h0 >= G$h0.SE | (h0 >= G$u & G$is.local.min))
        {
          clothesLine.ok <- h0 <= G$h0.SE
          lb <- line.break(g, G, h0)
          
          choose.lb <- T
          
          if (clothesLine.ok)
          {
            cl <- clothesLine(G, h0)
                        
            out.cl <- roughness(g, trimmed.body, cl)
            out.lb <- roughness(g, trimmed.body, lb)
            
            choose.lb <- better.result(out.lb, out.cl)
            
            
            if (!choose.lb)
            {
              if (G$is.last.riskband)  x95b <- g$R5b1
              else                     x95b <- g$R1b1
              
              return(list(h=cl$h, s=cl$s, x95b=x95b, is.last.riskband=G$is.last.riskband))
            }
          }
   
          
          # choose.lb = TRUE
          
          return(list(h=lb$h, s=lb$s, x95b=lb$x95b, is.last.riskband=G$is.last.riskband))
        }
        else if (!G$is.local.min)
        {
          if (foldable)  tmp <- slanted.tipi(g, G, h0, trimmed.body)
          else           tmp <- hill.top(G, h0)
          
          if (G$is.last.riskband)  x95b <- g$R5b2
          else                     x95b <- g$R1b2
          
          return(list(h=tmp$h, s=tmp$s, x95b=x95b, is.last.riskband=G$is.last.riskband))
        }
        else
        {
          lb <- line.break(g, G, h0)
          ht <- hill.top(G, h0)
              
            if (G$is.last.riskband)  x95b <- g$R5b2
            else                     x95b <- g$R1b2


          out.lb <- roughness(g, trimmed.body, lb)
          out.ht <- roughness(g, trimmed.body, ht)                            
          
          if (better.result(out.lb, out.ht))  return(list(h=lb$h, s=lb$s, x95b=x95b, is.last.riskband=G$is.last.riskband))
          else                                return(list(h=ht$h, s=ht$s, x95b=x95b, is.last.riskband=G$is.last.riskband))
        }
      } # end of closing.riskband
      
                       
      closing.riskband.flexible <- function(g, G, h0, trimmed.rd)
      {
        # G: closing riskband geometrics
        # trimmed.rd:  a list with dimensions (h, s)
        #
        # NOTE: This fct is called for regions that are 'endregion.is.local.min=T' ONLY

        convex <- h0 < G$h0.SE
        s1.b <- s1b(G, h0) 
        
        
        # ----------------------------------------------------------------------
        # We are searching for optimal (s1, s2, t)
        
        # Define range on s1 (and its initial value)
        
        if (!convex)
        {
          truncating.angle <- truncation(G, h0, g$epsilon)$angle
          
          rb5.range <- c(-pi/2, truncating.angle)
          rb5 <- bs.init(rb5.range, g$epsilon, margin=0.015, start.hi=T, untouchable.bounds=T, force.1pass=T) 
        }
        else
        {
          s.clothesLine <- clothesLine(G, h0)$s
          if (!G$is.last.riskband) s.clothesLine <- -s.clothesLine
          
          rb5.range <- c(atan(s.clothesLine), pi/2) 
          rb5 <- bs.init(rb5.range, g$epsilon, margin=0.015, start.lo=T, untouchable.bounds=T, force.1pass=T) 
        }
        
   

        while (rb5$continue)
        {
          s1 <- tan(rb5$value)  # stands for s5.a, really     
       
          # Find (s2, t) such that the integral \int_0^D h(x)l(x)dx = region prob
          #              [with  theta2 = atan(s2)]
          
          tmp <- NewtonRaphson.lastriskband.ts2(s1, h0, G, g$epsilon, convex=convex, s1.b=s1.b)     

          s2 <- tmp$s2
          t  <- tmp$t
          
          s5 <- c(s1, s2)
          h5 <- c(tmp$h.intersect, 0)
          
          closing.rb <- closing.riskband.parms(h5, s5, t, G$is.last.riskband)
          
          out <- roughness(g, trimmed.rd, closing.rb)
                
          rb5 <- bs.update(rb5, out, closing.rb)
        } # end of while-5 
        
        
        # Closing riskband parms
        
        parms <- rb5$monitor
    
        if (G$is.last.riskband)  x95b <- c(g$R5b1[1], g$R5b1[1] + parms$t/g$cos.theta, g$R5b1[2]) 
        else                     x95b <- c(g$R1b1[1], g$R1b1[2] - parms$t/g$cos.theta, g$R1b1[2])
        
        parms$x95b <- x95b
        parms
      } # end of closing.riskband.flexible
      
      
      closing.riskband.geo <- function(g, is.last.riskband)
      {
        r <- ifelse(is.last.riskband, g$R, 1)
        
        b <- g$b[r]
        p <- g$p[r]
        d <- g$d
        D <- d + b
        
        b2.2 <- b^2/2
        
        # South-East & South-South-East heights in closing riskbands
      
        h0.SSE <- 2 * p / b
        h0.SE  <- p * D / (b2.2 + b*d + d^2/3)
        
        endregion.index <- ifelse(is.last.riskband, 2, 1)
        is.local.min <- g$endregion.is.local.min[endregion.index]
        
        
        list(b=b, d=d, D=D, L=g$L, p=p, prob=g$region.prior.prob[r], u=g$u[r],
             is.local.min=is.local.min, h0.SE=h0.SE, h0.SSE=h0.SSE,
             b2.2 = b2.2, b3.6 = b^3/6, pd = p*d, d2.6 = d^2/6, 
             is.last.riskband=is.last.riskband)
      } # end of closing.riskband.geo
      
      
      closing.riskband.parms <- function(h, s, t, is.last.riskband)
      {        
        if (!is.last.riskband)
        {
          h <-  rev(h)
          s <- -rev(s)
        }
        
        list(h=h, s=s, t=t, is.last.riskband=is.last.riskband)
      } # end of closing.riskband.parms
      
           
      clothesLine <- function(G, h0)
      {
        # G: closing riskband geometrics
        
        b <- G$b
        D <- G$D
      
        a <- G$b2.2 + D*(D+b) / 2 - (D^3 - b^3) / (3*G$d)
        delta <- h0 * (D+b) / 2 
        s <- (G$p - delta) / a
        
        h <- h0 + s*D
        if (!G$is.last.riskband) s <- -s
        
        list(h=h, s=s, is.last.riskband=G$is.last.riskband)
      } # end of clothesLine
      
      
      drop.closing.riskband <- function(rd, drop.last.riskband=T)
      {
        n2drop <- ifelse(drop.last.riskband, length(rd$R5b)-1, length(rd$R1b)-1)
        
        if (drop.last.riskband)
        {  
          j2drop <- seq(to=length(rd$s), length=n2drop)
          rd$s <- rd$s[-j2drop]
          
          j2drop <- seq(to=length(rd$h), length=n2drop)
          rd$h <- rd$h[-j2drop]
        }
        else
        {
          j2drop <- seq(n2drop)
          rd$s <- rd$s[-j2drop]
          rd$h <- rd$h[-j2drop]
        }
        
        rd
      } # end of drop.closing.riskband
      
      
      drop.elements <- function(x, n2drop, tail=T)
      {
        if (tail) w <- seq(to=length(x), length=n2drop)
        else      w <- seq(n2drop)
        
        x[-w]
      } # end of drop.elements
      
      
      hill.top <- function(G, h.bottom)
      { 
        # G: closing riskband geometrics
          
        b <- G$b
        L <- G$L
        a <- L * G$d / 2
          
        k <- h.bottom * (L*b + a)
        s.coeff <- L* b^2 / 2 + a*b
        
        s <- (G$prob - k) / s.coeff
        sb <- s*b
        
        if (G$is.last.riskband)  s <- c(s,  0)
        else                     s <- c(0, -s)

        h <- rep(h.bottom + sb, 2)
        
        list(h=h, s=s, is.last.riskband=G$is.last.riskband)
      } # end of hill.top
                
        
      line.break <- function(g, G, h0, closer2vertical.fraction=0.25)
      {
        # Draw a broken line (into two linear pieces) from the next-to-last region height to 0 at the upper/outer value for x95;
        # the line break occurs at 'b' or before if necessary.
        #
        # G: closing riskband geometrics
        
        b <- G$b
        d <- G$d
        D <- G$D

    
        if (h0 < G$h0.SSE)
        {
          # Fix occurs in the upper triangle => we will break the line @ the x95 line going through
          #                                                                      inner corner (b)
          
          t <- c(b, D)
          a2 <- diff(t^2)/2
          a3 <- diff(t^3)/3
    
          a <- d*D^2 -2*D*a2 + a3
            
          s.coeff <- -b*d/2 - a/d
          k <- h0*b/2
          s <- (G$p - k) / s.coeff
            
          h.star <- -s*d
          s.star <- (h.star - h0) / b
          x <- b
        }
        else
        {
          # Fix occurs in the parallelepiped => we will break the line on a secant that
          # is 25% closer to the vertical (dropping to 0 on x=max(A)) than the truncating line
            
          truncating.angle <- truncation(G, h0, g$epsilon)$angle
            # slope of the truncating line  (we want to avoid truncating, hence the new slope 
            #                                s.star calculated below)
                       
          angle  <-  -pi/2 *closer2vertical.fraction + (1-closer2vertical.fraction)*truncating.angle
          s.star <-  tan(angle) 
          
          
          # Find where the line passing through (0, h0) with slope s.star and
          # line passing through (D,0) with slope s meet (call it x)
          #
          # Find x & s using a Newton-Raphson algorithm
          
          tmp <- NewtonRaphson.lastriskband.ts2(s.star, h0, G, g$epsilon, t.lt.b=T)
          
          s <- tmp$s2
          x <- tmp$t
          h.star <- tmp$h.intersect          
        }
        
        
        # Wrap up        
    
        h <- c(h.star, 0)
        s <- c(s.star, s)
        
        
        if (G$is.last.riskband)
        {
          lb <- g$R5b1[1] + x/g$cos.theta
          x95b <- c(g$R5b1[1], lb, g$R5b1[2])
        }
        else
        {
          h <-  rev(h)
          s <- -rev(s)
          
          lb <- g$R1b1[2] - x/g$cos.theta
          x95b <- c(g$R1b1[1], lb, g$R1b1[2])
        } 
        

        list(h=h, s=s, x95b=x95b, is.last.riskband=G$is.last.riskband)
      } # end of line.break
      
      
      NewtonRaphson.lastriskband.ts2 <- function(s1, h0, G, epsilon, convex=F, s1.b=0,
                                                 t.lt.b = h0 > G$h0.SSE | ifelse(convex, s1 > s1.b, s1 < s1.b))
      { 
        # G: closing riskband geometrics
        
        D <- G$D
          
          
        if (t.lt.b)
        {      
          # --- Find t < b -----------------------------------------------------
               
          t <- 0
          continue <- T
             
          while (continue)
          {
            s2 <- (h0 + s1*t) / (t-D)
            s2.prime <- -(h0 + s1*D) / (t-D)^2
            
            v <- (D-t)^2 / 2
            v.prime <- t - D
        
            f <- h0*t + s1*t^2/2 - s2*v + s2*G$d2.6  
            f.prime <- h0 + t*s1 - s2.prime*v - s2*v.prime + s2.prime*G$d2.6
            
            change <- (G$p - f) / f.prime
            t <- t + change
            continue <- abs(change) >= epsilon
          }
          
          s2 <- (h0 + s1*t) / (t-D)
          h.intersect <- -s2 * (D-t)
          
          return(list(s2=s2, t=t, h.intersect=h.intersect))            
        }
        else
        {                  
          # --- Find t > b -----------------------------------------------------
          
          target <- G$pd + h0*G$b2.2 + s1*G$b3.6
          a1 <- h0 * D
          a2 <- (s1*D - h0) / 2
          a3 <- -s1/3
           
          t <- G$b
          continue <- T
           
          while (continue)
          {
            u <- (h0 + s1*t) / 3
            v <- (D - t)^2
             
            f <- a1*t + a2*t^2 + a3*t^3 + u*v
            f.prime <- a1 + 2*t*a2 + 3*t^2*a3 + s1*v/3 - 2*u*(D-t)
             
            change <- (target - f) / f.prime
            t <- t + change
            continue <- abs(change) >= epsilon
          }
          
          h.intersect <- h0 + s1*t
          s2 <- h.intersect / (t-G$D) 
          
          return(list(s2=s2, t=t, h.intersect=h.intersect))
        }
      } # end of NewtonRaphson.lastriskband.ts2 --------------------------------
     
      
      parallelepiped <- function(r, g, entry.angle, h0)
      {
        b <- g$b[r]
          
        s0 <- tan(entry.angle)
        
        h0.middle <- h0 + s0*b/2
        if (h0.middle < 0 & h0.middle > -1e-10) h0.middle <- 0
        
        s1 <- g$s1.intercept[r] + g$s1.h0[r] * h0 - 3*s0
        h0.end <- h0.middle + s1*b/2
          
        if (h0.end < 0 & h0.end > -1e-10)  h0.end <- 0
        
        list(h=c(h0, h0.middle, h0.end), s=c(s0,s1))
      } # end of parallelepiped
      
      
      parallelepiped.connexion <- function(r, g, h0, h1, left2right=T)
      {
        # h0, h1: height @ left/right ends of riskband r [in that order]
        # if left2right = T, return s0
        #    otherwise       return -s1
        
        b <- g$b[r]
        s0 <- (4*g$hbar[r] - 3*h0 - h1) / b
        
        if (left2right)
        {
          return(s0)
        }
        else
        {
          s1 <- 2/b * (h1-h0) + s0
          return(-s1)
        }  
      } # end of parallelepiped.connexion
      
      
      parallelepiped.connexion.t <- function(r, g, h0, h1, t)
      {
        delta.h <- h1 - h0
        v <- c(delta.h, g$p[r] - h0[1]*g$b[r])
        s.solve(t, g$b[r], v)
      } # end of parallelepiped.connexion.t
      
      
      parallelepiped.range <- function(r, g, h0)
      {
        hbar <- g$hbar[r]
        b <- g$b[r]
        
        s.max <- (4*hbar - 3*h0 - g$ground[r]) / b  # stalagmite
         
        reflexion <- parallelepiped.reflexion(r, g, h0)
        if (reflexion > g$hi[r+1]) s.min <- parallelepiped.connexion(r, g, h0, g$hi[r+1])
        else s.min <- 2 / b * (g$lo[r] - h0)  # stalactite
    
    
        s <- c(s.min, s.max)
    
        atan(s)
      } # end of parallelepiped.range
      
      
      parallelepiped.reflexion <- function(r, g, h)
      {
        4*g$hbar[r] - 2*g$lo[r] - h
      } # end of parallelepiped.reflexion
      
        
      reg4.init <- function()
      {
        list(n=0, x.sum=0, x2.sum=0, xy.sum=0, y.sum=0, y2.sum=0)
      } # end of reg4.init
      
      
      reg4.update <- function(reg4, s4_entry, h04)
      {
        # Prepare for the calculation of regression parms between x & y, where
        #   x: h04
        #   y: entry_angle
        
        entry_angle <- atan(s4_entry)
        
        reg4$n      <- reg4$n + 1
        reg4$x.sum  <- reg4$x.sum   + h04
        reg4$x2.sum <- reg4$x2.sum  + h04**2
        reg4$xy.sum <- reg4$xy.sum  + h04 * entry_angle
        reg4$y.sum  <- reg4$y.sum   + entry_angle
        reg4$y2.sum <- reg4$y2.sum  + entry_angle**2
    
        return(reg4)
      } # end of reg4.update
      
      
      revisit.closing.riskband <- function(g, G, rd)
      {
        # G: closing riskband geometrics
        # rd: revisited density
        #
        # NOTE: This fct is called for regions that are 'endregion.is.local.min=T' ONLY
        
        rd <- drop.closing.riskband(rd, G$is.last.riskband)
         
          # height when entering closing riskband
          h0 <- ifelse(G$is.last.riskband, rd$h[length(rd$h)], rd$h[1])         
           
        rb5 <- closing.riskband.flexible(g, G, h0, rd)

        rd <- augmented.body(rd, rb5, G$is.last.riskband)

        if (G$is.last.riskband)  rd$R5b <- rb5$x95b
        else                     rd$R1b <- rb5$x95b
       
        rd
      } # end of revisit.closing.riskband
      

      revisit.closing2riskbands <- function(g, G, rd)
      {
        # G: closing riskband geometrics
        # rd: revisited density
        #
        # NOTE: This fct is called for the last two riskbands ONLY when the last one is 'endregion.is.local.min=T' 
        
        is.last2riskbands <- G$is.last.riskband
        
        last2rb <- ifelse(is.last2riskbands, 2, 1)       
        
        
        # _____________________________________________________________________________________
        # Remove elements corresponding to revisited riskbands        
        
        rb5.nSegments <- length(rd$R5b) - 1
        
        nSegments.closing.rb <- ifelse(is.last2riskbands, rb5.nSegments, length(rd$R1b) - 1)
        
        
        # Drop h & s from last 2 riskbands (do NOT drop any elements from rd$inner.b)
        
        rd$h <- drop.elements(rd$h, nSegments.closing.rb+2, tail=is.last2riskbands)
        rd$s <- drop.elements(rd$s, nSegments.closing.rb+2, tail=is.last2riskbands)
                                
        # _____________________________________________________________________________________
        
        # h at entry in the riskband *preceding* the closing riskband
        h04 <- ifelse(G$is.last.riskband, rd$h[length(rd$h)], rd$h[1])   
        
        # Riskband 4
        
        n2lr <- ifelse(is.last2riskbands, g$R-1, 2) # next-to-last riskband
        
        b4 <- g$b[n2lr]
          h0b <- h04 * b4
          p4 <- g$p[n2lr]
          
          s4.straight <- (p4 - h0b) * 2 / b4**2
          theta4.straight <- atan(s4.straight)
          
             
        rb4.t <- bs.init(c(0, b4), g$epsilon, init=b4/2, margin=b4/20, untouchable.bounds=T, force.1pass=T)
    
      
        while (rb4.t$continue)
        {
          t4 <- rb4.t$value
          
          if (h04 > g$u[n2lr])
          {
            s4.max <- parallelepiped.connexion.t(n2lr, g, h04, G$u, t4)[1]
            s4.min <- (g$ground[n2lr] - h04) / t4  # s-ground
            
            if (s4.max < s4.min) s4.max <- s4.min # possible for the highest values of h04
            
            s4.range <- c(s4.min, s4.max)  
          }
          else
          {
            s4.max <- parallelepiped.connexion.t(n2lr, g, h04, 0, t4)[1]
            s4.range <- c(s4.straight, s4.max)
          }
          
          
          theta4.range <- atan(s4.range)
          rb4.theta <- bs.init(theta4.range, g$epsilon, theta4.straight, 0.015)
    
    
          while (rb4.theta$continue)
          {
            theta4 <- rb4.theta$value
                   
            
            s4a <- tan(theta4)
            s4b <- -2 / (b4 - t4)^2 * (h04*t4 + s4a/2*t4^2 + (h04+s4a*t4)*(b4-t4) - p4)
                        
            s4 <- c(s4a, s4b)
            
            h4.t <- h04 + t4*s4a
            
            h4 <- h4.t + c(0, (b4 - t4) * s4b)
             
            
            # --- Riskband # 5 -------------------------------------------------
                        
            rb4 <- closing.riskband.parms(h4, s4, 0, is.last2riskbands)
            
            body4 <- augmented.body(rd, rb4, is.last2riskbands)
                
            closing.rb <- closing.riskband(g, G, body4)
            
            # ------------------------------------------------------------------
            
            
            rb4 <- augmented.body(rb4, closing.rb, is.last2riskbands)
            
            if (is.last2riskbands) rb5.nSegments <- length(closing.rb$s)
            
            
            complete.dens <- augmented.body(rd, rb4, is.last2riskbands)
            
            out <- roughness(g, complete.dens, rb5.nSegments=rb5.nSegments, last2rb=last2rb)
            
            parms <- list(rb4=rb4, x95b=closing.rb$x95b, rb5.nSegments=rb5.nSegments)
    
            rb4.theta <- bs.update(rb4.theta, out, parms)  
          } # end of while-rb4.theta loop
          
          
          parms <- rb4.theta$monitor
    
          complete.dens <- augmented.body(rd, parms$rb4, is.last2riskbands)
    
          out <- roughness(g, complete.dens, rb5.nSegments=parms$rb5.nSegments, last2rb=last2rb)
            parms$t4 <- t4     
      
          rb4.t <- bs.update(rb4.t, out, parms)
        } # end of while-rb4.t loop
        
      
        parms <- rb4.t$monitor
        
        
        # Wrap-up
        
        rd <- augmented.body(rd, parms$rb4, is.last2riskbands)
        
      
        if (is.last2riskbands)
        {
          rd$R5b <- parms$x95b
          rd$rb5.nSegments <- parms$rb5.nSegments
          
          rd$inner.b[6] <- rd$inner.b[5] + parms$t4 / g$cos.theta
        }
        else
        {
          rd$R1b <-  parms$x95b
          
          rd$inner.b[2] <- rd$inner.b[3] - parms$t4 / g$cos.theta
        }

        
        rd
      } # end of revisit.closing2riskbands

      
      revisit.inner.riskbands <- function(g, rd)
      {
        # rd: revisited density
        
        # Revisit parallelepiped inner riskbands; that is, we keep h-values at both ends unchanged
        # but see if we can change the breakpoint location between the two lines (and their slopes accordingly)
        # to get a smoother curve.
        
      
        rb1.nSegments <- length(rd$R1b) - 1
        rb5.nSegments <- length(rd$R5b) - 1
        
        h.index <- rb1.nSegments + 2*seq(3) # indices for the middle-riskband height of riskbands 2,3,4 [length 3]
        s.indices <- rb1.nSegments + seq(6) # indices for slopes in riskbands 2-3-4 (2 by rb)  [length 6]
          s.index2 <- s.indices[1:2]
          s.index3 <- s.indices[3:4]
          s.index4 <- s.indices[5:6]
                
        h0 <- rd$h[h.index-1] # h as the entry [left-side] in each riskband [length 3]
      
        delta.h <- rd$h[h.index+1] - h0  # increase in h in each inner riskband [length 3]
      
        v <- matrix(c(delta.h, g$p[2:4]-h0*g$b[2:4]), ncol=2) # 3 x 2 matrix
        
        
        h4.exit <- rd$h[h.index[3]+1]
    
        t4.min <- ifelse(h0[3] < g$u[4] & h4.exit > g$u[4], ((h0[3] + rd$h[h.index[3]+1])*g$b[4] - 2*g$p[4]) / delta.h[3], 0)
    
        # inits
        s3 <- rd$s[s.index3]
        s4 <- rd$s[s.index4]
        h.middle <- rd$h[h.index]
        
        
        # rb # 2     
        
        rb2 <- bs.init(c(0, g$b[2]), g$epsilon, g$b[2]/2, 0.05, untouchable.bounds=T)
    
          
        while (rb2$continue)
        {
          t2 <- rb2$value

          M2 <- matrix(c(t2, g$b[2]-t2, t2*g$b[2]-t2**2/2, (g$b[2]-t2)**2/2), ncol=2, byrow=T)
          s2 <- as.vector(solve(M2) %*% v[1,])
          
            h.middle[1] <- h0[1] + t2*s2[1]
            
            rd$h[h.index[1]] <- h.middle[1]
            rd$s[s.index2]   <- s2 
            
            
          # rb # 3
          
          rb3 <- bs.init(c(0, g$b[3]), g$epsilon, g$b[3]/2, 0.05, untouchable.bounds=T)
  
            
          while (rb3$continue)
          {
            t3 <- rb3$value

            M3 <- matrix(c(t3, g$b[3]-t3, t3*g$b[3]-t3**2/2, (g$b[3]-t3)**2/2), ncol=2, byrow=T)
            s3 <- as.vector(solve(M3) %*% v[2,])
            
              h.middle[2] <- h0[2] + t3*s3[1]
              
              rd$h[h.index[2]] <- h.middle[2]
              rd$s[s.index3]   <- s3  
              
              
            # rb # 4
            
            rb4 <- bs.init(c(t4.min, g$b[4]), g$epsilon, g$b[4]/2, 0.05, untouchable.bounds=T)
              
               
            while (rb4$continue)
            {
              t4 <- rb4$value

              M4 <- matrix(c(t4, g$b[4]-t4, t4*g$b[4]-t4**2/2, (g$b[4]-t4)**2/2), ncol=2, byrow=T)
              s4 <- as.vector(solve(M4) %*% v[3,])
              
                h.middle[3] <- h0[3] + t4*s4[1]
                
                rd$h[h.index[3]] <- h.middle[3]
                rd$s[s.index4]   <- s4        


              s <- c(s2, s3, s4)
              parms <- list(s=s4, t=t4, h.middle=h.middle)
            
              out <- roughness(g, rd, rb5.nSegments=rb5.nSegments, inner.rb=4)
              
              rb4 <- bs.update(rb4, out, parms)                          
            }  # end of while-4
            
            
            parms <- rb4$monitor
            
            s4 <- parms$s
            t4 <- parms$t
            h.middle <- parms$h.middle 

            rd$s[s.indices] <- c(s2, s3, s4)
            rd$h[h.index] <- h.middle
            
            
            parms <- list(s3=s3, t3=t3, s4=s4, t4=t4, h.middle=h.middle)
          
            out <- roughness(g, rd, rb5.nSegments=rb5.nSegments, inner.rb=3)
            
            rb3 <- bs.update(rb3, out, parms)                          
          } # end of while-3
          
          
          parms <- rb3$monitor
          
          s3 <- parms$s3
          t3 <- parms$t3
          s4 <- parms$s4
          t4 <- parms$t4
          
          h.middle <- parms$h.middle
          
          s <- c(s2, s3, s4)
          rd$s[s.indices] <- s
          rd$h[h.index] <- h.middle
          
  
          parms <- list(s=s, t=c(t2,t3,t4), h.middle=h.middle)
         
          out <- roughness(g, rd, rb5.nSegments=rb5.nSegments, inner.rb=2)       
          
          rb2 <- bs.update(rb2, out, parms)
        } # end of while-2
        
        
        parms <- rb2$monitor
        
        rd$h[h.index]   <- parms$h.middle
        rd$s[s.indices] <- parms$s
        
        j <- 2*seq(3)
        rd$inner.b[j] <- rd$inner.b[j-1] + parms$t / g$cos.theta
          
        rd$out <- roughness(g, rd, rb5.nSegments=rb5.nSegments) # we want to recalculate roughness without local considerations
    
        rd
      } # end of revisit.inner.riskbands
          
    
      roughness <- function(g, parms, closing.rb=list(),
                            rb5.nSegments = ifelse(length(parms$h)%%2 == 1, 1, 2),
                            inner.rb=0, last2rb=0)
      {
        # parms: a list with dimensions(h, s)
        # closing.rb: provided if one of the closing riskbands is re-estimated 
        #             (that is, AFTER the curve was estimated throughout, to see if we can find a 
        #              better density closing).
        #             > When provided, closing.rb is a list with dimensions $s, $h & $is.last.riskband
        #
        # last2rb: a numeric code, = 2 indicates we are focusing on the  last two riskbands
        #                          = 1 indicates we are focusing on the first two riskbands
        
        
        s <- parms$s
        h <- parms$h
        
        
        closing.rb.provided <- length(closing.rb) > 0
         
        
        if (closing.rb.provided)
        {
          if (closing.rb$is.last.riskband)
          {  
            s <- c(s, closing.rb$s)
            h <- c(h, closing.rb$h)
          }
          else
          {                              
            s <- c(closing.rb$s, s)
            h <- c(closing.rb$h, h)
          }

          if (closing.rb$is.last.riskband)  rb5.nSegments <- length(closing.rb$s)
          # else we keep the value received as a fct argument
        }

        
        s.len <- length(s)
        h.len <- length(h)
               
    
        # Correct for (most likely) numerical imprecision
        if (abs(s[1])     < 1e-10) s[1]     <- 0
        if (abs(s[s.len]) < 1e-10) s[s.len] <- 0
        
    
        # Find the sharpest turn / slope change
    
        sharpest.slopeChange <- sharpest.slope.change(s)
    
    
        # Count the number of sign changes
        
        s.sign <- sign(s)
        s.sign <- s.sign[s.sign != 0]
        number.of.sSign.changes <- sum(s.sign[-1] != s.sign[-length(s.sign)]) 
        
        
        if (!g$endregion.is.local.min[2] && rb5.nSegments > 1 && s[s.len-1] > 0 && s[s.len] < 0) number.of.sSign.changes <- number.of.sSign.changes - 1
        
        
        # Determine whether distribution is fully described or not
        
        rb1.nSegments <- s.len - 6 - rb5.nSegments
                
        dens.fct.is.fully.described <- ifelse(closing.rb.provided, length(parms$s) > 6, rb1.nSegments > 0)
        
            
        # See if density fct is going too deep (when compared to the piecewise uniform prior)

        if (!g$impossible)
        {        
          h.last.index <- h.len - rb5.nSegments
    
          h.first.index <- h.last.index - 6
          if (h.first.index < 1) h.first.index <- 1
          h.indices <- seq(from=h.first.index, to=h.last.index)
    
          ulow.indices <- seq(to=7, by=1, length=length(h.indices))
          too.deep <- any(h[h.indices] / g$u.low[ulow.indices] < g$max.deepness)
          
          
          # Final slope in closing riskbands should not be below some value, as indicated by s.fade.out 
          
          if (dens.fct.is.fully.described && h[1]     == 0 && s[1]     < g$s.fade.out[1])  too.deep <- T
          if (                               h[h.len] == 0 && s[s.len] > g$s.fade.out[2])  too.deep <- T
          
        }
        else too.deep <- -1
    
    
        # Consider initial/final value of density fct
    
        if (dens.fct.is.fully.described)
        {
          h.endpoints <- h[c(1, h.len)]     
          rb1.roughness <- ifelse(rb1.nSegments == 1, 0, abs(s[2]-s[1]))
          
          if (!g$endregion.is.local.min[1] && rb1.nSegments > 1 && s[1] > 0 && s[2] < 0) number.of.sSign.changes <- number.of.sSign.changes - 1
        }
        else
        {
          h.endpoints <- c(-1, -1)          
          rb1.roughness <- -1
        }
        
        
        rb5.roughness <- ifelse(rb5.nSegments == 1, 0, abs(s[s.len]-s[s.len-1]))
               
        endregions.roughness <- max(rb1.roughness, rb5.roughness)
        
        
        # Consider the number of inflexion points
        
        diffs.sign <- sign(diff(s))
        diffs.sign <- diffs.sign[diffs.sign!=0]
              
        n.inflexion.points <- sum(diffs.sign[-1] != diffs.sign[-length(diffs.sign)])
    
        
        extra.scale <- 0
        
        
        if (last2rb==2)
        {
          j <- seq(to=length(s), length=rb5.nSegments + 3)
          extra.scale <- sharpest.slope.change(s[j])
          endregions.roughness <- -1
        }
        else if (last2rb==1)
        {
          j <- seq(rb1.nSegments + 3)
          extra.scale <- sharpest.slope.change(s[j])
          endregions.roughness <- -1
        }
        else if (closing.rb.provided)
        {
          if (closing.rb$is.last.riskband)  j <- seq(to=s.len, length=3)
          else                              j <- seq(3)
          
          extra.scale <- sharpest.slope.change(s[j])
          endregions.roughness <- -1
        }
        else if (inner.rb > 0)
        {
          j <- seq(from=rb1.nSegments+2*inner.rb-3, to=rb1.nSegments+6)
          extra.scale <- max(abs(s[j]))
        }
        else
        {
          extra.scale <- max(abs(s))
          #extra.scale <- sd(s)
        }
        
        
        c(number.of.sSign.changes, too.deep, sharpest.slopeChange, max(h.endpoints), 
          n.inflexion.points, endregions.roughness, extra.scale)
      } # end of roughness
          
      
      s1b <- function(G, h0)
      {
        # G: closing riskband geometrics
        
        b <- G$b
        d <- G$d
        
        s1b <- (G$p - h0*(b + d/3)) / (b^2/2 + b*d/3)
                
        s1b
      } # end of s1b
      
      
      slanted.tipi <- function(g, G, h0, body)
      {
        # G: closing riskband geometrics
        
        d <- G$d
        b <- G$b
        
        j <- ifelse(G$is.last.riskband, 1, 2)
        theta.min <- atan(hill.top(G, h0)$s[j])
        if (!G$is.last.riskband)  theta.min <- -theta.min
        
        theta.max <- southEast.connexion(G, h0)
        
        f <- 6/d**2
        a <- (G$p - d/2*h0 - h0*b) * f
        s1.mult <- (b**2/2 + b*d/2) * f
        
        
        rb <- bs.init(c(theta.min, theta.max), g$epsilon)
        
        while (rb$continue)
        {
          s1 <- tan(rb$value)
          s2 <- a - s1*s1.mult
          
          s <- c(s1, s2)
          hb <- h0 + b*s1
          h <- c(hb, hb + s2*d)
          
          closing.rb <- closing.riskband.parms(h, s, 0, G$is.last.riskband)
        
          out <- roughness(g, body, closing.rb)
          
          rb <- bs.update(rb, out, closing.rb)
        }
      
        
        rb$monitor
      } # end of slanted.tipi
  
           
      southEast.connexion <- function(G, h0)
      {
        # G: closing riskband geometrics
        
        b <- G$b
        d <- G$d
        
        s <- (G$p - h0*(d/3+b)) / (b**2/2 + b*d/3)
        atan(s)
      } # end of southEast.connexion      
      
      
      s.solve <- function(t, b, v)
      {
        M <- c(t, b-t, t*b-t**2/2, (b-t)**2/2)
        M.det <- M[1]*M[4] - M[2]*M[3]
        s1 <- (M[4]*v[1] - M[2]*v[2]) / M.det
        s2 <-  (v[1] - s1*M[1]) / M[2]
        
        c(s1, s2)
      } # end of s.solve
                    
      
      sharpest.slope.change <- function(s)
      {
        max(abs(angle.between.slope.changes(s)))
      } # end of sharpest.slope.change    
      
      
      theta4.limits <- function(g, h04)
      {
        r <- g$R - 1  # region r = 4
        
        s.max <- parallelepiped.connexion(r, g, h04, g$u[g$R])
                
        # Plateau
        s.plateau <- 8/3/g$b[r] * (g$hbar[r] - h04)
        
        s <- c(s.plateau, s.max)
        
        atan(s)
      } # end of theta4.limits
            
      
      truncation <- function(G, h0, epsilon)
      {
        # G: closing riskband geometrics
        # Note: this fct does NOT adjust the sign of 's' for G$is.last.riskband on purpose
        
           
        if (h0 > G$h0.SSE)
        {
          # truncation occurs on t < b
          t <- 2*G$p / h0
        }
        else
        {
          target <- G$pd + h0 * G$b2.2
          a0 <- h0 * G$b3.6
          a1 <- h0 * G$D / 2
          a2 <- -h0 / 6


          t <- G$b
          continue <- T

          while (continue)
          {
            f <- a0/t + a1*t + a2*t^2
            f.prime <- -a0/t^2 + a1 + 2*t*a2

            change <- (f - target) / f.prime
            t <- t - change
            continue <- abs(change) >= epsilon
          }
        }
        
        
        s <- -h0/t
      
        list(s=s, t=t, angle=atan(s))
      } # end of truncation
      
    
      # --------------------------------------------------------------------------
      # Start of smooth
      
      g$epsilon <- epsilon
      reg4 <- reg4.init()

      g.left  <- closing.riskband.geo(g, F)
      g.right <- closing.riskband.geo(g, T)
    
       
      rb2.h0 <- bs.init(c(g$ground[1], g$hi[2]), epsilon, g$u[2], abs(g$u[2]-g$u[1]), force.1pass=T)
    
      
      while (rb2.h0$continue)
      {
        h02 <- rb2.h0$value
        
        # Region 2
       
        tmpr2 <- parallelepiped.range(r=2, g, h02)
    
        rb2  <- bs.init(tmpr2, epsilon, force.1pass=T)
         
    
        while (rb2$continue)
        {
          tmp <- parallelepiped(r=2, g, rb2$value, h02) 
          s2 <- tmp$s
          h2 <- tmp$h[-1]
          body2 <- list(h=c(h02, h2), s=s2)
      
    
          # Region 3
    
          h03 <- h2[2]
          tmpr3 <- parallelepiped.range(r=3, g, h03)
    
          rb3  <- bs.init(tmpr3, epsilon, force.1pass=T)
          
    
          while (rb3$continue)
          {
            tmp <- parallelepiped(r=3, g, rb3$value, h03)
            s3 <- tmp$s
            h3 <- tmp$h[-1]
            body3 <- list(h=c(h03, h3), s=s3)
         
    
            # Region 4
    
            h04 <- h3[2]
            if (g$endregion.is.local.min[2])  tmpr4 <- theta4.limits(g, h04)
            else                              tmpr4 <- parallelepiped.range(r=4, g, h04)
            
    
            bg4  <- best.guess(h04, reg4) 
            rb4  <- bs.init(tmpr4, epsilon, start=bg4, force.1pass=T)
          
    
            while (rb4$continue)
            {    
              tmp <- parallelepiped(r=4, g, rb4$value, h04)
              s4 <- tmp$s
              h4 <- tmp$h[-1]
              if (h4[2] < 0 & h4[2] > -1e-8) h4[2] <- 0 # correction for numerical imprecision
              body4 <- list(h=c(h04, h4), s=s4)
    
          
              # Region 5
              rb5 <- closing.riskband(g, g.right, body4)
                rb5$rb4 <- list(h=h4, s=s4)

              out4 <- roughness(g, body4, rb5)
              
              rb4 <- bs.update(rb4, out4, rb5)
            } # end of while-4 -------------------------------------------------
                        
            
            parms <- rb4$monitor  
                         
              reg4 <- reg4.update(reg4, rb5$rb4$s[1], h04)
              
            rb4 <- augmented.body(rb5$rb4, rb5)  
            tmp <- augmented.body(body3, rb4)
                                 
            out3 <- roughness(g, tmp)
    
            rb3 <- bs.update(rb3, out3, tmp)
          } # end of while-3 ---------------------------------------------------
         
          
          parms <- rb3$monitor
          tmp <- augmented.body(body2, parms)
                   
          out2 <- roughness(g, tmp)
     
          rb2 <- bs.update(rb2, out2, tmp)
        } # end of while-rb2 ---------------------------------------------------
        
        
        parms <- rb2$monitor
          
            
        # Complete the density by closing the 1st riskband
        
        if (g.left$is.local.min)  rb1 <- closing.riskband.flexible(g, g.left, h02, parms)
        else                      rb1 <- closing.riskband(g, g.left, parms, foldable=T)        
        
        out <- roughness(g, parms, rb1)
        
        rb2.h0 <- bs.update(rb2.h0, out, list(parms=parms, rb1=rb1))
      } # end of while-rb2.h0
                        
        
      tmp <- rb2.h0$monitor
      
      best <- tmp$parms # parms for rb#2 and after
        best$rb5.nSegments <- ifelse(length(best$s)%%2 == 1, 1, 2)
        
        best <- augmented.body(best, tmp$rb1, F)
        best$inner.b <- g$inner.b
        best$R1b <- tmp$rb1$x95b
        
        best$R5b <- switch(best$rb5.nSegments, g$R5b1, g$R5b2)
               
        
      # Revisit the first/last two riskbands if necessary      
      
      # (first two riskbands)
      if (g.left$is.local.min)
      { 
        if (best$s[length(best$R1b)] < 0)  best <- revisit.closing2riskbands(g, g.left, best)
        best <- revisit.closing.riskband(g, g.left, best)
      }
      else
      {
        # Revisit first riskband parms for a better fit
        best <- drop.closing.riskband(best, F)
        rb1 <- closing.riskband(g, g.left, best, foldable=T)
        best <- augmented.body(best, rb1, F)
        best$R1b <- rb1$x95b
      }
      

      
      # (last two riskbands)
      if (g.right$is.local.min)
      {
        if (best$s[length(best$s)+1-length(best$R5b)] > 0)  best <- revisit.closing2riskbands(g, g.right, best)
        best <- revisit.closing.riskband(g, g.right, best)
      }
      else
      {
        # Revisit last riskband parms for a better fit
        best <- drop.closing.riskband(best)
        rb5 <- closing.riskband(g, g.right, best, foldable=T) 
        best <- augmented.body(best, rb5)
        best$R5b <- rb5$x95b
      }
      

      
      # Revisit the inner regions breakpoints & their corresponding heights
            
      best <- revisit.inner.riskbands(g, best) # out is calculated herein for the last time
      
      
      # Rearrange best's dimensions in sensible order
      
      if (g$reverse)
      {
        best$h <-  rev(best$h)
        best$s <- -rev(best$s)
            
        tmp      <-  best$R5b
        best$R5b <- -rev(best$R1b)
        best$R1b <- -rev(tmp)
        
        best$inner.b <- -rev(best$inner.b)
      }
              
      
      best$r <- rep(seq(5), c(length(best$R1b)-1, 2, 2, 2, length(best$R5b)-1))
      best$x95 <- c(best$R1b, best$inner.b[-1], best$R5b[-1])
      best$reverse <- g$reverse
      
      best    
    } # end of smooth
      
  
    # --------------------------------------------------------------------------
    # Start of CPL.prior
        
    
    # Prepare geometry
    
    g <- geometry(region.prior.prob, lim, z, A, outcome.is.logNormally.distributed)
      g.rev <- geometry.reverse(g)
    
      g     <- geometry.completeDefn(g)    
      g.rev <- geometry.completeDefn(g.rev)
                  
      
    # Actually estimate the best density curve
    
    best <- smooth(g)

    tmp <- smooth(g.rev)
      if (better.result(tmp$out, best$out)) best <- tmp
      
      
    # Wrap-up
  
    g$sin.theta <- sqrt(1-g$cos.theta^2)
    best$geometry <- g
    
    
    best
  } # end of CPL.prior
  
  
  invgamma.logpcum <- function(lim, alpha, beta)
  {
    if (length(lim) > 1)
    {
      invgamma.lim <- c(min(lim), max(lim))
      gamma.lim <- rev(1/invgamma.lim)
      p1 <- pgamma(gamma.lim[1], alpha, beta, log.p=T, lower.tail=T)
      p2 <- pgamma(gamma.lim[2], alpha, beta, log.p=T, lower.tail=F)
      lower.tail <- p1 < p2
      gamma.lim <- rev(1/lim)
    }
    else
    {
      gamma.lim <- 1/lim
      gamma.mode <- (alpha-1)/beta
      lower.tail <- gamma.lim < gamma.mode
    }

    log.pcum <- pgamma(gamma.lim, alpha, beta, log.p=T, lower.tail=lower.tail)
    list(log.pcum=rev(log.pcum), lower.tail=!lower.tail)
  } # end of invgamma.logpcum
  
  
  Location <- function(theta, lim, z)
  { 
    # modif_0.15 This fct's code was in the body of the main function in earlier versions

    mu.lower <- lim$mu[1]
    mu.upper <- lim$mu[2]

    sigma.lower <- lim$sigma[1]
    sigma.upper <- lim$sigma[2]
    
    l.intersects <- (theta - mu.lower)/z  # left-side  intersects
    r.intersects <- (theta - mu.upper)/z  # right-side intersects


    # For each point in l.intersects and r.intersects, 
    # indicate whether it is:
    # a) ABOVE sigma$upper                    (1)
    # b) between sigma$lower and sigma$upper  (0)
    # c) BELOW sigma$lower                   (-1)
  
    tmp.lloc <- as.numeric(l.intersects > sigma.upper) - as.numeric(l.intersects < sigma.lower)
    tmp.rloc <- as.numeric(r.intersects > sigma.upper) - as.numeric(r.intersects < sigma.lower)

    
    list(l.lower = c(-1, tmp.lloc),
         l.upper = c(tmp.lloc, 1),
         r.lower = c(-1, tmp.rloc),
         r.upper = c(tmp.rloc, 1))
  } # end of Location
  
  
  logp.sample <- function(log.w)
  {
    p <- exp(log.w - max(log.w))
    sample(seq(along=p), 1, prob=p)
  } # end of logp.sample
  
  
  region.areas <- function(R, lim, z, theta, location)
  {
    # modif_0.15 This fct's code was in the body of the main fct in earlier versions
    
    mu.lower <- lim$mu[1]
    mu.upper <- lim$mu[2]
      mu.range <- diff(lim$mu)

    sigma.lower <- lim$sigma[1]
    sigma.upper <- lim$sigma[2]
      sigma.range <- diff(lim$sigma)


    region.area <- rep(NA, R)
    dtheta <- diff(theta)
  
    #d <- dtheta/sqrt(1+z^2)  # distance between regions
    v <- dtheta/z            # vertical distance
    
    # Define shape for each region
    region.shape <- 41 - location$r.upper - 3*location$r.lower - 9*location$l.upper - 27*location$l.lower
  
  
    # Compute area for each region

    for (r in 1:R)
    {
      if (region.shape[r] == 4)
      {      
        H <- sigma.upper - (theta[r-1] - mu.upper)/z
        B <- z*H
        region.area[r] <- B*H/2 
      }
      else if (region.shape[r] == 5)
      {   
        b <- mu.upper + z*sigma.upper - theta[c(r, r-1)] # vector of length 2
        h <- b/z
        region.area[r] <- diff(b*h)/2
      }
      else if (region.shape[r] == 7)
      {  
        H <- sigma.range
        B <- mu.upper - theta[r-1] + z*sigma.upper
        b <- B - z*H
        region.area[r] <- (B+b)/2 * H
      }
      else if (region.shape[r] == 8)
      {     
        H <- sigma.range
        B <- dtheta[r-1]
        b.prime <- z*H
        h <- sigma.upper - (theta[r]-mu.upper)/z
        b <- z*h
        region.area[r] <- H*(B+b) - (H*b.prime + b*h)/2 
      }
      else if (region.shape[r] == 9)
      {
        H <- sigma.range
        B <- dtheta[r-1]
        region.area[r] <- B*H
      }
      else if (region.shape[r] == 31)
      {
        B <- mu.range
        h <- sigma.upper - (theta[r-1] - c(mu.lower, mu.upper))/z # vector of length 2
        region.area[r] <- B*mean(h)
      }
      else if (region.shape[r] == 32)
      {         
        B <- mu.range
        H <- v[r-1]
        b <- theta[r] - z*sigma.upper - mu.lower
        h <- b/z
        region.area[r] <- B*H - b*h/2
      }
      else if (region.shape[r] == 34)
      {
        B <- mu.range
        H <- sigma.range
        h <- (theta[r-1] - mu.lower)/z - sigma.lower
        b <- h*z
        region.area[r] <- B*H - b*h/2       
      }
      else if (region.shape[r] == 35)
      {     
        H <- sigma.range
        B <- mu.range
        sigma <- (theta[c(r-1, r)] - c(mu.lower, mu.upper))/z  # vector of length 2
        h <- c(-sigma.lower, sigma.upper) + c(1, -1)*sigma     # vector of length 2
        b <- h*z                                               # vector of length 2
        region.area[r] <- B*H - sum(b*h)/2     
      }
      else if (region.shape[r] == 36)
      {
        b <- theta[r] - z*sigma.upper - mu.lower
        h <- sigma.upper - (theta[r-1]-mu.lower)/z
        H <- sigma.range
        B <- dtheta[r-1]
        region.area[r] <- B*H - h*(B-b)/2
      }
      else if (region.shape[r] == 41)
      {       
        B <- mu.range
        H <- v[r-1]
        region.area[r] <- B*H
      }
      else if (region.shape[r] == 44)
      {  
        B <- mu.range
        H <- v[r-1] 
        b <- mu.upper - theta[r-1] + z*sigma.lower
        h <- b/z
        region.area[r] <- B*H - b*h/2    
      }
      else if (region.shape[r] == 45)
      {       
        h <- (theta[c(r-1,r)] - mu.lower)/z - sigma.lower # vector of length 2 # modif_0.13
        b <- h*z                                          # vector of length 2
        region.area[r] <- diff(b*h)/2
      }
      else if (region.shape[r] == 61)
      {       
        region.area[r] <- mu.range * sigma.range
      }
      else if (region.shape[r] == 62)
      {      
        B <- mu.range
        H <- sigma.range
        h <- sigma.upper - (theta[r] - mu.upper)/z
        b <- z*h
        region.area[r] <- B*H - b*h/2
      }
      else if (region.shape[r] == 63)
      {      
        b <- theta[r] - z*c(sigma.lower, sigma.upper) - mu.lower # vector of length 2
        H <- sigma.range
        region.area[r] <- H*mean(b)
      }
      else if (region.shape[r] == 71)
      {   
        h <- (theta[r] - c(mu.lower, mu.upper))/z - sigma.lower # vector of length 2
        B <- mu.range
        region.area[r] <- mean(h)*B 
      }
      else if (region.shape[r] == 72)
      {     
        H <- (theta[r]-mu.lower)/z - sigma.lower
        B <- z*H
        region.area[r] <- B*H/2
      }
    }
    
    
    region.area
  } # end of region.areas

  
  # new_0.15
  rnorm.cpld <- function(cond.moments, f)
  {
    # Sample one value from N(cond.moments$mean, cond.moments$sd) weighed by a piecewise linear fct f 
    # [with dimensions {m, b, lower, upper}]
    
    lower <- f$lower
    upper <- f$upper
    r     <- f$r
    
    sqrt.2pi <- sqrt(2*pi)
    s <- 1 / 2 / cond.moments$sd^2
    
    # Compute each bands weight
    
      # The piecewise linear weight fct is given by a+b*x, where a=f$b, b=f$m
      a <- f$b
      b <- f$m
    
    # Parameters of the (unweighed) posterior distrn
    mu    <- cond.moments$mean
    sigma <- cond.moments$sd
    
    Phi.upper <- pnorm(upper, mu, sigma)
    Phi.lower <- pnorm(lower, mu, sigma)
      Phi <- Phi.upper - Phi.lower
      
    v2.upper <- exp(-s*(upper-mu)^2)
    v2.lower <- exp(-s*(lower-mu)^2)
    
    riskband.weight <- (a + b*mu) * Phi - b * sigma * (v2.upper - v2.lower) / sqrt.2pi
    
      # Protection against numerical imprecision
      w <- which(riskband.weight < 0)
      if (length(w) > 0)
      {
        tmp <- which(riskband.weight[w] > -1e-10)
        if (length(tmp) > 0) riskband.weight[w[tmp]] <- 0
        
        w <- which(riskband.weight < 0)
        if (length(w) > 0) stop("Negative probability found in riskband.weight = ", riskband.weight, "\n")
      }
    
    
    # Sample an interval
    
    j <- sample(seq(along=riskband.weight), 1, prob=riskband.weight)
    
    
    # Sample a value from the region picked above (by a mixed Newton-Raphson/bisectional algorithm)
    
    vol <- riskband.weight[j]
    lower <- lower[j]
    upper <- upper[j]
    a <- a[j]
    b <- b[j]
    
    A <- a + b*mu
    B <- - b * sigma / sqrt.2pi
    
      Bp <- - B / sigma^2
    
    u <- runif(1)
    target <- u*vol + A*Phi.lower[j] + B*v2.lower[j]
       
    continue <- T               
    x <- (lower + upper) / 2 # starting value
    

    while (continue)
    {
      f <- A*pnorm(x, mu, sigma) + B*exp(-s*(x-mu)^2)
      f.prime <- A * dnorm(x, mu, sigma) + Bp * (x - mu) * exp(-s * (x-mu)^2)
        change <- (f - target) / f.prime
        new.x <- x - change
    
      if      (new.x >= upper) new.x <- (x + upper) / 2
      else if (new.x <= lower) new.x <- (x + lower) / 2
      
      if (f > target) upper <- x
      else lower <- x
      
      continue <- abs(new.x - x) > 1e-6
      x <- new.x
    }                     


    list(region=r[j], mu=x)
  } # end of rnorm.cpld
  
  
  # new_0.15
  rsigma.cpld <- function(N, beta, f)
  {
    # Sample one value from f(mu|sigma) weighed by a piecewise linear fct f [with dimensions {breaks, m, b}]
    # when N > 2                       [use rsigma.cpld.N.le.2 when N <= 2]
    #   
    
    lower <- f$lower
    upper <- f$upper
    
    shape.a <- (N-1) / 2
    shape.b <- (N-2) / 2
  
    
    # Compute each bands weight
    
      # The piecewise linear weight fct is given by a+b*x,  (a=f$b, b=f$m)
      a <- f$b
      b <- f$m
      
    lambda.lower <- 1/upper^2
    lambda.upper <- 1/lower^2
    
    
    pa.upper <- pgamma(lambda.upper, shape=shape.a, rate=beta)
    pa.lower <- pgamma(lambda.lower, shape=shape.a, rate=beta)
    
    pb.upper <- pgamma(lambda.upper, shape=shape.b, rate=beta)
    pb.lower <- pgamma(lambda.lower, shape=shape.b, rate=beta)
    
    B <- sqrt(beta) * gamma((N-2)/2) / gamma((N-1)/2)
    
    vol <- a*(pa.upper-pa.lower) + b*B*(pb.upper-pb.lower)
  
    
      # Protection against numerical imprecision
      w <- which(vol < 0)
      if (length(w) > 0)
      {
        tmp <- which(vol[w] > -1e-10)
        if (length(tmp) > 0) vol[w[tmp]] <- 0
        
        w <- which(vol < 0)
        if (length(w) > 0) stop("Negative probability found in vol = ", vol, "\n")
      }
    
    
    # Randomly pick an interval
    
    j <- sample(seq(along=vol), 1, prob=vol)

    
    # Sample a value from the riskband selected above (through a mixed Bisectional/Newton-Raphson algorithm)
    
    a <- a[j]
    b <- b[j]
    lower <- lower[j]
    upper <- upper[j]
    lambda.upper <- lambda.upper[j]
    
  
    B <- b*B
    u <- runif(1)
    target <- u * vol[j] - a*pa.upper[j] - B*pb.upper[j]
  
    x <- (lower + upper) / 2
    continue <- T
  

    while (continue)
    {
      y <- 1/x^2
      pa.lower <- pgamma(y, shape=shape.a, rate=beta)
      pb.lower <- pgamma(y, shape=shape.b, rate=beta)
  
      f <- -a*pa.lower - B*pb.lower
      f.prime <- (a*dgamma(y, shape=shape.a, rate=beta) + B*dgamma(y, shape=shape.b, rate=beta)) / 2 / y^1.5
        change <- (f - target) / f.prime
        new.x <- x - change
      
      if      (new.x >= upper) new.x <- (x + upper) / 2
      else if (new.x <= lower) new.x <- (x + lower) / 2
        
      if (f > target) upper <- x
      else            lower <- x
        
      continue <- abs(new.x - x) > 1e-6
      x <- new.x
    }
  
    x
  } # end of rsigma.cpld
  
  
  # new_0.15
  rsigma.cpld.N.le.2 <- function(N, beta, f)
  {
    # Sample one value from f(mu|sigma) weighed by a piecewise linear fct f [with dimensions {breaks, m, b}]
    # when N <= 2                       [use rsigma.cpld when N > 2]
    #  
    
    lower <- f$lower
    upper <- f$upper
      I <- length(lower)
    
    shape.a <- (N-1) / 2
    
    
    # Compute each bands weight
    
      # The piecewise linear weight fct is given by a+b*x,  (a=f$b, b=f$m)
      a <- f$b
      b <- f$m
      
    lambda.lower <- 1/upper^2
    lambda.upper <- 1/lower^2
    
  
    pa <- rep(NA, I)
    pb <- pa
    
    ga <- function(t){t^((N-3)/2)*exp(-beta*t)} # m = N
    gb <- function(t){t^((N-4)/2)*exp(-beta*t)} # m = N - 1
  
  
    if (N == 1)
    {
      for (i in 1:I) pa[i] <- 1/2 * integrate(ga, lambda.lower[i], lambda.upper[i])$value
      for (i in 1:I) pb[i] <- 1/2 * integrate(gb, lambda.lower[i], lambda.upper[i])$value
        
      vol <- a*pa + b*pb
    }
    else
    {
      # N == 2
        
      k <- 2 * sqrt(beta/pi)
      pa.upper <- 1/k * pgamma(lambda.upper, shape=shape.a, rate=beta) 
      pa.lower <- 1/k * pgamma(lambda.lower, shape=shape.a, rate=beta)
        
      for (i in 1:I) pb[i] <- 1/2 * integrate(gb, lambda.lower[i], lambda.upper[i])$value
  
      vol <- a*(pa.upper - pa.lower) + b*pb
    }
    
    
    # Randomly pick an interval
    
    j <- sample(seq(along=vol), 1, prob=vol)

    
    # Sample a value from the riskband selected above (through a mixed Bisectional/Newton-Raphson algorithm)
    
    a <- a[j]
    b <- b[j]
    lower <- lower[j]
    upper <- upper[j]
    lambda.upper <- lambda.upper[j]
    
    
    u <- runif(1)
    target <- u * vol[j]
    if (N == 2) target <- target - a*pa.upper[j]
    
    x <- (lower + upper) / 2
    continue <- T
      
      
    if (N == 1)
    {
      while (continue)
      {
        y <- 1/x^2 # lambda.lower
          pa <- 1/2 * integrate(ga, y, lambda.upper)$value
          pb <- 1/2 * integrate(gb, y, lambda.upper)$value
            
        f <- a*pa + b*pb
        f.prime <- (a/2*ga(y) + b/2*gb(y)) / 2 / y^1.5
          change <- (f - target) / f.prime
          new.x <- x - change
        
        if      (new.x >= upper) new.x <- (x + upper) / 2
        else if (new.x <= lower) new.x <- (x + lower) / 2
          
        if (f > target) upper <- x
        else lower <- x
          
        continue <- abs(new.x - x) > 1e-6
        x <- new.x
      } 
    }
    else
    {
      # N == 2
        
      while (continue)
      {
        y <- 1/x^2  # lambda.lower
          pa.lower <- 1/k * pgamma(y, shape=shape.a, rate=beta)
          pb <- 1/2 * integrate(gb, y, lambda.upper)$value
  
        f <- -a*pa.lower + b*pb       
        f.prime <- (a/k * dgamma(y, shape=shape.a, rate=beta) + b/2*gb(y)) / 2 / y^1.5
          change <- (f - target) / f.prime
          new.x <- x - change
        
        if      (new.x > upper) new.x <- (x + upper) / 2
        else if (new.x < lower) new.x <- (x + lower) / 2
          
        if (f > target) upper <- x
        else            lower <- x
          
        continue <- abs(new.x - x) > 1e-6
        x <- new.x
      }  
    }
  
  
    x
  } # end of rsigma.cpld.N.le.2

  
  ##############################################################################
  # Start function

  # new_0.12: added arguments save.RData and RData.dir
  RData <- list(save=save.RData, dir=RData.dir) # new_0.12

  # new_0.12
  if (RData$save)
  {
    mu.save.objects    <- paste(RData$dir, '_mu.RData', sep='/')
    sigma.save.objects <- paste(RData$dir, '_sigma.RData', sep='/')
  }
  else
  {
    mu.save.objects    <- character(0)
    sigma.save.objects <- character(0)
  }


  # A is on the same scale as data
  if (outcome.is.logNormally.distributed) theta <- log(A)
  else theta <- A
  
  # modif_0.15 Embed mu & sigma limits in a list
  lim <- list(mu=c(mu.lower, mu.upper), sigma=c(sigma.lower, sigma.upper))


  z <- qnorm(quantile)

  
  # Make sure that each region delimited by terms in 'A' 
  # do intersect with the domain (mu.lower, mu.upper) x (sigma.lower, sigma.upper)
  
  location <- Location(theta, lim, z) # modif_0.15 values now embedded in a list
  
  
  # Check that no region is empty (does not cross prior domain)
  
  not.empty <- location$l.upper > location$l.lower | location$l.upper == 0 | location$l.lower > location$r.lower
  empty <- !not.empty
    
  if (any(empty))
  {
    msg <- character(0) # define empty string error message
    f.open  <- ifelse(outcome.is.logNormally.distributed, 'log(', '')
    f.close <- ifelse(outcome.is.logNormally.distributed, ')', '')
    
    if (empty[1])
    {
      empty.set <- paste('{(mu,sigma): sigma <= (', f.open, A[1], f.close, '-mu)/', z, '}', sep='')
      msg <- c(msg, empty.set)
    }
    
    
    which.empty <- which(empty)
    which.empty <- setdiff(which.empty, c(1, length(empty))) # Drop first and last regions
    if (length(which.empty) > 0)
    {
      empty.sets <- paste('{(mu,sigma): ', f.open, A[which.empty-1], f.close, '-mu)/', z, 
        ' < sigma <= (', f.open, A[which.empty], f.close, '-mu)/', z, '}', sep='')
      msg <- c(msg, empty.sets)
    }
    

    if (empty[length(empty)])
    {
      empty.set <- paste('{(mu,sigma): sigma > (', f.open, A[length(A)], f.close, '-mu)/', z, '}', sep='')
      msg <- c(msg, empty.set)
    }
    
    msg <- c('The region(s) below do not intersect with domain defined by (mu.lower, mu.upper) X (sigma.lower, sigma.upper):', msg)
    stop(msg)
  }
  
  
  # Make sure that the number of regions and region.prior.prob length correspond
  
  R <- length(A) + 1
  if (R != length(region.prior.prob)) 
  {
    msg <- paste("Number of region prior probabilities found in region.prior.prob should be ", R, '[length(A)+1]')
    stop(msg)
  }
  
  
  # Make sure that the values in region.prior.prob make sense
  
  undefined.region.prior.prob <- all(is.na(region.prior.prob))
  
  
  if (use.continuous.piecewise.linear.prior)
  {
    # new_0.15
    
    if (undefined.region.prior.prob)
    {
      if (!equally.probable.regions) stop('region.prior.prob must be provided or equally.probable.regions set to TRUE when use.continuous.piecewise.linear.prior is TRUE.')
      region.prior.prob <- rep(1/R, R)
    }
    
    if (any(is.na(region.prior.prob)) | any(region.prior.prob < 0))
    {
      stop('Values in region.prior.prob should all be non-negative numbers.')
    }
    else
    {
      tmp <- sum(region.prior.prob)
      if (tmp != 1) stop('Values in region.prior.prob must sum to 1.')
    }
    

    CPL <- CPL.prior(region.prior.prob, lim, z, A, outcome.is.logNormally.distributed)
  }
  else
  {
    # modif_0.15 embedded in the else-block
  
    uniform.prior.on.mu.and.sigma <- F # new_0.13
    
    region.area <- region.areas(R, lim, z, theta, location) # new_0.15

  
    if (undefined.region.prior.prob)
    { 
      uniform.prior.on.mu.and.sigma <- !equally.probable.regions # new_0.13
      
      # modif_0.13
      if (uniform.prior.on.mu.and.sigma) region.prior.prob <- region.area/sum(region.area)
      else                               region.prior.prob <- rep(1/R, R)
    }
    else if (any(is.na(region.prior.prob)) | any(region.prior.prob < 0))
    {
      stop('Values in region.prior.prob should all be non-negative numbers.')
    }
    else
    {
      tmp <- sum(region.prior.prob)
      if (tmp != 1) stop('Values in region.prior.prob must sum to 1.')
    }
    
    
    # Compute density in each region
                               
    if (uniform.prior.on.mu.and.sigma)
    {
      # modif_0.13 The condition is now: if uniform.prior.on.mu.and.sigma
      # modif_0.13 The code in this block is also changed
      
      tmp <- 1 / ((mu.upper-mu.lower)*(sigma.upper-sigma.lower))
      region.prior.density <- rep(tmp, R)
      region.prior.logdens <- rep(0, R) # new_0.13 Not important to correspond to be = log(tmp)
    }
    else
    {
      region.prior.density <- as.numeric(region.prior.prob > 0) # 0/1 variable indicating non-null region probabilities
      
      regions <- which(region.prior.density > 0)
      r.ref <- regions[1]    # first region with a non-null prior probability
      regions <- regions[-1] # drop r.ref from list 
    
      for (r in regions)
      {
        region.prior.density[r] <- region.prior.prob[r] * region.area[r.ref] / (region.prior.prob[r.ref] * region.area[r])
      }
      
      
      region.prior.density <- region.prior.density/sum(region.prior.density) # Standardize values in region.prior.density
      
      region.prior.logdens <- log(region.prior.density) # new_0.13
    }
  }
  
  
  ######################################################################################################################
  # The rest of the function is very similar to what was done in previous versions
  
  me <- any.me(me.sd.range, cv.range) # see if measurement error is desired [new_0.10]
  
  data <- data.summary(y=y, lt=lt, gt=gt, interval.lower=interval.lower, interval.upper=interval.upper,
                       logNormal.distrn=outcome.is.logNormally.distributed, me.through.cv=me$through.cv)
                          
  gen.y <- list()
  true.values <- list()
  
  default.inits <- Default.inits(data, outcome.is.logNormally.distributed, c(mu.lower, mu.upper), c(sigma.lower, sigma.upper)) # new_0.11
  
  tau.alpha <- (data$N - 1)/2
  
  # Make sure that initial values for mu are in allowed range
  
  mu.range <- mu.upper - mu.lower
  
  w <- which(init.mu > mu.upper)
  if (length(w) > 0) init.mu[w] <- mu.upper - mu.range/20
  
  w <- which(init.mu < mu.lower)
  if (length(w) > 0) init.mu[w] <- mu.lower + mu.range/20
  
  
  # Prepare dens.gen.icdf objects
  
  
  if (me$any) o.tv <- truevalue.gen.object(me, outcome.is.logNormally.distributed)
      
  
  
  if (!use.continuous.piecewise.linear.prior)
  {
    # modif_0.15 This block was made conditional
    
    if (me$through.cv & !outcome.is.logNormally.distributed)
    {
      o.mu <- mu.truncatedData.gen.object(data$N)
      o.sigma <- sigma.truncatedData.gen.object(data$N)
    }
    else if (data$N <= 1)
    {            
      o.sigma <- sigma.gen.object(data$N)
    }
    else
    {
      o.sigma <- list()
    }
  }
  else if (me$any)
  {
    stop('The option use.continuous.piecewise.linear.prior is not available (yet) in presence of measurement error.\n')
    # ICI voir ce qui se passe lorsque use.continuous.piecewise.linear.prior = TRUE
    # en presence d'erreur de mesure
  }
  
  
  if (me$any & !me$known) o.me <- me.gen.object(me, outcome.is.logNormally.distributed, data$N)

  
  # Prepare output objects with appropriate dimensions
  
  region.sample <- matrix(NA, nrow=n.chains, ncol=n.iter)

  sample <- empty.matrices(n.iter, n.chains, me)

  burnin <- list()
  if (monitor.burnin) burnin <- empty.matrices(n.burnin, n.chains, me)

  
  
  M.iter <- n.burnin + n.iter * n.thin
  
  for (ch in seq(n.chains))
  {
    saved.iter <- 0
    
    # Set initial values
    
    mu <- init.mu[ch]

    if (is.na(init.sigma[ch]) | init.sigma[ch] == 0) sigma <- default.inits$sigma
    else sigma <- init.sigma[ch]
    
    # Initialize measured values for subjects with censored values [modif_0.10]
    if (data$any.censored$any) gen.y <- y.gen.inits(data, mu, sigma, me$through.cv, outcome.is.logNormally.distributed)
    
    if (me$any) me$parm <- me$init # initialize Measurement Error parameter value
  
    
    for (iter in 1:M.iter)
    {
      # Sample true values (in presence of measurement error) [new_0.10]
      if (me$any) true.values <- truevalues.gen(gen.y, data, mu, sigma, me, outcome.is.logNormally.distributed, o=o.tv, RData=RData)
                       
      # Sample y values for censored observations
      if (data$any.censored$any) gen.y <- y.gen(true.values, data, sigma, me, outcome.is.logNormally.distributed, mu=mu)
      
      
      # Compute data points sum and square sum
      moments <- out.logout.moments(me$any, outcome.is.logNormally.distributed, data, gen.y, true.values)
            
      
      # Sample from f(sigma | mu)
      
      tau.beta <- (moments$sum2 - 2*mu*moments$sum + data$N*mu^2) / 2
      
      if (use.continuous.piecewise.linear.prior)
      {
         # new_0.15
         f <- cpld.fsigma.given.mu(mu, CPL, lim, z)
         
         if (data$N > 2) sigma <- rsigma.cpld(data$N, tau.beta, f)
         else            sigma <- rsigma.cpld.N.le.2(data$N, tau.beta, f)
      }
      else
      {
        # modif_0.15 block below now embedded in the else-bock
        # modif_0.10
        # modif_0.13 The condition in the 'if' below is now: uniform.prior.on.mu.and.sigma
        
        if (uniform.prior.on.mu.and.sigma)
        {
          # The joint prior for (mu, sigma) is uniform over the whole rectangle domain, hence no region delimiting is necessary
          
          sigma.range <- c(sigma.lower, sigma.upper)
          
          if (me$through.cv & !outcome.is.logNormally.distributed)
          {
            sigma <- sigma.truncatedData.gen(o.sigma, sigma.range, tau.beta, mu, sigma, RData=RData)
          }
          else
          {
            sigma <- sqrt.invertedGamma.gen(data$N, tau.beta, sigma.range, o=o.sigma, RData=RData)
          }
        }
        else
        {   
          sigma.cutoffs <- (theta - mu)/z
          w <- which(sigma.cutoffs >= sigma.lower & sigma.cutoffs <= sigma.upper)
          sigma.cutoffs <- unique(c(sigma.lower, sigma.cutoffs[w], sigma.upper))
          
          multiple.regions <- length(sigma.cutoffs) > 2
          if (multiple.regions) 
          {
            possible.regions <- c(w, max(w) + 1) # regions in which (mu, sigma) can fall
            log.w <- region.prior.logdens[possible.regions]
          }
          
          
          if (me$through.cv & !outcome.is.logNormally.distributed)
          {
            A <- c(o.sigma$A, list(b=tau.beta, mu=mu))
  
            if (multiple.regions)
            {
              r <- seq(along=possible.regions)
              region.wt <- rep(NA, length(r))          
              for (j in r) region.wt[j] <- integrate(o.sigma$f, sigma.cutoffs[j], sigma.cutoffs[j+1], A=A)$value
              region.prob <- region.wt/sum(region.wt) * exp(region.prior.logdens[possible.regions])
              j <- sample(r, 1, prob=region.prob)
            }
            else
            {
              j <- 1
            }
  
            sigma <- dens.gen.icdf(o.sigma, A, range=sigma.cutoffs[j+c(0,1)], save.objects=sigma.save.objects)
          }
          else if (data$N <= 1)
          {            
            A <- c(o.sigma$A, list(b=tau.beta))          
              
            if (multiple.regions)
            {
              region.wt <- rep(NA, length(possible.regions)) 
              for (i in seq(along=region.wt)) region.wt[i] <- integrate(o.sigma$f, sigma.cutoffs[i], sigma.cutoffs[i+1], A=A)$value
   
              log.w <- log(region.wt) + log.w
              j <- logp.sample(log.w) 
            }
            else
            {
              j <- 1
            }
              
            range <- sigma.cutoffs[j+c(0,1)]
            start <- mean(range)
            o.sigma.start <- o.sigma$start(A)
            if (o.sigma.start >= range[1] & o.sigma.start <= range[2]) start <- c(start, o.sigma.start)
              
            sigma <- dens.gen.icdf(o.sigma, A, range=range, start=start, save.objects=sigma.save.objects)
          }
          else
          { 
            # size > 1
  
            if (tau.beta == 0)
            {
              if (multiple.regions)
              {
                sigma.logp <- diff(sigma.cutoffs^(data$N-1))/(data$N-1) + log.w
                j <- logp.sample(sigma.logp)    
                sigma.cutoffs <- sigma.cutoffs[j+c(0,1)] # limits of the interval for sigma in which we want to sample
              }
                
              sigma <- random.pow(data$N, sigma.cutoffs)
            }
            else
            {
              log.pcum <- invgamma.logpcum(sigma.cutoffs^2, tau.alpha, tau.beta)
                
              if (multiple.regions)
              {
                sigma.logp <- logp.from.logpcum(log.pcum) + log.w
                j <- logp.sample(sigma.logp)
              }
              else
              {
                j <- 1
              }
          
              sigma.logpcum <- log.pcum$log.pcum[j+c(0,1)]
              log.U <- runif.logp(sigma.logpcum, lower.tail=log.pcum$lower.tail)
              tau <- qgamma(log.U, tau.alpha, tau.beta, log.p=T, lower.tail=!log.pcum$lower.tail)
              sigma <- 1/sqrt(tau)
            }
          }
        }
      }
            

      # Sample from f(mu | sigma)
      
      mu.cond <- list(mean = moments$sum / data$N,
                      sd = sigma/sqrt(data$N))
      
      if (use.continuous.piecewise.linear.prior)
      {
        # new_0.15
        
        f <- cpld.fmu.given.sigma(sigma, CPL, lim, z)
        
        tmp <- rnorm.cpld(mu.cond, f)
          mu <- tmp$mu
          region <- tmp$region
      }
      else
      {
        # modif_0.15 This block was embedded in the else-block
        # modif_0.10 [block]
        # modif_0.13 The condition is now: if uniform.prior.on.mu.and.sigma
        
        if (uniform.prior.on.mu.and.sigma)
        {
          # The joint prior for (mu, sigma) is uniform over the whole rectangle domain, hence no region delimiting is necessary
          
          mu.range <- c(mu.lower, mu.upper)
          
          if (me$through.cv & !outcome.is.logNormally.distributed)
          {
            mu <- mu.truncatedData.gen(o.mu, mu.range, mu.cond$mean, sigma=sigma, RData=RData)
          }
          else
          {          
            phi <- pnorm.logpcum(mu.range, mean=mu.cond$mean, sd=mu.cond$sd)
            log.U <- runif.logp(phi$log.pcum, lower.tail=phi$lower.tail)
            mu <- qnorm(log.U, mean=mu.cond$mean, sd=mu.cond$sd, log.p=T, lower.tail=phi$lower.tail)
          }
          
          # Define region in which (mu, sigma) is sitting
          tmp <- mu + z*sigma
          region <- 1 + sum(tmp > theta)
        }
        else
        { 
          mu.cutoffs <- theta - z*sigma
          w <- which(mu.cutoffs >= mu.lower & mu.cutoffs <= mu.upper)
          mu.cutoffs <- unique(c(mu.lower, mu.cutoffs[w], mu.upper))
          
          if (length(w) > 0) region <- c(w, max(w) + 1) # regions in which (mu, sigma) can fall
          else region <- 1
          
          if (me$through.cv & !outcome.is.logNormally.distributed)
          {        
            A <- c(o.mu$A, list(mu.mean=mu.cond$mean, s=sigma, s2=sigma^2))
          
            if (length(region) > 1)
            {
              r <- seq(along=region)
              region.wt <- rep(NA, length(region))  
              for (j in r) region.wt[j] <- integrate(o.mu$f, mu.cutoffs[j], mu.cutoffs[j+1], A=A)$value
              region.prob <- region.wt/sum(region.wt) * exp(region.prior.logdens[region])
              j <- sample(r, 1, prob=region.prob)
            }
            else
            {
              j <- 1
            }
            
            mu <- dens.gen.icdf(o.mu, A, range=mu.cutoffs[j+c(0,1)], save.objects=mu.save.objects)
          }
          else
          {
            phi <- pnorm.logpcum(mu.cutoffs, mean=mu.cond$mean, sd=mu.cond$sd)
            
            if (length(region) > 1)
            {
              logp <- logp.from.logpcum(phi) + region.prior.logdens[region]
              j <- logp.sample(logp)
            }
            else
            {
              j <- 1
            }
            
            log.pcum <- phi$log.pcum[j+c(0,1)] # vector of length 2
            log.U <- runif.logp(log.pcum, lower.tail=phi$lower.tail)
            mu <- qnorm(log.U, mean=mu.cond$mean, sd=mu.cond$sd, log.p=T, lower.tail=phi$lower.tail)
          }
          
          region <- region[j] # define region in which (mu, sigma) is sitting
        }
      }
    
    
      # Sample Measurement Error from its posterior density
      
      if (me$any & !me$known) me$parm <- me.gen(o.me, me, data, gen.y, true.values, RData=RData)
           
   
      # Save values only when iter# modulo thinning = 0

      if (iter <= n.burnin)
      {
        if (monitor.burnin)
        {
          burnin$mu[ch, iter] <- mu
          burnin$sd[ch, iter] <- sigma
          
          if (me$any & !me$known) burnin$me.parm[ch, iter] <- me$parm
        }
      }
      else if ((iter-n.burnin)%%n.thin == 0)
      {
        saved.iter <- saved.iter + 1
        sample$mu[ch, saved.iter] <- mu
        sample$sd[ch, saved.iter] <- sigma
        region.sample[ch, saved.iter] <- region
        
        if (me$any & !me$known) sample$me.parm[ch, saved.iter] <- me$parm
      }
    }
  }
  
  
  # Estimate region posterior probability (from the frequency table of 'region' variable)
  
  region.posterior.prob <- matrix(0, nrow=n.chains, ncol=R)
  
  for (ch in n.chains)
  {
    for (r in 1:R)
    {
      tmp <- region.sample == r
      region.posterior.prob[ch, r] <- mean(tmp)
    }
  }
  
  region.posterior.prob <- drop(region.posterior.prob) # turn the matrix into a vector if n.chains == 1


  # Collect information to return in output  
  
  # modif_0.15
  out <- list(quantile=quantile, z.quantile=z,
              region.prior.prob=region.prior.prob, R=R,
              use.continuous.piecewise.linear.prior=use.continuous.piecewise.linear.prior)
              
  if (use.continuous.piecewise.linear.prior)
  {
    out$geo <- list(u=CPL$geometry$u,
                    area=CPL$geometry$area, 
                    x95b=CPL$geometry$x95b, 
                    cos.theta=CPL$geometry$cos.theta,
                    R=CPL$geometry$R)
    
    CPL$geometry <- NULL
    out$CPL.prior <- CPL
  }
  else
  {
    # modif_0.15 Now added to output conditional on use.continuous.piecewise.linear.prior = FALSE
    out$region.area <- region.area 
    out$region.prior.density <- region.prior.density
  }
        
              
  tmp <- out.sample(sample, burnin, n.chains, monitor.burnin, me, outcome.is.logNormally.distributed)
  out$sample <- tmp$sample
  if (monitor.burnin) out$burnin <- tmp$burnin
    
  out$mcmc  <- list(n.chains=n.chains, n.burnin=n.burnin, n.iter=n.iter, n.thin=n.thin, monitor.burnin=monitor.burnin)
  out$inits <- list(mu=init.mu, sigma=init.sigma)
  out$data  <- list(y=y, lt=lt, gt=gt, interval.lower=interval.lower, interval.upper=interval.upper)
  out$parms <- list(A=A, mu.lower=mu.lower, mu.upper=mu.upper, sigma.lower=sigma.lower, sigma.upper=sigma.upper,
                    outcome.is.logNormally.distributed=outcome.is.logNormally.distributed)
                    
  out$region.posterior.prob <- region.posterior.prob
                   
  # new_0.14, modif_0.16
  out$sample$x95 <- out$sample$mu + z * out$sample$sd
  

  out
} # end of SEG.riskband

         
SEG.riskband.plot.CPL.prior <- function(o, color='purple')
{
  # o: an output of SEG.riskband
  
  if (!o$use.continuous.piecewise.linear.prior) stop("This output/object was not obtained with use.continuous.piecewise.linear.prior=TRUE\n")
  
  g <- o$geo

  h <- o$CPL.prior$h
  s <- o$CPL.prior$s * g$cos.theta
  x95 <- o$CPL.prior$x95
  R <- g$R
  
  
  # Draw the piecewise uniform prior distrn
  
  u <- g$u
  bounds <- g$x95b
    
  x <- c(bounds[1], rep(bounds[-c(1, R+1)], rep(2, R-1)), bounds[R+1])
  f <- rep(u, rep(2, R))
  ylim <- c(0, max(c(f,h)))
  plot(x, f, type='l', col='gray', ylim=ylim, xlab=expression(x[95]), ylab=expression(f(x[95])))
    title('Prior distribution')
    str <- paste(o$region.prior.prob, collapse=', ')
    str <- paste(c('region prior prob = c(', str, ')'), collapse='')
    mtext(str, side=3, at=par('usr')[2], adj=1, line=0.5)


  # Overlay the curve as indicated by (h, s)

  h.len <- length(h)
  interval.len <- diff(x95)

  h1 <- h[-h.len] + s * interval.len
  
  i.seq <- seq(along=h)
  

  for (i in seq(along=i.seq)[-h.len])
  {
    j <- i.seq[i]
    next.j <- i.seq[i+1]
    points(x95[c(j,next.j)], c(h[j], h1[j]), type='l', col=color)
  }


  # Plot a legend
  
    legend.lty <- c(1,1)
    legend.txt <- c('continuous piecewise linear [CPL]', 'piecewise uniform [for reference only]') 
    legend.size <- legend('bottom', legend.txt, lty=legend.lty, bty='n', plot=F)$rect # do not plot -- just note legend size
  
    # Choose the best location for the legend -- pick a spot from which the graphic points are as far as possible
    
    par.usr <- par('usr')
    par.pin <- par('pin')
    
    # Points per inches, on both dimensions
    x.ppi <- diff(par.usr[c(1,2)]) / par.pin[1]
    y.ppi <- diff(par.usr[c(3,4)]) / par.pin[2]
    
    x <- c(par.usr[1] + legend.size$w/2, mean(par.usr[c(1,2)]), par.usr[2] - legend.size$w/2)
    x <- rep(x, 2)
    y <- c(par.usr[3] + legend.size$h, par.usr[4] - legend.size$h)
    y <- rep(y, rep(3,2))
    ref <- data.frame(x=x, y=y, loc=c('bottomleft', 'bottom', 'bottomright', 'topleft', 'top', 'topright'))
    n.ref <- nrow(ref)
    d <- rep(Inf, n.ref)
    d2pt <- rep(NA, 2)
    
    for (i in seq(along=s))
    {
      x <- c(x95[i], x95[i+1])
      y <- c(h[i], h[i] + s[i]*diff(x))
      
      m <- -1/s[i] * (y.ppi/x.ppi)^2 # visual perpendicular slope
      a <- h[i] - s[i] * x[1] # intercept
      
      
      for (j in 1:n.ref)
      {
        for (k in 1:2) d2pt[k] <- sqrt((x[k]-ref$x[j])^2/x.ppi^2 + (y[k]-ref$y[j])^2/y.ppi^2)
        
        a.perpendicular <- ref$y[j] - m * ref$x[j] # intercept of the (visually) perpendicular line
        
        # coordinates of the point on the line that is closest to point ref[j]
  
        line.x <- ifelse(s[i]==0, ref$x[j], (a.perpendicular - a) / (s[i] - m))
        
        if (x[1] < line.x & line.x < x[2])
        {
          line.y <- a + s[i] * line.x
          this.d <- sqrt((line.x-ref$x[j])^2/x.ppi^2 + (line.y-ref$y[j])^2/y.ppi^2)
        }
        else this.d <- min(d2pt)

        d[j] <- min(d[j], this.d)
      }
    }
  
    w <- which.max(d)
    loc <- ref$loc[w]
  
  
  legend(loc, legend.txt, col=c(color, 'gray'), lty=legend.lty, bty='n')
} # end of SEG.riskband.plot.CPL.prior


SEG.riskband.posterior.plot <- function(o, z=o$z.quantile)
{
  # o: an output from SEG.riskband
  
  plot(o$sample$mu, o$sample$sd, type='p', pch='.', col='navy', 
       xlab=expression(mu), ylab=expression(sigma))
  title(expression(paste("(", mu, ", ", sigma, "): scatter plot of a sample obtained from their joint posterior distribution"), sep=''))

  color <- 'darkgoldenrod'
  hist(o$sample$x95, nclass=50, xlab='P95', prob=T,
       main='P95: histogram of values from a sample obtained\nfrom its posterior distribution')
       
    abline(v=o$geo$x95b, col=color, lty=2)
    y.txt <- par('usr')[4] - par('cxy')[2]/2
    plot.xlim <- par('usr')[c(1,2)]
    Cat.lbl <- c('Cat 0', 'Cat 1', 'Cat 2', 'Cat 3', 'Cat 4')
    rb.midpoint <- (o$geo$x95b[-1] + o$geo$x95b[-6]) / 2 
     
    for (i in 1:5)
    {
      if (o$geo$x95b[i] >= plot.xlim[1] & o$geo$x95b[i+1] <= plot.xlim[2])  text(rb.midpoint[i],                     y.txt, Cat.lbl[i], adj=0.5, col=color)
      else if (o$geo$x95b[i]   >= plot.xlim[1])                             text((o$geo$x95b[i] + plot.xlim[2])/2,   y.txt, Cat.lbl[i], adj=0.5, col=color)
      else if (o$geo$x95b[i+1] <= plot.xlim[2])                             text((o$geo$x95b[i+1] + plot.xlim[1])/2, y.txt, Cat.lbl[i], adj=0.5, col=color)
    }  
    
    legend.lbls <- paste(paste(c('Cat 0', 'Cat 1', 'Cat 2', 'Cat 3', 'Cat 4'), 
                         round(100*o$region.posterior.prob, digits=1), sep=': '), '%', sep='')
    legend('right', legend=legend.lbls, text.col=color, xjust=0.5, title='Posterior probabilities')
  
  
  hist(o$sample$mu, nclass=50, xlab=expression(mu), prob=T, 
       main=expression(paste(mu, ': histogram of values obtained from (', mu, ', ', sigma, ') joint posterior distribution')))
  
  
  hist(o$sample$sd, nclass=50, xlab=expression(sigma), prob=T, 
       main=expression(paste(sigma, ': histogram of values obtained from (', mu, ', ', sigma, ') joint posterior distribution')))
} # end of SEG.riskband.posterior.plot
