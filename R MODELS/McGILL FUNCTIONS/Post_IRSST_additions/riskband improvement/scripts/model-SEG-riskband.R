source("c:/users/patri/home/bin/R/fcts/dens.gen.icdf/dens-gen-icdf.R")
source("c:/users/patri/home/consultation/L/Lavoue/webexpo/R/fcts.R")
source("c:/users/patri/home/consultation/L/Lavoue/webexpo/riskband/R/CPL_prior.R")

source("c:/users/patri/home/bin/R/fcts/bisectionalSearch.R")
                                                                                    
                                                                                         
#  Version 0.26 (Feb 2022)
                                                                                                                              
                                                                    
# Change Log *******************************************************************          
#
#
# Version 0.26 (Feb 2022)
# -----------------------
#   The function used to calculate a Continuous-Piecewise Linear (CPL)
#   prior distribution --- namely CPL.prior --- is now a stand-alone function;
#   its code can be found in file CPL_prior.R
#
#   The function SEG.riskband.plot.CPL.prior was also moved to that file
#   AND renamed CPL.prior.plot
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
  riskband.prior.prob=rep(NA, R),
  gm.min=exp(mu.lower), gm.max=exp(mu.upper),
    mu.lower=-3, mu.upper=6.2,
  gsd.min=1.05, gsd.max=4, 
    sigma.lower=log(gsd.min), sigma.upper=log(gsd.max), 
  init.mu=rep(default.inits$mu, n.chains),
  init.sigma=rep(default.inits$sigma, n.chains),
  outcome.is.logNormally.distributed=TRUE,
  equally.probable.regions=TRUE,  # ignored when region.prior.prob is provided
  use.continuous.piecewise.linear.prior=length(cpl.prior)>0,
  cpl.prior=list(),
  quantile=0.95,
  me.sd.range=numeric(0), cv.range=numeric(0),
  save.RData=FALSE, RData.dir='c:/users/jerome')
{
  # If you want to use the Continuous-Piecewise Linear prior distribution,
  # you can either: 
  #    a) specify  use.continuous.piecewise.linear.prior = TRUE and let the code calculate the CPL prior distribution 
  # or b) provide the CPL prior distribution through the argument cpl.prior=my.CPL.prior 
  #       (assuming you have calculated the CPL prior beforehand and saved the result under the name my.CPL.prior)
  #       IMPORTANT/WARNING: if you provide a value for cpl.prior (the list returned by CPL.prior), 
  #                          the values (region.prior.prob,
  #                                      gm.min, gm.max, mu.lower, mu.upper, 
  #                                      gsd.min, gsd.max, sigma.lower, sigma.upper
  #                                      quantile, A)
  #                          given to SEG.riskband's arguments will be ignored.

  # IMPORTANT:
  # Selon la valeur de outcome.is.logNormally.distributed (TRUE/FALSE)
  # nous devrions demander à l'usager 
  #  i) seulement (gm.min, gm.max, gsd.min, gsd.max)               si outcome.is.logNormally.distributed = TRUE
  # ii) seulement (mu.lower, mu.upper, sigma.lower, sigma.upper)   sinon


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
  
  
  if (length(cpl.prior) > 0)
  {
    region.prior.prob <- cpl.prior$geometry$riskband.prior.prob
    
    lim <- list(mu   =cpl.prior$geometry$mu.lim,
                sigma=cpl.prior$geometry$sigma.lim)

      mu.lower <- lim$mu[1]
      mu.upper <- lim$mu[2]
    
    A <- cpl.prior$geometry$A
    z <- cpl.prior$geometry$z
  }
  else
  {
    lim <- list(mu=c(mu.lower, mu.upper), sigma=c(sigma.lower, sigma.upper))
    
    z <- qnorm(quantile)
    
    # A is on the same scale as data
    if (outcome.is.logNormally.distributed) theta <- log(A)
    else theta <- A
  }



  

  
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
    

    if (length(cpl.prior) > 0)  CPL <- cpl.prior
    else                        CPL <- CPL.prior(region.prior.prob, lim, A, outcome.is.logNormally.distributed, z)
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
