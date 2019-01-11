# Version 0.13 (September 2018)

## update J lavoué sept 14 : élimination of gsd and gm (confusion)


  # Change Log *****************************************************************
  #
  #   When updated, look for comments with new_* and modif_*
  #   to rapidly identify new/modified code.
  #
  # Version 0.13
  # ------------
  #   Previous versions were incorrect in the 
  #   conditional posterior distributions f(mu|sigma) and f(sigma|mu)
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
  A=c(0.1,0.5,1,5),
  region.prior.prob=rep(NA, n.regions),
  mu.lower= log(0.001), mu.upper=log(100),
  sigma.lower=log(1.05), sigma.upper=log(10), 
  init.mu=rep(default.inits$mu, n.chains),
  init.sigma=rep(default.inits$sigma, n.chains),
  outcome.is.logNormally.distributed=TRUE,
  quantile=0.95,
  me.sd.range=numeric(0), cv.range=numeric(0),
  save.RData=F, RData.dir='c:/users/jerome')
{
  # IMPORTANT:
  # mu and sigma are the log transformed gm and gsds if is.lognormal=TRUE


  # Useful function(s)
  
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
  
  
  logp.sample <- function(log.w)
  {
    p <- exp(log.w - max(log.w))
    sample(seq(along=p), 1, prob=p)
  } # end of logp.sample
  

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




  if (outcome.is.logNormally.distributed)
  {
    theta <- log(A)
  }
  else
  {
    theta <- A
  }
  
  z <- qnorm(quantile)
  
  
  # Make sure that each region delimited by terms in 'A' 
  # do intersect with the domain (mu.lower, mu.upper) x (sigma.lower, sigma.upper)
  
  l.intersects <- (theta - mu.lower)/z  # left-side  intersects
  r.intersects <- (theta - mu.upper)/z  # right-side intersects


  # For each point in l.intersects and r.intersects, 
  # indicate whether it is:
  # a) ABOVE sigma.upper                    (1)
  # b) between sigma.lower and sigma.upper  (0)
  # c) BELOW sigma.lower                   (-1)
  
  tmp.lloc <- as.numeric(l.intersects > sigma.upper) - as.numeric(l.intersects < sigma.lower)
  tmp.rloc <- as.numeric(r.intersects > sigma.upper) - as.numeric(r.intersects < sigma.lower)

  
  location.l.lower <- c(-1, tmp.lloc)
  location.l.upper <- c(tmp.lloc, 1)
  location.r.lower <- c(-1, tmp.rloc)
  location.r.upper <- c(tmp.rloc, 1)
  
  
  # Check that no region is empty (does not cross prior domain)
  
  not.empty <- location.l.upper > location.l.lower | location.l.upper == 0 | location.l.lower > location.r.lower
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
  
  n.regions <- length(A) + 1
  if (n.regions != length(region.prior.prob)) 
  {
    msg <- paste("Number of region prior probabilities found in region.prior.prob should be ", n.regions, '[length(A)+1]')
    stop(msg)
  }
  
  
  # Define shape for each region
  
  region.shape <- 41 - location.r.upper - 3*location.r.lower - 9*location.l.upper - 27*location.l.lower

  
  # Compute area for each region
  
  region.area <- rep(NA, n.regions)
  dtheta <- diff(theta)
  
  d <- dtheta/sqrt(1+z^2)  # distance between regions
  v <- dtheta/z            # vertical distance
  
  mu.range <- mu.upper - mu.lower
  sigma.range <- sigma.upper - sigma.lower
  
  for (r in 1:n.regions)
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
      h <- (theta[c(r-1,r)] - mu.lower)/z # vector of length 2
      b <- h*z                            # vector of length 2
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
  
  
  # Make sure that values in region.prior.prob make sense
  
  undefined.region.prior.prob <- all(is.na(region.prior.prob))
  
  if (undefined.region.prior.prob)
  { 
    region.prior.prob <- region.area/sum(region.area)
  }
  else if (any(is.na(region.prior.prob)) | any(region.prior.prob < 0))
  {
    stop('Values in region.prior.prob should all be non-negative numbers.')
  }
  else
  {
    tmp <- sum(region.prior.prob)
    if (tmp != 1) stop('Values in region.prior.prob must sum to 1.')
    
    ## region.prior.logprob <- log(region.prior.prob) # modif_0.13  Dropped term (not used anymore)
  }
  
  
  # Compute density in each region
  
  # modif_0.13
  #if (undefined.region.prior.prob)
  #{
  #  region.prior.density <- rep(1/n.regions, n.regions)
  #}
  #else
  #{
  
  # modif_0.13 The whole block of code below is not conditional anymore (it was in an else-block before)
     
  region.prior.density <- as.numeric(region.prior.prob > 0) # 0/1 variable indicating non-null region probabilities
  
  regions <- which(region.prior.density > 0)
  r.ref <- regions[1]    # first region with a non-null prior probability
  regions <- regions[-1] # drop r.ref from list 
  
  for (r in regions)
  {
    region.prior.density[r] <- region.prior.prob[r] * region.area[r.ref] / (region.prior.prob[r.ref] * region.area[r])
  }
    
  # modif_0.13 -> * region.area appears in the sum() in the denominator in line below
  region.prior.density <- region.prior.density / sum(region.prior.density * region.area) # Standardize values in region.prior.density
  region.prior.log.density <- log(region.prior.density) # new_0.13
    
  # modif_0.13  End of block above, now made non-conditional
  #}
  
  
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
  
  w <- which(init.mu > mu.upper)
  if (length(w) > 0) init.mu[w] <- mu.upper - mu.range/20
  
  w <- which(init.mu < mu.lower)
  if (length(w) > 0) init.mu[w] <- mu.lower + mu.range/20
  
  
  # Prepare dens.gen.icdf objects
  
  
  if (me$any) o.tv <- truevalue.gen.object(me, outcome.is.logNormally.distributed)
      
  
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
      
      # modif_0.10

      # modif_0.13 Block commented out
      
      #if (undefined.region.prior.prob)
      #{
      #  # The joint prior for (mu, sigma) is uniform over the whole rectangle domain, hence no region delimiting is necessary
      #  
      #  sigma.range <- c(sigma.lower, sigma.upper)
      #  
      #  if (me$through.cv & !outcome.is.logNormally.distributed)
      #  {
      #    sigma <- sigma.truncatedData.gen(o.sigma, sigma.range, tau.beta, mu, sigma, RData=RData)
      #  }
      #  else
      #  {
      #    sigma <- sqrt.invertedGamma.gen(data$N, tau.beta, sigma.range, o=o.sigma, RData=RData)
      #  }
      #}
      #else
      #{
      
      # modif_0.13 The whole block of code below is not conditional anymore (it was in an else-block before)
         
      sigma.cutoffs <- (theta - mu)/z
      w <- which(sigma.cutoffs >= sigma.lower & sigma.cutoffs <= sigma.upper)
      sigma.cutoffs <- unique(c(sigma.lower, sigma.cutoffs[w], sigma.upper))

        
      multiple.regions <- length(sigma.cutoffs) > 2
      if (multiple.regions) 
      {
        possible.regions <- c(w, max(w) + 1) # regions in which (mu, sigma) can fall
        log.w <- region.prior.log.density[possible.regions]
          # modif_0.13 Line above is now using region.prior.log.density
      }
        
        
      if (me$through.cv & !outcome.is.logNormally.distributed)
      {
        A <- c(o.sigma$A, list(b=tau.beta, mu=mu))

        if (multiple.regions)
        {
          r <- seq(along=possible.regions)
          region.wt <- rep(NA, length(r))          
          for (j in r) region.wt[j] <- integrate(o.sigma$f, sigma.cutoffs[j], sigma.cutoffs[j+1], A=A)$value
          region.prob <- region.wt * region.prior.density[possible.regions]
            # modif_0.13 Line above now uses region.prior.density 
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
        
        
      # modif_0.13  End of block above, now made non-conditional  
      #}
            

      # Sample from f(mu | sigma)
      
      mu.cond.mean <- moments$sum / data$N
      mu.cond.sd <- sigma/sqrt(data$N)
      
      # modif_0.10 [block]
      
      # modif_0.13 Block commented out
      
      #if (undefined.region.prior.prob)
      #{
      #  # The joint prior for (mu, sigma) is uniform over the whole rectangle domain, hence no region delimiting is necessary
      #  
      #  mu.range <- c(mu.lower, mu.upper)
      #  
      #  if (me$through.cv & !outcome.is.logNormally.distributed)
      #  {
      #    mu <- mu.truncatedData.gen(o.mu, mu.range, mu.cond.mean, sigma=sigma, RData=RData)
      #  }
      #  else
      #  {          
      #    phi <- pnorm.logpcum(mu.range, mean=mu.cond.mean, sd=mu.cond.sd)
      #    log.U <- runif.logp(phi$log.pcum, lower.tail=phi$lower.tail)
      #    mu <- qnorm(log.U, mean=mu.cond.mean, sd=mu.cond.sd, log.p=T, lower.tail=phi$lower.tail)
      #  }
      #  
      #  # Define region in which (mu, sigma) is sitting
      #  tmp <- mu + z*sigma
      #  region <- 1 + sum(tmp > theta)
      #}
      #else
      #{ 
      
      # modif_0.13 The whole block of code below is not conditional anymore (it was in an else-block before)
      
      mu.cutoffs <- theta - z*sigma
      w <- which(mu.cutoffs >= mu.lower & mu.cutoffs <= mu.upper)
      mu.cutoffs <- unique(c(mu.lower, mu.cutoffs[w], mu.upper))
        
      if (length(w) > 0) region <- c(w, max(w) + 1) # regions in which (mu, sigma) can fall
      else region <- 1
        
      if (me$through.cv & !outcome.is.logNormally.distributed)
      {        
        A <- c(o.mu$A, list(mu.mean=mu.cond.mean, s=sigma, s2=sigma^2))
        
        if (length(region) > 1)
        {
          r <- seq(along=region)
          region.wt <- rep(NA, length(region))  
          for (j in r) region.wt[j] <- integrate(o.mu$f, mu.cutoffs[j], mu.cutoffs[j+1], A=A)$value
          region.prob <- region.wt * region.prior.density[region] # modif_0.13
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
        phi <- pnorm.logpcum(mu.cutoffs, mean=mu.cond.mean, sd=mu.cond.sd)
          
        if (length(region) > 1)
        {
          logp <- logp.from.logpcum(phi) + region.prior.log.density[region] # modif_0.13
          j <- logp.sample(logp)
        }
        else
        {
          j <- 1
        }
          
        log.pcum <- phi$log.pcum[j+c(0,1)] # vector of length 2
        log.U <- runif.logp(log.pcum, lower.tail=phi$lower.tail)
        mu <- qnorm(log.U, mean=mu.cond.mean, sd=mu.cond.sd, log.p=T, lower.tail=phi$lower.tail)
      }
        
      region <- region[j] # define region in which (mu, sigma) is sitting
        
      # modif_0.13  End of block above, now made non-conditional
      #}
    
    
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
  
  region.count <- apply(region.sample, 1, table)
    
  if (n.chains == 1)
  {
    region.posterior.prob <- rep(0, n.regions)
    w <- as.numeric(rownames(region.count)) # regions for which count is > 0
    region.posterior.prob[w] <- region.count
  }
  else
  {
    region.posterior.prob <- matrix(0, nrow=n.chains, ncol=n.regions)
    
    for (ch in 1:n.chains)
    {
      tmp <- region.count[[ch]]
      w <- as.numeric(names(tmp)) # regions for which count is > 0
      region.posterior.prob[ch, w] <- tmp
    }
  }
    
  region.posterior.prob <- region.posterior.prob/n.iter


  # Collect information to return in output  
  
  out <- list(quantile=quantile, z.quantile=z,
              n.regions=n.regions, region.area=region.area, 
              region.prior.prob=region.prior.prob, region.prior.density=region.prior.density, region.posterior.prob=region.posterior.prob)
        
              
  tmp <- out.sample(sample, burnin, n.chains, monitor.burnin, me, outcome.is.logNormally.distributed)
  out$sample <- tmp$sample
  if (monitor.burnin) out$burnin <- tmp$burnin
    
  out$mcmc  <- list(n.chains=n.chains, n.burnin=n.burnin, n.iter=n.iter, n.thin=n.thin, monitor.burnin=monitor.burnin)
  out$inits <- list(mu=init.mu, sigma=init.sigma)
  out$data  <- list(y=y, lt=lt, gt=gt, interval.lower=interval.lower, interval.upper=interval.upper)
  out$parms <- list(A=A, mu.lower=mu.lower, mu.upper=mu.upper, sigma.lower=sigma.lower, sigma.upper=sigma.upper,
                    outcome.is.logNormally.distributed=outcome.is.logNormally.distributed)
  
  out
} # end of SEG.riskband