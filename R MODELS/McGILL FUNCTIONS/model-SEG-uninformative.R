# Version 0.12 (July 2017)


  # Change Log **************************************************************************
  #
  #   When updated, look for comments with new_* and modif_*
  #   to rapidly identify new/modified code.
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
  #  - mu.lim was renamed mu.range
  #
  #
  # Version 0.10 [replaces model.Uninformative]:
  # ------------
  # - me.sd argument was dropped and replaced by me.sd.range
  # - cv.range was added 
  #     only one of me.sd.range or cv.range must be specified [of length 2, if any]
  # - data$n was replaced data$size$y
  # - using empty.matrices
  # - using mu.truncatedData.gen and sigma.truncatedData.gen 
  # - using sigma.gen.object
  # - using sqrt.invertedGamma.gen
  # - using truevalues.gen 
  # - using y.gen and y.gen.inits
  # - dropped use.alt.sd.posterior
  # - dropped tau.range (and everything called tau.*, in favor of sigma.*)
  # - a few corrections were made (look for new_ & modif_*)
  # - objects mu.sample & sd.sample were regrouped in list 'sample'
  #     while burnin.mu & burnin.sd were regrouped in list 'burnin'
  #     > both are wrapped through function out.sample (defined in fcts.R) at the end of function
  #
  # - use.uniform.prior.on.sd, tau.alpha, tau.beta and tau.cond.alpha were dropped


SEG.uninformative <- function(y=numeric(0), lt=numeric(0), gt=numeric(0), 
  interval.lower=numeric(0), interval.upper=numeric(0),
  n.iter=15000, n.thin=1, n.chains=1, n.burnin=500, monitor.burnin=F,
  mu.lower=-100, mu.upper=100,
  init.mu=rep(default.inits$mu, n.chains), init.sd=rep(default.inits$sigma, n.chains),
  outcome.is.logNormally.distributed=TRUE, sd.range=c(0, 100), 
  me.sd.range=numeric(0), cv.range=numeric(0),
  save.RData=F, RData.dir='c:/users/jerome')
{
  # Notes:  
  # - me.sd.range is the range of Measurement Error SD (optional) [new_0.10]
  # - cv.range    is the range of Measurement Error Coefficient of Variation (optional)
  #   IMPORTANT: only one of me.sd.range or cv.range can be entered  

  # new_0.12: added arguments save.RData and RData.dir
  RData <- list(save=save.RData, dir=RData.dir) # new_0.12
    
  me <- any.me(me.sd.range, cv.range) # see if measurement error is desired [new_0.10]


  # Remove y values that are NAs and count number of censored values (y = NA and lt and/or gt not NA)
  data <- data.summary(y=y, lt=lt, gt=gt, interval.lower=interval.lower, interval.upper=interval.upper, 
                       logNormal.distrn=outcome.is.logNormally.distributed, me.through.cv=me$through.cv) 
                       # me$through.cv argument was added above [modif_0.10]
     
  mu.range <- c(mu.lower, mu.upper)
     
  # [block] new_0.10                 
  
  gen.y <- list()
  true.values <- list()
  
  default.inits <- Default.inits(data, outcome.is.logNormally.distributed, mu.range <- c(mu.lower, mu.upper), sd.range) # new_0.11
  
  
  # Initialize mu [modif_0.11]
  
  if (length(init.mu) > 0 & length(init.mu) != n.chains) stop("init.mu must be of length 'n.chains', that is: ", n.chains)
  
  if (length(init.mu) == 0)
  {
    init.mu <- rep(default.inits$mu, n.chains)
  }
  else if (length(init.mu) == 1 & n.chains > 1)
  {
    init.mu <- rep(init.mu, n.chains)
  }
      
      
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

  sample <- empty.matrices(n.iter, n.chains, me)

  burnin <- list()
  if (monitor.burnin) burnin <- empty.matrices(n.burnin, n.chains, me)
  
  
  M.iter <- n.burnin + n.iter * n.thin
  
  for (ch in seq(n.chains))
  {
    # Initial values for mu and sigma (modif_0.11)
    
    mu <- init.mu[ch]

    if (is.na(init.sd[ch]) | init.sd[ch] == 0) sigma <- default.inits$sigma
    else sigma <- init.sd[ch]
        

    # Initialize measured values for subjects with censored values [modif_0.10]
    if (data$any.censored$any) gen.y <- y.gen.inits(data, mu, sigma, me$through.cv, outcome.is.logNormally.distributed)
       
    if (me$any) me$parm <- me$init # initialize Measurement Error parameter value
       
    
    # Start MCMC

    saved.iter <- 0
  
    for (iter in 1:M.iter)
    {    
      # Sample true values (in presence of measurement error) [new_0.10]
      if (me$any) true.values <- truevalues.gen(gen.y, data, mu, sigma, me, outcome.is.logNormally.distributed, o=o.tv, RData=RData)
      
    
      # Sample y values for censored observations
      if (data$any.censored$any) gen.y <- y.gen(true.values, data, sigma, me, outcome.is.logNormally.distributed, mu=mu)
      
      
      # Compute data points sum and square sum
      moments <- out.logout.moments(me$any, outcome.is.logNormally.distributed, data, gen.y, true.values)
      

      # Sample from f(sigma | mu)
      # modif_0.10
      
      sigma.beta <- (moments$sum2 - 2*mu*moments$sum + data$N*mu^2) / 2
      if (sigma.beta < 1e-6) sigma.beta <- 0 # protection against numeric imprecision

      if (me$through.cv & !outcome.is.logNormally.distributed)
      {
        sigma <- sigma.truncatedData.gen(o.sigma, sd.range, sigma.beta, mu, current.sigma=sigma, RData=RData)
      }
      else
      {     
        sigma <- sqrt.invertedGamma.gen(data$N, sigma.beta, sd.range, o=o.sigma, RData=RData)
      }
      
    
      # Sample from f(mu | sigma)
      # modif_0.10
      
      
      mu.cond.mean <- moments$sum / data$N

      if (me$through.cv & !outcome.is.logNormally.distributed)
      {  
        mu <- mu.truncatedData.gen(o.mu, mu.range, mu.cond.mean, sigma, RData=RData)
      }
      else
      {
        mu.cond.sd <- sigma/sqrt(data$N)
        mu <- rnorm.censored(mu.cond.mean, mu.cond.sd, lower=mu.range[1], upper=mu.range[2])
      }
      
      
      # Sample Measurement Error from its posterior density
      
      if (me$any & !me$known) me$parm <- me.gen(o.me, me, data, gen.y, true.values, RData=RData)

      
      # Save values only when iter# modulo thinning = 0
      
      if (iter <= n.burnin)
      {
        if (monitor.burnin)
        {
          burnin$mu[ch, iter]  <- mu
          burnin$sd[ch, iter] <- sigma
          
          if (me$any & !me$known) burnin$me.parm[ch, iter] <- me$parm
        }
      }
      else if ((iter-n.burnin)%%n.thin == 0)
      {
        saved.iter <- saved.iter + 1
        sample$mu[ch, saved.iter] <- mu
        sample$sd[ch, saved.iter] <- sigma
        
        if (me$any & !me$known) sample$me.parm[ch, saved.iter] <- me$parm
      }
    }
  }
  
  
  out <- out.sample(sample, burnin, n.chains, monitor.burnin, me, outcome.is.logNormally.distributed)
  
  out$mcmc  <- list(n.chains=n.chains, n.burnin=n.burnin, n.iter=n.iter, n.thin=n.thin, monitor.burnin=monitor.burnin)
  out$inits <- list(mu=init.mu)
  out$data  <- list(y=y, lt=lt, gt=gt, interval.lower=interval.lower, interval.upper=interval.upper)
  
  
  parms.list <- list(mu.lower=mu.lower, mu.upper=mu.upper,
                    outcome.is.logNormally.distributed=outcome.is.logNormally.distributed, 
                    sd.range=sd.range)
                    
                    
  out$parms <- parms.list
  
  me[c("known", "parm")] <- NULL # Drop these items from 'me' list
  out$me <- me        
                    
  out
} # end of SEG.uninformative