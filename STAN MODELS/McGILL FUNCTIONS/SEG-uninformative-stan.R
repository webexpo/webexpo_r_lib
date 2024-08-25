############# WEBEXPO STAN LIBRARY #######################

##### SCRIPT FOR the UNINFFORMEDVAR FUNCTION

# Requires the RSTAN library to be active

# Requires sourcing the stan-fcts.R script

# Requires the STAN MODEL to be used to be pre-compiled with the correct name in a folder selected by the user.

# Version 0.7 (Aug 2024)


# ------------------------------------------------------------------------------
# New in
# Version 0.7 (Aug 2024)
#
#
#
#                                                            (end of Change Log)



SEG.uninformative.stan <- function(y=numeric(0), lt=numeric(0), gt=numeric(0),
                                   interval.lower=numeric(0), interval.upper=numeric(0),
                                   n.iter=15000, n.burnin=500,
                                   mu.lower=-100, mu.upper=100, sd.range=c(0, 100),
                                   init.mu=NULL, init.sd=NULL,
                                   outcome.is.logNormally.distributed=TRUE,
                                   me.sd.range=numeric(0), cv.range=numeric(0),
                                   models.folder=paste(stan.folder, '/models', sep=''),
                                   silent=TRUE)
{
  # Notes:
  # - me.sd.range is the range of Measurement Error SD (optional)
  # - cv.range    is the range of Measurement Error Coefficient of Variation (optional)
  #   IMPORTANT: only one of me.sd.range or cv.range can be entered
  
  # y  -> Data points known exactly
  # lt -> Right-censored Data points  (x < some value)
  # gt -> Left-censored Data points  (x > some value)
  # interval.lower & interval.upper -> Interval-censored Data points  (interval.lower[i] < x < interval.upper[i], 1 <= i <= length(interval.lower))
  
  
  if (length(interval.lower) != length(interval.upper))  stop('interval.lower & interval.upper must be of same length.')
  if (length(sd.range) != 2)  stop('sd.range must be of length 2.')
  
  
  # Prepare data & inits
  
  o <- webexpo.stan.inits(y, lt, gt, interval.lower, interval.upper,
                          init.mu, init.sd,
                          outcome.is.logNormally.distributed, me.sd.range, cv.range, models.folder, 'Uninformative',
                          mu.lower, mu.upper, sd.range[1], sd.range[2])
  
  
  # Add parameter values to data
  
  o$data$MU_LO <- mu.lower
  o$data$MU_HI <- mu.upper
  
  o$data$SIGMA_LO <- sd.range[1]
  o$data$SIGMA_HI <- sd.range[2]
  
  
  # Submit Stan model
  
  stan.out <- sampling(o$model, data=o$data, init=list(o$inits), pars=o$monitor,
                       chains=1, warmup=n.burnin, iter=n.burnin+n.iter, cores=1, show_messages=!silent, verbose=!silent)
  
  
  # Return MCMC sampled values
  
  out <- extracted.nodes(stan.out, o$monitor)
  
  return(out)
} # end of SEG.uninformative.stan
