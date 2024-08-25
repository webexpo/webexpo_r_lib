############# WEBEXPO STAN LIBRARY #######################

##### SCRIPT FOR the INFORMEDVAR FUNCTION

# Requires the RSTAN library to be active

# Requires sourcing the stan-fcts.R script

# Requires the STAN MODEL to be used to be pre-compiled with the correct name in a folder selected by the user.



# Version 0.8 (Aug 2024)


# ------------------------------------------------------------------------------
# New in
# Version 0.8 (Aug 2024)
#
#  Added argument past.data
#
# 



SEG.informedvar.stan <- function(y=numeric(0), lt=numeric(0), gt=numeric(0),
                                 interval.lower=numeric(0), interval.upper=numeric(0),
                                 n.iter=15000, n.burnin=500,
                                 mu.lower=-100, mu.upper=100,
                                 log.sigma.mu=-0.1744,
                                 log.sigma.sd=1/sqrt(log.sigma.prec), log.sigma.prec=2.5523,
                                 init.mu=NULL, init.sd=NULL,
                                 outcome.is.logNormally.distributed=TRUE,
                                 me.sd.range=numeric(0), cv.range=numeric(0),
                                 past.data=list(mean=numeric(0), sd=numeric(0), n=numeric(0)),
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
  
  
  
  # Prepare data & inits
  
  o <- webexpo.stan.inits(y, lt, gt, interval.lower, interval.upper,
                          init.mu, init.sd,
                          outcome.is.logNormally.distributed, me.sd.range, cv.range, models.folder, 'informedVar',
                          mu.lower, mu.upper, past.data=past.data)
  
  
  # Add parameter values to data
  
  o$data$MU_LO <- mu.lower
  o$data$MU_HI <- mu.upper
  
  o$data$LOGSIGMA_MEAN <- log.sigma.mu
  o$data$LOGSIGMA_SD   <- log.sigma.sd
  
  
  # Submit Stan model
  
  stan.out <- sampling(o$model, data=o$data, init=list(o$inits), pars=o$monitor,
                       chains=1, warmup=n.burnin, iter=n.burnin+n.iter, cores=1, show_messages=!silent, verbose=!silent)
  
  
  # Return MCMC sampled values
  
  out <- extracted.nodes(stan.out, o$monitor)
  
  return(out)
} # end of SEG.informedvar.stan