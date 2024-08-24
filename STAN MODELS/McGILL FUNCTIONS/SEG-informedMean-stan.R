
# Version 0.5 (Aug 2024)


# ------------------------------------------------------------------------------
# New in
# Version 0.5 (Aug 2024)
#
#  Split the Stan models into individual .RDS files
#
#                                                            (end of Change Log)



library(rstan)

stan.folder <- 'C:/Users/patri/home/consultation/L/Lavoue/webexpo/stan' # ICI modifier ce path


setwd(stan.folder)
source('stan-fcts.R')


SEG.informedmean.stan <- function(y=numeric(0), lt=numeric(0), gt=numeric(0),
  interval.lower=numeric(0), interval.upper=numeric(0),
  n.iter=15000, n.burnin=500,
  mu.mean=NULL, mu.sd=NULL,
  log.sigma.mu=-0.1744,
  log.sigma.sd=1/sqrt(log.sigma.prec), log.sigma.prec=2.5523,
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
  if (is.null(mu.mean) || is.null(mu.sd)) stop('mu.mean & mu.sd must be both defined.')


  # Prepare data & inits
  o <- webexpo.stan.inits(y, lt, gt, interval.lower, interval.upper,
                          init.mu, init.sd,
                          outcome.is.logNormally.distributed, me.sd.range, cv.range, models.folder, 'informedMean')


  # Add parameter values to data

  o$data$MU_MEAN <- mu.mean
  o$data$MU_SD   <- mu.sd

  o$data$LOGSIGMA_MEAN <- log.sigma.mu
  o$data$LOGSIGMA_SD   <- log.sigma.sd


  # Submit Stan model

  stan.out <- sampling(o$model, data=o$data, init=list(o$inits), pars=o$monitor,
                   chains=1, warmup=n.burnin, iter=n.burnin+n.iter, cores=1, show_messages=!silent, verbose=!silent)


  # Return MCMC sampled values

  out <- extracted.nodes(stan.out, o$monitor)

  return(out)
} # end of SEG.informedmean.stan
