
# Version 0.6 (Aug 2024)


# ------------------------------------------------------------------------------
# New in
# Version 0.6 (Aug 2024)
#
#  Added arg past.data to fct webexpo.stan.inits
#
#                                                            (end of Change Log)


any.me <- function(sd.minmax, cv.minmax)
{
  me <- list(any=FALSE, through.sd=FALSE, through.cv=FALSE, known=FALSE)
  
  if (length(sd.minmax) > 0)
  {
    if (length(sd.minmax) != 2) stop("me.sd.range must be of length 2.")
    me$any <- TRUE
    me$through.sd <- TRUE # indicates that measurement error was specified through a constant sd
    me$range <- sort(sd.minmax)
  }
  
  if (length(cv.minmax) > 0)
  {
    if (me$any) stop("Only me.sd.range or cv.range must be specified [if any measurement error is present].")
    if (length(cv.minmax) != 2) stop("cv.range must be of length 2.")
    me$any <- TRUE
    me$through.cv <- TRUE
    me$through.sd <- FALSE
    me$range <- sort(cv.minmax)
  }
  
  if (me$any)
  {
    me$init <- mean(me$range)
    me$known <- diff(me$range) == 0
  }
  
  me
} # end if any.me


extracted.nodes <- function(stan.out, monitor)
{
  out <- list()
  
  for (node in monitor)
  {
    tmp <- rstan::extract(stan.out, node, permuted=FALSE, inc_warmup=FALSE)
    out[[node]] <- tmp
  }
  
  return(out)
} # end of extracted.nodes


webexpo.stan.inits <- function(y, lt, gt, interval.lower, interval.upper,
                               mu.init, sigma.init, outcome.is.logNormally.distributed,
                               me.sd.range, cv.range, models.folder, priors.label,
                               mu.lower=-Inf, mu.upper=Inf, sigma.lower=0, sigma.upper=Inf,
                               past.data=list(mean=numeric(0), sd=numeric(0), n=numeric(0)))
{
  within.range <- function(tentative.theta, theta.lower, theta.upper, f=0.10)
  {
    x <- tentative.theta
    
    if (x < theta.lower)
    {
      if (is.infinite(theta.upper))  x <- theta.lower + 1
      else                           x <- theta.lower + f*(theta.upper-theta.lower)
    }
    else if (x > theta.upper)
    {
      if (is.infinite(theta.lower))  x <- theta.upper - 1
      else                           x <- theta.upper - f*(theta.upper-theta.lower)
    }
    
    return(x)
  } # end of within.range
  
  
  # Verify that past.data is correctly used (when used)
  
  mean.len <- length(past.data$mean)
  
  if (mean.len > 1 | mean.len != length(past.data$sd) | mean.len != length(past.data$n))
  {
    stop("Elements mean, sd and n in past.data must be of same length (0 or 1).")
  }
  
  if (mean.len == 1 && past.data$n == 1)  stop('Sample size (n) in past data must be > 1. (You may want to add them to current collected data.)')
  
  past.data.used <- mean.len == 1
  
  
  # Assess for presence of Measurement Error
  
  me <- any.me(me.sd.range, cv.range) # see if measurement error is desired
  
  
  # Prepare data
  
  interval <- matrix(c(interval.lower, interval.upper), ncol=2)
  
  
  data <- list(y=y, lt=lt, gt=gt, interval=interval)
  
  
  if (me$through.cv && any(c(data$y, data$lt, data$gt, data$interval[,1]) < 0))  stop('Measurement Error parametrized with CV while having negative data values makes little sense. Please reconsider.')
  
  
  if (is.null(mu.init) || is.null(sigma.init))
  {
    tmp <- data
    # exclude from inits calculation
    tmp$lt <- NULL
    tmp$gt <- NULL
    
    if (outcome.is.logNormally.distributed)  tmp <- lapply(tmp, log)
    
    tmp$interval <- apply(tmp$interval, 1, mean)
    tmp <- unlist(tmp, use.names=FALSE)
    
    
    if (is.null(mu.init))  mu.init <- mean(tmp)
    
    if (is.null(sigma.init))
    {
      sigma.init <- sqrt(var(tmp))
      if (is.na(sigma.init))  sigma.init <- -Inf # will be corrected by call to within.range() below
    }
  }
  
  
  mu.init    <- within.range(mu.init, mu.lower, mu.upper)
  sigma.init <- within.range(sigma.init, sigma.lower, sigma.upper)
  
  
  inits <- list(mu=mu.init, sigma=sigma.init)
  
  
  data$N <- length(y)
  data$L <- length(lt)
  data$G <- length(gt)
  data$I <- nrow(interval)
  
  
  # Complete data & inits with Measurement Error-related parameters & initial latent true values
  
  if (me$any)
  {
    if (me$through.cv)
    {
      data$CV_LO <- cv.range[1]
      data$CV_HI <- cv.range[2]
      
      inits$cv <- (data$CV_LO + data$CV_HI) / 2
      
      # Generate inits for (latent) true values
      
      inits$true_y   <- rnorm(data$N, data$y,   inits$cv*data$y)
      inits$true_lt  <- rnorm(data$L, data$lt,  inits$cv*data$lt)
      inits$true_gt  <- rnorm(data$G, data$gt,  inits$cv*data$gt)
    }
    else
    {
      data$ME_SD_LO <- me.sd.range[1]
      data$ME_SD_HI <- me.sd.range[2]
      
      inits$me_sd <- (data$ME_SD_LO + data$ME_SD_HI) / 2
      
      # Generate inits for (latent) true values
      
      inits$true_y   <- rnorm(data$N, data$y,   inits$me_sd)
      inits$true_lt  <- rnorm(data$L, data$lt,  inits$me_sd)
      inits$true_gt  <- rnorm(data$G, data$gt,  inits$me_sd)
    }
    
    
    interval <- data$interval
    if (outcome.is.logNormally.distributed)  interval <- log(interval)
    
    p0 <- pnorm(interval[,1], inits$mu, inits$sigma)
    p1 <- pnorm(interval[,2], inits$mu, inits$sigma)
    u <- runif(data$I, p0, p1)
    
    true_int <- qnorm(u, inits$mu, inits$sigma)
    if (outcome.is.logNormally.distributed)  true_int <- exp(true_int)
    inits$true_int <- true_int
  }
  
  
  if (outcome.is.logNormally.distributed && !me$any)
  {
    data$y <- log(data$y)
    data$lt <- log(data$lt)
    data$gt <- log(data$gt)
    data$interval <- log(data$interval)
  }
  
  
  # List of monitored nodes
  
  monitor <- c('mu', 'sigma')
  
  if (me$any)
  {
    tmp <- ifelse(me$through.cv, 'cv', 'me_sd')
    monitor <- c(monitor, tmp)
  }
  
  
  # Pick the appropriate Stan model to submit
  
  tmp <- c('SEG','stan', 'model', priors.label)  ##### "SEG" added by Jerome Lavoué August 25th for clarity in model terminology
  
  if (me$any)
  {
    distrn <- ifelse(outcome.is.logNormally.distributed, 'logNormal', 'Normal')
    me <- ifelse(me$through.cv, 'mecv', 'mesd')
    tmp <- c(tmp, distrn, me)
  }
  
  if (past.data.used)
  {
    tmp <- c(tmp, 'pastData')
    
    data$pastData_n    <- past.data$n
    data$pastData_mean <- past.data$mean
    data$pastData_sd   <- past.data$sd
  }
  
  
  # and load it!
  
  tmp <- paste(tmp, collapse='_')
  f <- paste(models.folder, '/', tmp, '.RDS', sep='')
  
  if (!file.exists(f)) stop('File not found -> ', f, '\n\tPlease compile the corresponding stan model first and resubmit.')
  model <- readRDS(f)
  
  return(list(data=data, inits=inits, monitor=monitor, model=model))
} # end of webexpo.stan.inits
