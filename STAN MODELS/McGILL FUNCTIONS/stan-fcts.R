
# Version 0.10 (Aug 2024)


# ------------------------------------------------------------------------------
# New in
# Version 0.10 (Aug 2024)
#
# Added argument 'recompile' to function augment.stan.models.list
# Added as.array protection in fct webexpo.stan.inits
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


augment.stan.models.list <- function(stan.models.list, stan.file, recompile=FALSE)
{
  # Call with recompile=TRUE when you want to force compilation of model
  #                          if it is already present in stan.models.list
  
  if (!is.list(stan.models.list))                              stop("Object stan.models.list is not a list. Please make it an empty list and resubmit.")
  if (!grepl('https:', stan.file) && !file.exists(stan.file))  stop("Stan file not found: ", stan.file)
  
  
  code <- readLines(stan.file)
  tmp <- grep('label', code, value=TRUE)
  model.label <- rev(unlist(strsplit(tmp, ' ')))[1]
  
  m <- match(model.label, names(stan.models.list), nomatch=-1)
  
  
  if (m > 0)
  {
    if (!recompile)  stop('Model read in ', stan.file, ' is already part of your list.')
    stan.models.list <- stan.models.list[-m]
  }
  
  
  cat('Compiling model; please be patient... ')
  t0 <- Sys.time()
  model <- stan_model(model_code=code)
  t1 <- Sys.time()
  t <- round(as.numeric(difftime(t1, t0)), 1)
  cat('Done (compiled in', t, 'seconds).\n')
  
  
  n <- length(stan.models.list)
  if (n == 0)  stan.models.list <- list(model)
  else         stan.models.list <- append(stan.models.list, model)
  
  names(stan.models.list)[n+1] <- model.label
  
  return(stan.models.list)
} # end of augment.stan.models.list


compiled.models.list <- function(stan.models.list)
{
  return(sort(names(stan.models.list)))
} # end of compiled.models.list


drop.model.from.list <- function(stan.models.list, model2drop.label)
{
  # model2drop.label: can be of length > 1
  
  if (!is.list(stan.models.list))  stop('stan.models.list is not a list.')
  model.names <- names(stan.models.list)
  
  if (is.null(model.names))  stop('stan.models.list is empty.')
  
  m <- match(model2drop.label, model.names, nomatch=-1)
  if (all(m < 0))  stop('Model(s) ', model2drop.label, ' is(are) absent from stan.models.list')
  
  m <- m[m > 0]
  stan.models.list <- stan.models.list[-m]
  
  return(stan.models.list)
} # end of drop.model.from.list


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
                               me.sd.range, cv.range, models.list, model.label,
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
  
  
  if (length(models.list) == 0)  stop('models.list is empty. Please submit a list with at list one compiled model in.')
  if (!is.list(models.list))     stop('models.list is not a R list.')
  
  
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
      my.cv <- mean(cv.range)
      
      if (me$known)
      {
        data$cv <- my.cv
      }
      else
      {
        data$CV_LO <- cv.range[1]
        data$CV_HI <- cv.range[2]
        
        inits$cv <- my.cv
      }
      
      # Generate inits for (latent) true values
      
      inits$true_y   <- rnorm(data$N, data$y,  my.cv*data$y)
      inits$true_lt  <- rnorm(data$L, data$lt, my.cv*data$lt)
      inits$true_gt  <- rnorm(data$G, data$gt, my.cv*data$gt)
    }
    else
    {
      my.me_sd <- mean(me.sd.range)
      
      if (me$known)
      {
        data$me_sd <- my.me_sd
      }
      else
      {
        data$ME_SD_LO <- me.sd.range[1]
        data$ME_SD_HI <- me.sd.range[2]
        
        inits$me_sd <- my.me_sd
      }
      
      # Generate inits for (latent) true values
      
      inits$true_y   <- rnorm(data$N, data$y,  my.me_sd)
      inits$true_lt  <- rnorm(data$L, data$lt, my.me_sd)
      inits$true_gt  <- rnorm(data$G, data$gt, my.me_sd)
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
    data$y  <- log(data$y)
    data$lt <- log(data$lt)
    data$gt <- log(data$gt)
    data$interval <- log(data$interval)
  }
  
  
  # List of monitored nodes
  
  monitor <- c('mu', 'sigma')
  
  if (me$any && !me$known)
  {
    tmp <- ifelse(me$through.cv, 'cv', 'me_sd')
    monitor <- c(monitor, tmp)
  }
  
  
  if (past.data.used)
  {
    data$pastData_n    <- past.data$n
    data$pastData_mean <- past.data$mean
    data$pastData_sd   <- past.data$sd
  }
  
  
  # Pick the appropriate Stan model to submit
  
  model <- webexpo.stan.model(model.label, me, outcome.is.logNormally.distributed, models.list, past.data.used)
  
  
  # Make sure a few objects are vectors/arrays
  
  data$y  <- as.array(data$y)
  data$lt <- as.array(data$lt)
  data$gt <- as.array(data$gt)
  
  inits$true_y  <- as.array(inits$true_y)
  inits$true_lt <- as.array(inits$true_lt)
  inits$true_gt <- as.array(inits$true_gt)
  
  
  return(list(data=data, inits=inits, monitor=monitor, model=model))
} # end of webexpo.stan.inits


webexpo.stan.model <- function(model.label, me, outcome.is.logNormally.distributed, models.list, past.data.used=FALSE, use.uniform.prior.on.sds=FALSE)
{
  tmp <- c('SEG', model.label)
  
  
  if (model.label == 'BetweenWorkers' || me$any)
  {
    distrn <- ifelse(outcome.is.logNormally.distributed, 'logNormal', 'Normal')
    tmp <- c(tmp, distrn)
  }
  
  
  if (model.label == 'BetweenWorkers')
  {
    sigma.prior <- ifelse(use.uniform.prior.on.sds, 'Uniform', 'logNormal')
    tmp <- c(tmp, sigma.prior)
  }
  
  
  if (me$any)
  {
    ME <- ifelse(me$through.cv, 'mecv', 'mesd')
    if (me$known)  ME <- paste(ME, 'known', sep='')
    tmp <- c(tmp, ME)
  }
  
  
  if (past.data.used)  tmp <- c(tmp, 'pastData')
  
  
  # Find the appropriate model in models.list
  
  tmp <- tolower(paste(tmp, collapse='_'))
  
  model.no <- match(tmp, names(models.list))
  if (is.na(model.no))  stop('Expected model (', tmp, ') not found in models.list; sorry.')
  
  
  model <- models.list[[model.no]]
  
  return(model)
} # end of webexpo.stan.model
