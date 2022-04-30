# WebExpo version 0.5h (April 2020)
      


# Change Log
# ----------

# Version 0.5h (April 2020)
# -------------------------
#   Added argument old.McNally=T to call to inits.BetweenWorker
#
# Version 0.5g (April 2020)
# ------------------------
#   Use the new fct inits.McNally (see file initsMcNally.R)
#     to better initialize the MCMC chain
#   Added the option use.old.inits
#
# Version 0.5f (April 2020)
# ------------------------
#   Modified path for RData file where current values are saved in case of error
#
# Version 0.5e (April 2020)
# -------------------------
#  Corrected the creation of the burnin output // Jerome Lavoue
#
# Version 0.5d (March 2020)
# -------------------------
#  Changed the value of b.min
#
# Version 0.5c (March 2020)
# -------------------------
#   Removed sigma.(within|between).range in calls to  
#
# Version 0.5b
# ------------
#   Corrected the first argument in two calls to the function rgamma.truncated
#
# Other
# -----
#   Coquille corrigée le 25 avril 2018


model.McNally <- function(y=numeric(0), worker=numeric(0), 
  lt=numeric(0), worker.lt=numeric(0),
  gt=numeric(0), worker.gt=numeric(0),
  interval.lower=numeric(0), interval.upper=numeric(0), worker.interval=numeric(0),
  n.chains=1, n.iter=15000, n.burnin=500, n.thin=1, monitor.burnin=F,
  log.sigma.between.mu=-0.8786, log.sigma.between.prec=1.634,
  log.sigma.within.mu=-0.4106, log.sigma.within.prec=1.9002,
  mu.overall.lower=-100, mu.overall.upper=100,
  outcome.is.logNormally.distributed=T,
  use.uniform.prior.on.sds=F, sigma.between.range=c(0,100), sigma.within.range=c(0,100),
  use.old.inits=T)
{
  # Notes:
  # - sigma.between.range & sigma.within.range           are relevant only when  use.uniform.prior.on.sds = TRUE
  # -   log.sigma.between.mu,   log.sigma.between.prec,
  #     log.sigma.within.mu and log.sigma.within.prec    are relevant only when  use.uniform.prior.on.sds = FALSE

  if (length(worker) != length(y))     stop("y and worker must be of same length.")
  if (length(worker.lt) != length(lt)) stop("lt and worker.lt must be of same length.")
  if (length(worker.gt) != length(gt)) stop("gt and worker.gt must be of same length.")
  if (length(worker.interval) != length(interval.lower) || length(worker.interval) != length(interval.upper)) stop("interval.lower, interval.upper and worker.interval must be of same length.")
  
  # b.min <- 1e-8 # modif_0.5d
  b.min <- 1e-12 # new_0.5d

  # new_0.5f
  error.file <- paste(tempdir(), "webExpo-McNally_error.RData", sep="/")
  error.file <- gsub("\\\\", "/", error.file)

  data <- data.summary(y=y, lt=lt, gt=gt , interval.lower=interval.lower, interval.upper=interval.upper, 
                       worker=worker, worker.lt=worker.lt, worker.gt=worker.gt, worker.interval=worker.interval,
                       take.log=outcome.is.logNormally.distributed)
  
  
  # Build data frame and identify lines correspond to left-, right- or interval-censored entries
  
  df <- data.frame(y=data$y, id=data$id$y)
  
  if (data$any.censored$any)
  {
    which.censored <- list()
    
    if (data$any.censored$l)
    {
      tmp.df <- data.frame(y=data$cens$l, id=data$id$l)
      df <- rbind.data.frame(df, tmp.df)
      which.censored$l <- seq(to=nrow(df), length=length(data$id$l))   
    }
    
    if (data$any.censored$r)
    {
      tmp.df <- data.frame(y=data$cens$r, id=data$id$r)
      df <- rbind.data.frame(df, tmp.df)
      which.censored$r <- seq(to=nrow(df), length=length(data$id$r))   
    }
    
    if (data$any.censored$i)
    {
      tmp.df <- data.frame(y=(data$cens$i$l + data$cens$i$r)/2, id=data$id$i)
      df <- rbind.data.frame(df, tmp.df)
      which.censored$i <- seq(to=nrow(df), length=length(data$id$i))   
    }
  }
  
  df$id <- factor(df$id)
  

  log.sigma.within.sd <- 1/sqrt(log.sigma.within.prec)
  log.sigma.between.sd <- 1/sqrt(log.sigma.between.prec)

  
  if (!data$any.censored$any)
  {
    y.worker.avg <- aggregate(y~id, df, mean)$y
    y.avg <- mean(df$y)
    
    if (use.uniform.prior.on.sds)
    {
      # If there is no variation observed within subjects, 
      # then we require a non-null lower bound for sigma.within
      
      problem <- all(data$worker$count==1)
    
      if (!problem)
      {
        y.worker.var <- aggregate(y~id, df, var)$y
        y.worker.var <- y.worker.var[!is.na(y.worker.var)]
        problem <- all(y.worker.var==0)
      } 
    
      if (problem && sigma.within.range[1] == 0) stop("Lower bound for sigma.within must be > 0 due to observed null within-subjects variance.")
    }
  }
  
  
  if (use.uniform.prior.on.sds)
  {    
    tau.within.range  <- 1/rev(sort(sigma.within.range))^2
    tau.between.range <- 1/rev(sort(sigma.between.range))^2
  }
  

  # Prepare output objects with appropriate dimensions

  mu.overall.sample  <- matrix(NA, nrow=n.chains, ncol=n.iter)
  mu.worker.sample <- matrix(NA, nrow=data$n.workers*n.chains, ncol=n.iter)
  sigma.within.sample <- matrix(NA, nrow=n.chains, ncol=n.iter)
  sigma.between.sample <- matrix(NA, nrow=n.chains, ncol=n.iter)


  if (monitor.burnin)
  {
    burnin.mu.overall  <- matrix(NA, nrow=n.chains, ncol=n.burnin)
    burnin.mu.worker <- matrix(NA, nrow=data$n.workers*n.chains, ncol=n.burnin)
    burnin.sigma.within <- matrix(NA, nrow=n.chains, ncol=n.burnin)
    burnin.sigma.between <- matrix(NA, nrow=n.chains, ncol=n.burnin)
  }


  sample.id <- data.frame(worker=rep(data$worker$patkey, n.chains), chain=rep(seq(n.chains), rep(data$n.workers, n.chains)))
  
  
  # Initial values
  
  # new_0.5g

  if (use.old.inits) 
  {
    inits <- list(converged=F)
  }
  else
  { 
    # modif_0.5h inverted order of arguments sigma.within.range & sigma.between.range
    # modif_0.5h and added argument old.McNally=T

    inits <- inits.BetweenWorker(data, use.uniform.prior.on.sds, sigma.within.range, sigma.between.range,
                                 log.sigma.within.mu, log.sigma.within.prec, 
                                 log.sigma.between.mu, log.sigma.between.prec, old.McNally=T)
  }
  

  # modif_0.5g (in earlier version, the block below was not conditional)

  if (!inits$converged)
  {
    # For mu.worker, take observed mean by worker
    mu.worker <- aggregate(y~id, data=df, mean)$y
        
    mu.overall <- mean(mu.worker)
    mu.worker <- mu.worker - mu.overall # center mu.worker#

    # For sigma.within, use sd of initial values around worker predicted means
    predicted <- mu.overall + mu.worker[df$id]
    sigma.within <- sd(df$y-predicted) 
    
    inits <- list(mu.overall=mu.overall, mu.worker=mu.worker, sigma.within=sigma.within)
  }

  M.iter <- n.burnin + n.iter * n.thin

  for (ch in seq(n.chains))
  {
    # Set initial values
    mu.overall <- inits$mu.overall
    mu.worker <- inits$mu.worker
    sigma.within <- inits$sigma.within
  
    saved.iter <- 0


    for (iter in 1:M.iter)
    {
      # Sample y latent values for censored measurements
      
      if (data$any.censored$any)
      {
        if (data$any.censored$l)
        {
          tmp.mean <- mu.overall + mu.worker[data$id$l]
          tmp <- rnorm.censored(tmp.mean, sigma.within, lower=data$cens$l)
          df[which.censored$l, 1] <- tmp
        }
     
        if (data$any.censored$r)
        {
          tmp.mean <- mu.overall + mu.worker[data$id$r]
          tmp <- rnorm.censored(tmp.mean, sigma.within, upper=data$cens$r)
          df[which.censored$r, 1] <- tmp
        }        

        if (data$any.censored$i)
        {
          tmp.mean <- mu.overall + mu.worker[data$id$i]
          tmp <- rnorm.censored(tmp.mean, sigma.within, lower=data$cens$i$l, upper=data$cens$i$r)
          df[which.censored$i, 1] <- tmp
        }   
         
        
        # Compute subject means including newly imputed values
        
        y.worker.avg <- aggregate(y~id, df, mean)$y
        y.avg <- mean(df$y)
      }
    
    
      # Sample from f(sigma.within | other parms)
      
      residuals <- df$y - mu.overall - mu.worker[df$id]
      b <- sum(residuals^2)/2
      
      if (use.uniform.prior.on.sds)
      {
        if (b < b.min)
        {
          a <- data$N
          u <- runif(1)
          sigma.within <- sum((sigma.within.range^(1-a))*c(1-u,u))^(1/(1-a))
        }
        else if (data$N <= 2)
        {
          sigma.within <- sigma.gen.icdf(data$N, b, Banerjee=T, range=sigma.within.range)
            # This is not Banerjee, but the posterior density for sigma.within 
            # is the same as in Banerjee model, hence the option Banerjee=T above
        }
        else
        {
          tau <- rgamma.truncated((data$N - 1) / 2, b, tau.within.range) # modif_0.5b corrected first argument
          sigma.within <- 1/sqrt(tau)
        }
      }
      else
      {
        # sigma.within <- sigma.gen.icdf(data$N, b, log.sigma.within.mu, log.sigma.within.sd, range=sigma.within.range)  # modif_0.5c (see line below)
        sigma.within <- sigma.gen.icdf(data$N, b, log.sigma.within.mu, log.sigma.within.sd) # new_0.5c
      }
      
      
      # Sample from f(sigma.between | other parms)
      
      b <- sum(mu.worker^2)/2
      
      if (use.uniform.prior.on.sds)
      {
        if (b < b.min)
        {
          a <- data$n.workers
          u <- runif(1)
          sigma.between <- sum((sigma.between.range^(1-a))*c(1-u,u))^(1/(1-a))
        }
        else if (data$n.workers <= 2)
        {
          sigma.between <- sigma.gen.icdf(data$n.workers, b, Banerjee=T, range=sigma.between.range)
            # This is not Banerjee, but the posterior density for sigma.between 
            # is the same as in Banerjee model, hence the option Banerjee=T above
        }
        else
        {
          tau <- rgamma.truncated((data$n.workers - 1) / 2, b, tau.between.range) # modif_0.5b corrected first argument
          sigma.between <- 1/sqrt(tau)
        }
        
        # ICI bloc tempo
        if (sigma.between == 0)
        {
          # modif_0.5f
          save('data', 'b', 'sigma.between.range', 'residuals', 'df',
               'mu.overall', 'mu.worker', 'sigma.within', file=error.file)
          msg <- c("fin pcq sigma.between = 0", "Current objects were saved to file", error.file)
          msg <- paste(msg, collapse='\n')
          stop(msg)
        }
      }
      else
      {
        # sigma.between <- sigma.gen.icdf(data$n.workers, b, log.sigma.between.mu, log.sigma.between.sd, range=sigma.between.range) # modif_0.5c (replaced by line below)
        sigma.between <- sigma.gen.icdf(data$n.workers, b, log.sigma.between.mu, log.sigma.between.sd) # new_0.5c
      }
      
      
      # Sample from from f(mu.overall | other parms)
      
      tmp.mean <- y.avg - sum(data$worker$count*mu.worker)/data$N
      tmp.sd <- sigma.within/sqrt(data$N)
    
      mu.overall <- rnorm.censored(tmp.mean, tmp.sd, lower=mu.overall.lower, upper=mu.overall.upper)

      
      # Sample from f(mu.worker's | other parms)
      
      sigma2.A <- (sigma.within^2)/data$worker$count
      sigma2.B <- sigma.between^2
      
      mu.A <- y.worker.avg - mu.overall
      
      tmp.mean <- mu.A*sigma2.B / (sigma2.A + sigma2.B)
      tmp.var <- sigma2.A * sigma2.B / (sigma2.A + sigma2.B)
      mu.worker <- rnorm(data$n.workers, mean=tmp.mean, sd=sqrt(tmp.var))

      
      # Save values only when iter# modulo thinning = 0
    
      if (iter <= n.burnin)
      {
        if (monitor.burnin)
        {
          burnin.mu.overall[ch, iter]    <- mu.overall
          burnin.sigma.within[ch, iter]  <- sigma.within
          burnin.sigma.between[ch, iter] <- sigma.between
          
          w <- seq(data$n.workers) + (ch-1) * data$n.workers + (iter-1) * n.chains * data$n.workers
          burnin.mu.worker[w] <- mu.worker
        }
      }
      else if ((iter-n.burnin)%%n.thin == 0)
      {
        saved.iter <- saved.iter + 1
        
        mu.overall.sample[ch, saved.iter]    <- mu.overall
        sigma.within.sample[ch, saved.iter]  <- sigma.within
        sigma.between.sample[ch, saved.iter] <- sigma.between
        
        w <- seq(data$n.workers) + (ch-1) * data$n.workers + (saved.iter-1) * n.chains * data$n.workers
        mu.worker.sample[w] <- mu.worker
      }
    }
  }
  
  
  out <- list(sample=list(which=sample.id, mu.overall=mu.overall.sample, mu.worker=mu.worker.sample,
                          sigma.within=sigma.within.sample, sigma.between=sigma.between.sample))


  
    
  if (monitor.burnin)
  {
    # modif_0.5e : burnin.log.sigma.* replaced by burnin.sigma.* in line below 
    out$burnin <- list(mu.overall=burnin.mu.overall, mu.worker=burnin.mu.worker,
                       sigma.within=burnin.sigma.within, sigma.between=burnin.sigma.between)
    
    if (n.chains == 1) out$burnin <- lapply(out$burnin, as.vector)
  }
  
                           
  
  out$mcmc <- list(n.chains=n.chains, n.burnin=n.burnin, n.iter=n.iter, n.thin=n.thin, monitor.burnin=monitor.burnin)
  out$inits <- inits
  out$data <- list(y=y, lt=lt, gt=gt, interval.lower=interval.lower, interval.upper=interval.upper, 
                       worker=worker, worker.lt=worker.lt, worker.gt=worker.gt, worker.interval=worker.interval)
     
                       
  if (use.uniform.prior.on.sds)
  {
    parms.list <- list(sigma.within.range=sigma.within.range, sigma.between.range=sigma.between.range)
  }
  else
  {
    parms.list <- list(log.sigma.between.mu=log.sigma.between.mu, log.sigma.between.prec=log.sigma.between.prec,
                      log.sigma.within.mu=log.sigma.within.mu, log.sigma.within.prec=log.sigma.within.prec)
  }
              
  parms.list$mu.overall.lower <- mu.overall.lower
  parms.list$mu.overall.upper <- mu.overall.upper
  parms.list$outcome.is.logNormally.distributed <- outcome.is.logNormally.distributed
  parms.list$use.uniform.prior.on.sds <- use.uniform.prior.on.sds
  
  out$parms <- parms.list         
                      
  out
} # end of model.McNally
