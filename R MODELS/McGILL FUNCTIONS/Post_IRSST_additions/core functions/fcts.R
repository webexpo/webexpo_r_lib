source("C:\\Users\\blitzrcr\\patrickb.stat\\Lavoue\\pkg.cpfiles\\rgen_icdf.R")


# Version 0.26 (Apr 2021)
# [distributed]

                                                                               
  # Change Log *****************************************************************
  #                                                      
  #   When updated, look for comments with new_* and modif_*
  #   to rapidly identify new/modified code.                                      
  #
  # Version 0.26 (Apr 2021)
  # -----------------------
  #   The function dens.gen.icdf was moved to file dens_gen_icdf.R
  #     AND replaced by a more user/programmer-friendly version,
  #     called rgen.icdf, saved to file rgen_icdf.R
  #
  #   The following fcts were also added:
  #     Dens.approx
  #     Quadratic.solution
  #  
  # Version 0.24.3 (Feb 2020)
  # --------------
  #   Added a comment in the fct sqrt.invertedGamma.gen
  # 
  # Version 0.24.2
  # --------------
  #   Defined (and used) inestimable.lower.limit in mu.truncatedData.gen
  #
  # Versions 0.24 & 0.24.1
  # ----------------------
  #   Added a line to out.logout.moments
  #   Modified mu.truncatedData.* functions
  #
  #
  # Version 0.19
  # ------------
  #   -> Modified log.f.inv.remote in truevalue.gen.object
  #   -> Modified dens.gen.icdf
  #   -> look for comments with _0.19 to viec each new/modified piece of code
  #
  # Version 0.17
  # ------------
  #   Added yet another alternative starting point to truevalue.gen.object
  #
  # Version 0.16
  # ------------
  #   Added argument RData to function truevalues.gen and truevalue.gen
  #     -> useful at development stage only, for monitoring parameters values when dens.gen.icdf is crashing
  #
  # Version 0.15
  # ------------
  #   Added an alternative starting point to truevalue.gen.object
  #
  # Versions 0.12, 0.13 & 0.14
  # --------------------------
  #   Modified function Default.inits to deal with data sets where only one *y* point is not censored
  #                                                                         
  # Version 0.11
  # ------------
  #   New function(s):
  #     - Default.inits
  #     - out.logout.moments
  #     - polynomial.product.coeff
  #     - take.log
  #
  #   All the functions log.f.inv.remote.left & log.f.inv.remote.right were 
  #     combined into log.f.inv.remote
  #     - the new function split.starting.values.leftRight was written to split the 'remote'
  #       solutions into left- or right- side (relatively to mode) starting points
  #
  #   The following functions were modified:
  #     - logPhi.quadratic.approx.coeff (the argument 'degree' was added)
  #     - me.gen
  #     - me.gen.object
  #         > the condition if (me$through.cv & !logNormal.distrn) 
  #           was changed to if (me$through.cv)                                           
  #           and the "else" block was split in two
  #     - real.cubic.roots -> dropped argument accept.when.abs.soln.lt, added arguments l & u
  #     - rnorm.censored (argument take.log was dropped)
  #     - secured.NR.search (more changes than what is indicated through new_* and modif_*; sorry...)
  #     - smoothed.area.estimate
  #     - truevalue.gen
  #     - truevalue.gen.object
  #     - truevalues.gen
  #     - y.gen
  #
  #   A$mesd and A$mesd2 were changed for A$ksi and A$ksi2, respectively
  #
  #   The function(s)
  #     datalist.sum
  #     was (were) dropped
  #
  #   We added remote.right argument when calling smoothed.area.estimate on a few occasions
  #   In every function, math.limits (vector of length 2) was changed to math.lower.limit (scalar)
  #
  #
  # Version 0.10 
  # ------------                             
  #   The following functions were added:                              
  #     - any.me                                                    
  #     - cond.values                                               
  #     - datalist.sum (dropped in version 0.11)
  #     - dens.gen.icdf
  #     - empty.matrices
  #     - logPhi.quadratic.approx
  #     - logPhi.quadratic.approx.coeff
  #     - lphi
  #     - me.gen
  #     - me.gen.object
  #     - mu.truncatedData.gen
  #     - mu.truncatedData.gen.object
  #     - out.sample
  #     - quadratic.solution
  #     - real.cubic.roots 
  #     - renamed.me.parm
  #     - sigma.gen.object
  #     - sigma.truncatedData.gen
  #     - sigma.truncatedData.gen.object 
  #     - sqrt.invertedGamma.gen             
  #     - truevalue.gen
  #     - truevalue.gen.object
  #     - truevalues.gen
  #     - y.gen (largely inspired from code present in model-Uninformative.R version 0.8)
  #     - y.gen.inits
  #
  #   - functions logp.from.logpcum and pnorm.logpcum were moved from sigmaGen.R to here
  #
  #   - rgamma.truncated was modified
  #   - rnorm.censored   was modified
  #   - runif.logp       was modified


any.me <- function(sd.minmax, cv.minmax)
{
  me <- list(any=F, through.sd=F, through.cv=F, known=F)
  
  if (length(sd.minmax) > 0)
  {
    if (length(sd.minmax) != 2) stop("me.sd.range must be of length 2.")
    me$any <- T
    me$through.sd <- T # indicates that measurement error was specified through a constant sd
    me$range <- sort(sd.minmax)
  }
  
  if (length(cv.minmax) > 0)
  {
    if (me$any) stop("Only me.sd.range or cv.range must be specified [if any measurement error is present].")
    if (length(cv.minmax) != 2) stop("cv.range must be of length 2.")
    me$any <- T
    me$through.cv <- T
    me$through.sd <- F
    me$range <- sort(cv.minmax)
  }
  
  if (me$any)
  { 
    me$init <- mean(me$range)
    me$known <- diff(me$range) == 0
  }
  
  me
} # end if any.me


cond.values <- function(cond, T.values, F.values)
{
  if (cond)
  {
    T.values
  }
  else
  {
    F.values
  }
} # end of cond.values



Default.inits <- function(data, logNormal.distrn, mu.range, sigma.range, include.censored.data=F)
{
  # modif_0.13
  # added argument include.censored.data to this fct

  sigma.lower <- min(sigma.range)
  sigma.upper <- max(sigma.range)

  # new_0.12
  normalized.y <- numeric(0)

  if (data$size$y > 0)
  {
    if (logNormal.distrn) normalized.y <- data$log$y
    else normalized.y <- data$y
  }


  # new_0.12

  if (include.censored.data)
  {
    if (data$size$gt > 0)
    {
      if (logNormal.distrn) normalized.y <- c(normalized.y, data$log$gt)
      else normalized.y <- c(normalized.y, data$gt)
    }

    if (data$size$lt > 0)
    {
      if (logNormal.distrn) normalized.y <- c(normalized.y, data$log$lt)
      else normalized.y <- c(normalized.y, data$lt)
    }

    if (data$size$i > 0)
    {
      if (logNormal.distrn) normalized.y <- c(normalized.y, data$log$i)
      else normalized.y <- c(normalized.y, data$i)
    }
  }

    
  if (length(normalized.y) > 0)
  {
    mu <- mean(normalized.y)
      
    mu.lower <- min(mu.range)
    if (mu < mu.lower) mu <- mu.lower
    else
    {
      mu.upper <- max(mu.range)
      if (mu > mu.upper) mu <- mu.upper
    }
    
    
    # modif_0.13

    if (length(normalized.y) > 1)
    {
      sigma <- sqrt(var(normalized.y))
    }
    else if (!include.censored.data)
    {
      # new_0.12
      # Call this fct again, this time using censored data
      tmp <- Default.inits(data, logNormal.distrn, mu.range, sigma.range, include.censored.data=T)
      sigma <- tmp$sigma
    }
    else
    {
      # new_0.13
      sigma <- 0 # will be corrected below
    }

    if (sigma < sigma.lower) sigma <- sigma.lower
    if (sigma == 0) sigma <- sigma.upper/10 # modif_0.14
  }
  else 
  {
    mu <- mean(mu.range)
    sigma <- ifelse(sigma.lower > 0, sigma.lower, sigma.upper/10)     
  }

  list(mu=mu, sigma=sigma)
} # end of Default.inits


Dens.approx <- function(x, o, return.remote.start=F, log.f.ratio.remote=numeric(0))
{
  remote.start <- c(NA, NA)
  approx <- list()
  mode <- NA
  error <- F
  
  logfp <- o$log.f.prime(x, o$A)
  logfs <- o$log.f.second(x, o$A)
      

  if (o$dens.approx == 1)
  {
    # Normal distrn approx
    
    error <- logfs > 0
    
    if (!error)
    {
      sigma <- sqrt(-1/logfs)
      mu <- x + logfp*sigma**2
      mode <- mu
      
      approx <- list(mu=mu, sigma=sigma)
      
      if (return.remote.start) remote.start <- mu + c(-1,1)*sigma*sqrt(-2*log.f.ratio.remote)
    }
  }
  else if (o$dens.approx == 2)
  {
    # log-Normal distrn approx
    
    lambda <- -logfs*x**2 - x*logfp
    error <- lambda <= 0 | logfs > 0  # when logfs > 0, the approx is so bad that we prefer to ignore it
    
    if (!error)
    {
      mu <- log(x) + (x*logfp + 1)/lambda
      sigma <- 1/sqrt(lambda)
      mode <- exp(mu - sigma**2)
      
      approx <- list(mu=mu, sigma=sigma)
      
      #if (return.remote.start) remote.start <- qlnorm(c(0.001, 0.999), mu, sigma)
      
      if (return.remote.start)
      {
        logf.mode <- o$log.f(mode, o$A)
        logfs.mode <- o$log.f.second(mode, o$A)
        lambda <- -logfs.mode*mode**2
        A <- -lambda/2
        B <- lambda*mu-1
        C <- -lambda/2*mu**2 - logf.mode - log.f.ratio.remote
        delta <- B**2 - 4*A*C
        remote.start <- exp((-B+c(1,-1)*sqrt(delta))/2/A)
      }
    }
  }
  else if (o$dens.approx == 3)
  {
    # Beta distrn approx
    
    error <- x == 0 | x == 1
    
    if (!error)
    {
      B <- -1/(1-x)**2 - 1/x/(1-x)
      b <- (logfs + logfp/x) / B
      a <- x * (logfp + b/(1-x))
    
      error <- a < 0 | b < 0
    }

    if (!error)
    {
      alpha <- a + 1
      beta  <- b + 1
      mode <- (alpha - 1) / (alpha + beta - 2)

      approx <- list(alpha=alpha, beta=beta)
    
      #if (return.remote.start) remote.start <- qbeta(c(0.001, 0.999), alpha, beta)
      
      if (return.remote.start)
      {
        logf.mode <- a*log(a) + b*log(b) - (a+b)*log(a+b)
        target <- logf.mode + log.f.ratio.remote
        remote.start <- c(exp(target/a), 1-exp(target/b))
        
        if (remote.start[1] == 0)
        {
          tmp <- target
          logfp <- o$log.f.prime(0, o$A)
          cont <- is.infinite(logfp)
          
          while (cont)
          {
            tmp <- tmp + 1
            remote.start[1] <- exp(tmp/a)
            logfp <- o$log.f.prime(remote.start[1], o$A)
            cont <- is.infinite(logfp)
          }
        }
        
        if (remote.start[2] == 1)
        {
          tmp <- target
          logfp <- o$log.f.prime(1, o$A)
          cont <- is.infinite(logfp)
          
          while (cont)
          {
            tmp <- tmp + 1
            remote.start[2] <- 1-exp(tmp/b)
            logfp <- o$log.f.prime(remote.start[2], o$A)
            cont <- is.infinite(logfp)
          }
        }          
      }
    }
  }
  else if (o$dens.approx == 4)
  {
    # square-root inverted gamma
    
    beta <- -logfp/4*x**3 -logfs/4*x**4
    a <- 2*beta/x**2 - x*logfp
    
    error <- a < 1 | beta < 0
    
    if (!error)
    {
      alpha <- (a-1)/2
      mode <- sqrt(beta/(alpha+1/2))
      approx <- list(alpha=alpha, beta=beta)

      if (return.remote.start) 
      {
         logf.mode <- o$log.f(mode, o$A)
         target <- logf.mode + log.f.ratio.remote
         remote.start <- c(sqrt(-beta/target), exp(-target/a))       
      }
    }
  } 
  
  
  if (length(approx) == 0)
  {
    approx <- list(NA, NA)
    if      (o$dens.approx == 1 | o$dens.approx == 2) approx.names <- c('mu', 'sigma')
    else if (o$dens.approx == 3 | o$dens.approx == 4) approx.names <- c('alpha', 'beta')
    
    names(approx) <- approx.names
  }
  
  list(remote.start=remote.start, approx=approx, mode=mode, error=error)
} # end of Dens.approx
  

empty.matrices <- function(n.iter, n.chains, me)
{
  # n.iter & n.chains: scalars
  # me: list
  
  template <- matrix(NA, nrow=n.chains, ncol=n.iter)
  
  out <- list(mu=template, sd=template)
  if (me$any & !me$known) out$me.parm <- template
    
  out
} # end of empty.matrices


logp.from.logpcum <- function(log.pcum)
{
  # log.pcum: vector of length 2

  if (log.pcum$lower.tail)
  {
    j <- 1
    s <- -1
  }
  else
  {
    j <- length(log.pcum$log.pcum)
    s <- 1
  }

  log.pcum$log.pcum[-j] + log(1-exp(s*diff(log.pcum$log.pcum))) # scalar output

} # end of logp.from.logpcum


logPhi.quadratic.approx.coeff <- function(degree=2)
{
  theta <- c(-log(2), sqrt(2/pi), -1/pi)
  if (degree == 3) theta <- c(theta, (4/pi-1)*sqrt(2/pi)/6)
  theta
} # end of logPhi.quadratic.approx.coeff


lphi <- function(z)
{
  # z: can be scalar or vector
  
  log.phi <- dnorm(z, log=T)      # log(phi(z))
  log.Phi <- pnorm(z, log.p=T)    # log(Phi(z))
  
  r  <- exp(log.phi-log.Phi)      # phi(z)/Phi(z)
  r2 <- exp(2*(log.phi-log.Phi))  # phi^2(z) / Phi^2 (z)
  
  list(r=r, r2=r2) # r and r2 are of same length as z
} # end of lphi


me.gen <- function(o, me, data, gen.y, true.values, RData=list(save=F))
{
  done <- F 

  # new_0.18
  if (RData$save) save.objects <- paste(RData$dir, '_me.RData', sep='/')
  else save.objects <- character(0)
  
    
  if (me$through.cv)
  {
    if (o$logNormal.distrn)
    {
      log.y <- c(data$log$y, gen.y$log$gt, gen.y$log$lt, gen.y$log$i)
      log.t <- c(true.values$log$y, true.values$log$gt, true.values$log$lt, true.values$log$i)
      b <- sum((exp(log.y-log.t) - 1)^2) / 2 
    }
    else
    {
      y <- c(data$y, gen.y$gt, gen.y$lt, gen.y$i)
      t <- c(true.values$y, true.values$gt, true.values$lt, true.values$i)
      b <- sum((y/t-1)^2) / 2 
    }
  
    A <- c(o$A, list(b=b, range=me$range))
    me.parm <- dens.gen.icdf(o, A, range=me$range, start=o$start(A), inestimable.lower.limit=me$range[1]==0, save.objects=save.objects)
    done <- T
  }
  else 
  {
    y <- c(data$y, gen.y$gt, gen.y$lt, gen.y$i)
    t <- c(true.values$y, true.values$gt, true.values$lt, true.values$i)
    b <- sum((y-t)^2) / 2   
  
    if (o$logNormal.distrn)
    {
      A <- c(o$A, list(b=b, true.values=t, range=me$range))
      me.parm <- dens.gen.icdf(o, A, range=me$range, inestimable.lower.limit=me$range[1]==0, save.objects=save.objects)
      done <- T
    }
  }


  if (!done) me.parm <- sqrt.invertedGamma.gen(data$N, b, me$range, o=o, RData=RData)
  
  me.parm
} # end of me.gen


me.gen.object <- function(me, logNormal.distrn, N)
{
  if (me$through.cv)
  {  
    A <- list(N=N, M=0)
    f <- function(v, A){exp(-A$N*log(v) - A$b/v^2 - A$N*pnorm(1/v, log.p=T) - A$M)}
    log.f <- function(v, A){-A$N*log(v) - A$b/v^2 - A$N*pnorm(1/v, log.p=T)}
    log.f.prime  <- function(v, A){l=lphi(1/v); -A$N/v + 2*A$b/v^3 + A$N*l$r/v^2}     
    log.f.second <- function(v, A){l=lphi(1/v); A$N/v^2 - 6*A$b/v^4 + A$N/v^4 * (l$r*(1/v - 2*v) + l$r2)}
    
      
    # new_0.11
    log.f.inv.remote <- function(target, A)
    {
      D <- -A$b
      B <- -A$N*(log(A$mode)-1) - target
      qA <- - A$N/A$mode
      roots <- real.cubic.roots(c(D, 0, B, qA), l=0, u=A$mode) # modif_0.11

      roots <- c(roots, exp(log(2)-target/A$N))

      roots
    } # end of log.f.inv.remote
    
    
    start <- function(A)
    {
      start <- sqrt(2*A$b/A$N)
      if (start < A$range[1] | start > A$range[2]) start <- mean(A$range)
      start
    } # end of start
    
        
    # modif_0.11
    o <- list(A=A, f=f, log.f=log.f, log.f.prime=log.f.prime, log.f.second=log.f.second,
              log.f.inv.remote=log.f.inv.remote,
              potentially.bimodal=F, math.lower.limit=0,
              start=start, logNormal.distrn=logNormal.distrn)
  }
  else if (logNormal.distrn)
  {
    A <- list(N=N, M=0)
    f <- function(ksi, A){exp(-A$N*log(ksi) - as.vector(pnorm(matrix(A$true.values, byrow=T, nrow=length(ksi), ncol=A$N)/ksi, log.p=T)%*%matrix(1,ncol=1,nrow=A$N)) - A$b/ksi^2 - A$M)}
    log.f <- function(ksi, A){-A$N*log(ksi) - sum(pnorm(A$true.values/ksi, log.p=T)) - A$b/ksi^2}
    log.f.prime <- function(ksi, A){l=lphi(A$true.values/ksi); -A$N/ksi + sum(l$r*A$true.values/ksi^2) + 2*A$b/ksi^3}  # modif_0.11
    log.f.second <- function(ksi, A){l=lphi(A$true.values/ksi); A$N/ksi^2 + sum(A$true.values*l$r*(A$true.values^2 + ksi*l$r*A$true.values - 2*ksi^2))/ksi^5  - 6*A$b/ksi^4} # modif_0.11
    

    # new_0.11
    log.f.inv.remote <- function(target, A)
    {      
      theta <- c(-A$b, 0, 3*A$N/2 - target, -2*A$N)
      roots <- real.cubic.roots(theta, l=0)
      
      roots <- c(roots, 2*exp(-target/A$N))
      roots
    } # end of log.f.inv.remote

    # modif_0.11
    o <- list(A=A, f=f, log.f=log.f, log.f.prime=log.f.prime, log.f.second=log.f.second,
              log.f.inv.remote=log.f.inv.remote,
              potentially.bimodal=F, math.lower.limit=0,
              logNormal.distrn=logNormal.distrn)
  }
  else
  {
    o <- sigma.gen.object(N)
    o$logNormal.distrn <- logNormal.distrn
  }

  o
} # end of me.gen.object


mu.truncatedData.gen <- function(o, range, mu.mean, sigma, prec=o$A$N/sigma^2, RData=list(save=F), start=mean(range))
{
  # o: list
  # range: vector of length 2
  # mu.mean & sigma: scalars
  
  # new_0.24.1 Added argument prec  to the function
  # new_0.24.2 Added argument start to the function

  # new_0.18
  if (RData$save) save.objects <- paste(RData$dir, '_muTruncatedData.RData', sep='/')
  else save.objects <- character(0)
  
  A <- c(o$A, list(mu.mean=mu.mean, s=sigma, s2=sigma^2, prec=prec)) # modif_0.24.1
  inestimable.lower.limit <- is.infinite(range[1]) # new_0.24.2
      
  mu <- dens.gen.icdf(o, A, range=range, 
                      start=start, inestimable.lower.limit=inestimable.lower.limit, # modif_0.24.2 (added start and inestimable.lower.limit)
                      save.objects=save.objects) 
  mu
}  # end of mu.truncatedData.gen


mu.truncatedData.gen.object <- function(N)
{
  A <- list(N=N, M=0) # list of arguments of f
  
  # modif_0.24 The four lines below were modified
  # [A$N/A$s^2 was changed to A$prec]
  # See next 4-lines block
  #f <- function(mu, A){exp(-A$N/2*((mu-A$mu.mean)/A$s)^2 - A$N*pnorm(mu/A$s, log.p=T) - A$M)}
  #log.f <- function(mu, A){-A$N/2*((mu-A$mu.mean)/A$s)^2 - A$N*pnorm(mu/A$s, log.p=T)}
  #log.f.prime  <- function(mu, A){z=mu/A$s; l=lphi(z); -A$N*(mu-A$mu.mean)/A$s2 - A$N*l$r/A$s}
  #log.f.second <- function(mu, A){z=mu/A$s; l=lphi(z); A$N/A$s2 * (z*l$r + l$r2 - 1)}
  
  # new_0.24
  f <- function(mu, A){exp(-A$prec/2*(mu-A$mu.mean)^2 - A$N*pnorm(mu/A$s, log.p=T) - A$M)}
  log.f <- function(mu, A){-A$prec/2*(mu-A$mu.mean)^2 - A$N*pnorm(mu/A$s, log.p=T)}
  log.f.prime  <- function(mu, A){z=mu/A$s; l=lphi(z); -A$prec*(mu-A$mu.mean) - A$N*l$r/A$s}
  log.f.second <- function(mu, A){z=mu/A$s; l=lphi(z); -A$prec + A$N/A$s2 * (z*l$r + l$r2)} # modif_0.24 (this line was more modified than the above three...)
  
  # new_0.11
  log.f.inv.remote <- function(target, A)
  {
    k <- -A$prec/2  # modif_0.24
    
    # right-side roots
    B <- -2*A$mu.mean
    C <- A$mu.mean^2
    coeff <- k * c(C, B, 1)
    roots <- quadratic.solution(coeff, target=target)  
    
    
    # left-side roots
    coeff <- coeff - A$N * logPhi.quadratic.approx.coeff() / c(1, A$s, A$s2)
    tmp <- quadratic.solution(coeff, target=target)
    roots <- c(roots, tmp)
     
    roots
  } # end of log.f.inv.remote
  
  
  o <- list(A=A, f=f, log.f=log.f, log.f.prime=log.f.prime, log.f.second=log.f.second,
            log.f.inv.remote=log.f.inv.remote,
            logNormal.distrn=F, # new_0.24.1
            potentially.bimodal=F, math.lower.limit=-Inf)
            
  o
} # end of mu.truncatedData.gen.object


out.logout.moments <- function(any.me, logNormal.distrn, data, gen.y, true.values)
{
  if (any.me)
  {
    if (logNormal.distrn)
    {
      y <- c(true.values$log$y, true.values$log$gt, true.values$log$lt, true.values$log$i)
    }
    else
    {
      y <- c(true.values$y, true.values$gt, true.values$lt, true.values$i)
    }
    
    out <- list(sum=sum(y), sum2=sum(y^2))
  }
  else 
  {
    if (logNormal.distrn)
    {
      out <- list(sum=data$log$uncensored.sum, sum2=data$log$uncensored.sum2)
      y <- c(gen.y$log$gt, gen.y$log$lt, gen.y$log$i)
    }
    else
    {
      out <- list(sum=data$uncensored.sum, sum2=data$uncensored.sum2)
      y <- c(gen.y$gt, gen.y$lt, gen.y$i)
    }
    
    out <- list(sum=out$sum + sum(y), sum2=out$sum2 + sum(y^2))
  }

  out$mean <- out$sum / data$N # new_0.24
  out   
} # end of out.logout.moments


out.sample <- function(sample, burnin, n.chains, monitor.burnin, me, logNormal.distrn)
{    
  if (me$any & !me$known)
  {
    sample <- renamed.me.parm(sample, me$through.cv)
    if (monitor.burnin) burnin <- renamed.me.parm(burnin, me$through.cv)
  }


  if (logNormal.distrn)
  {
    sample$gm  <- exp(sample$mu)
    sample$gsd <- exp(sample$sd)
  }
  
  if (n.chains == 1) sample <- lapply(sample, as.vector)
  out <- list(sample=sample)
  
  
  if (monitor.burnin)
  {    
    if (logNormal.distrn)
    {
      burnin$gm  <- exp(burnin$mu)
      burnin$gsd <- exp(burnin$sd)
    }
    
    if (n.chains == 1) burnin <- lapply(burnin, as.vector)
    out$burnin <- burnin
  }
  
  
  out
} # end of out.sample


pnorm.logpcum <- function(x, mean=0, sd=1)
{
  # x: can be vector or scalar
 
  if (length(x) > 1)
  {
    z <- c(min(x), max(x))
    z <- abs(z - mean)
    lower.tail <- z[1] > z[2]
  }
  else
  {
    lower.tail <- x < mean
  }

  log.pcum <- pnorm(x, mean=mean, sd=sd, log.p=T, lower.tail=lower.tail)
  list(log.pcum=log.pcum, lower.tail=lower.tail)
} # end of pnorm.logpcum


polynomial.product.coeff <- function(theta1, theta2)
{
  len1 <- length(theta1)
  len2 <- length(theta2)
  len <- len1 + len2 - 1
  coeff <- rep(0, len)
  
  for (i in seq(len1))
  {
    tmp <- theta1[i] * theta2
    w <- seq(from=i, to=i+len2-1)
    coeff[w] <- coeff[w] + tmp
  }
  
  coeff
} # end of polynomial.product.coeff


Quadratic.solution <- function(A, B, C)
{
  # Note: there is another version with a similar name -- quadratic.solution -- below
  delta <- B**2 - 4*A*C
  negative.delta <- delta < 0
  symmetry.axis <- -B/(2*A)
  
  if (negative.delta)
  {
    left <- NA
    right <- NA
  }
  else
  {
    m <- sqrt(delta)/(2*abs(A))
    left  <- symmetry.axis - m
    right <- symmetry.axis + m
  }
  
  list(left=left, right=right, negative.delta=negative.delta, symmetry.axis=symmetry.axis)
} # end of Quadratic.solution


# modif_0.11 (added target)
quadratic.solution <- function(theta, target=0, l=-Inf, u=Inf)
{
  # Note: there is another version with a similar name -- Quadratic.solution -- above
  
  C <- theta[1] - target
  B <- theta[2]
  A <- theta[3]
  
  delta <- B^2 - 4*A*C
  
  if (delta < 0)
  {
    soln <- NA
  }
  else
  {
    soln <- (-B + c(-1,1)*sqrt(delta))/2/A
    soln <- soln[soln > l & soln < u]
    if (length(soln) == 0) soln <- NA
  }
  
  soln
} # end of quadratic.solution


random.pow <- function(a, range)
{
  # sample a value from f(x) = 1/x^a on the range specified
  
  U <- runif(1)
  
  if (a==1)
  {
    if (is.infinite(range[2])) range[2] <- 1e8 
    log.z <- diff(c(U-1,U)*log(range))
    z <- exp(log.z)
  }
  else
  {
    fcum <- range^(1-a)
    tmp <- sum(c(1-U,U)*fcum) 
    z <- tmp^(1/(1-a))
  }
  
  z
} # end of random.pow


real.cubic.roots <- function(theta, target=0, l=numeric(0), u=numeric(0))
{
  # Find real cubic roots of   theta[1] + theta[2]*x + theta[3]*x^2 + theta[4]*x^3  = target
  
  # Returns an empty vector rather than NA when no solution is found [modif_0.11]
  
  # *** IMPORTANT ***
  # this function is based on the built-in function polyroot, which solves a polynomial of degree n (here of third degree):
  # if such a function is not accessible in your programming language, you could use the
  # Newton-Raphson-based algorithm alt.real.cubic.roots, found in alt-realCubicRoots.R
  
  theta[1] <- theta[1] - target
  solns <- polyroot(theta)
  w <- which(abs(Im(solns)) < 1e-6) # indicate which solns are real (with negligeable imaginary part)
  
  if (length(w) > 0)
  {
    solns <- Re(solns[w])
    
    # reject out-of-bounds solutions

    if (length(l) > 0)
    {
      w <- which(solns > l)
      solns <- solns[w]
    }    
    
    if (length(u) > 0)
    {
      w <- which(solns < u)
      solns <- solns[w]
    }     
  }
  else
  {
    solns <- numeric(0)
  }
  
  solns
} # end of real.cubic.roots


renamed.me.parm <- function(obj, through.cv)
{
  # rename me.parm either me.sd or me.cv
  obj.names <- names(obj)
  w <- which(obj.names=="me.parm")
  new.name <- ifelse(through.cv, "me.cv", "me.sd") 
  names(obj)[w] <- new.name
  obj
} # end of renamed.me.parm


rgamma.truncated <- function(alpha, beta, range)
{
  # We first look on which side of the distribution (comparing to its mode)
  # the range lies, in order to decide on which side we will compute the log.prob
  # (choosing the correct side gives more accuracy in the sampling in the eventuality of a remote range)
  
  if (alpha < 1)
  {
    lower.tail <- F
  }
  else
  {
    mode <- (alpha-1)/beta
    
    min.range <- min(range)
    max.range <- max(range)
    
    if (min.range > mode)
    {
      lower.tail <- F
    }
    else if (max.range < mode)
    {
      lower.tail <- T
    }
    else
    {
      logp.left  <- pgamma(min.range, alpha, beta, log.p=T, lower.tail=T)
      logp.right <- pgamma(max.range, alpha, beta, log.p=T, lower.tail=F)
      lower.tail <- logp.left < logp.right
    }
  }

  logp.lim <- pgamma(range, alpha, beta, log.p=T, lower.tail=lower.tail)
  logp <- runif.logp(logp.lim, lower.tail=lower.tail)
  qgamma(logp, alpha, beta, log.p=T, lower.tail=lower.tail)
} # end of rgamma.truncated


rnorm.censored <- function(mu, sd, lower=numeric(0), upper=numeric(0), negative.values.disallowed=F)
{
  # Useful functions
  
  rnorm.interval.censored <- function(mu, sd, lower, upper)
  {
    # mu, sd: same length (1, or same as length(lower))
    # lower and upper: same length, but length can be > 1
    
    log.p.lower <- pnorm(lower, mean=mu, sd=sd, log.p=T)
    log.p.upper <- pnorm(upper, mean=mu, sd=sd, log.p=T)
    
    l <- length(lower)
    log.p <- numeric(l)
    
    for (i in seq(l)) log.p[i] <- runif.logp(c(log.p.lower[i], log.p.upper[i]))
    y <- qnorm(log.p, mean=mu, sd=sd, log.p=T)
    y
  } # end of rnorm.interval.censored
  
  
  rnorm.left.censored <- function(mu, sd, lower)
  {
    -rnorm.right.censored(-mu, sd, -lower)
  } # end of rnorm.left.censored
  

  rnorm.right.censored <- function(mu, sd, upper)
  {
    # This function was incorrectly named rnorm.left.censored in previous versions
    
    # mu, sd: same length (1, or same as length(upper))
    
    log.p <- log(runif(length(upper))) + pnorm(upper, mean=mu, sd=sd, log.p=T)
    y <- qnorm(log.p, mean=mu, sd=sd, log.p=T)
    y
  } # end of rnorm.right.censored

  
  # Beginning of calculations
  
  
  # new_0.4
  # (note that negative.values.disallowed was added to fct arguments)

  if (negative.values.disallowed)
  {
    if (length(lower) == 0)
    {
      lower.len <- max(1, length(upper)) 
      lower <- rep(0, lower.len)
    }
  }
  
  left.censored <- length(lower) > 0 
  
  
  if (left.censored)
  {
    right.censored <- length(upper) > 0
    
    if (right.censored)
    {
      # it is interval-censored
      z <- rnorm.interval.censored(mu, sd, lower, upper)
    }
    else
    {
      # is is left-censored
      z <- rnorm.left.censored(mu, sd, lower)
    }
  }
  else
  {
    # then it is right-censored
    z <- rnorm.right.censored(mu, sd, upper)
  }
  
  z
} # end of rnorm.censored


runif.logp <- function(logp.lim, lower.tail=T, size=1, u=runif(size))
{
  # To sample 'size' values uniformly in c(exp(logp.lim[1]), exp(logp.lim[2]))
  w <- diff(logp.lim)
  if (!lower.tail) w <- abs(w)
  log.p <- log(u) + pexp(w, log.p=T)
  x <- qexp(log.p, log.p=T)
  max(logp.lim) - x
} # end of runif.logp



sigma.gen.object <- function(N, lnorm.mu=numeric(0), lnorm.sd=numeric(0))
{
  if (length(lnorm.mu) > 0)
  {
    A <- list(N=N, lm=lnorm.mu, ls2=lnorm.sd^2, M=0)
    f <- function(sigma, A){exp(-(A$N+1)*log(sigma) - A$b/sigma^2 - (log(sigma)-A$lm)^2/(2*A$ls2) - A$M)}
    log.f <- function(sigma, A){-(A$N+1)*log(sigma) - A$b/sigma^2 - (log(sigma)-A$lm)^2/(2*A$ls2)}
    log.f.prime  <- function(sigma, A){-(A$N+1)/sigma + 2*A$b/(sigma^3) - (log(sigma)-A$lm)/(sigma*A$ls2)}
    log.f.second <- function(sigma, A){(A$N+1)/sigma^2 - 6*A$b/sigma^4 + (log(sigma)-A$lm-1)/(A$ls2*sigma^2)}
    
    # new_0.11
    log.f.inv.remote  <- function(target, A)
    {
      # small values
      c3 <- - (A$N + 1) + A$lm/A$ls2
      c2 <- A$N + 1 - A$lm^2/2/A$ls2 - target
      c0 <- -A$b
      roots <- real.cubic.roots(c(c0, 0, c2, c3), l=0)

      # large values
      C <- -A$lm^2/2/A$ls2 - target
      B <- A$lm/A$ls2 - (A$N + 1)
      A <- -1/2/A$ls2
      tmp <- quadratic.solution(c(C, B, A), l=0)
      roots <- c(roots, exp(tmp))

      roots
    } # end of log.f.inv.remote
    
    
    start <- function(A)
    {
      start <- numeric(0)
      tmp <- 2*A$b/(A$N+1-A$lm/A$ls2)
      if (tmp > 0) start <- sqrt(tmp)
      start <- c(exp(A$lm), start)
      start
    } # end of start
  }
  else
  {
    A <- list(N=N, M=0)
    f <- function(sigma, A){exp(-A$N*log(sigma) - A$b/sigma^2 - A$M)}
    log.f <- function(sigma, A){-A$N*log(sigma) - A$b/sigma^2}
    log.f.prime  <- function(sigma, A){-A$N/sigma + 2*A$b/sigma^3}
    log.f.second <- function(sigma, A){A$N/sigma^2 - 6*A$b/sigma^4}
           
    # new_.11
    log.f.inv.remote  <- function(target, A)
    {
      c3 <- -A$N
      c2 <- A$N - target
      c0 <- -A$b 
      roots <- real.cubic.roots(c(c0, 0, c2, c3), l=0)
      
      roots <- c(roots, -target/A$N)
      roots
    } # end of log.f.inv.remote
    
    
    start <- function(A){sqrt(2*A$b/A$N)}
  }
  
  o <- list(f=f, log.f=log.f, log.f.prime=log.f.prime, log.f.second=log.f.second, A=A,
            log.f.inv.remote=log.f.inv.remote,
            potentially.bimodal=F, start=start, math.lower.limit=0)
            
  o
} # end of sigma.gen.object


sigma.truncatedData.gen <- function(o, range, b, mu, current.sigma, RData=list(save=F))
{  
  # new_0.18
  if (RData$save) save.objects <- paste(RData$dir, '_sigmaTruncatedData.RData', sep='/')
  else save.objects <- character(0)

  A <- c(o$A, list(b=b, mu=mu))
  start <- o$start(A)
  if (length(start) == 0) start <- current.sigma   
  sigma <- dens.gen.icdf(o, A, range=range, start=start, inestimable.lower.limit=range[1]==0, save.objects=save.objects)
  sigma
} # end of sigma.truncatedData.gen 


sigma.truncatedData.gen.object <- function(N)
{
  A <- list(N=N, M=0)
  f <- function(s, A){exp(-A$N*log(s) - A$b/s^2 - A$N*pnorm(A$mu/s,log.p=T) - A$M)}
  log.f <- function(s, A){-A$N*log(s) - A$b/s^2 - A$N*pnorm(A$mu/s,log.p=T)}
  log.f.prime  <- function(s, A){z=A$mu/s; l=lphi(z); -A$N/s + 2*A$b/s^3 + A$N*z*l$r/s}    
  log.f.second <- function(s, A){z=A$mu/s; l=lphi(z); A$N/s^2 - 6*A$b/s^4 + A$N*z*((z^2-2)*l$r+z*l$r2)/s^2}
  
  
  # new_.11
  log.f.inv.remote <- function(target, A)
  {
    c3 <- -A$N 
    c2 <- A$N - target
    c0 <- -A$b
    
    t <- logPhi.quadratic.approx.coeff(degree=3)
    t <- -A$N * t * c(1, A$mu, A$mu^2, A$mu^3)
    theta <- c(c0, 0, c2, c3) + t    
    roots <- real.cubic.roots(theta, l=0)

    roots <- c(roots, 2*exp(-target/A$N))
    roots
  } # end of log.f.inv.remote
  

  
  # modif_0.11
  start <- function(A)
  {
    s0 <- sqrt(2*A$b/A$N)
    g0 <- pnorm(A$mu/s0, log.p=T)
    gp0 <- -dnorm(A$mu/s0) *A$mu/s0^2 / pnorm(A$mu/s0)
    c3 <- A$N * (gp0 - 1)
    c2 <- A$N * (1 - g0 + gp0*s0)
    c0 <- -A$b
    theta <- c(c0, 0, c2, c3)
    tmp <- real.cubic.roots(theta, l=0)
    if (length(tmp) == 0) tmp <- s0
    tmp
  } # end of start
      
      
  o <- list(A=A, f=f, log.f=log.f, log.f.prime=log.f.prime, log.f.second=log.f.second,
            log.f.inv.remote=log.f.inv.remote,
            start=start, potentially.bimodal=F, math.lower.limit=0)
  o
} # end of sigma.truncatedData.gen.object


sqrt.invertedGamma.gen <- function(n, beta, xrange, o=list(), beta.min=1e-8, RData=list(save=F))
{
  # Sample a value from the distrn 
  #
  # f(x) =     1
  #         ------- . exp(-beta/x^2)
  #           x^n

  # new_0.18
  if (RData$save) save.objects <- paste(RData$dir, '_sqrtInvGamma.RData', sep='/')
  else save.objects <- character(0)

  if (n <= 1)
  {
    A <- c(o$A, list(b=beta))
    x <- dens.gen.icdf(o, A, range=xrange, start=sqrt(2*beta), inestimable.lower.limit=xrange[1]==0, save.objects=save.objects)
  }
  else if (beta < beta.min)
  {
    if (xrange[1] == 0) xrange[1] <- 0.0001
    x <- random.pow(n, xrange)
  }
  else
  {
    # If X ~ f, then Y = X^{-2} ~ Gamma((N-1)/2, b)
    #
    
    alpha <- (n-1)/2
    invx2.range <- 1/rev(xrange^2)
    invx2 <- rgamma.truncated(alpha, beta, invx2.range) 
    x <- 1/sqrt(invx2)
  }
        
  x
} # end of sqrt.invertedGamma.gen


take.log <- function(obj, include.ydim=F)
{
  out <- list()
  
  if (length(obj$gt) > 0) out$gt <- log(obj$gt)
  if (length(obj$lt) > 0) out$lt <- log(obj$lt)
  if (length(obj$i) > 0)  out$i  <- log(obj$i)

  if (include.ydim & length(obj$y) > 0) out$y <- log(obj$y)
  
  out
} # end of take.log


truevalue.gen <- function(o, me, y, mu, u=runif(1), log.y=numeric(0), RData=list(save=F))
{ 
  # o:     list
  # me:    list
  # mu:    scalar
  # y:     scalar

  # modif_0.16 added argument RData
  
  A <- o$A
  A$mu <- mu
  
  if (o$logNormal.distrn) A$log.y <- log.y

  A$y  <- y
  A$y2 <- y^2
  
  
  if (o$through.sd)
  {
    A$ksi  <- me$parm
    A$ksi2 <- me$parm^2
  }
  else
  {
    A$cv2 <- me$parm^2
  }
  
  
  start <- o$start(A) # can be a scalar or a vector of length 2
    
  # new_0.16
  if (RData$save) save.objects <- paste(RData$dir, '_tv.RData', sep='/')
  else save.objects <- character(0)

  x <- dens.gen.icdf(o, A, range=o$range, u=u, start=start, inestimable.lower.limit=!o$logNormal.distrn, save.objects=save.objects) # modif_0.11, modif_0.16
  
  if (o$logNormal.distrn) x <- exp(x) # new_0.11
  
  x
} # end of truevalue.gen 


truevalue.gen.object <- function(me, logNormal.distrn)
{
  A <- list(M=0)
  
     
  if (me$through.sd)
  {
    # modif_0.11
    # ME through SD and T ~ logNormal
    
      # We will sample a value for s = log(t)  [rather than for t]
      # and the functions below will therefore be expressed in terms of s, whose density function is given by:
      # dens(s) = f(t(s)) . |Jacobean|
      #         = f(exp(s)) exp(s) 
    
    range <- c(-Inf, Inf)
    
    f <- function(s, A){exp(- pnorm(exp(s)/A$ksi, log.p=T) - (A$y-exp(s))^2/2/A$ksi2 - (s-A$mu)^2/2/A$sigma2 - A$M)}
    log.f <- function(s, A){- pnorm(exp(s)/A$ksi, log.p=T) - (A$y-exp(s))^2/2/A$ksi2 - (s-A$mu)^2/2/A$sigma2} 
    log.f.prime <- function(s, A){z <- exp(s)/A$ksi; l <- lphi(z); -l$r*z + A$y*z/A$ksi - z^2 - (s-A$mu)/A$sigma2}
    log.f.second <- function(s, A){z <- exp(s)/A$ksi; l <- lphi(z); l$r*(z^3-z) + (l$r*z)^2 + (A$y*exp(s)-2*exp(2*s))/A$ksi2 - 1/A$sigma2}
      

    log.f.inv.remote <- function(target, A)
    {
      # solutions for large s
      tmp <- A$y + c(-1,1)*sqrt(2*abs(target)*A$ksi2) # vector of length 2
      tmp <- tmp[tmp>0]
      roots <- log(tmp)
      
      # solutions for small s (large negative values)
      C <- log(2) -A$y2/2/A$ksi2 - A$mu^2/2/A$sigma2
      B <- A$mu/A$sigma2
      qA <- -1/2/A$sigma2
      tmp <- quadratic.solution(c(C, B, qA), target=target)
      roots <- c(roots, tmp)
      
      # solutions around 0
      C <- log(1/A$ksi) - (A$y-1)^2/2/A$ksi2 - A$mu^2/2/A$sigma2
      B <- A$mu/A$sigma2
      qA <- -1/2/A$sigma2
      tmp <- quadratic.solution(c(C, B, qA), target=target)
      roots <- c(roots, tmp)
    
      roots
    } # end of log.f.inv.remote
      

    start <- function(A){c(A$log.y, A$mu)}    
  }
  else if (logNormal.distrn)
  {
    # modif_0.11
    # ME through CV and T ~ logNormal
    
      # We will sample a value for s = log(t)  [rather than for t]
      # and the functions below will therefore be expressed in terms of s, whose density function is given by:
      # dens(s) = f(t(s)) . |Jacobean|
      #         = f(exp(s)) exp(s) 

    range <- c(-Inf, Inf) 
    
    f <- function(s, A){exp(-s - (A$y*exp(-s)-1)^2/2/A$cv2 - (s-A$mu)^2/2/A$sigma2 - A$M)}
    log.f <- function(s, A){-s - (A$y*exp(-s)-1)^2/2/A$cv2 - (s-A$mu)^2/2/A$sigma2}
    log.f.prime <- function(s, A){-1 + (A$y2*exp(-2*s) - A$y*exp(-s))/A$cv2 - (s-A$mu)/A$sigma2}
    log.f.second <- function(s, A){-2*A$y2*exp(-2*s)/A$cv2 + A$y*exp(-s)/A$cv2 - 1/A$sigma2}
    

    start <- function(A)
    {
      C <- -1 + A$y2/A$cv2 - A$y/A$cv2 + A$mu/A$sigma2
      B <- -2*A$y2/A$cv2 + A$y/A$cv2 - 1/A$sigma2
      qA <- 2*A$y2/A$cv2 - A$y/2/A$cv2
      start <- quadratic.solution(c(C, B, qA))

      c3 <- (-4/3*A$y2 + A$y/6) / A$cv2 # new_0.15
      start <- c(start, real.cubic.roots(c(C, B, qA, c3))) # modif_0.17 --- dropped argument l=0

      start <- c(start, A$mu, A$mu-A$sigma2, -B/2/qA, log(A$y)) # new_0.17
         
      start
    } # end of start


    log.f.inv.remote <- function(target, A)
    {
      tmp <- 1 + c(-1,1)*sqrt(2*A$cv2*abs(target))
      tmp <- tmp[tmp>0]
      roots <- log(tmp/A$y)
    
      C <- -(A$y-1)^2/2/A$cv2 - A$mu^2/2/A$sigma2
      B <- -1 + (A$y-1)*A$y/A$cv2 + A$mu/A$sigma2
      qA <- -A$y2/2/A$cv2 - 1/2/A$sigma2
      tmp <- quadratic.solution(c(C, B, qA), target=target)
      
      roots <- c(roots, tmp)

      # new_0.19
      D <- (-A$y2+2*A$y-1)/2/A$cv2 - A$mu^2/2/A$sigma2
      C <- -1 + (A$y2-A$y)/A$cv2 + A$mu/A$sigma2
      B <- (-2*A$y2+A$y)/2/A$cv2 -1/2/A$sigma2
      qA <- (4*A$y2-A$y)/6/A$cv2
      theta <- c(D, C, B, qA)
      tmp <- real.cubic.roots(theta, target)
      roots <- c(roots, tmp)
            
      roots
    } # end of log.f.inv.remote
  }
  else
  {
    # ME through CV and T ~ Normal

    range <- c(0, Inf) 

    f <- function(t, A){exp(-log(t) - (A$y-t)^2/2/A$cv2/t^2 - (t-A$mu)^2/2/A$sigma2 - A$M)}
    log.f <- function(t, A){-log(t) - (A$y-t)^2/2/A$cv2/t^2 - (t-A$mu)^2/2/A$sigma2}       
    log.f.prime <- function(t, A){-1/t - (A$y/t^2-A$y2/t^3)/A$cv2 + (A$mu-t)/A$sigma2}   
    log.f.second <- function(t, A){1/t^2 + (2*A$y/t^3 - 3*A$y2/t^4)/A$cv2 -1/A$sigma2}


    # new_.11
    log.f.inv.remote <- function(target, A)
    {
      roots <- A$y / (1 + sqrt(-2*A$cv2*target))


      # We try a solution with a small t
      D <- - A$y2/2/A$cv2
      C <- A$y/A$cv2
      B <- 3/2 - 1/2/A$cv2 - A$mu^2/2/A$sigma2 - target
      qA <- -2 + A$mu/A$sigma2
      theta <- c(D, C, B, qA)

      tmp <- real.cubic.roots(theta)
      roots <- c(roots, tmp)


      # Look also for a solution with a large value for t
      tmp <- -2*A$sigma2*(target + 1/2/A$cv2)
      if (tmp > 0)
      { 
        root <- A$mu + sqrt(tmp)
        if (root > 0) roots <- c(roots, root)
      }

      roots
    } # end of log.f.inv.remote


    # modif_0.11
    start <- function(A)
    { 
      # solutions in small values for t
      C <- A$y2/A$cv2
      B <- -A$y/A$cv2
      roots.small <- quadratic.solution(c(C, B, -1), l=0)

      # solutions in large values for t
      D <- A$y/A$cv2
      B <- -A$mu/A$sigma2
      qA <- 1/A$sigma2
      theta <- c(D, 1, B, qA)
      roots.large <- real.cubic.roots(theta, l=0)

      roots <- c(roots.small, roots.large)
      roots[!is.na(roots)]
    } # end of start
  }

  
  # modif_0.11
  o <- list(f=f, log.f=log.f, log.f.prime=log.f.prime, log.f.second=log.f.second,
            A=A, range=range, logNormal.distrn=logNormal.distrn, 
            through.sd=me$through.sd, through.cv=me$through.cv,
            log.f.inv.remote=log.f.inv.remote,
            start=start, potentially.bimodal=T, math.lower.limit=range[1])    
  
  o
} # end of truevalue.gen.object        


truevalues.gen <- function(gen.y, data, mu, sigma, me, logNormal.distrn, o=list(), RData=list(save=F))
{
  # gen.y:             list
  # data:              list
  # mu:                scalar
  # sigma:             scalar
  # me:                list
  # logNormal.distrn:  logical scalar

  # new_0.16 added argument RData


  new.truevalues <- list()


  if (me$through.sd & !logNormal.distrn)
  {
    me.sd <- me$parm
    tau.star <- 1/sigma^2 + 1/me.sd^2
    sd.star <- 1/sqrt(tau.star)

    if (data$size$y > 0)
    {
      tmp.mean <- (data$y/me.sd^2 + mu/sigma^2) / tau.star
      new.truevalues$y <- rnorm(data$size$y, tmp.mean, sd.star)
    }

    if (data$size$gt > 0)
    {
      tmp.mean <- (gen.y$gt/me.sd^2 + mu/sigma^2) / tau.star
      new.truevalues$gt <- rnorm(data$size$gt, tmp.mean, sd.star)
    }

    if (data$size$lt > 0)
    {
      tmp.mean <- (gen.y$lt/me.sd^2 + mu/sigma^2) / tau.star
      new.truevalues$lt <- rnorm(data$size$lt, tmp.mean, sd.star)
    }

    if (data$size$i > 0)
    { 
      tmp.mean <- (gen.y$i/me.sd^2 + mu/sigma^2) / tau.star
      new.truevalues$i <- rnorm(data$size$i, tmp.mean, sd.star)
    }          
  }
  else
  {
    o$A$sigma2 <- sigma^2
    
    if (data$size$y > 0)  new.truevalues$y  <- rep(0, data$size$y)
    if (data$size$gt > 0) new.truevalues$gt <- rep(0, data$size$gt)
    if (data$size$lt > 0) new.truevalues$lt <- rep(0, data$size$lt)
    if (data$size$i > 0)  new.truevalues$i  <- rep(0, data$size$i)
      
  
    if (data$size$y > 0)
    {
      for (j in 1:data$size$y)
      {
        new.truevalues$y[j] <- truevalue.gen(o, me, data$y[j], mu, log.y=data$log$y[j], RData=RData) # modif_0.16
      }
    }
    

    if (data$size$gt > 0)
    {
      for (j in 1:data$size$gt)
      {
        new.truevalues$gt[j] <- truevalue.gen(o, me, gen.y$gt[j], mu, log.y=gen.y$log$gt[j], RData=RData) # modif_0.16
      }
    }


    if (data$size$lt > 0)
    {    
      for (j in 1:data$size$lt)
      {
        new.truevalues$lt[j] <- truevalue.gen(o, me, gen.y$lt[j], mu, log.y=gen.y$log$lt[j], RData=RData) # modif_0.16
      }
    }


    if (data$size$i > 0)
    { 
      for (j in 1:data$size$i)
      { 
        new.truevalues$i[j] <- truevalue.gen(o, me, gen.y$i[j], mu, log.y=gen.y$log$i[j], RData=RData) # modif_0.16
      }
    }
    
    
    if (logNormal.distrn) new.truevalues$log <- take.log(new.truevalues, include.ydim=T)
  }

  
  new.truevalues
} # end of truevalues.gen


Xinverse.y <- function(X, y)
{
  R <- nrow(X)
  z <- cbind(X, y)

  for (r in 1:(R-1))
  {
    for (s in (r+1):R) z[s,] <- z[s,] - z[s,r]/z[r,r]*z[r,]
  }
  
  for (r in 1:(R-1))
  {
    for (s in (r+1):R) z[r,] <- z[r,] - z[r,s]/z[s,s]*z[s,]
  }
  
  y <- z[,R+1]
  z <- diag(z)
  log.abs.theta <- log(abs(y))-log(abs(z))
  theta <- exp(log.abs.theta) * sign(y) * sign(z)
  theta  
} # end of Xinverse.y


y.gen <- function(true.values, data, sigma, me, logNormal.distrn, mu=numeric(0))
{
  new.gen.y <- list()
  

  if (me$any)
  {
    # Sample y (censored) values | true values

    if (data$any.censored$gt)
    {
      tmp.mean <- true.values$gt
      tmp.sd <- cond.values(me$through.cv, me$parm*tmp.mean, me$parm)
      new.gen.y$gt <- rnorm.censored(tmp.mean, tmp.sd, lower=data$cens$gt)
    }

    if (data$any.censored$lt)
    {
      tmp.mean <- true.values$lt
      tmp.sd <- cond.values(me$through.cv, me$parm*tmp.mean, me$parm)
      new.gen.y$lt <- rnorm.censored(tmp.mean, tmp.sd, upper=data$cens$lt, negative.values.disallowed=logNormal.distrn | me$through.cv)
    }

    if (data$any.censored$i)
    {
      tmp.mean <- true.values$i
      tmp.sd <- cond.values(me$through.cv, me$parm*tmp.mean, me$parm)
      new.gen.y$i <- rnorm.censored(tmp.mean, tmp.sd, lower=data$cens$i$gt, upper=data$cens$i$lt)
    }
    
    if (logNormal.distrn) new.gen.y$log <- take.log(new.gen.y)    
  }
  else if (logNormal.distrn)
  {    
    if (data$any.censored$gt) 
    {
      new.gen.y$log$gt <- rnorm.censored(mu, sigma, lower=data$log.cens$gt)
      new.gen.y$gt <- exp(new.gen.y$log$gt)
    } 
        
    if (data$any.censored$lt) 
    {
      new.gen.y$log$lt <- rnorm.censored(mu, sigma, upper=data$log.cens$lt)
      new.gen.y$lt <- exp(new.gen.y$log$lt) 
    }
    
    if (data$any.censored$i)
    {  
      new.gen.y$log$i <- rnorm.censored(mu, sigma, lower=data$log.cens$i$gt, upper=data$log.cens$i$lt)
      new.gen.y$i <- exp(new.gen.y$log$i)
    }
  }
  else
  {
    if (data$any.censored$gt) new.gen.y$gt <- rnorm.censored(mu, sigma, lower=data$cens$gt)   
    if (data$any.censored$lt) new.gen.y$lt <- rnorm.censored(mu, sigma, upper=data$cens$lt) 
    if (data$any.censored$i)  new.gen.y$i  <- rnorm.censored(mu, sigma, lower=data$cens$i$gt, upper=data$cens$i$lt)
  }
  
  new.gen.y
} # end of y.gen


y.gen.inits <- function(data, mu, sigma, me.through.cv, logNormal.distrn)
{
  gen.y <- list()
  
  if (logNormal.distrn)
  {
    if (data$any.censored$gt) 
    {
      gen.y$log$gt <- rnorm.censored(mu, sigma, lower=data$log.cens$gt)
      gen.y$gt <- exp(gen.y$log$gt)
    }
    
    if (data$any.censored$lt) 
    {
      gen.y$log$lt <- rnorm.censored(mu, sigma, upper=data$log.cens$lt, negative.values.disallowed=me.through.cv & !logNormal.distrn)
      gen.y$lt <- exp(gen.y$log$lt)
    }
    
    if (data$any.censored$i)  
    {
      gen.y$log$i  <- rnorm.censored(mu, sigma, lower=data$log.cens$i$gt, upper=data$log.cens$i$lt)
      gen.y$i <- exp(gen.y$log$i)
    }
  }
  else
  {
    if (data$any.censored$gt) gen.y$gt <- rnorm.censored(mu, sigma, lower=data$cens$gt)
    if (data$any.censored$lt) gen.y$lt <- rnorm.censored(mu, sigma, upper=data$cens$lt, negative.values.disallowed=me.through.cv & !logNormal.distrn)
    if (data$any.censored$i)  gen.y$i  <- rnorm.censored(mu, sigma, lower=data$cens$i$gt, upper=data$cens$i$lt)
  }
  
  gen.y
} # end of y.gen.inits
 
 
############################################################################################################
#
# Fcts of local interest only
#


date.fhead <- function()
{
  tmp <- as.character(Sys.time())
  tmp <- sub(" ", "_", tmp)
  tmp <- gsub("[-\\:]", "", tmp) 
  tmp
} # end of date.fhead


temp <- function(fhead, dir='C:\\Users\\blitzrcr\\patrickb.stat\\tmp', ext='RData', time.stamp=F)
{
  dot <- grep("\\.", fhead)
  fname <- fhead
  
  if (time.stamp)
  {
    tmp <- date.fhead()
    fname <- paste(fname, tmp, sep="-")
  }
  
  if (length(dot) == 0) fname <- paste(fname, ext, sep='.')
  else fname <- fhead
  
  prefix <- '_'
  continue <- T
  
  while (continue)
  {
    path <- paste(dir, fname, sep='/')
    continue <- file.exists(path)
    fname <- paste(prefix, fname, sep='')
  }
  
  path
} # end of temp


