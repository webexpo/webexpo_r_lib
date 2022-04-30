# Version 0.4 (Apr 2020) 
#


# _______________________________________________________________
# Change Log
#
#
# Version 0.4 (Apr 2020)
# ----------------------
#  - The code for the calculation of f1 = dl.dmui and f2 = dl.dmu
#    was unchanged but moved to function f12
#  - A protection against non-convergence was also added to function one.subject.estimate
#
# Version 0.3 (Apr 2020)
# ----------------------
#   Added argument old.McNally
#   Changed order of arguments sw.range & sb.range in fct call
#   Renamed function f3 -> dl.dsigmaw
#
# Version 0.2 (Apr 2020)
# ----------------------
#   Changed name of function dl.dsigmaw -> f3



inits.BetweenWorker <- function(data, 
                         unif.sds=use.uniform.prior.on.sds,
                         sw.range=sigma.within.range, 
                         sb.range=sigma.between.range,
                         log.sw.mu=log.sigma.within.mu,  log.sw.prec=log.sigma.within.prec,
                         log.sb.mu=log.sigma.between.mu, log.sb.prec=log.sigma.between.prec,
                         max.niter=500, epsilon=1e-6, old.McNally=F
                         )
{
  # As init values for theta = (mu.worker, mu.overall, sigma.within, sigma.between),
  # we will use the values that maximize the posterior distrn


  # --- Useful functions -----------------------------------------------------------------


  dg1.dmu <- function(a, sigma, mu.sign=-1)
  {
    # modif_0.3 no default value for sigma

    z <- a/sigma
    vphi <- varphi(z)
    mu.sign * sum(vphi$fp) / sigma^2
  } # end of dg1.dmu


  dg1.dsigma <- function(a, sigma)
  {
    # modif_0.3 no default value for sigma

    z <- a/sigma
    vphi <- varphi(z)

    - (sum(z*vphi$fp) + sum(vphi$f)) / sigma^2
  } # end of dg1.dsigma


  dg2.dmu <- function(a, sigma, mu.sign=-1)
  {
    # modif_0.3 no default value for sigma

    z <- a/sigma
    vphi <- varphi(z)

    mu.sign * (sum(vphi$f) + sum(z*vphi$fp)) / sigma^2
  } # end of dg2.dmu


  dg2.dsigma <- function(a, sigma)
  {
    # modif_0.3 no default value for sigma

    z <- a/sigma
    vphi <- varphi(z)

    - (sum(vphi$fp*z^2) + 2*sum(z*vphi$f)) / sigma^2
  } # end of dg2.dsigma


  dgd.dmu <- function(l, r, mu, sigma)
  {
    # modif_0.3 no default value for sigma

    z2 <- (r-mu) / sigma
    z1 <- (l-mu) / sigma
    log.sigma <- log(sigma)

    log.u <- log.dphi(z1, z2)
    log.v <- log.dPhi(z1, z2)

    log.u1 <- dnorm(z1, log=T)
    log.u2 <- dnorm(z2, log=T) 

    log.up1 <- list(x=log.u1 + log(abs(z1)) - log.sigma, sign=sign(z1))
    log.up2 <- list(x=log.u2 + log(abs(z2)) - log.sigma, sign=sign(z2))
    
    log.vp1 <- dPhi.dmu(l-mu, sigma, log=T)
    log.vp2 <- dPhi.dmu(r-mu, sigma, log=T)

    log.up <- log.diff(log.up1, log.up2)
    log.vp <- log.diff(log.vp1, log.vp2)
      
    log.d1 <- list(x=log.up$x + log.v, sign=log.up$sign)
    log.d2 <- list(x=log.u$x + log.vp$x, sign=log.u$sign * log.vp$sign)
    log.d <- log.diff(log.d2, log.d1)
    log.d$x <- log.d$x - 2 * log.v
    d <- exp(log.d$x) * log.d$sign

    sum(d) / sigma
  } # end of dgd.dmu


  dgm.dmu <- function(l, r, mu, sigma)
  {
    # modif_0.3 no default value for sigma

    u2 <- gm1(r-mu, sigma, log=T)
    u1 <- gm1(l-mu, sigma, log=T)
    u <- log.diff(u1, u2)

    v <- log.dPhi((l-mu)/sigma, (r-mu)/sigma)

    up2 <- dgm1.dmu(r-mu, sigma, log=T)
    up1 <- dgm1.dmu(l-mu, sigma, log=T)
    u.prime <- log.diff(up1, up2)

    vp2 <- dPhi.dmu(r-mu, sigma, log=T)
    vp1 <- dPhi.dmu(l-mu, sigma, log=T)
    v.prime <- log.diff(vp1, vp2)
      
    d1 <- log.mult(u, v.prime)
    d2 <- list(x=u.prime$x + v, sign=u.prime$sign)
    d <- log.diff(d1, d2)
    d <- exp(d$x - 2*v) * d$sign

    sum(d)
  } # end of dgm.dmu


  dgm.dsigma <- function(l, r, mu, sigma)
  {
    # modif_0.3 no default value for sigma

    log.u2 <- gm1(r-mu, sigma, log=T)
    log.u1 <- gm1(l-mu, sigma, log=T)
    log.u <- log.diff(log.u1, log.u2)

    log.v <- log.dPhi((l-mu)/sigma, (r-mu)/sigma)

    log.up1 <- dgm1.dsigma(l-mu, sigma, log=T)
    log.up2 <- dgm1.dsigma(r-mu, sigma, log=T)
    log.up <- log.diff(log.up1, log.up2)

    log.vp1 <- dPhi.dsigma(l-mu, sigma, log=T)
    log.vp2 <- dPhi.dsigma(r-mu, sigma, log=T)
    log.vp <- log.diff(log.vp1, log.vp2)
   
    log.d1 <- list(x=log.up$x + log.v, sign=log.up$sign)
    log.d2 <- log.mult(log.u, log.vp)
    log.d <- log.diff(log.d2, log.d1)
    log.d$x <- log.d$x - 2*log.v
   
    d <- exp(log.d$x) * log.d$sign
    sum(d)
  } # end of dgm.dsigma


  dgm1.dmu <- function(a, sigma, log=F)
  {
    # a: vector
    # sigma: scalar

    z <- a / sigma

    if (log)
    {
      b <- z^2 - 1
      x <- dnorm(z, log=T) + log(abs(b)) - 2 * log(sigma)
      s <- sign(b)
      out <- list(x=x, sign=s)
    }
    else
    {
      phi <- dnorm(z)
      out <- phi * (z^2 - 1) / sigma^2
    }

    out
  } # end of dgm1.dmu


  dgm1.dsigma <- function(a, sigma, log=F)
  {
    # a: vector
    # sigma: scalar

    z <- a / sigma

    if (log)
    {
      log.phi <- dnorm(z, log=T)
      x <- log.phi + log(abs(z)) + log(abs(z^2-2)) - 2*log(sigma)
      out <- list(x=x, sign=sign(z) * sign(z^2 - 2))
    }
    else
    {
      phi <- dnorm(z)
      out <- phi * z *  (z^2 - 2) / sigma^2
    }

    out
  } # end of dgm1.dsigma


  dl.dsigmaw <- function(sigma.within, mu.worker, mu.overall,
                         y, l, r, d, nrows,
                         unif.sds, log.sw.moments)
  {
    # modif_0.3 changed name of function f3 -> dl.dsigmaw
    # modif_0.2 changed name of function dl.dsigmaw -> f3

    f <-  -nrows$y / sigma.within
    f.prime <- nrows$y / sigma.within^2


    if (nrows$y > 0)
    {
      tmp <- y$y - mu.worker[y$id] - mu.overall  # vector of length = nrows$y
      tmp <- sum(tmp^2)

      f <- f + tmp/sigma.within^3
      f.prime <- f.prime - 3*tmp/sigma.within^4
    }


    if (nrows$r > 0)
    {
      mu <- r$r - mu.worker[r$id] - mu.overall

      f <- f - g2(mu, sigma.within) # modif_0.3 added argument sigma.within
      f.prime <- f.prime - dg2.dsigma(mu, sigma.within)  # modif_0.3 added sigma.within as argument
    }


    if (nrows$l > 0)
    {
      mu <- mu.worker[l$id] + mu.overall - l$l 

      f <- f - g2(mu, sigma.within)  # modif_0.3 added argument sigma.within
      f.prime <- f.prime - dg2.dsigma(mu, sigma.within)  # modif_0.3 added sigma.within as argument
    }


    if (nrows$d > 0)
    {
      mu <- mu.worker[d$id] + mu.overall

      f <- f - gm(d$l, d$r, mu, sigma.within)  # modif_0.3 added argument sigma.within
      f.prime <- f.prime - dgm.dsigma(d$l, d$r, mu, sigma.within)  # modif_0.3 added sigma.within as argument
    }   


    if (!unif.sds)
    {
      logn <- logN(sigma.within, log.sw.moments)
      f <- f + logn$f
      f.prime <- f.prime + logn$fp
    }


    list(f=f, fp=f.prime)
  } # end of dl.dsigmaw


  dPhi.dmu <- function(a, sigma, log=F)
  {
    z <- a/sigma

    if (log)
    {
      x <- dnorm(z, log=T) - log(sigma)
      out <- list(x=x, sign=rep(-1, length(x)))
    }
    else
    {
      phi <- dnorm(z)
      out <- - phi / sigma
    }
  
    out
  } # end of dPhi.dmu


  dPhi.dsigma <- function(a, sigma, log=F)
  {
    z <- a/sigma

    if (log)
    {
      x <- dnorm(z, log=T) + log(abs(a)) - 2*log(sigma)
      out <- list(x=x, sign=-sign(a))
    }
    else
    {
      phi <- dnorm(z)
      out <- - phi * z / sigma
    }

    out
  } # end of dPhi.dsigma


  # new_0.4
  f12 <- function(mu.worker, mu.overall, sigma.within, sigma.between,
                  y.count, l.count, r.count, d.count, y.sum)
  {
    workers <- seq(along=mu.worker) # new_0.4

    tmp <- y.sum - (mu.worker + mu.overall) * y.count  # vector of length = n.workers
    f1 <- tmp / sigma.within^2
    df1.dmui <- -y.count / sigma.within^2
    df1.dsigmaw <- -2*tmp / sigma.within^3


    for (i in workers)
    {
      my.mu <- mu.worker[i] + mu.overall # modif_0.4: changed theta[i] for mu.worker[i]

      if (r.count[i] > 0)
      {
        tmp <- subset(r, id==i) 

        f1[i] <- f1[i] - g1(tmp$r - my.mu, sigma.within)  # modif_0.3 added sigma.within as argument
        df1.dmui[i] <- df1.dmui[i] - dg1.dmu(tmp$r - my.mu, sigma.within)          # modif_0.3 specified sigma.within as argument
        df1.dsigmaw[i] <- df1.dsigmaw[i] - dg1.dsigma(tmp$r - my.mu, sigma.within) # modif_0.3 specified sigma.within as argument
      }


      if (l.count[i] > 0)
      {
        tmp <- subset(l, id==i)

        f1[i] <- f1[i] + g1(my.mu - tmp$l, sigma.within)  # modif_0.3 added sigma.within as argument
        df1.dmui[i] <- df1.dmui[i] + dg1.dmu(my.mu - tmp$l, sigma.within, mu.sign=1) # modif_0.3 specified sigma.within as argument
        df1.dsigmaw[i] <- df1.dsigmaw[i] + dg1.dsigma(my.mu - tmp$l, sigma.within)   # modif_0.3 specified sigma.within as argument
      }


      if (d.count[i] > 0)
      {
        tmp <- subset(d, id==i)

        f1[i] <- f1[i] - gd(tmp$l, tmp$r, my.mu, sigma.within)  # added argument sigma.within
        df1.dmui[i] <- df1.dmui[i] - dgd.dmu(tmp$l, tmp$r, my.mu, sigma.within)  # modif_0.3 specified sigma.within as argument
        df1.dsigmaw[i] <- df1.dsigmaw[i] - dgm.dmu(tmp$l, tmp$r, my.mu, sigma.within) # dgd.dsigma = dgm.dmu  # modif_0.3 specified sigma.within as argument
      }
    }


    df1.dmu <- df1.dmui

    f2 <- sum(f1)
    df2.dmu <- sum(df1.dmui)
    df2.dsigmaw <- sum(df1.dsigmaw)

    f1 <- f1 - mu.worker/sigma.between^2
    df1.dmui <- df1.dmui - 1/sigma.between^2


    list(f1=f1, df1.dmui=df1.dmui, df1.dmu=df1.dmu, df1.dsigmaw=df1.dsigmaw,
         f2=f2, df2.dmu=df2.dmu, df2.dsigmaw=df2.dsigmaw)
  } # end of f12


  g1 <- function(a, sigma)
  {
    # modif_0.3 no default value for sigma

    z <- a/sigma
    sum(varphi(z)$f) / sigma
  } # end of g1


  g2 <- function(a, sigma)
  {
    # modif_0.3 no default value for sigma

    z <- a/sigma
    sum(z * varphi(z)$f) / sigma
  } # end of g2


  gd <- function(l, r, mu, sigma)
  {
    # modif_0.3 no default value for sigma

    # l, r:       vectors of same length
    # mu, sigma:  scalars 

    z2 <- (r-mu) / sigma
    z1 <- (l-mu) / sigma

    num2.log <- dnorm(z2, log=T)
    num1.log <- dnorm(z1, log=T)
    num.log <- log.diff.exp(num1.log, num2.log)

    denom.log <- log.dPhi((l-mu)/sigma, (r-mu)/sigma)

    ratios <- exp(num.log$x - denom.log) * num.log$sign

    sum(ratios) / sigma
  } # end of gd


  gm <- function(l, r, mu, sigma)
  {
    # modif_0.3 no default value for sigma

    # l, r:       vectors of same length
    # mu, sigma:  scalars 

    u2 <- gm1(r-mu, sigma, log=T)
    u1 <- gm1(l-mu, sigma, log=T)
    u <- log.diff(u1, u2)

    log.v <- log.dPhi((l-mu)/sigma, (r-mu)/sigma)

    log.ratios <- u$x - log.v
    ratios <- u$sign * exp(log.ratios)

    sum(ratios)
  } # end of gm


  gm1 <- function(a, sigma, log=F)
  {
    # a: vector
    # sigma: scalar

    z <- a / sigma

    if (log)
    {
      x <- dnorm(z, log=T) + log(abs(a)) - 2 * log(sigma)
      out <- list(x=x, sign=sign(a))
    }
    else
    {
      out <- dnorm(z) * z / sigma
    }

    # if log = TRUE, return a list with dimensions (x, sign)
    # if log = FALSE, return a scalar

    out
  } # end of gm1


  log.diff <- function(l1, l2)
  {
    l1$sign <- - l1$sign
    log.sum(l1, l2)
  } # end of log.diff


  log.diff.exp <- function(l1, l2)
  {
    # l1, l2: vectors of same length
    # return log(abs(exp(l2) - exp(l1)) along with its sign

    s <- sign(l2 - l1)
    x <- rep(NA, length(s))

    w <- which(s > 0)
    if (length(w) > 0) x[w] <- l2[w] + log(1 - exp(l1[w]-l2[w]))
    
    w <- which(s < 0)
    if (length(w) > 0) x[w] <- l1[w] + log(1 - exp(l2[w]-l1[w]))
 
    w <- which(s == 0)
    if (length(w) > 0) x[w] <- -Inf

    list(x=x, sign=s)    
  } # end of log.diff.exp


  log.dphi <- function(z1, z2)
  {
    log1 <- dnorm(z1, log=T)
    log2 <- dnorm(z2, log=T)

    log.diff.exp(log1, log2)
  } # end of log.dphi


  log.dPhi <- function(z1, z2)
  {
    denom2.log <- pnorm(z2, log.p=T)
    denom1.log <- pnorm(z1, log.p=T)
    denom.log <- log.diff.exp(denom1.log, denom2.log)$x

    w <- which(is.infinite(denom.log))

    if (length(w) > 0)
    {
      tmp.denom2.log <- pnorm(z2[w], log.p=T, lower.tail=F)
      tmp.denom1.log <- pnorm(z1[w], log.p=T, lower.tail=F)
      denom.log[w] <- log.diff.exp(tmp.denom1.log, tmp.denom2.log)$x
    }

    denom.log
  } # end of log.dPhi


  log.mult <- function(l1, l2)
  {
    list(x=l1$x + l2$x, sign=l1$sign * l2$sign)
  } # end of log.mult


  log.sum <- function(l1, l2)
  {
    # l1 & l2 are two lists with dimension ($x, $sign)

    out.x    <- rep(NA, length(l1$x))
    out.sign <- rep(NA, length(l1$x))


    w <- which(l1$sign == l2$sign)

    if (length(w) > 0)
    {
      e <- data.frame(l1=l1$x[w], l2=l2$x[w])
      e.max <- pmax(e$l1, e$l2)
      e <- e - e.max
      x <- e.max + log(exp(e$l1) + exp(e$l2))
      out.x[w] <- x
      out.sign[w] <- l1$sign[w]

      w2 <- which(is.infinite(e.max))
      if (length(w2) > 0) out.x[w[w2]] <- -Inf # prevent NaN's
    }


    w <- which(l1$sign != l2$sign)

    if (length(w) > 0)
    {
      tmp <- log.diff.exp(l1$x[w], l2$x[w])
      out.x[w] <- tmp$x
      out.sign[w] <- tmp$sign * l2$sign[w]
    }

    out <- list(x=out.x, sign=out.sign)
    out
  } # end of log.sum


  logN <- function(sigma, m)
  {
    f <- -1/sigma * (1 + m$prec * (log(sigma) - m$mu))
    fp <- (1 - m$prec*(m$mu + 1 - log(sigma))) / sigma^2
  
    list(f=f, fp=fp)
  } # end of logN


  one.subject.estimate <- function(my.y, my.l, my.r, my.d, epsilon=1e-6, max.niter=250)
  {
    tmp.mu    <- numeric(0)
    tmp.sigma <- numeric(0)

    count.y <- nrow(my.y)
    count.l <- nrow(my.l)
    count.r <- nrow(my.r)
    count.d <- nrow(my.d)

    # First find initial values for Newton-Raphson algorithm

    if (count.y > 0)
    {
      mu <- mean(my.y$y)
      sigma <- sd(my.y$y)
      if (is.na(sigma)) sigma <- numeric(0)
      tmp.mu <- c(tmp.mu, mu)
      tmp.sigma <- c(tmp.sigma, sigma)
    }

    if (count.d > 0)
    {
      u <- median(my.d$r)
      l <- median(my.d$l)
      mu <- (u+l)/2
      sigma <- (u-l)/4
      tmp.mu <- c(tmp.mu, mu)
      tmp.sigma <- c(tmp.sigma, sigma)
    }

    if (length(tmp.mu) > 0)
    {
      mu <- median(tmp.mu)

      if (count.l > 0)
      {
        l <- median(my.l$l)
        sigma <- abs(mu-l)/2
        tmp.sigma <- c(tmp.sigma, sigma)
      }

      if (count.r > 0)
      {
        r <- median(my.r$r)
        sigma <- abs(mu-r)/2
        tmp.sigma <- c(tmp.sigma, sigma)
      }
    }
    else if (count.l > 0 & count.r > 0)
    {
      l <- median(my.l$l)
      r <- median(my.r$r)
      sigma <- abs(r-l)/4
      mu <- (l+r)/2
      tmp.mu <- c(tmp.mu, mu)
      tmp.sigma <- c(tmp.sigma, sigma)
    }

    mu <- median(tmp.mu)
    sigma <- 0.001


    # Prepare objects

    y <- my.y$y
    l <- my.l$l
    r <- my.r$r
    dl <- my.d$l
    dr <- my.d$r

    y.len <- length(y)
      my.ysum <- sum(y)
    r.len <- length(r)
    l.len <- length(l)
    d.len <- length(dl)

    continue <- T
    converged <- F
    counter <- 0

    J <- matrix(NA, 2, 2)


      # Make sure sigma is not too small to start with

      if (d.len > 0)
      {
        tmp <- gd(dl, dr, mu, sigma)

        while (is.nan(tmp) | is.infinite(tmp))
        {
          # tmp is not-a-number if the denominator in gd was null
          sigma <- sigma * 2
          tmp <- gd(dl, dr, mu, sigma)
        }
      }


    theta <- c(mu, sigma)


    # Run Newton-Raphson algorithm


    while (continue)
    {
      counter <- counter + 1

      a <- y - mu                  # vector of length = y.len
      a.sum <- my.ysum - y.len*mu  # scalar

      f1 <- a.sum / sigma^2                  # scalar
      f2 <- -y.len/sigma + sum(a^2)/sigma^3  # scalar


      # Define Jacobian (first part)

      df1.dmu <- -y.len / sigma^2
      df2.dmu <- -2*a.sum/sigma^3
      df2.dsigma <- y.len/sigma^2
        if (y.len > 0) df2.dsigma <- df2.dsigma - 3*sum(a^2)/sigma^4


      if (r.len > 0)
      {
        f1 <- f1 - g1(r-mu, sigma)
        f2 <- f2 - g2(r-mu, sigma)

        # Complete defn of Jacobian

        df1.dmu    <- df1.dmu    - dg1.dmu(r-mu, sigma)
        df2.dmu    <- df2.dmu    - dg2.dmu(r-mu, sigma)
        df2.dsigma <- df2.dsigma - dg2.dsigma(r-mu, sigma)
      }

      if (l.len > 0)
      {
        f1 <- f1 + g1(mu-l, sigma)
        f2 <- f2 - g2(mu-l, sigma)

        # Complete defn of Jacobian

        df1.dmu    <- df1.dmu    + dg1.dmu(mu-l, sigma, mu.sign=1)
        df2.dmu    <- df2.dmu    - dg2.dmu(mu-l, sigma, mu.sign=1)
        df2.dsigma <- df2.dsigma - dg2.dsigma(mu-l, sigma)   
      }

      if (d.len > 0)
      {
        f1 <- f1 - gd(dl, dr, mu, sigma)
        f2 <- f2 - gm(dl, dr, mu, sigma) 

        # Complete defn of Jacobian

        df1.dmu    <- df1.dmu    - dgd.dmu(dl, dr, mu, sigma) 
        df2.dmu    <- df2.dmu    - dgm.dmu(dl, dr, mu, sigma) 
        df2.dsigma <- df2.dsigma - dgm.dsigma(dl, dr, mu, sigma)
      }


      f <- matrix(c(f1, f2), nrow=2) # matrix of dimensions 2 x 1

      J[1,1] <- df1.dmu
      J[1,2] <- df2.dmu
      J[2,1] <- df2.dmu
      J[2,2] <- df2.dsigma


      previous.sigma <- sigma

      
      J.logdet <- log(abs(det(J))) # new_0.4

      # modif_0.4

      if (any(is.na(J)) | J.logdet < -80)
      {
        continue <- F
        converged <- F
      }
      else
      {
        change <- as.vector(solve(J) %*% f) # vector of length 2
        theta <- theta - change
        continue <- any(abs(change) > epsilon)
        converged <- !continue
      }
 

      if (theta[2] < 0)
      {
        theta[2] <- previous.sigma / 2
        continue <- T
      }

      if (continue) continue <- counter < max.niter

      mu <- theta[1]
      sigma <- theta[2]
    }


    list(converged=converged, mu=mu, sigma=sigma)
  } # end of one.subject.estimate


  varphi <- function(z)
  {
    log.phi <- dnorm(z, log=T)
    log.Phi <- pnorm(z, log.p=T)

    vphi <- exp(log.phi - log.Phi)
    vphi.prime <- - (z*vphi + vphi^2)

    list(f=vphi, fp=vphi.prime)
  } # end of varphi


  # -----------------------------------------------------------------------------------------
 
  # Start function

  log.sw.moments <- list(mu=log.sw.mu, prec=log.sw.prec)
  log.sb.moments <- list(mu=log.sb.mu, prec=log.sb.prec)

  if (unif.sds)
  {
    sw.lim <- list(l=min(sw.range), u=max(sw.range))
    sb.lim <- list(l=min(sb.range), u=max(sb.range))
  }
  else
  {
    sw.lim <- list(l=0, u=Inf)
    sb.lim <- list(l=0, u=Inf)
  }


  y <- data.frame(y=data$y, id=data$id$y)  # Known data points

  # Censored data

  if (old.McNally)
  {
    # modif_0.3 block made conditional

    r <- data.frame(r=data$cens$r, id=data$id$r)                     # of shape y < r
    l <- data.frame(l=data$cens$l, id=data$id$l)                     # of shape y > l
    d <- data.frame(l=data$cens$i$l, r=data$cens$i$r, id=data$id$i)  # of shape l < y < r
  }
  else
  {
    # new_0.3
   
    r <- data.frame(r=data$cens$lt, id=data$i$lt)                      # of shape y < lt
    l <- data.frame(l=data$cens$gt, id=data$i$gt)                      # of shape y > gt
    d <- data.frame(l=data$cens$i$gt, r=data$cens$i$lt, id=data$id$i)  # of shape gt < y < lt
  }


  nrows <- list(y=nrow(y), r=nrow(r), l=nrow(l), d=nrow(d))

  # Prepare a few objects

  n.workers <- data$n.workers
  
  tmp <- rep(NA, n.workers)
  f1 <- tmp
  df1.dmui <- tmp
  df1.dmu <- tmp
  df1.dsigmaw <- tmp

  # Count number of values in objects y, r, l & i for each worker

  tmp <- rep(0, n.workers)
  y.count <- tmp
  r.count <- tmp
  l.count <- tmp
  d.count <- tmp

  y.sum <- tmp


  tmp <- aggregate(y~id, data=y, length)
  y.count[tmp$id] <- tmp$y

  tmp <- aggregate(r~id, data=r, length)
  r.count[tmp$id] <- tmp$r
 
  tmp <- aggregate(l~id, data=l, length)
  l.count[tmp$id] <- tmp$l

  tmp <- aggregate(r~id, data=d, length)
  d.count[tmp$id] <- tmp$r


  tmp <- aggregate(y~id, data=y, sum)
  y.sum[tmp$id] <- tmp$y
 

  #######################################################################
  # Define initial values for the quest for optimal theta
  # (by first getting individual estimates first)


  estimable <- y.count > 0 | d.count > 0 | (r.count > 0 & l.count > 0) # logical vector of length = n.workers
  mu.estimates <- rep(NA, length(estimable))
  sw.estimates <- rep(NA, length(estimable))

  for (i in seq(along=estimable))
  {
    if (estimable[i])
    {
      my.y <- subset(y, id==i)
      my.l <- subset(l, id==i)
      my.r <- subset(r, id==i)
      my.d <- subset(d, id==i)

      run.one.subject.estimate <- T
      if (length(my.y$y) == 1) run.one.subject.estimate <- F
      if (length(my.y$y) > 1) run.one.subject.estimate <- sd(my.y$y) > 0


      if (!run.one.subject.estimate) 
      {
        my.ybar <- mean(my.y$y)
        run.one.subject.estimate <- any(my.l$l > my.ybar) | any(my.d$l > my.ybar) | any(my.r$r < my.ybar) | any(my.d$r < my.ybar)
      }
      
      if (!run.one.subject.estimate)
      {
        tmp.l <- c(my.d$l, my.r$l)
        tmp.r <- c(my.d$r, my.r$r)
        tmp.l <- median(tmp.l)
        tmp.r <- median(tmp.r)
        tmp <- list(converged=T, mu=my.ybar, sigma=abs(tmp.r-tmp.l)/6) 
      }
      else
      {
        tmp <- one.subject.estimate(my.y, my.l, my.r, my.d)
      }

      if (tmp$converged)
      {
        mu.estimates[i] <- tmp$mu
        sw.estimates[i] <- tmp$sigma
      }
    }
  }


  # Initial values for the multivariate Newron-Raphson algorithm

  if (any(!is.na(mu.estimates)))
  {
    sigma.between <- sd(mu.estimates, na.rm=T)
    mu.overall <- mean(mu.estimates, na.rm=T)
    mu.worker <- mu.estimates - mu.overall
    mu.worker[is.na(mu.worker)] <- 0
    sigma.within <- median(sw.estimates, na.rm=T)
  }


# Ce bloc est inutile -- je le garde seulement au cas ou on aurait besoin d'y revenir...
#    # Prepare data used to derive initial values
#
#    df <- data.frame(y=numeric(0), id=numeric(0))
#
#    if (length(data$id$y) > 0)
#    {
#      df.tmp <- data.frame(y=data$y, id=data$id$y)
#      df <- rbind.data.frame(df, df.tmp)
#    }
#
#    if (length(data$id$l) > 0)
#    {
#      df.tmp <- data.frame(y=data$cens$l, id=data$id$l)
#      df <- rbind.data.frame(df, df.tmp)
#    }
#
#    if (length(data$id$r) > 0)
#    {
#      df.tmp <- data.frame(y=data$cens$r, id=data$id$r)
#      df <- rbind.data.frame(df, df.tmp)
#    }
#
#    if (length(data$id$i) > 0)
#    {
#      df.tmp <- data.frame(y=(data$cens$i$l + data$cens$i$r)/2, id=data$id$i)
#      df <- rbind.data.frame(df, df.tmp)
#    }
#    
#
#    # For mu.worker, take the observed mean by worker
#    mu.worker <- aggregate(y~id, data=df, mean)$y
#        
#    mu.overall <- mean(mu.worker)
#    mu.worker <- mu.worker - mu.overall # center mu.worker
#
#    # For sigma.within, use sd of initial values around worker predicted means
#    predicted <- mu.overall + mu.worker[df$id]
#    sigma.within <- sd(df$y-predicted)
#
#    sigma.between <- sd(mu.worker)


  # .................................................................
  # Run a univariate Newton-Raphson algorith to find an initial value
  # for sigma.within, with fixed values for mu.overall & mu.worker
  # as calculated above

  continue <- T
  converged <- F
  counter <- 0

  while (continue)
  {
    # modif_0.3 changed name f3 -> dl.dsigmaw
    # modif_0.2 changed name dl.dsigmaw -> f3

    tmp <- dl.dsigmaw(sigma.within, mu.worker, mu.overall,
                      y, l, r, d, nrows, unif.sds, log.sw.moments)

    change <- tmp$f / tmp$fp
    sigma.within <- sigma.within - change
      
    continue <- abs(change) > epsilon
    converged <- !continue
    counter <- counter + 1

    if (continue) continue <- counter < max.niter
  }

  # .................................................................


  # Make sure that sigma.within is on
  # the right side of its optimal value (that is, shows a negative slope)

  continue <- T
    
  if (sigma.within  > sw.lim$u) sigma.within  <- sw.lim$u
  if (sigma.between > sb.lim$u) sigma.between <- sb.lim$u

  while (continue)
  {
    # modif_0.3 changed name f3 -> dl.dsigmaw
    # modif_0.2 changed name dl.dsigmaw -> f3
    tmp <- dl.dsigmaw(sigma.within, mu.worker, mu.overall,
                      y, l, r, d, nrows, unif.sds, log.sw.moments)

    continue <- tmp$fp > 0
    if (continue) sigma.within <- (sw.lim$l + sigma.within) / 2
  }

  sigma.within <- (sigma.within + sw.lim$l) / 2


  # *** Newton-Raphson algorithm to find optimal theta ********************


  theta0 <- c(mu.worker, mu.overall, sigma.within)
  theta <- theta0
  workers <- seq(n.workers)

  continue <- T
  converged <- F
  counter <- 0





  while (continue)
  {
    # ----------------------------------------------------------------------- 
    # f1 = dl/dmu_i, i = 1, 2, ..., n.workers -------------------------------
    # f2 = dl/dmu


    tmp <- f12(mu.worker, mu.overall, sigma.within, sigma.between, y.count, l.count, r.count, d.count, y.sum)

      f1 <- tmp$f1
      f2 <- tmp$f2

      df1.dmui    <- tmp$df1.dmui
      df1.dmu     <- tmp$df1.dmu
      df1.dsigmaw <- tmp$df1.dsigmaw

      f2          <- tmp$f2
      df2.dmu     <- tmp$df2.dmu
      df2.dsigmaw <- tmp$df2.dsigmaw


    # -------------------------------------------------------------------------------------------------
    # Start building Jacobian J

    J.top.left <- diag(df1.dmui) # diagonal matrix of dimension n.workers x n.workers

    # df1.dsigmab <- 2*mu.worker/sigma.between^3 # vector of length = n.workers

    #J.bottom.left <- c(df1.dmu, df1.dsigmaw, df1.dsigmab)   # vector of length = 3 * n.workers

    J.bottom.left <- c(df1.dmu, df1.dsigmaw)                        # vector of length = 2 * n.workers
    J.bottom.left <- matrix(J.bottom.left, byrow=T, ncol=n.workers) # matrix of dim = 2 x n.workers
    J.top.right <- t(J.bottom.left)

    J.left <- rbind.data.frame(J.top.left, J.bottom.left)


    # --- f3 = dl/dsigmaw -----------------------------------------------------

    # modif_0.3 changed name f3 -> dl.dsigmaw
    # modif_0.2 changed name dl.dsigmaw -> f3
    tmp <- dl.dsigmaw(sigma.within, mu.worker, mu.overall,
                      y, l, r, d, nrows, unif.sds, log.sw.moments)

    f3 <- tmp$f
    df3.dsigmaw <- tmp$fp
 

    # --- Complete Jacobian ----------------------------------------------------------

    J.bottom.right <- c(df2.dmu, df2.dsigmaw, df2.dsigmaw, df3.dsigmaw) # vector of length 4

    J.bottom.right <- matrix(J.bottom.right, nrow=2) # matrix of dimensions 2 x 2

    J.right <- rbind.data.frame(J.top.right, J.bottom.right)
    J <- cbind.data.frame(J.left, J.right)


    f <- matrix(c(f1, f2, f3), ncol=1)

    change <- as.vector(solve(J) %*% f) # vector of length = n.workers + 2
    theta <- theta - change

    continue <- any(abs(change) > epsilon)
    converged <- !continue
    counter <- counter + 1


    previous.sw <- sigma.within

    mu.worker     <- theta[workers]
    mu.overall    <- theta[n.workers + 1]
    sigma.within  <- theta[n.workers + 2]


      # safety against sigma's out of range

      if (sigma.within < sw.lim$l)
      {
        sigma.within <- (sw.lim$l + previous.sw)/2
        theta[n.workers + 2] <- sigma.within
        continue <- T
      }
      else if (sigma.within > sw.lim$u)
      {
        sigma.within <- (previous.sw + sw.lim$u)/2
        theta[n.workers + 2] <- sigma.within
        continue <- T
      } 


    if (continue) continue <- counter < max.niter
  }


  sigma.between <- sd(mu.worker)

  list(converged=converged, mu.worker=mu.worker, mu.overall=mu.overall, sigma.within=sigma.within, sigma.between=sigma.between)
} # end of inits.BetweenWorker
