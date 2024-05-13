# this.folder=this.dir()    on pourrait s'en servir (si on intalle le pkg 'this.path')
# puis faire setwd(this.folder)
source('c:/users/patri/home/bin/R/fcts/algebra.R')
source('C:/Users/patri/home/bin/R/fcts/bisectionalSearch.R')
source("c:/users/patri/home/bin/R/fcts/dens.gen.icdf/dens-gen-icdf.R")
source('c:/users/patri/home/consultation/l/Lavoue/webexpo/R/data-summary.R')
source("c:/users/patri/home/consultation/L/Lavoue/webexpo/R/fcts.R")
source("c:/users/patri/home/consultation/L/Lavoue/webexpo/riskband/R/CPL_prior.R")


# Version 1.3 (May 2024)


# -----------------------------------


# ------------------------------------------------------------------------------
# New in
# Version 1.3 (May 2024)
#
#
#
#                                                            (end of Change Log)
               
                                       
SEG.riskband <- function(y=numeric(0), lt=numeric(0), gt=numeric(0), 
  interval.gt=numeric(0), interval.lt=numeric(0),
  n.iter=15000, n.burnin=500, n.thin=1, n.chains=1, monitor.burnin=FALSE,
  A=c(0.01, 0.1, 0.5, 1), 
  riskband.prior.prob=rep(NA, R), 
  gm.min=NA, gm.max=NA, 

#    mu.lower=ifelse(is.na(gm.min), -3.0, log(gm.min)),
#    mu.upper=ifelse(is.na(gm.max),  6.2, log(gm.max)),
  mu.lower=ifelse(is.na(gm.min), -7.600902, log(gm.min)),
  mu.upper=ifelse(is.na(gm.max),  1.609438, log(gm.max)),

  gsd.min=1.05, gsd.max=4,  
    sigma.lower=log(gsd.min), sigma.upper=log(gsd.max),  
  init.mu=rep(default.inits$mu, n.chains), 
  init.sigma=rep(default.inits$sigma, n.chains), 
  outcome.is.logNormally.distributed=TRUE, 
  equally.probable.riskbands=TRUE,  # ignored when riskband.prior.prob is provided
  use.continuous.piecewise.linear.prior=length(cpl.prior)>0, 
  cpl.prior=list(), 
  quantile=0.95, 
  me.sd.range=numeric(0), cv.range=numeric(0)) 
{ 
  # If you want to use the Continuous-Piecewise Linear prior distribution,
  # you can either: 
  #    a) specify  use.continuous.piecewise.linear.prior = TRUE and let the code calculate the CPL prior distribution 
  # or b) provide the CPL prior distribution through the argument cpl.prior=my.CPL.prior 
  #       (assuming you have calculated the CPL prior beforehand and saved the result under the name my.CPL.prior)
  #       IMPORTANT/WARNING: if you provide a value for cpl.prior (the list returned by CPL.prior), 
  #                          the values (riskband.prior.prob,
  #                                      gm.min, gm.max, mu.lower, mu.upper, 
  #                                      gsd.min, gsd.max, sigma.lower, sigma.upper
  #                                      quantile, A)
  #                          given to SEG.riskband's arguments will be ignored.

  # IMPORTANT:
  # Selon la valeur de outcome.is.logNormally.distributed (TRUE/FALSE)
  # nous devrions demander à l'usager 
  #  i) seulement (gm.min, gm.max, gsd.min, gsd.max)               si outcome.is.logNormally.distributed = TRUE
  # ii) seulement (mu.lower, mu.upper, sigma.lower, sigma.upper)   sinon


  # ___ Useful function(s) ____________________________________________________________________________
  

  CPL.rnd.mu <- function(sigma, prior, mu.moments, NR.epsilon=1e-6)
  {
    exp.minusz2_2 <- function(x, mu, sd)
    {
      z <- (x - mu) / sd
      return(exp(-z^2/2))
    } # end of exp.minusz2



    ybar    <- mu.moments$mean
    cond.sd <- mu.moments$sd

    mu0 <- prior$geometry$mu.lim[1]
    mu1 <- prior$geometry$mu.lim[2]

    z <- prior$geometry$z


    # Sampling from f(mu|data,sigma)

    # mu cut-offs
    mu.lim <- prior$x95 - z*sigma
    mu.lim <- pmax(mu.lim, mu0)
    mu.lim <- pmin(mu.lim, mu1)

    x95.0 <- prior$x95[-length(prior$x95)]

    mulim.left  <- mu.lim[-length(mu.lim)]
    mulim.right <- mu.lim[-1]

    s <- prior$s * prior$geometry$cos.theta
    h0 <- prior$h[-length(prior$h)]
    a <- h0 - s*x95.0 + s*z*sigma

    Phi <- pnorm(mu.lim, ybar, cond.sd)
    F1 <- (a + s*ybar) * diff(Phi)
    E <- exp(-(mu.lim-ybar)^2/2/cond.sd^2)
    E <- -diff(E)
    F2 <- s*cond.sd/sqrt(2*pi)*E
    F <- F1 + F2
      if (any(F < -1e-10))  stop('error in F//1')

      F <- pmax(F, 0) # correction for numeric imprecision

    # -----------------------------------------------------------------------------------
    # Find in which segment the median is found
    Fcum <- cumsum(F)
    Fsum <- sum(F)
    j.median <- min(which(Fcum/Fsum >= 0.5))

    # -----------------------------------------------------------------------------------

    j <- sample(seq(along=F), 1, prob=F)  # sampled segment
    riskband <- prior$r[j]

    u <- runif(1)*F[j]
    my.a <- a[j]
    my.s <- s[j]
    mu0 <- mulim.left[j]
    range <- c(mu0, mulim.right[j])

    A <- my.a + my.s*ybar
    B <- my.s*cond.sd/sqrt(2*pi)
    F2.mu0 <- exp(-(mu0-ybar)^2/2/cond.sd^2)
    pnorm.mu0 <- pnorm(mu0, ybar, cond.sd)
    b <- B/cond.sd


      # Starting point

      if (j == j.median)
      {
        extrapolation <- function(start.side, range, ybar, cond.sd, pnorm.mu0, F2.mu0, A, B, u)
        {
          mu <- range[start.side]
          F2.mu <- exp.minusz2_2(mu, ybar, cond.sd)
          z <- (mu - ybar) / cond.sd

          fcum <- A * (pnorm(mu, ybar, cond.sd) - pnorm.mu0) + B*(F2.mu0 - F2.mu)
            F2p <- - F2.mu * z / cond.sd
            F2s <- - 1/cond.sd * (F2p*z + F2.mu/cond.sd)
          f <- A * dnorm(mu, ybar, cond.sd) - B*F2p # f = F'  (where F = fcum)
          fp <- -dnorm((mu-ybar)/cond.sd) * z / cond.sd^2 - B * F2s

          # Extrapolation via the quadratic approximation around mu (to find a starting point closer to the solution [with a better slope])

          qA <- fp/2
          qB <- f - 2*qA*mu
          qC <- fcum - qA*mu^2-qB*mu

          if (qA == 0 && ybar >= range[1] && ybar <= range[2])  return(ybar)

          delta <- qB^2 - 4*qA*(qC-u)

          if (delta < 0) return(NA)

          direction <- ifelse(start.side==1, 1, -1)

          mu <- -qB/(2*qA) + direction*sqrt(delta)/2/abs(qA)

          if (mu > range[1] && mu < range[2])  return(mu)
          return(NA)
        } # end of extrapolation


        mu                <- extrapolation(1, range, ybar, cond.sd, pnorm.mu0, F2.mu0, A, B, u)
        if (is.na(mu)) mu <- extrapolation(2, range, ybar, cond.sd, pnorm.mu0, F2.mu0, A, B, u)
        if (is.na(mu)) mu <- mean(range)
      }
      else mu <- ifelse(j > j.median, range[1], range[2])


      continue <- TRUE

      while (continue)
      {
        F2.mu <- exp.minusz2_2(mu, ybar, cond.sd)
        z <- (mu - ybar) / cond.sd
        fcum <- A * (pnorm(mu, ybar, cond.sd) - pnorm.mu0) + B * (F2.mu0 - F2.mu)
        f <- A * dnorm(mu, ybar, cond.sd) + b * F2.mu * z
        change <- -(fcum - u) / abs(f)
        i <- ifelse(change > 0, 1, 2)
        range[i] <- mu
        mu <- mu + change
        if (mu > range[2] || mu < range[1])  mu <- mean(range)
        converged <- abs(change) < NR.epsilon
        continue <- !converged

        if (continue)
        {
          squeezed <- diff(range) < NR.epsilon
          continue <- !squeezed
        }
      }


    return(list(mu=mu, riskband=riskband))
  } # end of CPL.rnd.mu


  CPL.rnd.sigma <- function(mu, beta, prior, data, NR.epsilon=1e-6)
  {
    # alpha, gamma.alpha: vectors of length 2
    # prior: the CPL prior

    dSqrtInvGamma <- function(x, alpha, beta)
    {
      J <- 2/x^3
      return(J * dgamma(1/x^2, alpha, beta))
    } # end of dSqrtInvGamma


    pSqrtInvGamma <- function(x, alpha, beta)
    {
      pgamma(1/x^2, alpha, beta, lower.tail=FALSE)
    } # end of pSqrtInvGamma


    qSqrtInvGamma <- function(p, alpha, beta)
    {
      tau <- qgamma(p, alpha, beta, lower.tail=FALSE)
      return(1/sqrt(tau))
    } # end of qSqrtInvGamma


    SqrtInvertedGamma.constant <- function(beta, data, j)
    {
      2 * beta^data$alpha[j] / data$gamma.alpha[j]
    } # end of SqrtInvertedGamma.constant


    # Compute sigma-intervals and their respective weights

    alpha <- data$alpha
    N     <- data$N

    z <- prior$geometry$z

    x95 <- prior$x95
    x95 <- x95[-c(1, length(x95))]

    sigma <- (x95 - mu) / z
    w <- which(sigma > prior$geometry$sigma.lim[1] & sigma < prior$geometry$sigma.lim[2])


    h0 <- prior$h[-length(prior$h)]
    s  <- prior$s * prior$geometry$cos.theta


    if (length(w) == 0)
    {
      m <- sum(sigma < prior$geometry$sigma.lim[1])
      w <- m + 1
      sigma <- prior$geometry$sigma.lim
    }
    else
    {
      sigma <- sigma[w]
      w <- c(w, max(w)+1)
      sigma <- c(prior$geometry$sigma.lim[1], sigma, prior$geometry$sigma.lim[2])
    }


    x95 <- prior$x95[-length(prior$x95)]
    a <- h0 - s*x95 + s*mu
    b <- z*s

    a <- a[w]
    b <- b[w]

    if (data$N >= 2)  c1 <- SqrtInvertedGamma.constant(beta, data, 1)


    # Compute intervals probabilities and sample from them
    # [N-dependent algorithm]


    if (data$N > 2)
    {
      # Compute integrals

      c2 <- SqrtInvertedGamma.constant(beta, data, 2)

      pr1 <- pSqrtInvGamma(sigma, (N-1)/2, beta)
      I1 <- diff(pr1)

      pr2 <- pSqrtInvGamma(sigma, (N-2)/2, beta)
      I2 <- diff(pr2)


      F <- a/c1*I1 + b/c2*I2
      if (any(F < -1e-10))  stop('error in F//2')
        F <- pmax(F, 0) # correction for numeric imprecision

      j <- sample(seq(along=F), 1, prob=F) # Sample sigma segment

      a <- a[j]
      b <- b[j]
      u <- runif(1)

        k1 <- a / c1
        k2 <- b / c2


      if (a > 0 & b > 0)
      {
        # Sample from a mixture of two Sqrt-Inverted-Gamma distrns

        wt <- c(k1*I1[j], k2*I2[j])
        pr <- wt[2] / sum(wt)
        i <- rbinom(1, 1, pr)

        if (i == 0)  pr <- pr1[c(j, j+1)]
        else         pr <- pr2[c(j, j+1)]

        u <- pr[1] + u*diff(pr) # value between pr[1] & pr[2]
        sigma <- qSqrtInvGamma(u, alpha[i+1], beta)
        return(sigma)
      }


      # a < 0 || b < 0

      sigma0 <- sigma[j]
      sigma  <- sigma[j+1] # upper limit value, used as a starting point

      p01 <- pSqrtInvGamma(sigma0, alpha[1], beta)
      p02 <- pSqrtInvGamma(sigma0, alpha[2], beta)

      range <- c(sigma0, sigma)

      fcum <- k1 * (pSqrtInvGamma(sigma,alpha[1],beta) - p01) + k2 * (pSqrtInvGamma(sigma,alpha[2],beta) - p02)
      u <- u*fcum


      continue <- TRUE

      while (continue)
      {
        fcum <- k1 * (pSqrtInvGamma(sigma,alpha[1],beta) - p01) + k2 * (pSqrtInvGamma(sigma,alpha[2],beta) - p02)
        f <- k1 * dSqrtInvGamma(sigma,alpha[1],beta) + k2 * dSqrtInvGamma(sigma,alpha[2],beta)
        change <- - (fcum-u) / f

        i <- ifelse(change>0, 1, 2)
        range[i] <- sigma

        sigma <- sigma + change
        if (sigma > range[1] || sigma < range[2])  sigma <- mean(range)

        converged <- abs(change) < NR.epsilon
        continue <- !converged

        if (continue)
        {
          squeezed <- diff(range) < NR.epsilon
          continue <- !squeezed
        }
      }

      return(sigma)
    }


    # --- N <= 2 ---------------------------------------------------------------


    e1 <- function(beta, sigma){exp(-beta/sigma^2)/sigma}

    I <- expInt_E1(beta/sigma^2)/2
    k1 <- a


    if (data$N == 2)
    {
      k1 <- k1/c1
      I1 <- pSqrtInvGamma(sigma, 0.5, beta)
      I2 <- I
    }
    else
    {
      # N == 1

      G <- function(sigma, beta){sigma*exp(-beta/sigma^2) + sqrt(4*beta*pi)*pnorm(sqrt(2*beta)/sigma)}

      I1 <- I
      I2 <- G(sigma, beta)
    }


    F1 <- k1 * diff(I1)
    F2 <-  b * diff(I2)
    F <- F1 + F2
      if (any(F < -1e-10))  stop('error in F//3')
      F <- pmax(F, 0) # correction for numeric imprecision


    j <- sample(seq(along=F), 1, prob=F)
      p01 <- I1[j]
      p02 <- I2[j]
      k1 <- k1[j]
      b <- b[j]

    target <- runif(1) * F[j]
    range <- sigma[j + c(0, 1)]

    sigma <- range[2] # starting point for N-R algorithm below


    continue <- TRUE

    while (continue)
    {
      if (data$N == 2)
      {
        F1 <- pSqrtInvGamma(sigma, 0.5, beta)
        f1 <- dSqrtInvGamma(sigma, 0.5, beta)

        F2 <- expInt_E1(beta/sigma^2) / 2
        f2 <- e1(beta, sigma)
      }
      else
      {
        F1 <- expInt_E1(beta/sigma^2) / 2
        f1 <- e1(beta, sigma)
        F2 <- G(sigma, beta)
        f2 <- exp(-beta/sigma^2)
      }


      F <- k1*(F1-p01) + b*(F2-p02)
      f <- k1*f1 + b*f2

      change <- - (F - target) / f

      i <- ifelse(change > 0, 1, 2)
      range[i] <- sigma

      sigma <- sigma + change
      if (sigma > range[1] || sigma < range[2])  sigma <- mean(range)

      converged <- abs(change) < NR.epsilon
      continue <- !converged

      if (continue)
      {
        squeezed <- diff(range) < NR.epsilon
        continue <- !squeezed
      }
    }


    return(sigma)
  } # end of CPL.rnd.sigma


  expInt_E1 <- function(x)
  {
    # source: Abramowitz & Stegun (1972), Dover

    E1.polynomial <- function(x, a)
    {
      j <- length(a)
      powers <- seq(from=0, to=j-1)
      powers <- rep(powers, rep(length(x), j))
      x <- matrix(x^powers, ncol=j)
      z <- x %*% matrix(a, ncol=1)

      return(as.vector(z))
    } # end of E1.polynomial


    E1 <- rep(NA, length(x))


    w <- which(x <= 1)

    if (length(w) > 0)
    {
       # For x <= 1
      a <- c(-0.57721566, 0.99999193, -0.24991055, 0.05519968, -0.00976004, 0.00107857)
      E1[w] <- E1.polynomial(x[w], a) - log(x[w])
    }


    w <- which(x > 1)

    if (length(w) > 0)
    {
      # For x > 1
      a <- c(0.2677737343,  8.6347608925, 18.0590169730, 8.5733287401, 1)
      b <- c(3.9584969228, 21.0996530287, 25.6329561486, 9.5733223454, 1)

      E1[w] <- E1.polynomial(x[w], a) / E1.polynomial(x[w], b) / x[w] * exp(-x[w])
    }


    return(E1)
  } # end of expInt_E1
  

  invgamma.logpcum <- function(lim, alpha, beta) 
  { 
    if (length(lim) > 1) 
    { 
      invgamma.lim <- c(min(lim), max(lim)) 
      gamma.lim <- rev(1/invgamma.lim) 
      p1 <- pgamma(gamma.lim[1], alpha, beta, log.p=TRUE, lower.tail=TRUE)
      p2 <- pgamma(gamma.lim[2], alpha, beta, log.p=TRUE, lower.tail=FALSE)
      lower.tail <- p1 < p2 
      gamma.lim <- rev(1/lim) 
    } 
    else 
    { 
      gamma.lim <- 1/lim 
      gamma.mode <- (alpha-1)/beta 
      lower.tail <- gamma.lim < gamma.mode 
    } 

    log.pcum <- pgamma(gamma.lim, alpha, beta, log.p=TRUE, lower.tail=lower.tail)
    list(log.pcum=rev(log.pcum), lower.tail=!lower.tail) 
  } # end of invgamma.logpcum 
  
  
  Location <- function(theta, lim, z) 
  {  
    mu.lower <- lim$mu[1] 
    mu.upper <- lim$mu[2] 

    sigma.lower <- lim$sigma[1] 
    sigma.upper <- lim$sigma[2] 
    
    l.intersects <- (theta - mu.lower)/z  # left-side  intersects 
    r.intersects <- (theta - mu.upper)/z  # right-side intersects 


    # For each point in l.intersects and r.intersects, 
    # indicate whether it is:
    # a) ABOVE sigma$upper                    (1)
    # b) between sigma$lower and sigma$upper  (0)
    # c) BELOW sigma$lower                   (-1)
  
    tmp.lloc <- as.numeric(l.intersects > sigma.upper) - as.numeric(l.intersects < sigma.lower) 
    tmp.rloc <- as.numeric(r.intersects > sigma.upper) - as.numeric(r.intersects < sigma.lower) 

    
    list(l.lower = c(-1, tmp.lloc), 
         l.upper = c(tmp.lloc, 1), 
         r.lower = c(-1, tmp.rloc), 
         r.upper = c(tmp.rloc, 1)) 
  } # end of Location 
  
  
  logp.sample <- function(log.w) 
  { 
    p <- exp(log.w - max(log.w)) 
    sample(seq(along=p), 1, prob=p) 
  } # end of logp.sample 
  
  
  riskband.areas <- function(R, lim, z, theta, location)
  { 
    mu.lower <- lim$mu[1] 
    mu.upper <- lim$mu[2] 
      mu.range <- diff(lim$mu) 

    sigma.lower <- lim$sigma[1] 
    sigma.upper <- lim$sigma[2] 
      sigma.range <- diff(lim$sigma) 


    riskband.area <- rep(NA, R)
    dtheta <- diff(theta) 
  
    #d <- dtheta/sqrt(1+z^2)  # distance between regions
    v <- dtheta/z            # vertical distance 
    
    # Define shape for each riskband
    riskband.shape <- 41 - location$r.upper - 3*location$r.lower - 9*location$l.upper - 27*location$l.lower
  
  
    # Compute area for each riskband

    for (r in 1:R) 
    { 
      if (riskband.shape[r] == 4)
      {       
        H <- sigma.upper - (theta[r-1] - mu.upper)/z 
        B <- z*H 
        riskband.area[r] <- B*H/2
      } 
      else if (riskband.shape[r] == 5)
      {    
        b <- mu.upper + z*sigma.upper - theta[c(r, r-1)] # vector of length 2 
        h <- b/z 
        riskband.area[r] <- diff(b*h)/2
      } 
      else if (riskband.shape[r] == 7)
      {   
        H <- sigma.range 
        B <- mu.upper - theta[r-1] + z*sigma.upper 
        b <- B - z*H 
        riskband.area[r] <- (B+b)/2 * H
      } 
      else if (riskband.shape[r] == 8)
      {      
        H <- sigma.range 
        B <- dtheta[r-1] 
        b.prime <- z*H 
        h <- sigma.upper - (theta[r]-mu.upper)/z 
        b <- z*h 
        riskband.area[r] <- H*(B+b) - (H*b.prime + b*h)/2
      } 
      else if (riskband.shape[r] == 9)
      { 
        H <- sigma.range 
        B <- dtheta[r-1] 
        riskband.area[r] <- B*H
      } 
      else if (riskband.shape[r] == 31)
      { 
        B <- mu.range 
        h <- sigma.upper - (theta[r-1] - c(mu.lower, mu.upper))/z # vector of length 2 
        riskband.area[r] <- B*mean(h)
      } 
      else if (riskband.shape[r] == 32)
      {          
        B <- mu.range 
        H <- v[r-1] 
        b <- theta[r] - z*sigma.upper - mu.lower 
        h <- b/z 
        riskband.area[r] <- B*H - b*h/2
      } 
      else if (riskband.shape[r] == 34)
      { 
        B <- mu.range 
        H <- sigma.range 
        h <- (theta[r-1] - mu.lower)/z - sigma.lower 
        b <- h*z 
        riskband.area[r] <- B*H - b*h/2
      } 
      else if (riskband.shape[r] == 35)
      {      
        H <- sigma.range 
        B <- mu.range 
        sigma <- (theta[c(r-1, r)] - c(mu.lower, mu.upper))/z  # vector of length 2 
        h <- c(-sigma.lower, sigma.upper) + c(1, -1)*sigma     # vector of length 2 
        b <- h*z                                               # vector of length 2 
        riskband.area[r] <- B*H - sum(b*h)/2
      } 
      else if (riskband.shape[r] == 36)
      { 
        b <- theta[r] - z*sigma.upper - mu.lower 
        h <- sigma.upper - (theta[r-1]-mu.lower)/z 
        H <- sigma.range 
        B <- dtheta[r-1] 
        riskband.area[r] <- B*H - h*(B-b)/2
      } 
      else if (riskband.shape[r] == 41)
      {        
        B <- mu.range 
        H <- v[r-1] 
        riskband.area[r] <- B*H
      } 
      else if (riskband.shape[r] == 44)
      {   
        B <- mu.range 
        H <- v[r-1]  
        b <- mu.upper - theta[r-1] + z*sigma.lower 
        h <- b/z 
        riskband.area[r] <- B*H - b*h/2
      } 
      else if (riskband.shape[r] == 45)
      {        
        h <- (theta[c(r-1,r)] - mu.lower)/z - sigma.lower # vector of length 2
        b <- h*z                                          # vector of length 2 
        riskband.area[r] <- diff(b*h)/2
      } 
      else if (riskband.shape[r] == 61)
      {        
        riskband.area[r] <- mu.range * sigma.range
      } 
      else if (riskband.shape[r] == 62)
      {       
        B <- mu.range 
        H <- sigma.range 
        h <- sigma.upper - (theta[r] - mu.upper)/z 
        b <- z*h 
        riskband.area[r] <- B*H - b*h/2
      } 
      else if (riskband.shape[r] == 63)
      {       
        b <- theta[r] - z*c(sigma.lower, sigma.upper) - mu.lower # vector of length 2 
        H <- sigma.range 
        riskband.area[r] <- H*mean(b)
      } 
      else if (riskband.shape[r] == 71)
      {    
        h <- (theta[r] - c(mu.lower, mu.upper))/z - sigma.lower # vector of length 2 
        B <- mu.range 
        riskband.area[r] <- mean(h)*B
      } 
      else if (riskband.shape[r] == 72)
      {      
        H <- (theta[r]-mu.lower)/z - sigma.lower 
        B <- z*H 
        riskband.area[r] <- B*H/2
      } 
    } 
    
    
    riskband.area
  } # end of riskband.areas
  
  
  ##############################################################################
  # Start function
  
  
  if (length(cpl.prior) > 0) 
  { 
    riskband.prior.prob <- cpl.prior$geometry$riskband.prior.prob 
    
    lim <- list(mu   =cpl.prior$geometry$mu.lim, 
                sigma=cpl.prior$geometry$sigma.lim) 

      mu.lower <- lim$mu[1] 
      mu.upper <- lim$mu[2] 
    
    A <- cpl.prior$geometry$A 
    z <- cpl.prior$geometry$z 
  } 
  else 
  { 
    lim <- list(mu=c(mu.lower, mu.upper), sigma=c(sigma.lower, sigma.upper)) 
    
    z <- qnorm(quantile) 
  } 


  # A is on the same scale as data
  if (outcome.is.logNormally.distributed) theta <- log(A)
  else theta <- A

  
  # Make sure that each riskband delimited by terms in 'A'
  # do intersect with the domain (mu.lower, mu.upper) x (sigma.lower, sigma.upper)
  
  location <- Location(theta, lim, z)
  
  
  # Check that no riskband is empty (does not cross prior domain)
  
  not.empty <- location$l.upper > location$l.lower | location$l.upper == 0 | location$l.lower > location$r.lower 
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
    which.empty <- setdiff(which.empty, c(1, length(empty))) # Drop first and last riskbands
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
    
    msg <- c('The riskband(s) below do not intersect with domain defined by (mu.lower, mu.upper) X (sigma.lower, sigma.upper):', msg)
    stop(msg) 
  } 
  
                                   
  # Make sure that the number of riskbands and riskband.prior.prob length correspond
  
  R <- length(A) + 1 
  if (R != length(riskband.prior.prob))  
  { 
    msg <- paste("Number of riskband prior probabilities found in riskband.prior.prob should be ", R, '[length(A)+1]')
    stop(msg) 
  } 
  
  
  # Make sure that the values in riskband.prior.prob make sense
  
  undefined.riskband.prior.prob <- all(is.na(riskband.prior.prob)) 
  
  
  if (use.continuous.piecewise.linear.prior) 
  { 
    if (undefined.riskband.prior.prob) 
    { 
      if (!equally.probable.riskbands) stop('riskband.prior.prob must be provided or equally.probable.riskbands set to TRUE when use.continuous.piecewise.linear.prior is TRUE.')
      riskband.prior.prob <- rep(1/R, R) 
    } 
    
    if (any(is.na(riskband.prior.prob)) | any(riskband.prior.prob < 0)) 
    { 
      stop('Values in riskband.prior.prob should all be non-negative numbers.') 
    } 
    else 
    { 
      tmp <- sum(riskband.prior.prob) 
      if (tmp != 1) stop('Values in riskband.prior.prob must sum to 1.') 
    } 
    

    if (length(cpl.prior) > 0)  CPL <- cpl.prior 
    else                        CPL <- CPL.prior(riskband.prior.prob, lim, A, outcome.is.logNormally.distributed, z) 
  } 
  else 
  { 
    uniform.prior.on.mu.and.sigma <- FALSE
    
    riskband.area <- riskband.areas(R, lim, z, theta, location)

  
    if (undefined.riskband.prior.prob) 
    {  
      uniform.prior.on.mu.and.sigma <- !equally.probable.riskbands
      
      if (uniform.prior.on.mu.and.sigma) riskband.prior.prob <- riskband.area/sum(riskband.area)
      else                               riskband.prior.prob <- rep(1/R, R) 
    } 
    else if (any(is.na(riskband.prior.prob)) | any(riskband.prior.prob < 0)) 
    { 
      stop('Values in riskband.prior.prob should all be non-negative numbers.') 
    } 
    else 
    { 
      tmp <- sum(riskband.prior.prob) 
      if (tmp != 1) stop('Values in riskband.prior.prob must sum to 1.') 
    } 
    
    
    # Compute density in each riskband
                               
    if (uniform.prior.on.mu.and.sigma) 
    { 
      tmp <- 1 / ((mu.upper-mu.lower)*(sigma.upper-sigma.lower)) 
      riskband.prior.density <- rep(tmp, R)
      riskband.prior.logdens <- rep(0, R)
    } 
    else 
    { 
      riskband.prior.density <- as.numeric(riskband.prior.prob > 0) # 0/1 variable indicating non-null riskband probabilities
      
      riskbands <- which(riskband.prior.density > 0)
      r.ref <- riskbands[1]       # first rb with a non-null prior probability
      riskbands <- riskbands[-1]  # drop r.ref from list
    
      for (r in riskbands)
      { 
        riskband.prior.density[r] <- riskband.prior.prob[r] * riskband.area[r.ref] / (riskband.prior.prob[r.ref] * riskband.area[r])
      } 
      
      
      riskband.prior.density <- riskband.prior.density/sum(riskband.prior.density) # Standardize values in riskband.prior.density
      
      riskband.prior.logdens <- log(riskband.prior.density)
    } 
  } 
  
  
  ######################################################################################################################
  # The rest of the function is very similar to what was done in previous versions
  
  me <- any.me(me.sd.range, cv.range) # see if measurement error is desired

  data <- data.summary(y=y, lt=lt, gt=gt, interval.lower=interval.gt, interval.upper=interval.lt,
                       logNormal.distrn=outcome.is.logNormally.distributed, me.through.cv=me$through.cv) 


    if (use.continuous.piecewise.linear.prior)
    {
      data$alpha <- (data$N - c(1,2)) / 2
      alpha <- data$alpha
      if (data$N <= 2)  alpha[alpha<=0] <- NA # to prevent warnings when trying to compute gamma(0)
      data$gamma.alpha <- gamma(alpha)
    }


                          
  gen.y <- list() 
  true.values <- list() 
  
  default.inits <- Default.inits(data, outcome.is.logNormally.distributed, c(mu.lower, mu.upper), c(sigma.lower, sigma.upper))


  tau.alpha <- (data$N - 1)/2 
  
  # Make sure that initial values for mu are in allowed range
  
  mu.range <- mu.upper - mu.lower 
  
  w <- which(init.mu > mu.upper) 
  if (length(w) > 0) init.mu[w] <- mu.upper - mu.range/20 
  
  w <- which(init.mu < mu.lower) 
  if (length(w) > 0) init.mu[w] <- mu.lower + mu.range/20 
  
  
  # Prepare dens.gen.icdf objects
  
  
  if (me$any) o.tv <- truevalue.gen.object(me, outcome.is.logNormally.distributed) 
      
  
  
  if (!use.continuous.piecewise.linear.prior) 
  { 
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
  } 
  else if (me$any) 
  { 
    stop('The option use.continuous.piecewise.linear.prior is not available (yet) in presence of measurement error.\n') 
    # ICI voir ce qui se passe lorsque use.continuous.piecewise.linear.prior = TRUE
    # en presence d'erreur de mesure
  } 
  
  
  if (me$any & !me$known) o.me <- me.gen.object(me, outcome.is.logNormally.distributed, data$N) 

  
  # Prepare output objects with appropriate dimensions

  sample <- empty.matrices(n.iter, n.chains, me, c('mu', 'sd', 'riskband'))

  burnin <- list() 
  if (monitor.burnin) burnin <- empty.matrices(n.burnin, n.chains, me, c('mu', 'sd', 'riskband'))

  
  
  M.iter <- n.burnin + n.iter * n.thin 
  
  for (ch in seq(n.chains)) 
  { 
    saved.iter <- 0 
    
    # Set initial values
    
    mu <- init.mu[ch] 

    if (is.na(init.sigma[ch]) | init.sigma[ch] == 0) sigma <- default.inits$sigma 
    else sigma <- init.sigma[ch] 


    if (data$N == 1)
    {
      my.y <- ifelse(outcome.is.logNormally.distributed, data$log$y, data$y)

      if (my.y == mu)
      {
        if (use.continuous.piecewise.linear.prior)
        {
          # Find prior mode for x95, take sigma = midpoint of its range, compute mu on mode-line given sigma
          w <- which.max(cpl.prior$h)
          x95 <- cpl.prior$x95[w]

          sigma <- (x95 - c(mu.upper, mu.lower)) / z
            if (sigma[1] < sigma.lower)  sigma[1] <- sigma.lower
            if (sigma[2] > sigma.upper)  sigma[2] <- sigma.upper

          sigma <- runif(1, sigma[1], sigma[2])

          mu <- x95 - z*sigma
        }
        else
        {
          w <- which.max(riskband.prior.logdens) # Riskband with higher prior density

          # Take a value for x95 in that riskband, and then randomly pick a value for mu corresponding to that x95
          x95.lower <- mu.lower + z*sigma.lower
          x95.upper <- mu.upper + z*sigma.upper
          tmp <- c(x95.lower, theta, x95.upper)
          tmp <- tmp[w+c(0,1)]
          x95 <- runif(1, tmp[1], tmp[2])

          # Range of possible values of sigma for the above selected value of x95
          sigma <- (x95 - c(mu.upper, mu.lower)) / z
            if (sigma[1] < sigma.lower)  sigma[1] <- sigma.lower
            if (sigma[2] > sigma.upper)  sigma[2] <- sigma.upper

          sigma <- runif(1, sigma[1], sigma[2])

          # Range for mu
          mu <- tmp - z*sigma
            if (mu[1] < mu.lower)  mu[1] <- mu.lower
            if (mu[2] > mu.upper)  mu[2] <- mu.upper

          mu <- runif(1, mu[1], mu[2])
        }
      }
    }

    
    # Initialize measured values for subjects with censored values
    if (data$any.censored$any) gen.y <- y.gen.inits(data, mu, sigma, me$through.cv, outcome.is.logNormally.distributed) 
    
    if (me$any) me$parm <- me$init # initialize Measurement Error parameter value 
  
    
    for (iter in 1:M.iter) 
    { 
      # Sample true values (in presence of measurement error)
      if (me$any) true.values <- truevalues.gen(gen.y, data, mu, sigma, me, outcome.is.logNormally.distributed, o=o.tv) 
                       
      # Sample y values for censored observations
      if (data$any.censored$any) gen.y <- y.gen(true.values, data, sigma, me, outcome.is.logNormally.distributed, mu=mu) 
      
      
      # Compute data points sum and square sum
      moments <- out.logout.moments(me$any, outcome.is.logNormally.distributed, data, gen.y, true.values) 
            
      
      # Sample from f(sigma | mu)
      
      tau.beta <- (moments$sum2 - 2*mu*moments$sum + data$N*mu^2) / 2 
      
      if (use.continuous.piecewise.linear.prior) 
      { 
        sigma <- CPL.rnd.sigma(mu, tau.beta, CPL, data)
      } 
      else 
      { 
        if (uniform.prior.on.mu.and.sigma) 
        { 
          # The joint prior for (mu, sigma) is uniform over the whole rectangle domain, hence no riskband delimiting is necessary
          
          sigma.range <- c(sigma.lower, sigma.upper) 
          
          if (me$through.cv & !outcome.is.logNormally.distributed) 
          { 
            sigma <- sigma.truncatedData.gen(o.sigma, sigma.range, tau.beta, mu, sigma) 
          } 
          else 
          { 
            sigma <- sqrt.invertedGamma.gen(data$N, tau.beta, sigma.range, o=o.sigma) 
          } 
        } 
        else 
        {    
          sigma.cutoffs <- (theta - mu)/z 
          w <- which(sigma.cutoffs >= sigma.lower & sigma.cutoffs <= sigma.upper) 
          sigma.cutoffs <- unique(c(sigma.lower, sigma.cutoffs[w], sigma.upper)) 
          
          multiple.riskbands <- length(sigma.cutoffs) > 2
          if (multiple.riskbands)
          { 
            possible.riskbands <- c(w, max(w) + 1) # riskbands in which (mu, sigma) can fall
            log.w <- riskband.prior.logdens[possible.riskbands]
          } 
          
          
          if (me$through.cv & !outcome.is.logNormally.distributed) 
          { 
            A <- c(o.sigma$A, list(b=tau.beta, mu=mu)) 
  
            if (multiple.riskbands)
            { 
              r <- seq(along=possible.riskbands)
              riskband.wt <- rep(NA, length(r))
              for (j in r) riskband.wt[j] <- integrate(o.sigma$f, sigma.cutoffs[j], sigma.cutoffs[j+1], A=A)$value
              riskband.prob <- riskband.wt/sum(riskband.wt) * exp(riskband.prior.logdens[possible.riskbands])
              j <- sample(r, 1, prob=riskband.prob)
            } 
            else 
            { 
              j <- 1 
            } 
  
            sigma <- dens.gen.icdf(o.sigma, A, range=sigma.cutoffs[j+c(0,1)]) 
          } 
          else if (data$N <= 1) 
          {             
            A <- c(o.sigma$A, list(b=tau.beta))           
              
            if (multiple.riskbands)
            { 
              riskband.wt <- rep(NA, length(possible.riskbands))
              for (i in seq(along=riskband.wt)) riskband.wt[i] <- integrate(o.sigma$f, sigma.cutoffs[i], sigma.cutoffs[i+1], A=A)$value
   
              log.w <- log(riskband.wt) + log.w
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
              
            sigma <- dens.gen.icdf(o.sigma, A, range=range, start=start) 
          } 
          else 
          {  
            # size > 1
  
            if (tau.beta == 0) 
            { 
              if (multiple.riskbands)
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
                
              if (multiple.riskbands)
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
              tau <- qgamma(log.U, tau.alpha, tau.beta, log.p=TRUE, lower.tail=!log.pcum$lower.tail)
              sigma <- 1/sqrt(tau) 
            } 
          } 
        } 
      } 
            

      # Sample from f(mu | sigma)
      
      mu.cond <- list(mean = moments$sum / data$N, sd = sigma/sqrt(data$N))
      
      if (use.continuous.piecewise.linear.prior) 
      { 
         tmp <- CPL.rnd.mu(sigma, CPL, mu.cond)
           mu <- tmp$mu
           riskband <- tmp$riskband
      } 
      else 
      { 
        if (uniform.prior.on.mu.and.sigma) 
        { 
          # The joint prior for (mu, sigma) is uniform over the whole rectangle domain, hence no riskband delimiting is necessary
          
          mu.range <- c(mu.lower, mu.upper) 
          
          if (me$through.cv & !outcome.is.logNormally.distributed) 
          { 
            mu <- mu.truncatedData.gen(o.mu, mu.range, mu.cond$mean, sigma=sigma) 
          } 
          else 
          {           
            phi <- pnorm.logpcum(mu.range, mean=mu.cond$mean, sd=mu.cond$sd) 
            log.U <- runif.logp(phi$log.pcum, lower.tail=phi$lower.tail) 
            mu <- qnorm(log.U, mean=mu.cond$mean, sd=mu.cond$sd, log.p=TRUE, lower.tail=phi$lower.tail)
          } 
          
          # Define riskband in which (mu, sigma) is sitting
          tmp <- mu + z*sigma 
          riskband <- 1 + sum(tmp > theta)
        } 
        else 
        {  
          mu.cutoffs <- theta - z*sigma 
          w <- which(mu.cutoffs >= mu.lower & mu.cutoffs <= mu.upper) 
          mu.cutoffs <- unique(c(mu.lower, mu.cutoffs[w], mu.upper)) 
          
          if (length(w) > 0) riskband <- c(w, max(w) + 1) # riskbands in which (mu, sigma) can fall
          else riskband <- 1
          
          if (me$through.cv & !outcome.is.logNormally.distributed) 
          {         
            A <- c(o.mu$A, list(mu.mean=mu.cond$mean, s=sigma, s2=sigma^2)) 
          
            if (length(riskband) > 1)
            { 
              r <- seq(along=riskband)
              riskband.wt <- rep(NA, length(riskband))
              for (j in r) riskband.wt[j] <- integrate(o.mu$f, mu.cutoffs[j], mu.cutoffs[j+1], A=A)$value
              riskband.prob <- riskband.wt/sum(riskband.wt) * exp(riskband.prior.logdens[riskband])
              j <- sample(r, 1, prob=riskband.prob)
            } 
            else 
            { 
              j <- 1 
            } 
            
            mu <- dens.gen.icdf(o.mu, A, range=mu.cutoffs[j+c(0,1)]) 
          } 
          else 
          { 
            phi <- pnorm.logpcum(mu.cutoffs, mean=mu.cond$mean, sd=mu.cond$sd) 
            
            if (length(riskband) > 1)
            { 
              logp <- logp.from.logpcum(phi) + riskband.prior.logdens[riskband]
              j <- logp.sample(logp) 
            } 
            else 
            { 
              j <- 1 
            } 
            
            log.pcum <- phi$log.pcum[j+c(0,1)] # vector of length 2 
            log.U <- runif.logp(log.pcum, lower.tail=phi$lower.tail) 
            mu <- qnorm(log.U, mean=mu.cond$mean, sd=mu.cond$sd, log.p=TRUE, lower.tail=phi$lower.tail)
          } 
          
          riskband <- riskband[j] # define riskband in which (mu, sigma) is sitting
        } 
      } 
    
    
      # Sample Measurement Error from its posterior density
      
      if (me$any & !me$known) me$parm <- me.gen(o.me, me, data, gen.y, true.values) 
           
   
      # Save values only when iter# modulo thinning = 0

      if (iter <= n.burnin) 
      { 
        if (monitor.burnin) 
        { 
          burnin$mu[ch, iter] <- mu 
          burnin$sd[ch, iter] <- sigma
          burnin$riskband[ch, iter] <- riskband
          
          if (me$any & !me$known) burnin$me.parm[ch, iter] <- me$parm 
        } 
      } 
      else if ((iter-n.burnin)%%n.thin == 0) 
      { 
        saved.iter <- saved.iter + 1 
        sample$mu[ch, saved.iter] <- mu 
        sample$sd[ch, saved.iter] <- sigma
        sample$riskband[ch, saved.iter] <- riskband
        
        if (me$any & !me$known) sample$me.parm[ch, saved.iter] <- me$parm 
      } 
    } 
  } 
  
  
  # Estimate riskband posterior probability (from the frequency table of 'riskband' variable)
  
  riskband.posterior.prob <- matrix(0, nrow=n.chains, ncol=R)
  
  for (ch in n.chains) 
  { 
    for (r in 1:R)  riskband.posterior.prob[ch, r] <- mean(sample$riskband==r)
  } 
  
  riskband.posterior.prob <- drop(riskband.posterior.prob) # turn the matrix into a vector if n.chains == 1


  # Collect information to return in output  
  
  out <- list(quantile=quantile, z.quantile=z, 
              riskband.prior.prob=riskband.prior.prob, R=R, 
              use.continuous.piecewise.linear.prior=use.continuous.piecewise.linear.prior) 
              
  if (use.continuous.piecewise.linear.prior) 
  { 
    out$geo <- list(u=CPL$geometry$u, 
                    area=CPL$geometry$area,  
                    x95b=CPL$geometry$x95b,  
                    cos.theta=CPL$geometry$cos.theta, 
                    R=CPL$geometry$R) 
    
    CPL$geometry <- NULL 
    out$CPL.prior <- CPL 
  } 
  else 
  { 
    out$riskband.area <- riskband.area
    out$riskband.prior.density <- riskband.prior.density
  } 
        
              
  tmp <- out.sample(sample, burnin, n.chains, monitor.burnin, me, outcome.is.logNormally.distributed) 
  out$sample <- tmp$sample 
  if (monitor.burnin) out$burnin <- tmp$burnin 
    
  out$mcmc  <- list(n.chains=n.chains, n.burnin=n.burnin, n.iter=n.iter, n.thin=n.thin, monitor.burnin=monitor.burnin) 
  out$inits <- list(mu=init.mu, sigma=init.sigma) 
  out$data  <- list(y=y, lt=lt, gt=gt, interval.gt=interval.gt, interval.lt=interval.lt)
  out$parms <- list(A=A, mu.lower=mu.lower, mu.upper=mu.upper, sigma.lower=sigma.lower, sigma.upper=sigma.upper, 
                    outcome.is.logNormally.distributed=outcome.is.logNormally.distributed) 
                    
  out$riskband.posterior.prob <- riskband.posterior.prob
  out$sample$x95 <- out$sample$mu + z * out$sample$sd

  return(out)
} # end of SEG.riskband 


SEG.riskband.posterior.plot <- function(o, one.plot.per.page=is.pdf,
                                           plot=list(scatter.plot=TRUE, x95.hist=TRUE, mu.hist=TRUE, sigma.hist=TRUE),
                                           title.cex=1)
{ 
  # o: an output from SEG.riskband or CPL.jags

  dev.mfrow <- par('mfrow')

  is.pdf <- names(dev.cur())[1] == 'pdf'

  if (one.plot.per.page)  par(mfrow=c(1,1))
  else
  {
    nplots <- sum(unlist(plot))

    if (nplots == 1)  par(mfrow=c(1,1))
    else              par(mfrow=c(2,2))
  }

  
  if (plot$scatter.plot)
  {
    plot(o$sample$mu, o$sample$sd, type='p', pch='.', col='navy',
         xlab=expression(mu), ylab=expression(sigma))
    title(expression(paste("(", mu, ", ", sigma, "): scatter plot of a sample obtained from their joint posterior distribution"), sep=''), cex.main=title.cex)
  }


  if (plot$x95.hist)
  {
    color <- 'darkgoldenrod'
    hist(o$sample$x95, nclass=50, xlab='x95' , prob=TRUE,
         main='x95: histogram of values from a sample obtained\nfrom its posterior distribution',
         cex=title.cex)


      A.lo <- o$parms$mu.lower + o$z.quantile*o$parms$sigma.lower
      A.hi <- o$parms$mu.upper + o$z.quantile*o$parms$sigma.upper
      A <- c(A.lo, log(o$parms$A), A.hi)

      abline(v=A, col=color, lty=2)
      y.txt <- par('usr')[4] - par('cxy')[2]/2
      plot.xlim <- par('usr')[c(1,2)]
      Cat.lbls <- c('Cat 1', 'Cat 2', 'Cat 3', 'Cat 4', 'Cat 5')
      rb.midpoint <- (A[-1] + A[-6]) / 2

      for (i in 1:5)
      {
        if (A[i] >= plot.xlim[1] & A[i+1] <= plot.xlim[2])  text(rb.midpoint[i],            y.txt, Cat.lbls[i], adj=0.5, col=color)
        else if (A[i]   >= plot.xlim[1])                    text((A[i] + plot.xlim[2])/2,   y.txt, Cat.lbls[i], adj=0.5, col=color)
        else if (A[i+1] <= plot.xlim[2])                    text((A[i+1] + plot.xlim[1])/2, y.txt, Cat.lbls[i], adj=0.5, col=color)
      }

      legend.lbls <- paste(paste(Cat.lbls,
                           round(100*o$riskband.posterior.prob, digits=1), sep=': '), '%', sep='')
      legend('right', legend=legend.lbls, text.col=color, xjust=0.5, title='Posterior probabilities')
  }
  
  
  if (plot$mu.hist)
  {
    hist(o$sample$mu, nclass=50, xlab=expression(mu), prob=TRUE,
         main=expression(paste(mu, ': histogram of values obtained from (', mu, ', ', sigma, ') joint posterior distribution')),
         cex=title.cex)
  }
  
  
  if (plot$sigma.hist)
  {
    hist(o$sample$sd, nclass=50, xlab=expression(sigma), prob=TRUE,
         main=expression(paste(sigma, ': histogram of values obtained from (', mu, ', ', sigma, ') joint posterior distribution')),
         cex=title.cex)
  }

  par(mfrow=dev.mfrow)
} # end of SEG.riskband.posterior.plot 
