# Version 0.1 (Feb 2022)


CPL.prior <- function(riskband.prior.prob, 
                      lim=list(mu=c(mu.lower, mu.upper), sigma=c(sigma.lower, sigma.upper)), 
                      A=c(0.01, 0.1, 0.5, 1), 
                      outcome.is.logNormally.distributed=T, 
                      z=qnorm(quantile), quantile=0.95,
                        gm.min=exp(mu.lower), gm.max=exp(mu.upper),
                        mu.lower=-3, mu.upper=6.2,
                        gsd.min=1.05, gsd.max=4, 
                        sigma.lower=log(gsd.min), sigma.upper=log(gsd.max),
                        precision=1e-4)
{     
  # This function is used to offer an alternative to the piecewise-uniform prior distribution
  # generally used in the Riskband model.
  #
  # Arguments:
  # ----------
  # riskband.prior.prob: a vector of length equal to the number of riskbands
  #                      giving the prior probability of each riskband
  #
  # lim: a list with dimensions $mu and $sigma, each of length 2, giving the range (domain) for mu and sigma respectively
  # A: a vector giving the cut-offs delimiting the riskbands


  # --- Useful functions -----------------------------------------------------  


  better.result <- function(result1, result2)
  {
    b1 <- c(result1 < result2, T)
    b2 <- c(result1 > result2, T)
    
    wb1 <- which(b1)[1]
    wb2 <- which(b2)[1]
    
    wb1 <= wb2 # include equality as being better -- even if not 'better' strictly speaking
  } # end of better.result

  
  geometry <- function(riskband.prior.prob, lim, z, A, outcome.is.logNormally.distributed=T, max.deepness=1/3)
  {
    if (outcome.is.logNormally.distributed) my.A <- log(A)
    else my.A <- A
    
    theta <- atan(z)
    cos.theta <- cos(theta)
    R <- length(my.A) + 1

    mu.lim    <- lim$mu
    sigma.lim <- lim$sigma

    H <- diff(sigma.lim)
    L <- H / cos.theta
    Bt <- H*z               # base of corner triangle
    d <- Bt * cos.theta     # length of chord going from (mu, sigma) domain outer corner to triangle diagonal

    m <- L / d

      mu0 <- mu.lim[1]
      mu1 <- mu.lim[2]
      sigma1 <- sigma.lim[2]
      mu.top <- c(mu0, my.A - z*sigma1, mu1 - Bt)

    B <- diff(mu.top)
    b <- B*cos.theta

    area <- B*H
      w <- c(1,5)
      area[w] <- area[w] + Bt * H / 2

    u <- riskband.prior.prob / area
    p <- riskband.prior.prob / L
    hbar <- p / b  # meaningless for the left-end and right-end riskbands
              
        
    # regression coefficients that permit to calculate s1 from s2 in paralleliped riskbands
    
    s1.intercept <- 8/b * hbar
    s1.h0 <- -8/b
       
       
    # Calculate the floor level to respect max deepness
    # and the highest height (on either side of the riskband)
    # to be able to respect max deepness
    
    lo <- u * max.deepness    # max deepness = 1/3 of the piecewise-uniform level
    hi <- 4 * hbar - 3 * lo
      hi[c(1,5)] <- Inf # meaningless for the left-end and right-end riskbands
      
      # (Note that even in the absence of floor/deepness consideration, the height
      #  of the distrn fct at the border of any riskband is still constrained
      #  by the fact that it cannot go in the negative values [at the riskband mid-section 
      #  or on the other side of the riskband] to compensate.)

    
    # At the border between two riskbands, the lower value allowed for the distrn height
    # is the lower of the two adjascent riskand 'lo' levels
    
    ground <- pmin(lo[-R], lo[-1])
      
    lo[c(1,5)] <- 0    


    # See whether or not it is possible to design a distrn curve that does not go 'deep' in at least one riskband
    
    impossible <- any(hi[-1] < lo[-R])
    impossible <- impossible | any(hi[-R] < lo[-1])   
        
    # inner-riskband lower points 
    #   (for the calculation of density curves appropriateness, 
    #   adding a cost to curves going too low when compared to piecewise-uniform prior inner-riskbands)
    
    u.low <- pmin(u[-R], u[-1])
    u.low <- c(u.low[1], u[2], u.low[2], 
                         u[3], u.low[3],
                         u[4], u.low[4])

    # Region bounds

    x95b <- c(mu.lim[1] + z*sigma.lim[1], my.A, mu.lim[2] + z*sigma.lim[2])
    
 
    R1b1 <- c(x95b[1], x95b[2])
    R1b2 <- c(x95b[1] + c(0, Bt), x95b[2])
    
    
    inner.b <- rep(NA, 2*R-3)
      tmp.inner.b <- x95b[-c(1,R+1)]
      inner.b[2*seq(R-1)-1] <- tmp.inner.b
      # inner riskbands midpoints
      inner.b[2*seq(R-2)] <- c(mean(tmp.inner.b[c(1,2)]),
                               mean(tmp.inner.b[c(2,3)]),
                               mean(tmp.inner.b[c(3,4)]))
    
    
    R5b1 <- c(x95b[R], x95b[R+1])
    R5b2 <- c(x95b[R], x95b[R+1] - c(Bt, 0))
    
    
    endriskband.is.local.min <- rep(NA, 2)
    uniq.u <- u[!duplicated(u)]
    endriskband.is.local.min[1] <- uniq.u[1] < uniq.u[2]
    
    uniq.u <- rev(u)
    uniq.u <- uniq.u[!duplicated(uniq.u)]
    endriskband.is.local.min[2] <- uniq.u[1] < uniq.u[2]
    
    j <- c(1,5)
    s.fade.out <- u[j] * max.deepness / (d + b[j]) * c(1, -1)
         

    list(z=z, A=A, H=H, L=L, Bt=Bt, d=d, m=m, cos.theta=cos.theta,
         mu.lim=mu.lim, sigma.lim=sigma.lim, R=R,
         B=B, b=b, area=area, u=u, u.low=u.low,
         R1b1=R1b1, R1b2=R1b2, inner.b=inner.b,
         R5b1=R5b1, R5b2=R5b2, 
         riskband.prior.prob=riskband.prior.prob, hbar=hbar, p=p,
         endriskband.is.local.min=endriskband.is.local.min, s.fade.out=s.fade.out,
         hi=hi, lo=lo, ground=ground, impossible=impossible,        
         outcome.is.logNormally.distributed=outcome.is.logNormally.distributed,
         x95b=x95b, # useful for call to *plot fct only
         s1.intercept=s1.intercept, s1.h0=s1.h0, reverse=F, 
         max.deepness=max.deepness)
  } # end of geometry
  
  
  geometry.completeDefn <- function(g)
  {
    if (g$endriskband.is.local.min[2]) 
    {
      # Correction to higher entry level in riskband R-1 when riskband R is a local min
      # (to ensure a density curve going down at entry of riskband R, at a higher point than u[R])
      
      r <- g$R - 1
      u <- g$u[g$R]
      hbar <- g$hbar[r]
      
      new.hi <- 4*hbar - 3*u
      if (g$hi[r] > new.hi) g$hi[r] <- new.hi
    } 
    
    
    if (!g$impossible & g$endriskband.is.local.min[1])
    {
      d <- g$d
      b <- g$b[1]
      p <- g$p[1]
      
      k1 <- d^2/3 + d*b
      
      #           __
      # plateau  /
      
      s <- p / k1
      plateau.h <- s * d
      
      if (plateau.h > g$hi[2]) g$impossible <- T
    }
    
    g
  } # end of geometry.completeDefn

  
  geometry.reverse <- function(g)
  {
    g$reverse <- !g$reverse
    
    g$B  <- rev(g$B)
    g$b  <- rev(g$b)
      
    g$hbar    <- rev(g$hbar)
    g$hi      <- rev(g$hi)
    g$lo      <- rev(g$lo)
    g$ground  <- rev(g$ground)
    g$p       <- rev(g$p)
      
    g$u                         <- rev(g$u)
    g$u.low                     <- rev(g$u.low)
    g$riskband.prior.prob       <- rev(g$riskband.prior.prob)
    g$endriskband.is.local.min  <- rev(g$endriskband.is.local.min)
      
    g$s1.intercept <- rev(g$s1.intercept)
    g$s1.h0        <- rev(g$s1.h0)
    
    g$inner.b <- -rev(g$inner.b)
    g$x95b    <- -rev(g$x95b)
    

          
    tmp.R5b1 <- g$R5b1
    tmp.R5b2 <- g$R5b2
      
    g$R5b1 <- -rev(g$R1b1)
    g$R5b2 <- -rev(g$R1b2)
      
    g$R1b1 <- -rev(tmp.R5b1)
    g$R1b2 <- -rev(tmp.R5b2)
    
      
    g
  } # end of geometry.reverse 
  

  smooth <- function(g, epsilon=precision)
  {
    # -- Useful fcts ---------------------------------------------------------
    
    angle.between.slope.changes <- function(s)
    {
      s.module <- sqrt(1+s^2)
        
      S <- length(s)
      scalar.product <- 1 + s[-1]*s[-S]
      tmp <- pmin(scalar.product / s.module[-1] / s.module[-S], 1) # protection against numerical imprecision
  
      acos(tmp)
    } # end of angle.between.slope.changes 
    
    
    augmented.body <- function(o, appendix, right.side=T) 
    {
      trim.appendix <- length(appendix$h) > length(appendix$s)

      if (right.side)
      {
        if (trim.appendix) appendix$h <- appendix$h[-1]
        
        o$h <- c(o$h, appendix$h)
        o$s <- c(o$s, appendix$s)
      }
      else
      {
        if (trim.appendix) o$h <- o$h[-1] # we trim o$h [same result as trimming new.h]

        o$h <- c(appendix$h, o$h)
        o$s <- c(appendix$s, o$s)
      }
      
      o
    } # end of augmented.body
  
  
    best.guess <- function(h0, reg4)
    {
      n <- reg4$n
      
      if (n > 2)
      {
        # Do a linear extrapolation
        
        b.denom <- n*reg4$x2.sum - reg4$x.sum**2
        
        if (b.denom == 0)
        {
          theta <- Inf
          margin <- 0 # irrelevant
        }
        else
        {
          b <- (n*reg4$xy.sum - reg4$x.sum*reg4$y.sum) / b.denom
          a <- (reg4$y.sum - b*reg4$x.sum) / n
          s <- a + b*h0 
          theta <- atan(s)
        
          sd2 <- (reg4$y2.sum - 2*a*reg4$y.sum - 2*b*reg4$xy.sum + n*a*a + 2*a*b*reg4$x.sum + reg4$x2.sum*b*b) / (n-1)
          if (sd2 < 0) sd2 <- 0 # protection against numerical imprecision
          
        
          if (sd2 > 0)
          {
            sd <- sqrt(sd2)
            s.upper <- s + sign(s)*sd
            margin <- abs(atan(s.upper) - theta)
          }
          else margin <- 0.01
        }
      }
      else
      {
        theta <- Inf
        margin <- 0 # irrelevant
      }
      
      
      list(init=theta, margin=margin)
    } # end of best.guess 
              
    
    closing.riskband <- function(g, G, trimmed.body, foldable=F)
    {        
      # G: closing riskband geometrics
      
      h0 <- ifelse(G$is.last.riskband, trimmed.body$h[length(trimmed.body$h)], trimmed.body$h[1]) # the height of the density on the inner side of
      # the closing riskband


      if (h0 >= G$h0.SE | (h0 >= G$u & G$is.local.min))
      {
        clothesLine.ok <- h0 <= G$h0.SE
        lb <- line.break(g, G, h0)
        
        choose.lb <- T
        
        if (clothesLine.ok)
        {
          cl <- clothesLine(G, h0)
                      
          out.cl <- roughness(g, trimmed.body, cl)
          out.lb <- roughness(g, trimmed.body, lb)
          
          choose.lb <- better.result(out.lb, out.cl)
          
          
          if (!choose.lb)
          {
            if (G$is.last.riskband)  x95b <- g$R5b1
            else                     x95b <- g$R1b1
            
            return(list(h=cl$h, s=cl$s, x95b=x95b, is.last.riskband=G$is.last.riskband))
          }
        }
 
        
        # choose.lb = TRUE
        
        return(list(h=lb$h, s=lb$s, x95b=lb$x95b, is.last.riskband=G$is.last.riskband))
      }
      else if (!G$is.local.min)
      {
        if (foldable)  tmp <- slanted.tipi(g, G, h0, trimmed.body)
        else           tmp <- hill.top(G, h0)
        
        if (G$is.last.riskband)  x95b <- g$R5b2
        else                     x95b <- g$R1b2
        
        return(list(h=tmp$h, s=tmp$s, x95b=x95b, is.last.riskband=G$is.last.riskband))
      }
      else
      {
        lb <- line.break(g, G, h0)
        ht <- hill.top(G, h0)
            
          if (G$is.last.riskband)  x95b <- g$R5b2
          else                     x95b <- g$R1b2


        out.lb <- roughness(g, trimmed.body, lb)
        out.ht <- roughness(g, trimmed.body, ht)                            
        
        if (better.result(out.lb, out.ht))  return(list(h=lb$h, s=lb$s, x95b=x95b, is.last.riskband=G$is.last.riskband))
        else                                return(list(h=ht$h, s=ht$s, x95b=x95b, is.last.riskband=G$is.last.riskband))
      }
    } # end of closing.riskband
    
                     
    closing.riskband.flexible <- function(g, G, h0, trimmed.rd)
    {
      # G: closing riskband geometrics
      # trimmed.rd:  a list with dimensions (h, s)
      #
      # NOTE: This fct is called for riskbands that are 'endriskband.is.local.min=T' ONLY

      convex <- h0 < G$h0.SE
      s1.b <- s1b(G, h0) 
      
      
      # ----------------------------------------------------------------------
      # We are searching for optimal (s1, s2, t)
      
      # Define range on s1 (and its initial value)
      
      if (!convex)
      {
        truncating.angle <- truncation(G, h0, g$epsilon)$angle
        
        rb5.range <- c(-pi/2, truncating.angle)
        rb5 <- bs.init(rb5.range, g$epsilon, margin=0.015, start.hi=T, untouchable.bounds=T, force.1pass=T) 
      }
      else
      {
        s.clothesLine <- clothesLine(G, h0)$s
        if (!G$is.last.riskband) s.clothesLine <- -s.clothesLine
        
        rb5.range <- c(atan(s.clothesLine), pi/2) 
        rb5 <- bs.init(rb5.range, g$epsilon, margin=0.015, start.lo=T, untouchable.bounds=T, force.1pass=T) 
      }
      
 

      while (rb5$continue)
      {
        s1 <- tan(rb5$value)  # stands for s5.a, really     
     
        # Find (s2, t) such that the integral \int_0^D h(x)l(x)dx = riskband prob
        #              [with  theta2 = atan(s2)]
        
        tmp <- NewtonRaphson.lastriskband.ts2(s1, h0, G, g$epsilon, convex=convex, s1.b=s1.b)     

        s2 <- tmp$s2
        t  <- tmp$t
        
        s5 <- c(s1, s2)
        h5 <- c(tmp$h.intersect, 0)
        
        closing.rb <- closing.riskband.parms(h5, s5, t, G$is.last.riskband)
        
        out <- roughness(g, trimmed.rd, closing.rb)
              
        rb5 <- bs.update(rb5, out, closing.rb)
      } # end of while-5 
      
      
      # Closing riskband parms
      
      parms <- rb5$monitor
  
      if (G$is.last.riskband)  x95b <- c(g$R5b1[1], g$R5b1[1] + parms$t/g$cos.theta, g$R5b1[2]) 
      else                     x95b <- c(g$R1b1[1], g$R1b1[2] - parms$t/g$cos.theta, g$R1b1[2])
      
      parms$x95b <- x95b
      parms
    } # end of closing.riskband.flexible
    
    
    closing.riskband.geo <- function(g, is.last.riskband)
    {
      r <- ifelse(is.last.riskband, g$R, 1)
      
      b <- g$b[r]
      p <- g$p[r]
      d <- g$d
      D <- d + b
      
      b2.2 <- b^2/2
      
      # South-East & South-South-East heights in closing riskbands
    
      h0.SSE <- 2 * p / b
      h0.SE  <- p * D / (b2.2 + b*d + d^2/3)
      
      endriskband.index <- ifelse(is.last.riskband, 2, 1)
      is.local.min <- g$endriskband.is.local.min[endriskband.index]
      
      
      list(b=b, d=d, D=D, L=g$L, p=p, prob=g$riskband.prior.prob[r], u=g$u[r],
           is.local.min=is.local.min, h0.SE=h0.SE, h0.SSE=h0.SSE,
           b2.2 = b2.2, b3.6 = b^3/6, pd = p*d, d2.6 = d^2/6, 
           is.last.riskband=is.last.riskband)
    } # end of closing.riskband.geo
    
    
    closing.riskband.parms <- function(h, s, t, is.last.riskband)
    {        
      if (!is.last.riskband)
      {
        h <-  rev(h)
        s <- -rev(s)
      }
      
      list(h=h, s=s, t=t, is.last.riskband=is.last.riskband)
    } # end of closing.riskband.parms
    
         
    clothesLine <- function(G, h0)
    {
      # G: closing riskband geometrics
      
      b <- G$b
      D <- G$D
    
      a <- G$b2.2 + D*(D+b) / 2 - (D^3 - b^3) / (3*G$d)
      delta <- h0 * (D+b) / 2 
      s <- (G$p - delta) / a
      
      h <- h0 + s*D
      if (!G$is.last.riskband) s <- -s
      
      list(h=h, s=s, is.last.riskband=G$is.last.riskband)
    } # end of clothesLine
    
    
    drop.closing.riskband <- function(rd, drop.last.riskband=T)
    {
      n2drop <- ifelse(drop.last.riskband, length(rd$R5b)-1, length(rd$R1b)-1)
      
      if (drop.last.riskband)
      {  
        j2drop <- seq(to=length(rd$s), length=n2drop)
        rd$s <- rd$s[-j2drop]
        
        j2drop <- seq(to=length(rd$h), length=n2drop)
        rd$h <- rd$h[-j2drop]
      }
      else
      {
        j2drop <- seq(n2drop)
        rd$s <- rd$s[-j2drop]
        rd$h <- rd$h[-j2drop]
      }
      
      rd
    } # end of drop.closing.riskband
    
    
    drop.elements <- function(x, n2drop, tail=T)
    {
      if (tail) w <- seq(to=length(x), length=n2drop)
      else      w <- seq(n2drop)
      
      x[-w]
    } # end of drop.elements
    
    
    hill.top <- function(G, h.bottom)
    { 
      # G: closing riskband geometrics
        
      b <- G$b
      L <- G$L
      a <- L * G$d / 2
        
      k <- h.bottom * (L*b + a)
      s.coeff <- L* b^2 / 2 + a*b
      
      s <- (G$prob - k) / s.coeff
      sb <- s*b
      
      if (G$is.last.riskband)  s <- c(s,  0)
      else                     s <- c(0, -s)

      h <- rep(h.bottom + sb, 2)
      
      list(h=h, s=s, is.last.riskband=G$is.last.riskband)
    } # end of hill.top
              
      
    line.break <- function(g, G, h0, closer2vertical.fraction=0.25)
    {
      # Draw a broken line (into two linear pieces) from the next-to-last riskband height to 0 at the upper/outer value for x95;
      # the line break occurs at 'b' or before if necessary.
      #
      # G: closing riskband geometrics
      
      b <- G$b
      d <- G$d
      D <- G$D

  
      if (h0 < G$h0.SSE)
      {
        # Fix occurs in the upper triangle => we will break the line @ the x95 line going through
        #                                                                      inner corner (b)
        
        t <- c(b, D)
        a2 <- diff(t^2)/2
        a3 <- diff(t^3)/3
  
        a <- d*D^2 -2*D*a2 + a3
          
        s.coeff <- -b*d/2 - a/d
        k <- h0*b/2
        s <- (G$p - k) / s.coeff
          
        h.star <- -s*d
        s.star <- (h.star - h0) / b
        x <- b
      }
      else
      {
        # Fix occurs in the parallelepiped => we will break the line on a secant that
        # is 25% closer to the vertical (dropping to 0 on x=max(A)) than the truncating line
          
        truncating.angle <- truncation(G, h0, g$epsilon)$angle
          # slope of the truncating line  (we want to avoid truncating, hence the new slope 
          #                                s.star calculated below)
                     
        angle  <-  -pi/2 *closer2vertical.fraction + (1-closer2vertical.fraction)*truncating.angle
        s.star <-  tan(angle) 
        
        
        # Find where the line passing through (0, h0) with slope s.star and
        # line passing through (D,0) with slope s meet (call it x)
        #
        # Find x & s using a Newton-Raphson algorithm
        
        tmp <- NewtonRaphson.lastriskband.ts2(s.star, h0, G, g$epsilon, t.lt.b=T)
        
        s <- tmp$s2
        x <- tmp$t
        h.star <- tmp$h.intersect          
      }
      
      
      # Wrap up        
  
      h <- c(h.star, 0)
      s <- c(s.star, s)
      
      
      if (G$is.last.riskband)
      {
        lb <- g$R5b1[1] + x/g$cos.theta
        x95b <- c(g$R5b1[1], lb, g$R5b1[2])
      }
      else
      {
        h <-  rev(h)
        s <- -rev(s)
        
        lb <- g$R1b1[2] - x/g$cos.theta
        x95b <- c(g$R1b1[1], lb, g$R1b1[2])
      } 
      

      list(h=h, s=s, x95b=x95b, is.last.riskband=G$is.last.riskband)
    } # end of line.break
    
    
    NewtonRaphson.lastriskband.ts2 <- function(s1, h0, G, epsilon, convex=F, s1.b=0,
                                               t.lt.b = h0 > G$h0.SSE | ifelse(convex, s1 > s1.b, s1 < s1.b))
    { 
      # G: closing riskband geometrics
      
      D <- G$D
        
        
      if (t.lt.b)
      {      
        # --- Find t < b -----------------------------------------------------
             
        t <- 0
        continue <- T
           
        while (continue)
        {
          s2 <- (h0 + s1*t) / (t-D)
          s2.prime <- -(h0 + s1*D) / (t-D)^2
          
          v <- (D-t)^2 / 2
          v.prime <- t - D
      
          f <- h0*t + s1*t^2/2 - s2*v + s2*G$d2.6  
          f.prime <- h0 + t*s1 - s2.prime*v - s2*v.prime + s2.prime*G$d2.6
          
          change <- (G$p - f) / f.prime
          t <- t + change
          continue <- abs(change) >= epsilon
        }
        
        s2 <- (h0 + s1*t) / (t-D)
        h.intersect <- -s2 * (D-t)
        
        return(list(s2=s2, t=t, h.intersect=h.intersect))            
      }
      else
      {                  
        # --- Find t > b -----------------------------------------------------
        
        target <- G$pd + h0*G$b2.2 + s1*G$b3.6
        a1 <- h0 * D
        a2 <- (s1*D - h0) / 2
        a3 <- -s1/3
         
        t <- G$b
        continue <- T
         
        while (continue)
        {
          u <- (h0 + s1*t) / 3
          v <- (D - t)^2
           
          f <- a1*t + a2*t^2 + a3*t^3 + u*v
          f.prime <- a1 + 2*t*a2 + 3*t^2*a3 + s1*v/3 - 2*u*(D-t)
           
          change <- (target - f) / f.prime
          t <- t + change
          continue <- abs(change) >= epsilon
        }
        
        h.intersect <- h0 + s1*t
        s2 <- h.intersect / (t-G$D) 
        
        return(list(s2=s2, t=t, h.intersect=h.intersect))
      }
    } # end of NewtonRaphson.lastriskband.ts2 --------------------------------
   
    
    parallelepiped <- function(r, g, entry.angle, h0)
    {
      b <- g$b[r]
        
      s0 <- tan(entry.angle)
      
      h0.middle <- h0 + s0*b/2
      if (h0.middle < 0 & h0.middle > -1e-10) h0.middle <- 0
      
      s1 <- g$s1.intercept[r] + g$s1.h0[r] * h0 - 3*s0
      h0.end <- h0.middle + s1*b/2
        
      if (h0.end < 0 & h0.end > -1e-10)  h0.end <- 0
      
      list(h=c(h0, h0.middle, h0.end), s=c(s0,s1))
    } # end of parallelepiped
    
    
    parallelepiped.connexion <- function(r, g, h0, h1, left2right=T)
    {
      # h0, h1: height @ left/right ends of riskband r [in that order]
      # if left2right = T, return s0
      #    otherwise       return -s1
      
      b <- g$b[r]
      s0 <- (4*g$hbar[r] - 3*h0 - h1) / b
      
      if (left2right)
      {
        return(s0)
      }
      else
      {
        s1 <- 2/b * (h1-h0) + s0
        return(-s1)
      }  
    } # end of parallelepiped.connexion
    
    
    parallelepiped.connexion.t <- function(r, g, h0, h1, t)
    {
      delta.h <- h1 - h0
      v <- c(delta.h, g$p[r] - h0[1]*g$b[r])
      s.solve(t, g$b[r], v)
    } # end of parallelepiped.connexion.t
    
    
    parallelepiped.range <- function(r, g, h0)
    {
      hbar <- g$hbar[r]
      b <- g$b[r]
      
      s.max <- (4*hbar - 3*h0 - g$ground[r]) / b  # stalagmite
       
      reflexion <- parallelepiped.reflexion(r, g, h0)
      if (reflexion > g$hi[r+1]) s.min <- parallelepiped.connexion(r, g, h0, g$hi[r+1])
      else s.min <- 2 / b * (g$lo[r] - h0)  # stalactite
  
  
      s <- c(s.min, s.max)
  
      atan(s)
    } # end of parallelepiped.range
    
    
    parallelepiped.reflexion <- function(r, g, h)
    {
      4*g$hbar[r] - 2*g$lo[r] - h
    } # end of parallelepiped.reflexion
    
      
    reg4.init <- function()
    {
      list(n=0, x.sum=0, x2.sum=0, xy.sum=0, y.sum=0, y2.sum=0)
    } # end of reg4.init
    
    
    reg4.update <- function(reg4, s4_entry, h04)
    {
      # Prepare for the calculation of regression parms between x & y, where
      #   x: h04
      #   y: entry_angle
      
      entry_angle <- atan(s4_entry)
      
      reg4$n      <- reg4$n + 1
      reg4$x.sum  <- reg4$x.sum   + h04
      reg4$x2.sum <- reg4$x2.sum  + h04**2
      reg4$xy.sum <- reg4$xy.sum  + h04 * entry_angle
      reg4$y.sum  <- reg4$y.sum   + entry_angle
      reg4$y2.sum <- reg4$y2.sum  + entry_angle**2
  
      return(reg4)
    } # end of reg4.update
    
    
    revisit.closing.riskband <- function(g, G, rd)
    {
      # G: closing riskband geometrics
      # rd: revisited density
      #
      # NOTE: This fct is called for riskbands that are 'endriskband.is.local.min=T' ONLY
      
      rd <- drop.closing.riskband(rd, G$is.last.riskband)
       
        # height when entering closing riskband
        h0 <- ifelse(G$is.last.riskband, rd$h[length(rd$h)], rd$h[1])         
         
      rb5 <- closing.riskband.flexible(g, G, h0, rd)

      rd <- augmented.body(rd, rb5, G$is.last.riskband)

      if (G$is.last.riskband)  rd$R5b <- rb5$x95b
      else                     rd$R1b <- rb5$x95b
     
      rd
    } # end of revisit.closing.riskband
    

    revisit.closing2riskbands <- function(g, G, rd)
    {
      # G: closing riskband geometrics
      # rd: revisited density
      #
      # NOTE: This fct is called for the last two riskbands ONLY when the last one is 'endriskband.is.local.min=T' 
      
      is.last2riskbands <- G$is.last.riskband
      
      last2rb <- ifelse(is.last2riskbands, 2, 1)       
      
      
      # _____________________________________________________________________________________
      # Remove elements corresponding to revisited riskbands        
      
      rb5.nSegments <- length(rd$R5b) - 1
      
      nSegments.closing.rb <- ifelse(is.last2riskbands, rb5.nSegments, length(rd$R1b) - 1)
      
      
      # Drop h & s from last 2 riskbands (do NOT drop any elements from rd$inner.b)
      
      rd$h <- drop.elements(rd$h, nSegments.closing.rb+2, tail=is.last2riskbands)
      rd$s <- drop.elements(rd$s, nSegments.closing.rb+2, tail=is.last2riskbands)
                              
      # _____________________________________________________________________________________
      
      # h at entry in the riskband *preceding* the closing riskband
      h04 <- ifelse(G$is.last.riskband, rd$h[length(rd$h)], rd$h[1])   
      
      # Riskband 4
      
      n2lr <- ifelse(is.last2riskbands, g$R-1, 2) # next-to-last riskband
      
      b4 <- g$b[n2lr]
        h0b <- h04 * b4
        p4 <- g$p[n2lr]
        
        s4.straight <- (p4 - h0b) * 2 / b4**2
        theta4.straight <- atan(s4.straight)
        
           
      rb4.t <- bs.init(c(0, b4), g$epsilon, init=b4/2, margin=b4/20, untouchable.bounds=T, force.1pass=T)
  
    
      while (rb4.t$continue)
      {
        t4 <- rb4.t$value
        
        if (h04 > g$u[n2lr])
        {
          s4.max <- parallelepiped.connexion.t(n2lr, g, h04, G$u, t4)[1]
          s4.min <- (g$ground[n2lr] - h04) / t4  # s-ground
          
          if (s4.max < s4.min) s4.max <- s4.min # possible for the highest values of h04
          
          s4.range <- c(s4.min, s4.max)  
        }
        else
        {
          s4.max <- parallelepiped.connexion.t(n2lr, g, h04, 0, t4)[1]
          s4.range <- c(s4.straight, s4.max)
        }
        
        
        theta4.range <- atan(s4.range)
        rb4.theta <- bs.init(theta4.range, g$epsilon, theta4.straight, 0.015)
  
  
        while (rb4.theta$continue)
        {
          theta4 <- rb4.theta$value
                 
          
          s4a <- tan(theta4)
          s4b <- -2 / (b4 - t4)^2 * (h04*t4 + s4a/2*t4^2 + (h04+s4a*t4)*(b4-t4) - p4)
                      
          s4 <- c(s4a, s4b)
          
          h4.t <- h04 + t4*s4a
          
          h4 <- h4.t + c(0, (b4 - t4) * s4b)
           
          
          # --- Riskband # 5 -------------------------------------------------
                      
          rb4 <- closing.riskband.parms(h4, s4, 0, is.last2riskbands)
          
          body4 <- augmented.body(rd, rb4, is.last2riskbands)
              
          closing.rb <- closing.riskband(g, G, body4)
          
          # ------------------------------------------------------------------
          
          
          rb4 <- augmented.body(rb4, closing.rb, is.last2riskbands)
          
          if (is.last2riskbands) rb5.nSegments <- length(closing.rb$s)
          
          
          complete.dens <- augmented.body(rd, rb4, is.last2riskbands)
          
          out <- roughness(g, complete.dens, rb5.nSegments=rb5.nSegments, last2rb=last2rb)
          
          parms <- list(rb4=rb4, x95b=closing.rb$x95b, rb5.nSegments=rb5.nSegments)
  
          rb4.theta <- bs.update(rb4.theta, out, parms)  
        } # end of while-rb4.theta loop
        
        
        parms <- rb4.theta$monitor
  
        complete.dens <- augmented.body(rd, parms$rb4, is.last2riskbands)
  
        out <- roughness(g, complete.dens, rb5.nSegments=parms$rb5.nSegments, last2rb=last2rb)
          parms$t4 <- t4     
    
        rb4.t <- bs.update(rb4.t, out, parms)
      } # end of while-rb4.t loop
      
    
      parms <- rb4.t$monitor
      
      
      # Wrap-up
      
      rd <- augmented.body(rd, parms$rb4, is.last2riskbands)
      
    
      if (is.last2riskbands)
      {
        rd$R5b <- parms$x95b
        rd$rb5.nSegments <- parms$rb5.nSegments
        
        rd$inner.b[6] <- rd$inner.b[5] + parms$t4 / g$cos.theta
      }
      else
      {
        rd$R1b <-  parms$x95b
        
        rd$inner.b[2] <- rd$inner.b[3] - parms$t4 / g$cos.theta
      }

      
      rd
    } # end of revisit.closing2riskbands

    
    revisit.inner.riskbands <- function(g, rd)
    {
      # rd: revisited density
      
      # Revisit parallelepiped inner riskbands; that is, we keep h-values at both ends unchanged
      # but see if we can change the breakpoint location between the two lines (and their slopes accordingly)
      # to get a smoother curve.
      
    
      rb1.nSegments <- length(rd$R1b) - 1
      rb5.nSegments <- length(rd$R5b) - 1
      
      h.index <- rb1.nSegments + 2*seq(3) # indices for the middle-riskband height of riskbands 2,3,4 [length 3]
      s.indices <- rb1.nSegments + seq(6) # indices for slopes in riskbands 2-3-4 (2 by rb)  [length 6]
        s.index2 <- s.indices[1:2]
        s.index3 <- s.indices[3:4]
        s.index4 <- s.indices[5:6]
              
      h0 <- rd$h[h.index-1] # h as the entry [left-side] in each riskband [length 3]
    
      delta.h <- rd$h[h.index+1] - h0  # increase in h in each inner riskband [length 3]
    
      v <- matrix(c(delta.h, g$p[2:4]-h0*g$b[2:4]), ncol=2) # 3 x 2 matrix
      
      
      h4.exit <- rd$h[h.index[3]+1]
  
      t4.min <- ifelse(h0[3] < g$u[4] & h4.exit > g$u[4], ((h0[3] + rd$h[h.index[3]+1])*g$b[4] - 2*g$p[4]) / delta.h[3], 0)
  
      # inits
      s3 <- rd$s[s.index3]
      s4 <- rd$s[s.index4]
      h.middle <- rd$h[h.index]
      
      
      # rb # 2     
      
      rb2 <- bs.init(c(0, g$b[2]), g$epsilon, g$b[2]/2, 0.05, untouchable.bounds=T)
  
        
      while (rb2$continue)
      {
        t2 <- rb2$value

        M2 <- matrix(c(t2, g$b[2]-t2, t2*g$b[2]-t2**2/2, (g$b[2]-t2)**2/2), ncol=2, byrow=T)
        s2 <- as.vector(solve(M2) %*% v[1,])
        
          h.middle[1] <- h0[1] + t2*s2[1]
          
          rd$h[h.index[1]] <- h.middle[1]
          rd$s[s.index2]   <- s2 
          
          
        # rb # 3
        
        rb3 <- bs.init(c(0, g$b[3]), g$epsilon, g$b[3]/2, 0.05, untouchable.bounds=T)

          
        while (rb3$continue)
        {
          t3 <- rb3$value

          M3 <- matrix(c(t3, g$b[3]-t3, t3*g$b[3]-t3**2/2, (g$b[3]-t3)**2/2), ncol=2, byrow=T)
          s3 <- as.vector(solve(M3) %*% v[2,])
          
            h.middle[2] <- h0[2] + t3*s3[1]
            
            rd$h[h.index[2]] <- h.middle[2]
            rd$s[s.index3]   <- s3  
            
            
          # rb # 4
          
          rb4 <- bs.init(c(t4.min, g$b[4]), g$epsilon, g$b[4]/2, 0.05, untouchable.bounds=T)
            
             
          while (rb4$continue)
          {
            t4 <- rb4$value

            M4 <- matrix(c(t4, g$b[4]-t4, t4*g$b[4]-t4**2/2, (g$b[4]-t4)**2/2), ncol=2, byrow=T)
            s4 <- as.vector(solve(M4) %*% v[3,])
            
              h.middle[3] <- h0[3] + t4*s4[1]
              
              rd$h[h.index[3]] <- h.middle[3]
              rd$s[s.index4]   <- s4        


            s <- c(s2, s3, s4)
            parms <- list(s=s4, t=t4, h.middle=h.middle)
          
            out <- roughness(g, rd, rb5.nSegments=rb5.nSegments, inner.rb=4)
            
            rb4 <- bs.update(rb4, out, parms)                          
          }  # end of while-4
          
          
          parms <- rb4$monitor
          
          s4 <- parms$s
          t4 <- parms$t
          h.middle <- parms$h.middle 

          rd$s[s.indices] <- c(s2, s3, s4)
          rd$h[h.index] <- h.middle
          
          
          parms <- list(s3=s3, t3=t3, s4=s4, t4=t4, h.middle=h.middle)
        
          out <- roughness(g, rd, rb5.nSegments=rb5.nSegments, inner.rb=3)
          
          rb3 <- bs.update(rb3, out, parms)                          
        } # end of while-3
        
        
        parms <- rb3$monitor
        
        s3 <- parms$s3
        t3 <- parms$t3
        s4 <- parms$s4
        t4 <- parms$t4
        
        h.middle <- parms$h.middle
        
        s <- c(s2, s3, s4)
        rd$s[s.indices] <- s
        rd$h[h.index] <- h.middle
        

        parms <- list(s=s, t=c(t2,t3,t4), h.middle=h.middle)
       
        out <- roughness(g, rd, rb5.nSegments=rb5.nSegments, inner.rb=2)       
        
        rb2 <- bs.update(rb2, out, parms)
      } # end of while-2
      
      
      parms <- rb2$monitor
      
      rd$h[h.index]   <- parms$h.middle
      rd$s[s.indices] <- parms$s
      
      j <- 2*seq(3)
      rd$inner.b[j] <- rd$inner.b[j-1] + parms$t / g$cos.theta
        
      rd$out <- roughness(g, rd, rb5.nSegments=rb5.nSegments) # we want to recalculate roughness without local considerations
  
      rd
    } # end of revisit.inner.riskbands
        
  
    roughness <- function(g, parms, closing.rb=list(),
                          rb5.nSegments = ifelse(length(parms$h)%%2 == 1, 1, 2),
                          inner.rb=0, last2rb=0)
    {
      # parms: a list with dimensions(h, s)
      # closing.rb: provided if one of the closing riskbands is re-estimated 
      #             (that is, AFTER the curve was estimated throughout, to see if we can find a 
      #              better density closing).
      #             > When provided, closing.rb is a list with dimensions $s, $h & $is.last.riskband
      #
      # last2rb: a numeric code, = 2 indicates we are focusing on the  last two riskbands
      #                          = 1 indicates we are focusing on the first two riskbands
      
      
      s <- parms$s
      h <- parms$h
      
      
      closing.rb.provided <- length(closing.rb) > 0
       
      
      if (closing.rb.provided)
      {
        if (closing.rb$is.last.riskband)
        {  
          s <- c(s, closing.rb$s)
          h <- c(h, closing.rb$h)
        }
        else
        {                              
          s <- c(closing.rb$s, s)
          h <- c(closing.rb$h, h)
        }

        if (closing.rb$is.last.riskband)  rb5.nSegments <- length(closing.rb$s)
        # else we keep the value received as a fct argument
      }

      
      s.len <- length(s)
      h.len <- length(h)
             
  
      # Correct for (most likely) numerical imprecision
      if (abs(s[1])     < 1e-10) s[1]     <- 0
      if (abs(s[s.len]) < 1e-10) s[s.len] <- 0
      
  
      # Find the sharpest turn / slope change
  
      sharpest.slopeChange <- sharpest.slope.change(s)
  
  
      # Count the number of sign changes
      
      s.sign <- sign(s)
      s.sign <- s.sign[s.sign != 0]
      number.of.sSign.changes <- sum(s.sign[-1] != s.sign[-length(s.sign)]) 
      
      
      if (!g$endriskband.is.local.min[2] && rb5.nSegments > 1 && s[s.len-1] > 0 && s[s.len] < 0) number.of.sSign.changes <- number.of.sSign.changes - 1
      
      
      # Determine whether distribution is fully described or not
      
      rb1.nSegments <- s.len - 6 - rb5.nSegments
              
      dens.fct.is.fully.described <- ifelse(closing.rb.provided, length(parms$s) > 6, rb1.nSegments > 0)
      
          
      # See if density fct is going too deep (when compared to the piecewise uniform prior)

      if (!g$impossible)
      {        
        h.last.index <- h.len - rb5.nSegments
  
        h.first.index <- h.last.index - 6
        if (h.first.index < 1) h.first.index <- 1
        h.indices <- seq(from=h.first.index, to=h.last.index)
  
        ulow.indices <- seq(to=7, by=1, length=length(h.indices))
        too.deep <- any(h[h.indices] / g$u.low[ulow.indices] < g$max.deepness)
        
        
        # Final slope in closing riskbands should not be below some value, as indicated by s.fade.out 
        
        if (dens.fct.is.fully.described && h[1]     == 0 && s[1]     < g$s.fade.out[1])  too.deep <- T
        if (                               h[h.len] == 0 && s[s.len] > g$s.fade.out[2])  too.deep <- T
        
      }
      else too.deep <- -1
  
  
      # Consider initial/final value of density fct
  
      if (dens.fct.is.fully.described)
      {
        h.endpoints <- h[c(1, h.len)]     
        rb1.roughness <- ifelse(rb1.nSegments == 1, 0, abs(s[2]-s[1]))
        
        if (!g$endriskband.is.local.min[1] && rb1.nSegments > 1 && s[1] > 0 && s[2] < 0) number.of.sSign.changes <- number.of.sSign.changes - 1
      }
      else
      {
        h.endpoints <- c(-1, -1)          
        rb1.roughness <- -1
      }
      
      
      rb5.roughness <- ifelse(rb5.nSegments == 1, 0, abs(s[s.len]-s[s.len-1]))
             
      endriskbands.roughness <- max(rb1.roughness, rb5.roughness)
      
      
      # Consider the number of inflexion points
      
      diffs.sign <- sign(diff(s))
      diffs.sign <- diffs.sign[diffs.sign!=0]
            
      n.inflexion.points <- sum(diffs.sign[-1] != diffs.sign[-length(diffs.sign)])
  
      
      extra.scale <- 0
      
      
      if (last2rb==2)
      {
        j <- seq(to=length(s), length=rb5.nSegments + 3)
        extra.scale <- sharpest.slope.change(s[j])
        endriskbands.roughness <- -1
      }
      else if (last2rb==1)
      {
        j <- seq(rb1.nSegments + 3)
        extra.scale <- sharpest.slope.change(s[j])
        endriskbands.roughness <- -1
      }
      else if (closing.rb.provided)
      {
        if (closing.rb$is.last.riskband)  j <- seq(to=s.len, length=3)
        else                              j <- seq(3)
        
        extra.scale <- sharpest.slope.change(s[j])
        endriskbands.roughness <- -1
      }
      else if (inner.rb > 0)
      {
        j <- seq(from=rb1.nSegments+2*inner.rb-3, to=rb1.nSegments+6)
        extra.scale <- max(abs(s[j]))
      }
      else
      {
        extra.scale <- max(abs(s))
        #extra.scale <- sd(s)
      }
      
      
      c(number.of.sSign.changes, too.deep, sharpest.slopeChange, max(h.endpoints), 
        n.inflexion.points, endriskbands.roughness, extra.scale)
    } # end of roughness
        
    
    s1b <- function(G, h0)
    {
      # G: closing riskband geometrics
      
      b <- G$b
      d <- G$d
      
      s1b <- (G$p - h0*(b + d/3)) / (b^2/2 + b*d/3)
              
      s1b
    } # end of s1b
    
    
    slanted.tipi <- function(g, G, h0, body)
    {
      # G: closing riskband geometrics
      
      d <- G$d
      b <- G$b
      
      j <- ifelse(G$is.last.riskband, 1, 2)
      theta.min <- atan(hill.top(G, h0)$s[j])
      if (!G$is.last.riskband)  theta.min <- -theta.min
      
      theta.max <- southEast.connexion(G, h0)
      
      f <- 6/d**2
      a <- (G$p - d/2*h0 - h0*b) * f
      s1.mult <- (b**2/2 + b*d/2) * f
      
      
      rb <- bs.init(c(theta.min, theta.max), g$epsilon)
      
      while (rb$continue)
      {
        s1 <- tan(rb$value)
        s2 <- a - s1*s1.mult
        
        s <- c(s1, s2)
        hb <- h0 + b*s1
        h <- c(hb, hb + s2*d)
        
        closing.rb <- closing.riskband.parms(h, s, 0, G$is.last.riskband)
      
        out <- roughness(g, body, closing.rb)
        
        rb <- bs.update(rb, out, closing.rb)
      }
    
      
      rb$monitor
    } # end of slanted.tipi

         
    southEast.connexion <- function(G, h0)
    {
      # G: closing riskband geometrics
      
      b <- G$b
      d <- G$d
      
      s <- (G$p - h0*(d/3+b)) / (b**2/2 + b*d/3)
      atan(s)
    } # end of southEast.connexion      
    
    
    s.solve <- function(t, b, v)
    {
      M <- c(t, b-t, t*b-t**2/2, (b-t)**2/2)
      M.det <- M[1]*M[4] - M[2]*M[3]
      s1 <- (M[4]*v[1] - M[2]*v[2]) / M.det
      s2 <-  (v[1] - s1*M[1]) / M[2]
      
      c(s1, s2)
    } # end of s.solve
                  
    
    sharpest.slope.change <- function(s)
    {
      max(abs(angle.between.slope.changes(s)))
    } # end of sharpest.slope.change    
    
    
    theta4.limits <- function(g, h04)
    {
      r <- g$R - 1  # riskband r = 4
      
      s.max <- parallelepiped.connexion(r, g, h04, g$u[g$R])
              
      # Plateau
      s.plateau <- 8/3/g$b[r] * (g$hbar[r] - h04)
      
      s <- c(s.plateau, s.max)
      
      atan(s)
    } # end of theta4.limits
          
    
    truncation <- function(G, h0, epsilon)
    {
      # G: closing riskband geometrics
      # Note: this fct does NOT adjust the sign of 's' for G$is.last.riskband on purpose
      
         
      if (h0 > G$h0.SSE)
      {
        # truncation occurs on t < b
        t <- 2*G$p / h0
      }
      else
      {
        target <- G$pd + h0 * G$b2.2
        a0 <- h0 * G$b3.6
        a1 <- h0 * G$D / 2
        a2 <- -h0 / 6


        t <- G$b
        continue <- T

        while (continue)
        {
          f <- a0/t + a1*t + a2*t^2
          f.prime <- -a0/t^2 + a1 + 2*t*a2

          change <- (f - target) / f.prime
          t <- t - change
          continue <- abs(change) >= epsilon
        }
      }
      
      
      s <- -h0/t
    
      list(s=s, t=t, angle=atan(s))
    } # end of truncation
    
  
    # --------------------------------------------------------------------------
    # Start of smooth
    
    g$epsilon <- epsilon
    reg4 <- reg4.init()

    g.left  <- closing.riskband.geo(g, F)
    g.right <- closing.riskband.geo(g, T)
  
     
    rb2.h0 <- bs.init(c(g$ground[1], g$hi[2]), epsilon, g$u[2], abs(g$u[2]-g$u[1]), force.1pass=T)
  
    
    while (rb2.h0$continue)
    {
      h02 <- rb2.h0$value
      
      # Region 2
     
      tmpr2 <- parallelepiped.range(r=2, g, h02)
  
      rb2  <- bs.init(tmpr2, epsilon, force.1pass=T)
       
  
      while (rb2$continue)
      {
        tmp <- parallelepiped(r=2, g, rb2$value, h02) 
        s2 <- tmp$s
        h2 <- tmp$h[-1]
        body2 <- list(h=c(h02, h2), s=s2)
    
  
        # Region 3
  
        h03 <- h2[2]
        tmpr3 <- parallelepiped.range(r=3, g, h03)
  
        rb3  <- bs.init(tmpr3, epsilon, force.1pass=T)
        
  
        while (rb3$continue)
        {
          tmp <- parallelepiped(r=3, g, rb3$value, h03)
          s3 <- tmp$s
          h3 <- tmp$h[-1]
          body3 <- list(h=c(h03, h3), s=s3)
       
  
          # Region 4
  
          h04 <- h3[2]
          if (g$endriskband.is.local.min[2])  tmpr4 <- theta4.limits(g, h04)
          else                                tmpr4 <- parallelepiped.range(r=4, g, h04)
          
  
          bg4  <- best.guess(h04, reg4) 
          rb4  <- bs.init(tmpr4, epsilon, start=bg4, force.1pass=T)
        
  
          while (rb4$continue)
          {    
            tmp <- parallelepiped(r=4, g, rb4$value, h04)
            s4 <- tmp$s
            h4 <- tmp$h[-1]
            if (h4[2] < 0 & h4[2] > -1e-8) h4[2] <- 0 # correction for numerical imprecision
            body4 <- list(h=c(h04, h4), s=s4)
  
        
            # Region 5
            rb5 <- closing.riskband(g, g.right, body4)
              rb5$rb4 <- list(h=h4, s=s4)

            out4 <- roughness(g, body4, rb5)
            
            rb4 <- bs.update(rb4, out4, rb5)
          } # end of while-4 -------------------------------------------------
                      
          
          parms <- rb4$monitor  
                       
            reg4 <- reg4.update(reg4, rb5$rb4$s[1], h04)
            
          rb4 <- augmented.body(rb5$rb4, rb5)  
          tmp <- augmented.body(body3, rb4)
                               
          out3 <- roughness(g, tmp)
  
          rb3 <- bs.update(rb3, out3, tmp)
        } # end of while-3 ---------------------------------------------------
       
        
        parms <- rb3$monitor
        tmp <- augmented.body(body2, parms)
                 
        out2 <- roughness(g, tmp)
   
        rb2 <- bs.update(rb2, out2, tmp)
      } # end of while-rb2 ---------------------------------------------------
      
      
      parms <- rb2$monitor
        
          
      # Complete the density by closing the 1st riskband
      
      if (g.left$is.local.min)  rb1 <- closing.riskband.flexible(g, g.left, h02, parms)
      else                      rb1 <- closing.riskband(g, g.left, parms, foldable=T)        
      
      out <- roughness(g, parms, rb1)
      
      rb2.h0 <- bs.update(rb2.h0, out, list(parms=parms, rb1=rb1))
    } # end of while-rb2.h0
                      
      
    tmp <- rb2.h0$monitor
    
    best <- tmp$parms # parms for rb#2 and after
      best$rb5.nSegments <- ifelse(length(best$s)%%2 == 1, 1, 2)
      
      best <- augmented.body(best, tmp$rb1, F)
      best$inner.b <- g$inner.b
      best$R1b <- tmp$rb1$x95b
      
      best$R5b <- switch(best$rb5.nSegments, g$R5b1, g$R5b2)
             
      
    # Revisit the first/last two riskbands if necessary      
    
    # (first two riskbands)
    if (g.left$is.local.min)
    { 
      if (best$s[length(best$R1b)] < 0)  best <- revisit.closing2riskbands(g, g.left, best)
      best <- revisit.closing.riskband(g, g.left, best)
    }
    else
    {
      # Revisit first riskband parms for a better fit
      best <- drop.closing.riskband(best, F)
      rb1 <- closing.riskband(g, g.left, best, foldable=T)
      best <- augmented.body(best, rb1, F)
      best$R1b <- rb1$x95b
    }
    

    
    # (last two riskbands)
    if (g.right$is.local.min)
    {
      if (best$s[length(best$s)+1-length(best$R5b)] > 0)  best <- revisit.closing2riskbands(g, g.right, best)
      best <- revisit.closing.riskband(g, g.right, best)
    }
    else
    {
      # Revisit last riskband parms for a better fit
      best <- drop.closing.riskband(best)
      rb5 <- closing.riskband(g, g.right, best, foldable=T) 
      best <- augmented.body(best, rb5)
      best$R5b <- rb5$x95b
    }
    

    
    # Revisit the inner riskbands breakpoints & their corresponding heights
          
    best <- revisit.inner.riskbands(g, best) # out is calculated herein for the last time
    
    
    # Rearrange best's dimensions in sensible order
    
    if (g$reverse)
    {
      best$h <-  rev(best$h)
      best$s <- -rev(best$s)
          
      tmp      <-  best$R5b
      best$R5b <- -rev(best$R1b)
      best$R1b <- -rev(tmp)
      
      best$inner.b <- -rev(best$inner.b)
    }
            
    
    best$r <- rep(seq(5), c(length(best$R1b)-1, 2, 2, 2, length(best$R5b)-1))
    best$x95 <- c(best$R1b, best$inner.b[-1], best$R5b[-1])
    best$reverse <- g$reverse
    
    best    
  } # end of smooth
    

  # --------------------------------------------------------------------------
  # Start of CPL.prior
      
  
  # Prepare geometry
  
  g <- geometry(riskband.prior.prob, lim, z, A, outcome.is.logNormally.distributed)
    g.rev <- geometry.reverse(g)
  
    g     <- geometry.completeDefn(g)    
    g.rev <- geometry.completeDefn(g.rev)
                
    
  # Actually estimate the best density curve
  
  best <- smooth(g)

  tmp <- smooth(g.rev)
    if (better.result(tmp$out, best$out)) best <- tmp
    
    
  # Wrap-up

  g$sin.theta <- sqrt(1-g$cos.theta^2)
  best$geometry <- g
  
  
  best
} # end of CPL.prior

       
CPL.prior.plot <- function(o, color='purple')
{
  # o: an output of SEG.riskband
  
  if (!o$use.continuous.piecewise.linear.prior) stop("This output/object was not obtained with use.continuous.piecewise.linear.prior=TRUE\n")
  
  g <- o$geo

  h <- o$CPL.prior$h
  s <- o$CPL.prior$s * g$cos.theta
  x95 <- o$CPL.prior$x95
  R <- g$R
  
  
  # Draw the piecewise uniform prior distrn
  
  u <- g$u
  bounds <- g$x95b
    
  x <- c(bounds[1], rep(bounds[-c(1, R+1)], rep(2, R-1)), bounds[R+1])
  f <- rep(u, rep(2, R))
  ylim <- c(0, max(c(f,h)))
  plot(x, f, type='l', col='gray', ylim=ylim, xlab=expression(x[95]), ylab=expression(f(x[95])))
    title('Prior distribution')
    str <- paste(o$riskband.prior.prob, collapse=', ')
    str <- paste(c('riskand prior prob = c(', str, ')'), collapse='')
    mtext(str, side=3, at=par('usr')[2], adj=1, line=0.5)


  # Overlay the curve as indicated by (h, s)

  h.len <- length(h)
  interval.len <- diff(x95)

  h1 <- h[-h.len] + s * interval.len
  
  i.seq <- seq(along=h)
  

  for (i in seq(along=i.seq)[-h.len])
  {
    j <- i.seq[i]
    next.j <- i.seq[i+1]
    points(x95[c(j,next.j)], c(h[j], h1[j]), type='l', col=color)
  }


  # Plot a legend
  
    legend.lty <- c(1,1)
    legend.txt <- c('continuous piecewise linear [CPL]', 'piecewise uniform [for reference only]') 
    legend.size <- legend('bottom', legend.txt, lty=legend.lty, bty='n', plot=F)$rect # do not plot -- just note legend size
  
    # Choose the best location for the legend -- pick a spot from which the graphic points are as far as possible
    
    par.usr <- par('usr')
    par.pin <- par('pin')
    
    # Points per inches, on both dimensions
    x.ppi <- diff(par.usr[c(1,2)]) / par.pin[1]
    y.ppi <- diff(par.usr[c(3,4)]) / par.pin[2]
    
    x <- c(par.usr[1] + legend.size$w/2, mean(par.usr[c(1,2)]), par.usr[2] - legend.size$w/2)
    x <- rep(x, 2)
    y <- c(par.usr[3] + legend.size$h, par.usr[4] - legend.size$h)
    y <- rep(y, rep(3,2))
    ref <- data.frame(x=x, y=y, loc=c('bottomleft', 'bottom', 'bottomright', 'topleft', 'top', 'topright'))
    n.ref <- nrow(ref)
    d <- rep(Inf, n.ref)
    d2pt <- rep(NA, 2)
    
    for (i in seq(along=s))
    {
      x <- c(x95[i], x95[i+1])
      y <- c(h[i], h[i] + s[i]*diff(x))
      
      m <- -1/s[i] * (y.ppi/x.ppi)^2 # visual perpendicular slope
      a <- h[i] - s[i] * x[1] # intercept
      
      
      for (j in 1:n.ref)
      {
        for (k in 1:2) d2pt[k] <- sqrt((x[k]-ref$x[j])^2/x.ppi^2 + (y[k]-ref$y[j])^2/y.ppi^2)
        
        a.perpendicular <- ref$y[j] - m * ref$x[j] # intercept of the (visually) perpendicular line
        
        # coordinates of the point on the line that is closest to point ref[j]
  
        line.x <- ifelse(s[i]==0, ref$x[j], (a.perpendicular - a) / (s[i] - m))
        
        if (x[1] < line.x & line.x < x[2])
        {
          line.y <- a + s[i] * line.x
          this.d <- sqrt((line.x-ref$x[j])^2/x.ppi^2 + (line.y-ref$y[j])^2/y.ppi^2)
        }
        else this.d <- min(d2pt)

        d[j] <- min(d[j], this.d)
      }
    }
  
    w <- which.max(d)
    loc <- ref$loc[w]
  
  
  legend(loc, legend.txt, col=c(color, 'gray'), lty=legend.lty, bty='n')
} # end of CPL.prior.plot
