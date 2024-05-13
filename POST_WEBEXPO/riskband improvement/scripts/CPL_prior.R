source("c:/users/patri/home/bin/R/fcts/bisectionalSearch.R")


# Version 1.1 (May 2024)


# -----------------------------------


CPL.prior <- function(riskband.prior.prob, 
                      lim=list(mu=c(mu.lower, mu.upper), sigma=c(sigma.lower, sigma.upper)), 
                      A=c(0.01, 0.1, 0.5, 1), 
                      outcome.is.logNormally.distributed=TRUE,
                      z=qnorm(quantile), quantile=0.95,
                      # gm.min=exp(mu.lower), gm.max=exp(mu.upper), # ICI voir avec Jerome si l'usager moyen prefererait
                            # entrer son prior en utilisant gm.min & gm.max plutot que (mu.lower, mu.upper) ou bien lim$mu

#                      mu.lower=-3, mu.upper=6.2, # ICI peut-etre changer ces valeurs par defaut pour avoir les memes que dans la version JavaScript
                      mu.lower=-7.600902, mu.upper=1.609438,

                      gsd.min=1.05, gsd.max=4, 
                      sigma.lower=log(gsd.min), sigma.upper=log(gsd.max),
                      precision=1e-6)
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


  angles <- function(s)
  {
    S <- length(s)
    s.module <- sqrt(1+s^2)

    s.j <- s[-S]
    s.module.j <- s.module[-S]
    # following slope
    s.k <- s[-1]
    s.module.k <- s.module[-1]

    scalar.product <- 1 + s.j*s.k
    tmp <- pmin(scalar.product / s.module.j / s.module.k, 1)  # protection against numerical imprecision
    tmp <- acos(tmp)

    tmp
  } # end of angles


  better.score <- function(result1, result2, negligeable.diff=1e-12)
  {    
    w <- which(abs(result1-result2) > negligeable.diff)
    
    if (length(w) > 0)
    {
      w <- w[1]
      return(result1[w] < result2[w])
    }


    b1 <- c(result1 < result2, TRUE)
    b2 <- c(result1 > result2, TRUE)
    
    wb1 <- which(b1)[1]
    wb2 <- which(b2)[1]
    
    wb1 <= wb2 # include equality as being better -- even if not 'better' strictly speaking
  } # end of better.score


  corrected.onePiece.innerRiskbands <- function(rd, g)
  {
    # See also simplified.onePiece.innerRiskbands --  which DOES remove superfluous segments (at the very end of CPL prior elicitation)
    #                                                 while THIS fct does NOT remove any slope/h/r term: it simply splits one-piece
    #                                                 riskbands into two equal-lengths equally-sloped segments
    #

    nSegments0 <- rd$endrb.nSegments[1]
    h.pivot.indices <- c(1,3,5,7) + nSegments0

    h.pivot <- rd$h[h.pivot.indices]
    hbar <- (h.pivot[-1] + h.pivot[-4]) / 2
    delta.h <- diff(h.pivot)

    u <- g$u[c(2,3,4)]
    one.piece <- abs(hbar-u) < 1e-5 # relaxed condition

    w <- which(one.piece)
    if (length(w) == 0)  return(list(corrected=FALSE))


    innerb.pivots <- rd$inner.b[c(1,3,5,7)]
    innerb.midpoints <- (innerb.pivots[-1] + innerb.pivots[-4]) / 2


    for (j in w)
    {
      r <- 1 + j  # r in {2,3,4}                              2

      j <- 2*r - 2
      rd$inner.b[j] <- innerb.midpoints[r-1]

      j <- c(0, 1) + 2*r + nSegments0 - 3
      straight.slope <- delta.h[r-1] / g$b[r]
      rd$s[j] <- rep(straight.slope, 2)


      j <- nSegments0 - 2 + 2*r
      rd$h[j] <- hbar[r-1]
    }


    return(list(corrected=TRUE, d=rd))
  } # end of corrected.onePiece.innerRiskbands


  corrected.onePiece.rb3 <- function(rb3, g)
  {
    h.pivot <- rb3$h[c(1,3)]
    hbar <- mean(h.pivot)
    one.piece <- abs(hbar-g$u[3]) < 1e-5 # relaxed condition

    if (!one.piece)  return(list(corrected=FALSE))

    delta.h <- diff(h.pivot)
    s.straight <- delta.h / g$b[3]

    rb3$h[2] <- hbar
    rb3$s[1:2] <- rep(s.straight, 2)

    return(list(corrected=TRUE, rb3=rb3))
  } # end of corrected.onePiece.rb3


  corrected.onePiece.rb4 <- function(rb4, G)
  {
    h.pivot <- rb4$h[c(1,3)]
    hbar <- mean(h.pivot)
    one.piece <- abs(hbar-G$n2l$u) < 1e-5 # relaxed condition

    if (!one.piece)  return(list(corrected=FALSE))

    delta.h <- diff(h.pivot)
    s.straight <- delta.h / G$n2l$b

    rb4$h[2] <- hbar
    rb4$s[1:2] <- rep(s.straight, 2)

    return(list(corrected=TRUE, rb4=rb4))
  } # end of corrected.onePiece.rb4


  endrb.x95b <- function(G, t)
  {
    if (t == G$D)                return(G$x95b)
    if (G$rightSide.closing.rb)  return(c(G$x95b[1], G$x95b[1] + t/G$cos.theta, G$x95b[2]))
    else                         return(c(G$x95b[1], G$x95b[2] - t/G$cos.theta, G$x95b[2]))
  } # end of endrb.x95b

  
  geometry <- function(riskband.prior.prob, lim, z, A, outcome.is.logNormally.distributed=TRUE, max.deepness=1/3)
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

    
    # At the border between two riskbands, the lower value allowed for the distrn height
    # is the lower of the two adjascent riskand 'lo' levels
    
    ground <- pmin(lo[-R], lo[-1])
      
    lo[c(1,5)] <- 0



    hi <- 4 * hbar - 3 * lo
      hi[c(1,5)] <- Inf # meaningless for the left-end and right-end riskbands

      alt.hi <- c(Inf, Inf, hi[-c(1,R)])
      hi <- pmin(hi, alt.hi)
      hi[5] <- Inf

      # (Note that even in the absence of floor/deepness consideration, the height
      #  of the distrn fct at the border of any riskband is still constrained
      #  by the fact that it cannot go in the negative values [at the riskband mid-section
      #  or on the other side of the riskband] to compensate.)


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

    # Riskband bounds

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
    
    
    endrb.is.local.min <- rep(NA, 2)
    uniq.u <- u[!duplicated(u)]
    endrb.is.local.min[1] <- uniq.u[1] < uniq.u[2]
    
    uniq.u <- rev(u)
    uniq.u <- uniq.u[!duplicated(uniq.u)]
    endrb.is.local.min[2] <- uniq.u[1] < uniq.u[2]
    
    j <- c(1,5)
    s.fade.out <- u[j] * max.deepness / (d + b[j]) * c(1, -1)


    gt.prec <- u[2:4] > u[1:3]
    gt.next <- u[2:4] > u[3:5]
    is.local.max <- c(NA, gt.prec & gt.next, NA)

    lt.prec <- u[2:4] < u[1:3]
    lt.next <- u[2:4] < u[3:5]
    is.local.min <- c(NA, lt.prec & lt.next, NA)


    going.up <- c(NA, lt.next, NA)
    going.dn <- c(NA, gt.next, NA)

      going.west.up <- c(NA, lt.prec, NA)
      going.west.dn <- c(NA, gt.prec, NA)


    middle.step.dn <- c(NA, gt.next & lt.prec, NA)
    middle.step.up <- c(NA, lt.next & gt.prec, NA)
      middle.step <- middle.step.dn | middle.step.up

    next.gt.next <- c(gt.next[-1], FALSE)
    next.lt.next <- c(lt.next[-1], FALSE)

    going.up.next2rb <- c(NA, lt.next & next.lt.next, NA)
    going.dn.next2rb <- c(NA, gt.next & next.gt.next, NA)

    prec.gt.prec <- c(FALSE, gt.prec[-3])
    prec.lt.prec <- c(FALSE, lt.prec[-3])

    going.dn.prec2rb <- c(NA, gt.prec & prec.gt.prec, NA) # going west
    going.up.prec2rb <- c(NA, lt.prec & prec.lt.prec, NA)



    list(z=z, A=A, H=H, L=L, Bt=Bt, d=d, m=m, cos.theta=cos.theta,
         mu.lim=mu.lim, sigma.lim=sigma.lim, R=R,
         B=B, b=b, area=area, u=u, u.low=u.low,

         going.up=going.up, going.dn=going.dn,
         middle.step.dn=middle.step.dn, middle.step.up=middle.step.up, middle.step=middle.step,
         going.up.next2rb=going.up.next2rb, going.dn.next2rb=going.dn.next2rb,
         going.dn.prec2rb=going.dn.prec2rb, going.up.prec2rb=going.up.prec2rb,
         going.west.up=going.west.up, going.west.dn=going.west.dn,

         R1b1=R1b1, R1b2=R1b2, inner.b=inner.b,
         R5b1=R5b1, R5b2=R5b2, 
         riskband.prior.prob=riskband.prior.prob, hbar=hbar, p=p,
         is.local.max=is.local.max, is.local.min=is.local.min,
         endrb.is.local.min=endrb.is.local.min, s.fade.out=s.fade.out,
         hi=hi, ground=ground, impossible=impossible,
         outcome.is.logNormally.distributed=outcome.is.logNormally.distributed,
         x95b=x95b, # useful for call to plot fct only
         s1.intercept=s1.intercept, s1.h0=s1.h0,
         forward=TRUE, max.deepness=max.deepness, my.A=my.A,
         side=0) # side = 0 identifies the 'central' nature of the returned geometrics (side = -1 | +1 for g.left | g.right)
  } # end of geometry
  
  
  geometry.closing.riskband <- function(g, rightSide.closing.rb=FALSE)
  {
    r    <- ifelse(rightSide.closing.rb, g$R,  1)
    side <- ifelse(rightSide.closing.rb,   1, -1)
    
    prob <- g$riskband.prior.prob[r]
    b <- g$b[r]
    p <- g$p[r]
    d <- g$d
    D <- d + b
    L <- g$L
    my.u <- g$u[r]

    
    b2.2 <- b^2/2
    
    b2.2d <- b2.2/d
    D.d <- D/d
    Db2.d <- D*(b^2)/d
    b3 <- b^3
    b3.6 <- b3 / 6
    b3.3d <- b3 / 3 / d

    D3mb3.3 <- (D^3 - b3) / 3
    
    d2 <- d^2
    d2.6 <- d2 / 6

    # South-East, South-South-East & East-South-South-East heights in closing riskbands
  
    h0.SSE <- 2 * p / b
    h0.SE  <- p * D / (b2.2 + b*d + d2/3)
    
    # (fadeOutConnexion meets the fading-out slope at x=b)
    j <- ifelse(rightSide.closing.rb, 2, 1)
    sfo <- g$s.fade.out[j]
    if (!rightSide.closing.rb) sfo <- -sfo
    h0.fadeOutConnexion <- sfo * (2/3*d2/b + d) + 2/b*p
    
    sfo.d <- sfo / d
    sfo.3d <- sfo.d / 3
    
    is.local.min <- g$endrb.is.local.min[j]
    

    # clothesLine parms

    a <- b2.2 + D*(D+b) / 2 - D3mb3.3/d
    cl <- list(a=p/a, b=(D+b) / 2 / a)
    

    # hill.top parms

    ht.a <- L * d / 2
    k <- L*b + ht.a
    s.coeff <- L* b2.2 + ht.a*b

    ht <- list(a=prob/s.coeff, b=k/s.coeff)


    # line.break parms

    a2 <- (D^2 - b^2) / 2
    lb.a <- d*D^2 -2*D*a2 + D3mb3.3

    s.coeff <- -b*d/2 - lb.a/d

    lb <- list(a=p/s.coeff, b=b/2/s.coeff)
    
    target.tltb <- p - sfo*d2.6


    # southEast connexion

    div <- b2.2 + b*d/3
    SE.connexion <- list(a=p/div, b=(b+d/3)/div)


    
    if (rightSide.closing.rb)
    { 
      x95b  <- g$R5b1
      x95b2 <- g$R5b2
    }
    else
    {
      x95b  <- g$R1b1
      x95b2 <- g$R1b2
    }


    # Neighbouring riskband parms (n2l = next-to-last)

    r  <- ifelse(rightSide.closing.rb, g$R-1, 2) # next-to-last riskband number
    gr <- ifelse(rightSide.closing.rb, g$R-1, 1)

    n2l <- list(b=g$b[r], p=g$p[r], u=g$u[r], ground=g$ground[gr],
                middle.step=g$middle.step[r],
                is.local.min=g$is.local.min[r],
                is.local.max=g$is.local.max[r],
                side=ifelse(rightSide.closing.rb, 1, -1),
                epsilon=g$epsilon)


      if (rightSide.closing.rb)
      {
        n2l$middle.step.dn <- g$middle.step.dn[r]
        n2l$middle.step.up <- g$middle.step.up[r]

        n2l$going.up <- g$going.up[r]
        n2l$going.dn <- g$going.dn[r]

        n2l$going.west.dn <- g$going.up[r-1]
        n2l$going.west.up <- g$going.dn[r-1]
      }
      else
      {
        n2l$middle.step.dn <- g$middle.step.up[r]
        n2l$middle.step.up <- g$middle.step.dn[r]

        # going up or down as seen when going from right to left, in this case

        n2l$going.west.dn <- g$going.dn[r]
        n2l$going.west.up <- g$going.up[r]

        n2l$going.up <- !g$endrb.is.local.min[1]
        n2l$going.dn <-  g$endrb.is.local.min[1]
      }

    
    list(R=g$R, b=b, d=d, D=D, L=L, m=g$m, p=p, my.u=my.u, prior.prob=prob,
         is.local.min=is.local.min, 
         h0.SE=h0.SE, h0.SSE=h0.SSE, h0.fadeOutConnexion=h0.fadeOutConnexion,
         b2.2 = b2.2, b3.6 = b3.6, pd = p*d, d2.6 = d2.6,
         b2.2d = b2.2d, D.d = D.d, Db2.d = Db2.d, b3.3d = b3.3d, 
         my.s.fade.out=sfo, sfo.3d = sfo.3d, sfo.d = sfo.d, target.tltb = target.tltb,
         cl=cl, ht=ht, lb=lb, SE.connexion=SE.connexion, n2l=n2l,
         rightSide.closing.rb=rightSide.closing.rb, x95b=x95b, x95b2=x95b2,
         u=g$u, s.fade.out=g$s.fade.out, u.low=g$u.low, endrb.is.local.min=g$endrb.is.local.min, 
         cos.theta=g$cos.theta, epsilon=g$epsilon, max.deepness=g$max.deepness,
         impossible=g$impossible, forward=g$forward,
         side=side)
  } # end of geometry.closing.riskband
  
  
  geometry.completeDefn <- function(g)
  {
    if (g$endrb.is.local.min[2]) 
    {
      # Correction to higher entry level in riskband R-1 when riskband R is a local min
      # (to ensure a density curve going down at entry of riskband R, at a higher point than u[R])
      
      r <- g$R - 1
      u <- g$u[g$R]
      hbar <- g$hbar[r]
      
      new.hi <- 4*hbar - 3*u
      if (g$hi[r] > new.hi) g$hi[r] <- new.hi
    } 
    
    
    if (!g$impossible & g$endrb.is.local.min[1])
    {
      d <- g$d
      b <- g$b[1]
      p <- g$p[1]
      
      k1 <- d^2/3 + d*b
      
      #           __
      # plateau  /
      
      s <- p / k1
      plateau.h <- s * d
      
      if (plateau.h > g$hi[2]) g$impossible <- TRUE
    }
    
    g
  } # end of geometry.completeDefn

  
  geometry.backward <- function(g)
  {
    g$forward <- !g$forward
    
    g$B  <- rev(g$B)
    g$b  <- rev(g$b)
      
    g$hbar    <- rev(g$hbar)
    g$hi      <- rev(g$hi)
    g$ground  <- rev(g$ground)
    g$p       <- rev(g$p)
      
    g$u                    <-  rev(g$u)
    g$u.low                <-  rev(g$u.low)
    g$riskband.prior.prob  <-  rev(g$riskband.prior.prob)
    g$endrb.is.local.min   <-  rev(g$endrb.is.local.min)
    g$s.fade.out           <- -rev(g$s.fade.out)
      
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

    g$is.local.max <- rev(g$is.local.max)
    g$is.local.min <- rev(g$is.local.min)


    g$going.up.next2rb <- rev(g$going.up.prec2rb)
    g$going.dn.next2rb <- rev(g$going.dn.prec2rb)

      g$going.up.prec2rb <- NULL
      g$going.dn.prec2rb <- NULL


    tmp.up <- g$going.up
    tmp.dn <- g$going.dn

    g$going.up <- rev(g$going.west.up)
    g$going.dn <- rev(g$going.west.dn)

      g$going.west.up <- rev(tmp.up)
      g$going.west.dn <- rev(tmp.dn)


    tmp <- g$middle.step.up
    g$middle.step.up <- rev(g$middle.step.dn)
    g$middle.step.dn <- rev(tmp)
    g$middle.step    <- rev(g$middle.step)


    g
  } # end of geometry.backward


  reg.init <- function(epsilon)
  {
    zero <- rep(0, 3) # we will use for 3 riskbands [nos. 2, 3 & 4]
    o <- list(n=zero, sum=list(x=zero, x2=zero, xy=zero, y=zero, y2=zero))

    reg <- o
    reg$forward <- TRUE
    reg$epsilon <- epsilon
    reg$backward <- o

    return(reg)
  } # end of reg.init


  riskband.plot.regions <- function(g)
  {
    A         <- g$my.A
    mu.lim    <- g$mu.lim
    sigma.lim <- g$sigma.lim
    z         <- g$z

    frame.x <- rep(mu.lim, rep(2, 2))
    frame.y <- c(sigma.lim, rev(sigma.lim))
      frame.x <- c(frame.x, frame.x[1])
      frame.y <- c(frame.y, frame.y[1])

    plot(frame.x, frame.y, col='white', xlab='mu', ylab='sigma')

    # Plot riskband delimitations

    X <- matrix(NA, nrow=2, ncol=2)
    Y <- X

    for (i in seq(along=A))
    {
      # bottom

      x.i <- A[i] - z*sigma.lim[1]
      x.i <- min(x.i, mu.lim[2])
      y.i <- ifelse(x.i < mu.lim[2], sigma.lim[1], (A[i]-x.i)/z)

      # top

      X.i <- A[i] - z*sigma.lim[2]
      X.i <- max(X.i, mu.lim[1])
      Y.i <- ifelse(X.i > mu.lim[1], sigma.lim[2], (A[i]-X.i)/z)

      x.coord <- c(x.i, X.i)
      y.coord <- c(y.i, Y.i)

      points(x.coord, y.coord, type='l', col='firebrick')

      if (i==1)
      {
        X[,1] <- x.coord
        Y[,1] <- y.coord
      }
      else if (i == length(A))
      {
        X[,2] <- x.coord
        Y[,2] <- y.coord
      }
    }


    # Paint the two end riskbands

    # [South-West corner]

    x <- X[,1]
    y <- Y[,1]

    if (x[2] > mu.lim[1])
    {
      # riskband is trapezoid-shaped
      x <- c(x, rep(mu.lim[1], 2))
      y <- c(y, rev(sigma.lim))
    }
    else
    {
      # riskband is triangular
      x <- c(x, mu.lim[1])
      y <- c(y, sigma.lim[1])
    }


    x <- c(x, x[1])
    y <- c(y, y[1])

    polygon(x, y, col='wheat', border=NA)


    # [North-East corner]

    x <- X[,2]
    y <- Y[,2]

    if (x[2] < mu.lim[2])
    {
      # riskband is trapezoid-shaped
      x <- c(x, rep(mu.lim[2], 2))
      y <- c(y, rev(sigma.lim))
    }
    else
    {
      # riskband is triangular
      x <- c(x, mu.lim[2])
      y <- c(y, sigma.lim[2])
    }


    x <- c(x, x[1])
    y <- c(y, y[1])

    polygon(x, y, col='wheat', border=NA)


    # Redraw the end riskbands delimitations

    for (i in 1:2) points(X[,i], Y[,i], type='l', col='firebrick')


    # Plot domain limits (frame the domain)

    polygon(frame.x, frame.y, border='gray')

  } # end of riskband.plot.regions


  simplified.onePiece.innerRiskbands <- function(rd, g)
  {
    # See also corrected.onePiece.innerRiskbands -- which does NOT remove any slope/h/r term: it simply splits one-piece
    #                                               riskbands into two equal-lengths equally-sloped segments,
    #                                               while THIS function DOES remove superfluous segments (at the very end of CPL prior elicitation)

    nSegments0 <- rd$endrb.nSegments[1]
    h.pivot.indices <- c(1,3,5,7) + nSegments0

    h.pivot <- rd$h[h.pivot.indices]
    hbar <- (h.pivot[-1] + h.pivot[-4]) / 2
    delta.h <- diff(h.pivot)

    u <- g$u[c(2,3,4)]
    one.piece <- abs(hbar-u) < 1e-5 # relaxed condition

    w <- rev(which(one.piece))
    if (length(w) == 0)  return(list(simplified=FALSE))

    for (j in w)
    {
      r <- 1 + j  # r in {2,3,4}

      j <- 2*r - 2
      rd$inner.b <- rd$inner.b[-j]

      j <- 2*r + nSegments0 - 3
      rd$s <- rd$s[-j]
      rd$s[j] <- delta.h[r-1] / g$b[r]

      rd$r <- rd$r[-j]

      j <- h.pivot.indices[r]-1
      rd$h <- rd$h[-j]
      rd$x95 <- rd$x95[-j]
    }


    return(list(simplified=TRUE, d=rd))
  } # end of simplified.onePiece.innerRiskbands
  

  smooth <- function(g, g.bwd, enabled, reg)
  {
    # --- Useful fcts ----------------------------------------------------------


    angle2 <- function(s)
    {
      # s: vector of length 2
      #
      # Return the angle between two slopes

      return(abs(diff(atan(s))))
    } # end of angle2


    ask.the.wizard <- function(key)
    {
      J <- ifelse(key$incoming.slope.included, 2, 1) # 1st riskband midpoint index

      # First riskband described/covered by key$h
      r <- 5 - (length(key$h) - 2) / 2
               # Note that key$h was right-trimmed in roughness() --- hence key$h is of even length


      if (key$w%%2 == J%%2)
      {
        rb <- r + (key$w - J) / 2
        j <- key$w + ifelse(key$incoming.slope.included, 0, 1)
        return(-sign(key$h[j] - key$u[rb]))
      }

      if (abs(key$w - J)%%4 == 3)  return(key$sign)

      return(-key$sign)
    } # end of ask.the.wizard


    ask.the.wizard.h03 <- function(key)
    {
      J <- 1  # first h index


      if (key$w%%2 == J)
      {
        # Left side of riskbands

        rb <- 2 + (key$w - J) / 2   # rb in {2, 3, 4, 5}

        if (rb == 2)
        {
          dir <- -key$sign
          return(dir)
        }

        if (rb == 3)  return(key$sign)
        return(0)
      }


      # Riskband midpoints

      rb <- 1 + (key$w - J + 1) / 2  # rb in {2, 3, 4}


      if (rb == 2)  return(0)

      if (rb == 3)
      {
        if (key$is.local.min3 && abs(key$h[5] - key$ground[3]) < 1e-10)  return(key$sSum.sign)
        if (key$middle.step[3] && key$null.s1)                           return(key$sign)
        return(-key$sign)
      }


      # rb == 4

      if (key$middle.step[4] && key$null.s1)  return(-key$sign)
      return(key$sign)
    } # end of ask.the.wizard.h03


    ask.the.wizard.t <- function(key, rb.t)
    {
      if (key$w == 1)  return( key$sign * rb.t$convexity)
      if (key$w == 3)  return(-key$sign * rb.t$convexity)
      if (key$w == 2)  return(sign(rb.t$soft - rb.t$x))

      return(0)
    } # end of ask.the.wizard.t

    
    augmented.body <- function(o, appendix, right.side=FALSE)
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
        if (trim.appendix) o$h <- o$h[-1] # we trim o$h [same result as trimming appendix$h]

        o$h <- c(appendix$h, o$h)
        o$s <- c(appendix$s, o$s)
      }
      
      
      if (appendix$nSegments > 0)
      {
        if (is.null(o$endrb.nSegments))  o$endrb.nSegments <- rep(0, 2)

        j <- ifelse(right.side, 2, 1)
        o$endrb.nSegments[j] <- appendix$nSegments

        if (right.side)  o$R5b <- appendix$x95b
        else             o$R1b <- appendix$x95b
      }
      
      return(o)
    } # end of augmented.body
      

    closing.riskband.flexible <- function(G, h0, s4, rb4.incoming.slope=NULL)
    {
      # G: closing riskband geometrics
      # s4: vector of length 2
      #
      # NOTE: This fct is called for riskbands that are 'endrb.is.local.min=T' ONLY
      

      if (h0 == G$h0.SE)
      {        
        s <- h0/G$D
        if (G$rightSide.closing.rb) s <- -s
        
        rb <- list(h=0, s=s, x95b=G$x95b, rightSide.closing.rb=G$rightSide.closing.rb, nSegments=1)
        
        return(rb)
      }


      convex <- h0 < G$h0.SE
      s1b <- southEast.connexion.slope(G, h0) 
      

      # ----------------------------------------------------------------------
      # We are searching for optimal (s1, s2, t)
      
      # Define range on s1 (and its initial value)
      
      if (!convex)
      {
        max.angle <- atan(NewtonRaphson.lastriskband.fade.out.s1(G, h0)$s1)
        rb5.range <- c(-pi/2, max.angle)
        
        rb5 <- my.bs.init(rb5.range, change=0.015, start.hi=TRUE, untouchable.bounds=rep(TRUE, 2), epsilon=G$epsilon)
      }
      else
      {
        untouchable.bounds <- rep(TRUE, 2)
        s.clothesLine <- clothesLine(G, h0)$s
        if (!G$rightSide.closing.rb) s.clothesLine <- -s.clothesLine
        
        rb5.range <- c(atan(s.clothesLine), pi/2)

        if (G$is.local.min)
        {
          rb5.range[2] <- 0
          untouchable.bounds[2] <- FALSE
        }

        rb5 <- my.bs.init(rb5.range, change=0.015, start.lo=TRUE, untouchable.bounds=untouchable.bounds, epsilon=G$epsilon)
      }

      incoming.slope.sign <- sign(s4[2])


      while (rb5$continue)
      {
        s1 <- tan(rb5$x)  # stands for s5.a, really
          if (abs(s1) < 1e-10) s1 <- 0
     
        # Find (s2, t) such that the integral \int_0^D h(x)l(x)dx = riskband prob
        #              [with  theta2 = atan(s2)]
        
        tmp <- NewtonRaphson.lastriskband.fixed.s1(s1, h0, G, convex=convex, s1b=s1b)     
        s2 <- tmp$s2
        t  <- tmp$t
        
        s5 <- c(s1, s2)
        h5 <- c(h0, tmp$h.intersect, 0)
        
        closing.rb <- list(h=h5, s=s5, t=t, rightSide.closing.rb=TRUE, endrb.nSegments=c(0, 2)) # Yes, rightSide.closing.rb = TRUE (always)


        rough <- roughness(G, closing.rb, incoming.slope.sign=incoming.slope.sign)
          score <- rough$score[seq(2)]


        s <- c(rb4.incoming.slope, s4, s5, 0)
        a <- angles(s)
        score <- c(score, sum(a^2))


        rb5 <- bs.update(rb5, score, closing.rb)
      } # end of while-5 
      
      
      # Closing riskband parms
      
      parms <- rb5$best.parms
      parms$x95b <- endrb.x95b(G, parms$t)
      parms$nSegments <- 2
      parms$endrb.nSegments <- NULL
      parms$rightSide.closing.rb <- G$rightSide.closing.rb
      parms$h <- parms$h[-1]

          # We do NOT reverse/flip $h and $s as the fct calling closing.riskband.flexible expects these parameters
          # to be in the order of a right-side closing riskband (and will reverse them when comes the time)

      
      return(parms)
    } # end of closing.riskband.flexible
    
         
    clothesLine <- function(G, h0)
    {
      # G: closing riskband geometrics

      if (G$my.u == h0)  return(list(h=h0, s=0, rightSide.closing.rb=G$rightSide.closing.rb, nSegments=1))

      s <- G$cl$a - h0*G$cl$b
      
      h <- h0 + s*G$D
      if (!G$rightSide.closing.rb) s <- -s
      
      list(h=h, s=s, rightSide.closing.rb=G$rightSide.closing.rb, nSegments=1)
    } # end of clothesLine


    dive.slope <- function(G, h0)
    {
      # G: must be g.left [or g.right]

      # Return the slope that dives near the cliff
      # and leading to the expected probability coverage
      # when merging with the fadeOutConnexion line

      s <- G$side * NewtonRaphson.lastriskband.fade.out.s1(G, h0)$s1
      return(s)
    } # end of dive.slope
    
    
    drop.elements <- function(x, n2drop, tail=TRUE)
    {
      if (tail) w <- seq(to=length(x), length=n2drop)
      else      w <- seq(n2drop)
      
      x[-w]
    } # end of drop.elements


    endrb.cost <- function(nSegments, s, h1, right.side=TRUE)
    {
      if (right.side)  j <- seq(to=length(s), length=nSegments + 3)
      else             j <- seq(nSegments + 3)

      s <- s[j]


      if (h1 == 0)
      {
        if (right.side)  s <- c(s, 0)
        else             s <- c(0, s)

        hidden.angle <- 0
      }
      else
      {
        hidden.angle <- ifelse(right.side, pi/2 + atan(s[length(s)]), pi/2 - atan(s[1]))
        hidden.angle <- hidden.angle^2 + (pi/2)^2
      }


      a <- angles(s)

      return(sum(a^2) + hidden.angle)
    } # end of endrb.cost


    end2rb.roughness <- function(s, h, nSegments, right.side=TRUE)
    {
      # h does not need to be complete --- as long as it contains the final h value

      len <- 3 + nSegments

      if (right.side)  j <- seq(to=length(s), length=len)
      else             j <- seq(len)

      s <- s[j]

      last.h <- ifelse(right.side, h[length(h)], h[1])


      if (last.h == 0)
      {
        if (right.side)  s <- c(s, 0)
        else             s <- c(0, s)

        hidden.angle <- 0
      }
      else
      {
        hidden.angle <- ifelse(right.side, pi/2 + atan(s[length(s)]), pi/2 - atan(s[1]))
        hidden.angle <- hidden.angle^2 + (pi/2)^2
      }


      a <- angles(s)

      return(sum(a^2) + hidden.angle)
    } # end of end2rb.roughness


    fadeOut.coords <- function(G, h0, s)
    {
      abs.s   <- abs(s)
      abs.sfo <- abs(G$my.s.fade.out)
      t <- (h0 - abs.sfo * G$D) / (abs.s - abs.sfo)

      h <- c(0, abs.sfo * (G$D - t))
      s <- c(abs.sfo, abs.s)

      if (G$side > 0)
      {
        h <-  rev(h)
        s <- -rev(s)
        x95b <- c(G$x95b[1] + c(0, t)/G$cos.theta, G$x95b[2])
      }
      else x95b <- c(G$x95b[1], G$x95b[2] - c(t, 0)/G$cos.theta)


      list(h=h, s=s, x95b=x95b, rightSide.closing.rb=G$rightSide.closing.rb, nSegments=2)
    } # end of fadeOut.coords


    first.draft <- function(g, g.right, g.bwd, g.bwd.right, reg, print=FALSE)
    {
      h03_range <- h03.range(g)
      h03 <- h03.inits(g)
        incoming.slope.sign3 <- h03$incoming.slope.sign


      rb3.h0 <- my.bs.init(h03_range, h03$init, h03$change, suggested.2nd.step=h03$suggested.2nd.step, epsilon=g$epsilon)

      if (print) catn('EG =', g$use.educated.guess, g.bwd$use.educated.guess)


      print3 <- FALSE
      print4 <- FALSE


      while (rb3.h0$continue)
      {
        h03 <- rb3.h0$x

        tmp <- riskband3.onward(h03, g, g.right, reg, incoming.slope.sign3, print1=print3)
          reg <- tmp$reg
          rb3 <- tmp$rb3
          rb4.incoming.slope <- -rb3$s[1] # it is actually rb2 incoming slope from its right side


        tmp <- riskband4.onward(h03, g.bwd, g.bwd.right, reg, rb4.incoming.slope, for.rb12=TRUE, print4=print4)
          reg <- tmp$reg
          rb12 <- tmp$rb4                       # At this point, parms in rb12 are in reverse order

          rb12 <- flipped.rb12(rb12)            # Now rb12 parms are in the right/expected order


        rough <- roughness(g, rb3, rb12)
          score <- rough$score

          # not combining rb3 & rb12 before the calculation of roughness
          # allows the use of local roughness for rb12 in the choice of the best combination


        where2go.next <-  ifelse(g$use.educated.guess,
                                 ask.the.wizard.h03(rough$key),
                                 0)

        parms <- list(rb12=rb12, rb3=rb3)


        if (print)
        {
          catn('---------------')
          catn('step#', rb3.h0$n.steps)
          catn('h03', h03)
          catn('score', score)
          catn('key', unlist(rough$key))
          catn('where2go.next', where2go.next)
        }



        rb3.h0 <- bs.update(rb3.h0, score, parms, where2go.next)

        if (print)
        {
          catn('continue', rb3.h0$continue)
          if (rb3.h0$last.score.was.the.best)  catn('<-- best h03 so far')
        }
      }


      best <- rb12.rb3.combo(rb3.h0$best.parms)
        best$score <- bs.best_score(rb3.h0)
        best$inner.b <- g$inner.b

      return(list(d=best, reg=reg))
    } # end of first.draft


    flat.start <- function(g, r, h0)
    {
      # Returns the height to reach after a flat start to get the right coverage
      #
      #        h(final)
      #       /
      # h0 __/

      return(4*g$p[r]/g$b[r] - 3*h0)
    } # end of flat.start


    flipped.rb12 <- function(rb12)
    {
      rb12$h <- rev(rb12$h[-1])
      rb12$s <- -rev(rb12$s)

      rb12$x95b <- -rev(rb12$x95b)
      rb12$rightSide.closing.rb <- FALSE

      return(rb12)
    } # end of flipped.rb12


    h03.inits <- function(g)
    {
      if      (g$is.local.min[2])  inits <- list(init=g$ground[2], suggested.2nd.step=g$u[2])
      else if (g$is.local.max[2])  inits <- list(init=g$u[2],      suggested.2nd.step=g$u[3])
      else if (g$middle.step[2])   inits <- list(init=g$u[2],      suggested.2nd.step=(g$u[2]+g$u[3])/2)
      else                         inits <- list(init=g$u[3],      suggested.2nd.step=(g$u[2]+g$u[3])/2)

      inits$change <- abs(inits$init - inits$suggested.2nd.step)

      incoming.slope.sign <- sign(g$u[3] - g$u[2])
        if (incoming.slope.sign == 0)  incoming.slope.sign <- ifelse(g$endrb.is.local.min[1], -1, 1)
        inits$incoming.slope.sign <- incoming.slope.sign

      return(inits)
    } # end of h03.inits


    h03.range <- function(g)
    {
      range <- c(g$ground[2], g$hi[3])

      if      (g$middle.step.dn[3])  range[1] <- g$u[3]
      else if (g$middle.step.up[3])  range[2] <- g$u[3]


      if (g$is.local.min[3])
      {
        tmp <- 4*g$p[3]/g$b[3] - 3*g$ground[3]
        if (tmp < range[2])  range[2] <- tmp
      }


      if (g$middle.step.dn[2])
      {
        tmp <- g$u[2]
        if (tmp < range[2])  range[2] <- tmp
      }
      else if (g$middle.step.up[2])
      {
        tmp <- g$u[2]
        if (tmp > range[1])  range[1] <- tmp

        tmp <- flat.start(g, 2, g$u[1])
        if (tmp < range[2])  range[2] <- tmp
      }


      if (g$going.dn.next2rb[3])
      {
        tmp <- flat.start(g, 3, g$u[4])
        if (tmp < range[2])  range[2] <- tmp
      }


      if (g$is.local.min[2])
      {
        tmp <- 4*g$p[2]/g$b[2] - 3*g$ground[1]
        if (tmp < range[2])  range[2] <- tmp
      }
      else if (g$middle.step.dn[3])
      {
        tmp <- plateau(3, g, g$hi[4])$h
        if (tmp > range[1])  range[1] <- tmp
      }


      return(range)
    } # end of h03.range
    
    
    line.break <- function(G, h0)
    {
      # Draw a broken line (into two linear pieces) from the next-to-last riskband height to 0 at the upper/outer value for x95;
      # the line break occurs at 'b' or before if necessary.
      #
      # G: closing riskband geometrics
      
      
      if (h0 < G$h0.fadeOutConnexion)
      {
        # Fix occurs above the intersection of s-fade-out slope and vertical line at x=b
        #  => we will break the line at x = b
        
        s2 <- G$lb$a - h0 * G$lb$b
        h.intersect <- -s2*G$d
        s1 <- (h.intersect - h0) / G$b
        
        s <- c(s1, s2)
        t <- G$b
      }
      else
      {
        # Fix occurs with a line of slope s1 (to be calculated) intersecting
        # with the line of slope G$s.fade.out going through south-east corner (D,0)

        tmp <- NewtonRaphson.lastriskband.fade.out.s1(G, h0)
        
        s <- c(tmp$s1, G$my.s.fade.out)
        t <- tmp$t
        h.intersect <- tmp$h.intersect 
      }
      
      
      # Wrap up        
  
      h <- c(h.intersect, 0)
      x95b <- endrb.x95b(G, t)
      
      if (!G$rightSide.closing.rb)
      {
        h <-  rev(h)
        s <- -rev(s)
      } 

      
      list(h=h, s=s, x95b=x95b, rightSide.closing.rb=G$rightSide.closing.rb, nSegments=2)
    } # end of line.break
                  

    my.bs.init <- function(range, init=start$init, change=start$change, untouchable.bounds=rep(FALSE, 2),
                           start.lo = FALSE, start.hi = FALSE, start = list(init = Inf, change = Inf), suggested.2nd.step=NULL,
                           epsilon=g$epsilon)
    {
      bs.init(range, epsilon, init=init, change=change, start.lo=start.lo, start.hi=start.hi,
              untouchable.bounds=untouchable.bounds, force.1pass=TRUE, better.score=better.score,
              suggested.2nd.step=suggested.2nd.step)
    } # end of my.bs.init
              

    NewtonRaphson.lastriskband.fade.out.s1 <- function(G, h0)
    { 
      # G: closing riskband geometrics
      
      # s2 [fixed] = s.fade.out
      
      D <- G$D
      s2 <- G$my.s.fade.out
        
      t.lt.b <- h0 > G$h0.fadeOutConnexion
      continue <- TRUE
                
      if (t.lt.b)
      {      
        # --- Find t < b -------------------------------------------------------
             
        t <- 0 
        target <- G$target.tltb
        s2.2 <- s2/2
        f.prime.a <- 3/2 * (h0 + D*s2)

        
        while (continue)
        {
          s1.t <- (t-D)*s2 - h0 
          
          f <- h0*t + s1.t * t/2 - s2.2 * (D-t)**2
          f.prime <- f.prime.a + s1.t - s2*t
          
          change <- (target - f) / f.prime
          t <- t + change
          continue <- abs(change) >= G$epsilon
        }
        

        h.intersect <- h0 + s1.t
        s1 <- ifelse(t==0, s2, s1.t/t)
        
        return(list(s1=s1, t=t, h.intersect=h.intersect))            
      }
      else
      {                  
        # --- Find t > b -------------------------------------------------------
        
        t <- G$b
        target <- G$p + h0*G$b2.2d
        a <- h0 * G$D.d
        s1.prime.numerator <- h0 + D*s2
        k <- 1/3/G$d 
        h0D.d <- h0 * G$D.d
        h0.d <- h0 / G$d
        
 
        while (continue)
        {
          s1 <- ((t-D)*s2 - h0) / t
          t2 <- t**2
          t2.2 <- t**2 / 2
          t3 <- t**3
          s1.prime <- s1.prime.numerator / t2
                 
          f <- s1*G$b2.2 + a*t + s1*G$D.d*(t2.2 - G$b2.2) - h0.d*t2.2 - s1*(k*t3 - G$b3.3d) - G$sfo.3d * (D-t)**3
          f.prime <- s1.prime * G$b2.2 +  h0D.d + s1.prime*G$D.d*(t2.2 - G$b2.2) + (s1*D-h0)*t/G$d - s1.prime*(k*t3 - G$b3.3d) - s1/G$d*t2 + G$sfo.d * (D-t)**2
           
          change <- (target - f) / f.prime
          t <- t + change
          continue <- abs(change) >= G$epsilon
        }
        
        
        h.intersect <- h0 + s1*t
        
        return(list(s1=s1, t=t, h.intersect=h.intersect))
      }
    } # end of NewtonRaphson.lastriskband.fade.out.s1
    
    
    NewtonRaphson.lastriskband.fixed.s1 <- function(s1, h0, G, convex=FALSE, s1b=0,
                                                    t.lt.b = h0 > G$h0.SSE | ifelse(convex, s1 > s1b, s1 < s1b))
    { 
      # G: closing riskband geometrics
      
      D <- G$D
        
        
      if (t.lt.b)
      {      
        # --- Find t < b -----------------------------------------------------
             
        t <- 0
        continue <- TRUE
           
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
          continue <- abs(change) >= G$epsilon
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
        continue <- TRUE
         
        while (continue)
        {
          u <- (h0 + s1*t) / 3
          v <- (D - t)^2
           
          f <- a1*t + a2*t^2 + a3*t^3 + u*v
          f.prime <- a1 + 2*t*a2 + 3*t^2*a3 + s1*v/3 - 2*u*(D-t)
           
          change <- (target - f) / f.prime
          t <- t + change
          continue <- abs(change) >= G$epsilon
        }
        
        h.intersect <- h0 + s1*t
        s2 <- h.intersect / (t-G$D) 
        
        return(list(s2=s2, t=t, h.intersect=h.intersect))
      }
    } # end of NewtonRaphson.lastriskband.fixed.s1 -----------------------------


    n.inflexion.points <- function(s)
    {
      s.diff <- diff(s)
      w0 <- which(s.diff==0)

      if (length(w0) > 0)
      {
        s.diff <- s.diff[-w0]
        s <- s[-w0]
      }

      s.diff.sign <- sign(s.diff)
      j <- 2 + which(s.diff.sign[-1] != s.diff.sign[-length(s.diff.sign)])
      if (length(j) == 0)  return(0)


      inflexion <- s[j]*s[j-1] >= 0 & s[j]*s[j-2] >= 0 & s[j-1]*s[j-2] >= 0

      n.inflex <- sum(inflexion)
      return(n.inflex)
    } # end of n.inflexion.points
    
    
    n.sign.changes <- function(x)
    {
      s <- sign(x)
      s <- s[s!=0]
            
      n <- sum(s[-1] != s[-length(s)])
      return(n)
    } # end of n.sign.changes
      
    
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
    
    
    parallelepiped.connexion <- function(r, g, h0, h1, left2right=TRUE)
    {
      # h0, h1: height @ left/right ends of riskband r [in that order]
      # if left2right = T, return s0
      #    otherwise       return -s1
      
      b <- g$b[r]
      s0 <- (4*g$hbar[r] - 3*h0 - h1) / b
      
      if (left2right)              return(s0)
      else
      {
        s1 <- 2/b * (h1-h0) + s0;  return(-s1)
      }  
    } # end of parallelepiped.connexion


    parallelepiped.fixedh.minTheta2Sum <- function(t.parms, t.start, s0l, s0r, NR.epsilon=1e-6)
    {
      continue <- TRUE
      range <- t.parms$range

      t <- t.start
        tmp <- parallelepiped.fixedh.theta2Sum(t, t.parms, s0l, s0r)
        if (is.nan(tmp$fp))  t <- mean(t.parms$range)


      while (continue)
      {
        tmp <- parallelepiped.fixedh.theta2Sum(t, t.parms, s0l, s0r)
        change <- - tmp$fp / tmp$fs

        j <- ifelse(change > 0, 1, 2)
        range[j] <- t

        t <- t + change

        if (t < range[1] || t > range[2])  t <- mean(range)

        converged <- abs(change) < NR.epsilon
        continue <- !converged

        if (continue)
        {
          squeezed <- diff(range) < NR.epsilon
          continue <- !squeezed
        }
      }

      tmp <- parallelepiped.fixedh.theta2Sum(t, t.parms, s0l, s0r)

      return(list(t=t, s=tmp$s, converged=converged))
    } # end of parallelepiped.fixedh.minTheta2Sum


    parallelepiped.fixedh.theta2Sum <- function(t, t.parms, s0l, s0r)
    {
      u <- function(theta)
      {
        f  <- theta$f^2
        fp <- 2 * theta$f*theta$fp
        fs <- 2 * (theta$fp^2 + theta$f*theta$fs)

        return(list(f=f, fp=fp, fs=fs))
      } # end of u


      k1 <- t.parms$k1
      k2 <- t.parms$k2
      k3 <- t.parms$k3
      alpha <- t.parms$alpha
      b <- t.parms$b


      s1 <- alpha + k1 / t
      s2 <- (k2 - alpha*t) / (b - t)

      s1p <- -k1 / t^2          # First derivatives
      s2p <-  k3 / (b - t)^2

      s1s <- 2 * k1 / t^3       # Second derivatives
      s2s <- 2 * k3 / (b-t)^3


      phi1 <- slope2angle(s1, s1p, s1s)
      phi2 <- slope2angle(s2, s2p, s2s)

      phi.l <- slope2angle(s0l)  # Constant (external) angles
      phi.r <- slope2angle(s0r)


      # Angles between consecutive slopes (s0l, s1, s2, s0r)

      theta1 <- slopes2angle(phi.l, phi1)
      theta2 <- slopes2angle(phi1, phi2)
      theta3 <- slopes2angle(phi2, phi.r)


      # We want to minimize f = theta1^2 + theta2^2 + theta3^2
      # that is, find t such that f' = 0.
      #
      # We will use Newton-Raphson
      # and hence need to compute the first and second derivatives of f

      u1 <- u(theta1)
      u2 <- u(theta2)
      u3 <- u(theta3)

      f  <- u1$f  + u2$f  + u3$f
      fp <- u1$fp + u2$fp + u3$fp  # First  derivative
      fs <- u1$fs + u2$fs + u3$fs  # Second derivative

      return(list(f=f, fp=fp, fs=fs, s=c(s1, s2)))
    } # end of parallelepiped.fixedh.theta2Sum
    
    
    parallelepiped.range <- function(r, g, h0)
    {
      hbar <- g$hbar[r]
      b <- g$b[r]
      
      s.max <- (4*hbar - 3*h0 - g$ground[r]) / b  # stalagmite
       
      reflexion <- parallelepiped.reflexion(r, g, h0)
      if (reflexion > g$hi[r+1]) s.min <- parallelepiped.connexion(r, g, h0, g$hi[r+1])
      else s.min <- 2 / b * (g$ground[r] - h0)  # stalactite
  
      s <- c(s.min, s.max)
      return(atan(s))
    } # end of parallelepiped.range
    
    
    parallelepiped.reflexion <- function(r, g, h)
    {
      4*g$hbar[r] - 2*g$ground[r] - h
    } # end of parallelepiped.reflexion


    plateau <- function(r, g, h0)
    {
      p <- g$p[r]
      b <- g$b[r]

      s <- 8/3 * (p/b^2 - h0/b)

      theta <- atan(s)
      h <- h0 + s*b/2

      return(list(theta=theta, h=h))
    } # end of plateau


    quick.close <- function(G, h0)
    {
      # This fct does return the parms of the closing region AS IF it were the right-hand side closing region
      # on purpose.

      # G: g.bwd.right

      if (G$is.local.min)
      {
        if (h0 <= G$h0.SE)
        {
          tmp <- clothesLine(G, h0)
          tmp$x95b <- G$x95b
          return(tmp)
        }


        if (h0 <= G$h0.fadeOutConnexion)  return(line.break(G, h0))

        s <- dive.slope(G, h0)
        return(fadeOut.coords(G, h0, s))
      }


      # is local max, then.

      if (h0 <= G$h0.SE)
      {
        tmp <- clothesLine(G, h0)
        tmp$x95b <- G$x95b
        return(tmp)
      }


      if (h0 <= G$h0.SSE)
      {
        s <- southEast.connexion.slope(G, h0)
        tmp <- southEast.connexion(G, h0, s)
        return(tmp)
      }


      if (h0 <= G$h0.fadeOutConnexion)  return(line.break(G, h0))


      s <- dive.slope(G, h0)
      return(fadeOut.coords(G, h0, s))
    } # end of quick.close


    rb12.rb3.combo <- function(o)
    {
      h <- c(o$rb12$h, o$rb3$h)
      s <- c(o$rb12$s, o$rb3$s)

      R1b <- o$rb12$x95b
      R5b <- o$rb3$x95b

      endrb.nSegments <- c(o$rb12$nSegments, o$rb3$endrb.nSegments[2])

      return(list(h=h, s=s, R1b=R1b, R5b=R5b, endrb.nSegments=endrb.nSegments))
    } # end of rb12.rb3.combo


    rb4x.range <- function(g, h04)
    {
      if (g$endrb.is.local.min[2])  rb4x_range <- theta4.limits(g, h04)
      else                          rb4x_range <- parallelepiped.range(r=4, g, h04)

      rb4x_range <- revisited.range(rb4x_range, 4, g, h04)

      if (!g$endrb.is.local.min[2])
      {
        h1.sup <- parallelepiped.reflexion(4, g, h04)
        h1.sup <- min(h1.sup, g$right$h0.SE)
        theta  <- atan(parallelepiped.connexion(4, g, h04, h1.sup))
        if (theta < rb4x_range[1])  rb4x_range[1] <- theta
      }

      return(rb4x_range)
    } # end of rb4x.range


    rb.t.inits <- function(my.g, h0, h1, r=NULL, incoming.slope=NULL)
    {
      # my.g: either g or g.left$n2l or g.right$n2l

      t.parms <- t.info(my.g, h0, h1, r)
      t.eq <- ifelse(diff(t.parms$range) == 0, t.parms$soft, t.equal_angles(t.parms, incoming.slope)$t)
        if (is.na(t.eq))  t.eq <- t.parms$b/2

      rb.t <- my.bs.init(t.parms$range, t.eq, diff(t.parms$range)/50,  untouchable.bounds=t.parms$untouchable.bounds,
                         suggested.2nd.step=t.parms$soft, epsilon=my.g$epsilon)

        rb.t$soft      <- t.parms$soft
        rb.t$delta.h   <- t.parms$delta.h
        rb.t$convexity <- t.parms$convexity

      return(rb.t)
    } # end of rb.t.inits


    reg.predict <- function(reg, r, h0)
    {
      # r: riskband no --- one of {2, 3, 4}
      #
      # Do a linear extrapolation to get an individual prediction at x = h0
      
      j <- r - 1
      n <- reg$n[j]
              
                        
      if (n <= 2)  return(list(pred=Inf, margin=Inf))

      
      # n > 2
      
      x.sum  <- reg$sum$x[j]
      x2.sum <- reg$sum$x2[j]
      
      x.bar <- x.sum / n
      SS <- x2.sum - n * x.bar^2  # sum of squares sum((xi-xbar)^2)  --- also b-hat denominator


      if (SS <= 0)  return(list(pred=Inf, margin=Inf))
        

      # b-hat denominator [SS] > 0  
      
      xy.sum <- reg$sum$xy[j]
      y.sum  <- reg$sum$y[j]
      y2.sum <- reg$sum$y2[j]
      y.bar  <- y.sum / n
      
      b <- (xy.sum - n * x.bar * y.bar) / SS
      a <- y.bar - b * x.bar
      pred <- a + b * h0
          
      sigma2.hat <- (y2.sum - 2*a*y.sum - 2*b*xy.sum + n*a*a + 2*a*b*x.sum + x2.sum*b*b) / (n - 2)
      sigma.hat  <- ifelse(sigma2.hat > 0, sqrt(sigma2.hat), 0)  # protection against numerical imprecision
    
    
      SE.individual <- sigma.hat * sqrt(1 + 1/n + (h0 - x.bar)^2 / SS)
          
      margin <- 3 * SE.individual # a conservative value
      margin <- max(margin, 100*reg$epsilon)
    
      return(list(pred=pred, margin=margin))
    } # end of reg.predict
    
    
    reg.update <- function(reg, r, rb)
    {
      # r: riskband no --- one of {2, 3, 4}, for each inner riskband, in order
      # rb: riskband object


      # Prepare for the calculation of regression parms between x & y, where
      #   x: h0
      #   y: entry_angle

      j <- r - 1

      h0 <- rb$h[1]
      s_entry <- rb$s[1]
      
      entry_angle <- atan(s_entry)
      
      reg$n[j]      <- reg$n[j]       + 1
      reg$sum$x[j]  <- reg$sum$x[j]   + h0
      reg$sum$x2[j] <- reg$sum$x2[j]  + h0^2
      reg$sum$xy[j] <- reg$sum$xy[j]  + h0 * entry_angle
      reg$sum$y[j]  <- reg$sum$y[j]   + entry_angle
      reg$sum$y2[j] <- reg$sum$y2[j]  + entry_angle^2


        if (reg$forward)
        {
          #  We register what's happening at the right-hand side of the riskband
          #  in order to use that information when we calculate the best CPL distribution
          #  going backward (where initial values based on regression predictions will take advantage
          #                  of the results obtained in the forward calculations).

          j <- 4 - j

          h0 <- rb$h[3]
          entry_angle <- -atan(rb$s[2])

          reg$backward$n[j]      <- reg$backward$n[j]       + 1
          reg$backward$sum$x[j]  <- reg$backward$sum$x[j]   + h0
          reg$backward$sum$x2[j] <- reg$backward$sum$x2[j]  + h0^2
          reg$backward$sum$xy[j] <- reg$backward$sum$xy[j]  + h0 * entry_angle
          reg$backward$sum$y[j]  <- reg$backward$sum$y[j]   + entry_angle
          reg$backward$sum$y2[j] <- reg$backward$sum$y2[j]  + entry_angle^2
        }
  
      return(reg)
    } # end of reg.update


    revisit.closing2riskbands <- function(G, g, rd, print5=FALSE)
    {
      # G: closing riskband geometry
      # rd: revisited density
      #
      # NOTE: This fct is called for the last two riskbands ONLY when the last one is 'endrb.is.local.min=T' 
      #  -> Nope! This condition is not necessary anymore.


      complete.with.rb45 <- function(rd, g, rightSide.closing2rb, rb4=parms$rb4, rb5=parms$rb5, parms)
      {
        if (!rightSide.closing2rb)
        {
          rb4$h <-  rev(rb4$h)
          rb4$s <- -rev(rb4$s)

          rb5$h    <-  rev(rb5$h)
          rb5$s    <- -rev(rb5$s)
          rb5$x95b <- -rev(rb5$x95b)
        }


        rd <- augmented.body(rd, rb4, rightSide.closing2rb)
        rd <- augmented.body(rd, rb5, rightSide.closing2rb)

        rd$score <- roughness(g, rd)$score


        if (rightSide.closing2rb)  rd$inner.b[6] <- rd$inner.b[5] + rb4$t / g$cos.theta
        else                       rd$inner.b[2] <- rd$inner.b[3] - rb4$t / g$cos.theta

        return(rd)
      } # end of complete.with.rb45


      h05.parms <- function(G, h04, s0)
      {
        # G: either g.left or g.right

        if (G$n2l$is.local.max)
        {
          range <- c(G$my.u, Inf)
          suggested.2nd.step <- G$n2l$u
        }
        else if (G$n2l$is.local.min)
        {
          tmp <- 2*G$n2l$p/G$n2l$b
          h1.1piece <- tmp - h04
          range <- c(min(max(h1.1piece, G$n2l$ground), G$n2l$u), tmp - G$n2l$ground)
            untouchable.bounds <- c(FALSE, TRUE)

          suggested.2nd.step <- G$n2l$u
        }
        else if (G$n2l$middle.step && abs(h04 - G$n2l$u) < 1e-10)
        {
          range <- rep(G$n2l$u, 2)
          suggested.2nd.step <- NULL
        }
        else if (G$n2l$going.up)
        {
          range <- c(G$n2l$u, Inf)
          suggested.2nd.step <- G$my.u
        }
        else
        {
          range <- c(G$my.u, G$n2l$u)
          suggested.2nd.step <- G$my.u
        }


        if (!G$n2l$is.local.min)
        {
          untouchable.bounds <- is.infinite(range)
            if      (G$n2l$going.dn)  untouchable.bounds[2] <- h04 != G$n2l$u
            else if (G$n2l$going.up)  untouchable.bounds[1] <- h04 != G$n2l$u
        }



        if (G$h05.1piece >= range[1] && G$h05.1piece <= range[2])
        {
          suggested.2nd.step <- G$h05.1piece
        }
        else if (G$middle.step && h04 != G$n2l$u)
        {
          tmp <- orig.slope.continuation(h04, G$n2l$b, G$n2l$p, s0)
          if (!is.na(tmp$t))  suggested.2nd.step <- tmp$h
        }
        else if (G$n2l$ground == range[1])
        {
          suggested.2nd.step <- G$n2l$ground
        }


        return(list(range=range, suggested.2nd.step=suggested.2nd.step, untouchable.bounds=untouchable.bounds))
      } # end of h05.parms


      orig.slope.continuation <- function(h04, b, p, s0)
      {
        A <- -s0/2
        B <- b*s0
        C <- h04*b - p
        delta <- B^2 - 4*A*C

        if (delta < 0 || A == 0)  return(list(t=NA))

        t.solns <- (-B + c(-1,1) * sqrt(delta)) / (2*A)
        t <- t.solns[t.solns >= 0 & t.solns <= b]

        if (length(t) == 0)  return(list(t=NA))

        h <- h04 + s0 * t
        return(list(h=h, t=t))
      } # end of orig.slope.continuation


      rb45 <- function(G, h04, h05, v2, rb4.incoming.slope, use.educated.guess, current.ssc, print.t)
      {
        b <- G$n2l$b


        if (h05 == G$h05.1piece)
        {
          s4 <- (h05 - h04) / b
          rb4 <- list(s=rep(s4, 2), h=c(h04, h05, h05), t=b, nSegments=-1)
          rb5 <- closing.riskband.flexible(G, h05, rb4$s, rb4.incoming.slope)

          parms <- list(rb4=rb4, rb5=rb5)
        }
        else if (abs(h05 - G$n2l$u) > 1e-10 || abs(h04 - G$n2l$u) > 1e-10)
        {
          h4.bar <- (h04 + h05) / 2
          one.piece.rb4 <- abs(h4.bar - G$n2l$u) < 1e-5


          if (one.piece.rb4)
          {
            s.straight <- (h05 - h04) / b

            s4  <- rep(s.straight, 2)
            rb4 <- list(s=s4, h=c(h04, h4.bar, h05), t=b/2, nSegments=-1)
            rb5 <- closing.riskband.flexible(G, h05, s4, rb4.incoming.slope)

            parms <- list(rb4=rb4, rb5=rb5)
          }
          else
          {
            rb4.t <- rb.t.inits(G$n2l, h04, h05, incoming.slope=rb4.incoming.slope)

            v <- matrix(c(rb4.t$delta.h, v2), ncol=1) # 2 x 1 matrix

            if (print.t) catn('starting @ soft', rb4.t$soft)


            while (rb4.t$continue)
            {
              t <- rb4.t$x
              if (print.t) catn('t', t)

              M4 <- matrix(c(t, b-t, t*b-t**2/2, (b-t)**2/2), ncol=2, byrow=TRUE)
              s4 <- as.vector(solve(M4) %*% v)

              # area.p(4, g, h04, s4[1], s4[2], t) -- oui, ca fonctionne
              # ou bien [lorsque !rightSide.closing2rb]: area.p(2, g, h04, s4[1], s4[2], t)

              ht <- h04 + t*s4[1]
              rb4 <- list(s=s4, h=c(h04, ht, h05), t=t, nSegments=-1)
              rb5 <- closing.riskband.flexible(G, h05, s4, rb4.incoming.slope)


              rough <- roughness(G, rb4, rb5, incoming.slope=rb4.incoming.slope, ssc.threshold=current.ssc)
                score <- rough$score
                parms <- list(rb4=rb4, rb5=rb5)
                where2go.next <- 0

                if (use.educated.guess)  where2go.next <- ask.the.wizard.t(rough$key, rb4.t)

                if (print.t) catn('[t] where2go.next', where2go.next)

              rb4.t <- bs.update(rb4.t, score, parms, where2go.next)
            } # end of loop on rb4.t


            parms <- rb4.t$best.parms
          }
        }
        else
        {
          rb4 <- list(s=rep(0, 2), h=c(rep(G$n2l$u, 2), h05), t=0, nSegments=-1)
          rb5 <- closing.riskband.flexible(G, h05, rb4$s, rb4.incoming.slope)

          parms <- list(rb4=rb4, rb5=rb5)
        }


        return(parms)
      } # end of rb45


      print.t <- FALSE

      if (G$n2l$middle.step.dn)  current.ssc <- Sharpest.slope.change(rd)
      else                       current.ssc <- -Inf

      rightSide.closing2rb <- G$rightSide.closing.rb
      
      j <- ifelse(rightSide.closing2rb, 2, 1)
      nSegments.closing.rb <- rd$endrb.nSegments[j]

      rd0.end2rb.roughness <- end2rb.roughness(rd$s, rd$h, nSegments.closing.rb, rightSide.closing2rb)
      
      
      # entry value in last riskband obtained in original CPL curve
      h05.init <- ifelse(rightSide.closing2rb, rd$h[length(rd$h)-nSegments.closing.rb],
                                               rd$h[1+nSegments.closing.rb])

      orig.entry.s <- ifelse(rightSide.closing2rb, rd$s[rd$endrb.nSegments[1]+5],
                                                  -rd$s[rd$endrb.nSegments[1]+2])
      
      
      # Drop h & s from last 2 riskbands (do NOT drop any elements from rd$inner.b)
      
      rd$h <- drop.elements(rd$h, nSegments.closing.rb+2, tail=rightSide.closing2rb)
      rd$s <- drop.elements(rd$s, nSegments.closing.rb+2, tail=rightSide.closing2rb)


      # h at entry in the riskband *preceding* the closing riskband
      h04 <- ifelse(rightSide.closing2rb, rd$h[length(rd$h)], rd$h[1])

      
      # Modify G
      #           [we do that as roughness() expects a curve that starts anywhere BUT ends
      #            at the rightmost endpoint of x95b]
      
      if (!rightSide.closing2rb)
      {
        G$u                     <-  rev(G$u)
        G$endrb.is.local.min    <-  rev(G$endrb.is.local.min)
        G$s.fade.out            <- -rev(G$s.fade.out)
        G$u.low                 <-  rev(G$u.low)
        G$rightSide.closing.rb  <-  TRUE
        G$x95b                  <- -rev(G$x95b)
        G$x95b2                 <- -rev(G$x95b2)
      }

      G$middle.step  <- G$n2l$middle.step
      G$h05.1piece   <- 2*G$n2l$u - h04 # straight line from h04->h05
  
                                
      # _____________________________________________________________________________________
    
      b <- G$n2l$b
      p <- G$n2l$p
      v2 <- p - h04 * b   

      rb4.incoming.slope <- ifelse(rightSide.closing2rb, rd$s[length(rd$s)], -rd$s[1])

      h05 <- h05.parms(G, h04, orig.entry.s)
      if (print5) catn('h05-range', h05$range)

      rb5.h0 <- my.bs.init(h05$range, h05.init, untouchable.bounds=h05$untouchable.bounds,
                           suggested.2nd.step=h05$suggested.2nd.step, epsilon=G$epsilon)
      

      while (rb5.h0$continue)
      {
        h05 <- rb5.h0$x
         
        parms <- rb45(G, h04, h05, v2, rb4.incoming.slope, g$use.educated.guess, current.ssc, print.t)

        rough <- roughness(g, parms$rb4, parms$rb5, incoming.slope=rb4.incoming.slope, ssc.threshold=current.ssc)
          score <- rough$score


        if (print5)
        {
          catn('-----------')
          catn('step # ', rb5.h0$n.steps)
          catn('h05', h05)
          catn('score', score)
          catn('key$w', parms$key$w)
          catn('t4', parms$rb4$t)
          catn('h4', parms$rb4$h)
        }



        rb5.h0 <- bs.update(rb5.h0, score, parms)


        if (print5 && rb5.h0$last.score.was.the.best)
        {
          catn('\\_ was the best so far')
          catn('----------');
        }
      } # end of loop on rb5.h0

      
      rb4 <- rb5.h0$best.parms$rb4    
      rb5 <- rb5.h0$best.parms$rb5

      
      # Wrap-up
       
      s <- c(rb4.incoming.slope, rb4$s, rb5$s)
      rd.end2rb.roughness <- end2rb.roughness(s, rb5$h, rb5$nSegments)

      if (rd0.end2rb.roughness < rd.end2rb.roughness)  return(list(better=FALSE))


      # The new distrn is accepted -- complete it properly

      d <- complete.with.rb45(rd, g, rightSide.closing2rb, rb4=rb4, rb5=rb5)


      tmp <- corrected.onePiece.innerRiskbands(d, g)
      if (tmp$corrected)  d <- tmp$d

      return(list(better=TRUE, d=d))
    } # end of revisit.closing2riskbands


    revisit.closing2riskbands.fixedh <- function(G, g, rd)
    {
      # G: either g.left or g.right


      revisited.innerb <- function(t, inner.b, innerb.index, cos.theta, right.side)
      {
        if (right.side)
        {
          ref.index <- innerb.index - 1
        }
        else
        {
          ref.index <- innerb.index + 1
          t <- -t
        }

        return(inner.b[ref.index] + t/cos.theta)
      } # end of revisited.innerb


      # ------------------------------------------------------------------------
      # Start of revisit.closing2riskbands.fixedh

      right.side <- G$side > 0
      outside.dir <- ifelse(right.side, 1, -1)

      side <- ifelse(right.side, 2, 1)
      nSegments <- rd$endrb.nSegments[side]

      s <- truncated.s(rd)


      # Find h0 & h1 for the last parallelepipedic riskband (# 4)

      h1.index <- ifelse(right.side, length(rd$h) - nSegments, 1 + nSegments)
      h1 <- rd$h[h1.index]

      h0.index <- h1.index - 2*outside.dir
      h0 <- rd$h[h0.index]


      s0l.index <- ifelse(right.side, 5, 4)
      s0l <- s[s0l.index] * outside.dir

      s0r.index <- s0l.index + 3*outside.dir
      s0r <- s[s0r.index] * outside.dir


      innerb.index <- ifelse(right.side, 6, 2)

      t4.start <- (rd$inner.b[innerb.index] - rd$inner.b[innerb.index-1]) * g$cos.theta

      # ------------------------------------------------------------------------
      t.parms <- t.info(G$n2l, h0, h1)

      if (!right.side)  t4.start <- t.parms$b - t4.start


      # If h.bar =~ u then there is no point searching for an optimal breakpoint t
      # in riskband 4 as (s1, s2) is very close to a straight line
      # (with coordinates stored below in straight4)

      if (t.parms$one.piece)  rb4.1piece <- t.parms$straight


      # ------------------------------------------------------------------------

      if (G$is.local.min && abs(h1 - G$my.u) < 1e-10)
      {
        # Special case, where the function h cannot go down at 0.
        #
        # Last riskband can only be flat at h = G$my.u
        # therefore there is no loop on t for the closing riskband,
        # only on t for the next-to-last riskband.

        if (t.parms$one.piece)  return(rd)


        tmp <- parallelepiped.fixedh.minTheta2Sum(t.parms, t4.start, s0l, 0)
        s <- tmp$s

        j <- ifelse(right.side, length(rd$h)-nSegments-1, 2+nSegments)
        rd$h[j] <- h0 + s[1]*tmp$t


        j <- ifelse(right.side, length(rd$s)-nSegments-1, nSegments+1)
        j <- seq(from=j, length=2)

        if (!right.side)  s <- -rev(s)
        rd$s[j] <- s


        rd$inner.b[innerb.index] <- revisited.innerb(tmp$t, rd$inner.b, innerb.index, g$cos.theta, right.side)

        rd$score <- roughness(g, rd)$score

        return(rd)
      } # end of special case


      # ----------------------------------------------------------------------------
      h05 <- h1

      convex5 <- h05 < G$h0.SE
      s1b <- southEast.connexion.slope(G, h05)


      # We will use the current slope/angle of first segment of closing riskband
      # as a starting value for the optimal value

      rb5x.start <- atan(s0r)


      # We are searching for optimal (s1, s2, t) --- on closing riskband

      # Define range on s1 (and its initial value)

      if (convex5)
      {
        untouchable.bounds <- rep(TRUE, 2)
        s.clothesLine <- clothesLine(G, h05)$s
        if (!G$rightSide.closing.rb) s.clothesLine <- -s.clothesLine


        rb5.range <- c(atan(s.clothesLine), pi/2)

        if (G$is.local.min)
        {
          rb5.range[2] <- 0
          untouchable.bounds[2] <- FALSE
        }

        rb5 <- my.bs.init(rb5.range, rb5x.start, 0.015, untouchable.bounds=untouchable.bounds, epsilon=G$epsilon)
      }
      else
      {
        max.angle <- atan(NewtonRaphson.lastriskband.fade.out.s1(G, h05)$s1)
        rb5.range <- c(-pi/2, max.angle)
        untouchable.bounds <- c(TRUE, FALSE)

        rb5 <- my.bs.init(rb5.range, rb5x.start, 0.015, untouchable.bounds=untouchable.bounds, epsilon=G$epsilon)
      }


      # ----------------------------------------------------------------------------


      while (rb5$continue)
      {
        s1 <- tan(rb5$x)  # stands for s5.a, really
          if (abs(s1) < 1e-10) s1 <- 0

        # Find (s2, t) such that the integral \int_0^D h(x)l(x)dx = riskband prob
        #              [with  theta2 = atan(s2)]

        tmp <- NewtonRaphson.lastriskband.fixed.s1(s1, h1, G, convex=convex5, s1b=s1b)
          parms <- list(t=tmp$t, s=c(s1, tmp$s2), h5m=tmp$h.intersect)


        # Run a NR search on t4

        if (t.parms$one.piece)  tmp <- rb4.1piece
        else                    tmp <- parallelepiped.fixedh.minTheta2Sum(t.parms, t4.start, s0l, s1)


        s <- c(s0l, tmp$s, parms$s)
        score <- endrb.cost(2, s, 0)

        parms$t4 <- tmp$t
        parms$s4 <- tmp$s

        rb5 <- bs.update(rb5, score, parms)
      }


      tmp <- rb5$best.parms

        s5  <- tmp$s
        t5  <- tmp$t
        h5m <- tmp$h5m

        t4 <- tmp$t4
        s4 <- tmp$s4


      # Wrap-up

      rd$h <- drop.elements(rd$h, nSegments + 2, right.side)
      rd$s <- drop.elements(rd$s, nSegments + 2, right.side)

      ht  <- h0 + t4*s4[1]
      h.etc <- c(ht, h1, h5m, 0)

      s.etc <- c(s4, s5)


      if (right.side)
      {
        rd$h <- c(rd$h, h.etc)
        rd$s <- c(rd$s, s.etc)

        rd$R5b <- g$R5b2
        rd$R5b[2] <- rd$R5b[1] + t5/g$cos.theta
      }
      else
      {
        h.etc <-  rev(h.etc)
        s.etc <- -rev(s.etc)

        rd$h <- c(h.etc, rd$h)
        rd$s <- c(s.etc, rd$s)

        rd$R1b <- g$R1b2
        rd$R1b[2] <- rd$R1b[3] - t5/g$cos.theta
      }


      rd$inner.b[innerb.index] <- revisited.innerb(t4, rd$inner.b, innerb.index, g$cos.theta, right.side)

      rd$endrb.nSegments[side] <- 2
      rd$score <- roughness(g, rd)$score

      return(rd)
    } # end of revisit.closing2riskbands.fixedh


    revisit.inner.riskbands <- function(g, rd)
    {
      # rd: revisited density

      # Revisit parallelepiped inner riskbands; that is, we keep h-values at both ends unchanged
      # but see if we can change the breakpoint location between the two lines (and their slopes accordingly)
      # to get a smoother curve.


      Avoid.entry.bump <- function(r, s)
      {
        if (r == 3)  return(FALSE)

        # r in {2, 4}

        sign2avoid <- ifelse(r == 2, -1, 1)
                 j <- ifelse(r == 2,  1, 8)

        return(sign(s[j]) == sign2avoid)
      } # end of Avoid.entry.bump


      best.side.parms <- function(tmp.l, tmp.r)
      {
        s.surround <- list(l=NA, r=NA)

        if (is.null(tmp.l$s) || !tmp.l$accepted)
        {
          angle.l <- Inf
        }
        else
        {
          s.surround$l <- c(tmp.l$s[2], tmp.l$external.slope)
          angle.l <- angle2(s.surround$l)
        }


        if (is.null(tmp.r$s) || !tmp.r$accepted)
        {
          angle.r <- Inf
        }
        else
        {
          s.surround$r <- c(tmp.r$external.slope, tmp.r$s[1])
          angle.r <- angle2(s.surround$r)
        }


        left.side <- angle.l < angle.r


        if (left.side)
        {
          s <- tmp.l$s
          t <- tmp.l$t
          r <- tmp.l$r
          angle.sides <- s.surround$l
        }
        else
        {
          s <- tmp.r$s
          t <- tmp.r$t
          r <- tmp.r$r
          angle.sides <- s.surround$r
        }


        list(s=s, t=t, r=r, angle.sides=angle.sides)
      } # end of best.side.parms


      t.mid_angle <- function(g, r, h.pivot, incoming.slope, outgoing.slope)
      {
        h0 <- h.pivot[r-1]
        h1 <- h.pivot[r]

        avoid.entry.bump <- g$middle.step[r]

        t.parms <- t.info(g, h0, h1, r, avoid.entry.bump)


        # If h.bar =~ u then there is no point searching for an optimal breakpoint t
        # in this riskband as (s1, s2) is very close to a straight line

        if (t.parms$one.piece)  return(t.parms$straight)



        t.start <- t.parms$soft
        if (t.start < t.parms$range[1] || t.start > t.parms$range[2])  t.start <- mean(t.parms$range)

        tmp <- parallelepiped.fixedh.minTheta2Sum(t.parms, t.start, incoming.slope, outgoing.slope)
        tmp$accepted <- TRUE

        return(tmp)
      } # end of t.mid_angle


      what.if.we.improve.angle <- function(g, w, s, h.pivot, a, ssc)
      {
        # Return a list with
        #
        #  i: the index of the riskband modified t score
        #  t: the corresponding modified t score
        #  s: the two slopes of the modified riskband
        #  j: the indices of the new segment slopes
        #  angle.change: the absolute change between angle[w] [original angle] and modified angle


        angle.w <- a[w]
          if (angle.w == 0) return(list(accepted=FALSE))


        if (w%%2 == 0)
        {
          r <- 1 + w/2  # r in {2, 3, 4}

          incoming.slope <- s[2*r-3]
          outgoing.slope <- s[2*r]


          tmp <- t.mid_angle(g, r, h.pivot, incoming.slope, outgoing.slope)
            if (!tmp$accepted) return(list(accepted=FALSE))
            tmp$r <- r
            angle.sides <- tmp$s
        }
        else if (w == 1)
        {
          r <- 2

          incoming.slope <- s[1]
          avoid.entry.bump <- Avoid.entry.bump(r, s)
          t.parms <- t.info(g, h.pivot[r-1], h.pivot[r], r, avoid.entry.bump)


          outgoing.slope <- s[4]
          tmp <- t.mid_angle(g, r, h.pivot, incoming.slope, outgoing.slope)
            if (!tmp$accepted) return(list(accepted=FALSE))
            tmp$r <- r
            angle.sides <- tmp$s
        }
        else if (w == 7)
        {
          r <- 4

          outgoing.slope <- s[8]
          avoid.entry.bump <- Avoid.entry.bump(r, s)
          t.parms <- t.info(g, h.pivot[r-1], h.pivot[r], r, avoid.entry.bump)


          incoming.slope <- s[5]
          tmp <- t.mid_angle(g, r, h.pivot, incoming.slope, outgoing.slope)
            if (!tmp$accepted) return(list(accepted=FALSE))
            tmp$r <- r
            angle.sides <- tmp$s
        }
        else
        {
          # w = 3 | 5

          r <- 1 + (w-1) / 2   # r in {2, 3} -- no of riskband to the left of 'w'


          # Try moving t in the riskband to the left of 'w'

          outgoing.slope <- s[2*r]
          avoid.entry.bump <- Avoid.entry.bump(r, s)
          t.parms <- t.info(g, h.pivot[r-1], h.pivot[r], r, avoid.entry.bump)
          tmp.l <- t.equal_angles(t.parms, outgoing.slope, TRUE)
            tmp.l$r <- r


          # Try moving t in the riskband to the right of 'w'

          r <- r + 1
          incoming.slope <- s[2*r-3]
          avoid.entry.bump <- Avoid.entry.bump(r, s)
          t.parms <- t.info(g, h.pivot[r-1], h.pivot[r], r, avoid.entry.bump)
          tmp.r <- t.equal_angles(t.parms, incoming.slope)
            tmp.r$r <- r


          # We pick the side with the best result (smaller modified angle [calculated in best.side.parms])

          if (!tmp.l$accepted && !tmp.r$accepted)  return(list(accepted=FALSE))

          tmp <- best.side.parms(tmp.l, tmp.r)
            angle.sides <- tmp$angle.sides
        }



        i <- tmp$r - 1
        j <- 2*i + c(0, 1)
        angle <- angle2(angle.sides)


        # Compute the new sharpest.slope.change$score: if it is decreased with the new above-computed slopes [tmp$s],
        # then the suggested slope is refused

        s[j] <- tmp$s
        tmp.ssc <- sharpest.slope.change(s)$value
        accepted <- tmp.ssc <= ssc


        return(list(s=tmp$s, t=tmp$t, i=i, j=j, angle.change=abs(angle.w-angle), accepted=accepted))
      } # end of what.if.we.improve.angle


      which.close2max <- function(x, epsilon=1e-9)
      {
        m <- max(x)
        which(abs(x-m) < epsilon)
      } # end of which.close2max


      # --- Start of revisit.inner.riskbands -----------------------------------

      any.change <- FALSE

      rb1.size <- rd$endrb.nSegments[1]
      h.index  <- rb1.size + 2*seq(3) # indices for the middle-riskband height of riskbands 2,3,4 [length 3]


      h05 <- rd$h[rb1.size + 7]

      h.pivot <- rd$h[h.index-1] # h at the entry [left-side] in each riskband [length 3]
      h.pivot <- c(h.pivot, h05)

      s <- truncated.s(rd)
      t <- (rd$inner.b[c(2,4,6)] - rd$inner.b[c(1,3,5)]) * g$cos.theta


      continue <- TRUE


      while (continue)
      {
        # Identify the angles 'w' that are (very) close to max angle,
        # loop over them and keep the one that was most modified by a change in t

        a <- angles(s)
        which.sharper <- which.close2max(a)
          larger.angle.change <- -Inf
          ssc <- sharpest.slope.change(s)$value


        for (w in which.sharper)
        {
          tmp <- what.if.we.improve.angle(g, w, s, h.pivot, a, ssc)

          if (tmp$accepted)
          {
            if (tmp$angle.change > larger.angle.change)
            {
              # Yes, perform the above two if-checks separately (angle.change may very well not be defined when !accepted)
              larger.angle.change <- tmp$angle.change
              best.wres <- tmp
            }
          }
        }


        if (is.infinite(larger.angle.change))  break


        # Accept the modified segment

        continue <- abs(t[best.wres$i] - best.wres$t) > 1e-5

        t[best.wres$i] <- best.wres$t
        s[best.wres$j] <- best.wres$s
        any.change <- TRUE
      }


      # See if we can still improve median angle for each inner riskband

      J <- c(2, 4, 6)

      a <- angles(s)
      ssc <- sharpest.slope.change(s)$value

      a.median_angles <- a[J]
      o <- order(a.median_angles, decreasing=TRUE)
      ordered.w <- J[o]


      for (w in ordered.w)
      {
        tmp <- what.if.we.improve.angle(g, w, s, h.pivot, a, ssc)

        if (tmp$accepted)
        {
          t[tmp$i] <- tmp$t
          s[tmp$j] <- tmp$s

          a <- angles(s)
          ssc <- sharpest.slope.change(s)$value

          any.change <- TRUE
        }
      }


      # Wrap-up (update rd's coordinates if any change was performed/accepted)

      if (!any.change)  return(list(better=FALSE))


      h.pivot <- h.pivot[-4]
      rd$h[h.index] <- h.pivot + t * s[c(2,4,6)]

      s <- s[-c(1,8)]
      j <- seq(from=rb1.size+1, length=6)
      rd$s[j] <- s

      j <- 2*seq(3)
      rd$inner.b[j] <- rd$inner.b[j-1] + t / g$cos.theta

      rd$score <- roughness(g, rd)$score


      return(list(better=TRUE, d=rd))
    } # end of revisit.inner.riskbands


    revisited.range <- function(range, r, g, h0)
    {
      if (g$is.local.max[r])
      {
        tmp <- plateau(r, g, h0)$theta
        if (tmp > range[1])  range[1] <- tmp

        tmp <- atan(parallelepiped.connexion(r, g, h0, g$u[r+1]))
        if (tmp < range[2])  range[2] <- tmp
      }
      else if (g$is.local.min[r])
      {
        tmp <- atan(parallelepiped.connexion(r, g, h0, g$u[r+1]))
        if (tmp > range[1])  range[1] <- tmp
      }
      else if (g$going.up[r])
      {
        alt.max <- plateau(r, g, h0)$theta
        if (alt.max < range[2])  range[2] <- alt.max

        if (h0 <= g$u[r] && g$middle.step.up[r] && range[1] < 0)  range[1] <- 0
      }
      else if (g$going.dn[r])
      {
        alt.min <- plateau(r, g, h0)$theta
        if (alt.min > range[1])  range[1] <- alt.min

        if (h0 >= g$u[r] && g$middle.step.dn[r] && range[2] > 0)  range[2] <- 0
      }


      if      (g$middle.step.dn[r])  range <- pmin(range, 0)
      else if (g$middle.step.up[r])  range <- pmax(range, 0)


      if (g$going.up.next2rb[r])
      {
        theta <- atan(parallelepiped.connexion(r, g, h0, g$u[r+1]))
        if (theta > range[1])  range[1] <- theta
      }
      else if (g$going.dn.next2rb[r])
      {
        theta <- atan(parallelepiped.connexion(r, g, h0, g$u[r+1]))
        if (theta < range[2])  range[2] <- theta
      }


      return(range)
    } # end of revisited.range
    

    riskband3.onward <- function(h03, g, g.right, reg, incoming.slope.sign3,
                                 rb3.suggested.2nd.step=NULL, print1=FALSE)
    {
      ask.the.wizard.3onward <- function(key)
      {
        if (key$w  > 2)  return(0)
        if (key$w == 2)  return(-key$sign)
        return(key$sign)
      } # end of ask.the.wizard.3onward


      rb3.inits.revisited <- function(g, rb3, h03)
      {
        rb4.suggested.2nd.step <- NULL

        # If the last three riskbands are going up, then force h04 to start at g$u[4]
        # to ensure a point in the search that will not have a local max in riskband # 4

        if (g$u[3] <= g$u[4] && g$u[4] <= g$u[5])
        {
          tentative.x <- atan(parallelepiped.connexion(r=3, g, h03, g$u[4]))

          if (tentative.x >= rb3$range[1] && tentative.x <= rb3$range[2])
          {
            rb3$suggested.2nd.step <- rb3$x # relegate the calculated initial step to 2nd step
            rb3$x <- tentative.x
            rb3$started.on.boundj <- which1(rb3$x == rb3$range)
            rb3$x.on.bound <- rb3$started.on.boundj > 0
            if (rb3$started.on.boundj == 0)  rb3$x <- rb3$x + 1e-8 # make sure to get h04 a wee bit below g$u[4]

            rb4.suggested.2nd.step <- 0
          }
        }

        return(list(rb3=rb3, rb4.suggested.2nd.step=rb4.suggested.2nd.step))
      } # end of rb3.inits.revisited


      # Riskband 3

      bg <- reg.predict(reg, 3, h03)

      rb3x_range <- parallelepiped.range(r=3, g, h03)
        rb3x_range <- revisited.range(rb3x_range, 3, g, h03)


      rb3 <- my.bs.init(rb3x_range, bg$pred, bg$margin, suggested.2nd.step=rb3.suggested.2nd.step, epsilon=g$epsilon)
        tmp <- rb3.inits.revisited(g, rb3, h03)
        rb3 <- tmp$rb3
        rb4.suggested.2nd.step <- tmp$rb4.suggested.2nd.step


      print4 <- FALSE

      if (print1)
      {
        catn('**** h03 *******************', h03)
        catn('rb3x_range', rb3x_range)
      }
        
  
      while (rb3$continue)
      {
        tmp <- parallelepiped(r=3, g, rb3$x, h03)
        s3 <- tmp$s
        h3 <- tmp$h
        body3 <- list(h=h3, s=s3, nSegments=-1) #, endrb.nSegments=c(0,0))


        tmp <- corrected.onePiece.rb3(body3, g)
        if (tmp$corrected)  body3 <- tmp$rb3
       
  
        # Riskband 4
 
        h04 <- h3[3]
        rb4.incoming.slope <- body3$s[2]

         
        o <- riskband4.onward(h04, g, g.right, reg, rb4.incoming.slope, rb4.suggested.2nd.step, print4=print4)
          reg <- o$reg

        rb4 <- o$rb4
        rb4$endrb.nSegments <- c(0, rb4$nSegments)

        tmp <- augmented.body(rb4, body3)
                     
        rough <- roughness(g, tmp, incoming.slope.sign=incoming.slope.sign3)
          out3 <- rough$score

          where2go.next <- ifelse(g$use.educated.guess,
                                  ask.the.wizard.3onward(rough$key),
                                  0)


        if (print1)
        {
          catn('#steps', rb3$n.steps)
          catn('rb3.x', rb3$x)
          catn('h04', h04)
          catn('out3', out3)
          catn('where2go.next', where2go.next)
        }

        rb3 <- bs.update(rb3, out3, tmp, where2go.next)
        rb4.suggested.2nd.step <- NULL


        if (print1 && rb3$last.score.was.the.best) catn('<-- best h04 so far')
      }
      
      rb3 <- rb3$best.parms
      reg <- reg.update(reg, 3, rb3)


      
      list(rb3=rb3, reg=reg)
    } # end of riskband3.onward


    riskband4.onward <- function(h04, g, G, reg, incoming.slope, suggested.2nd.step=NULL, for.rb12=FALSE, print4=FALSE)
    {
      # G: either g.right (for riskbands 4 & 5)
      #        or g.left  (for riskbands 1 & 2)


      # Riskbands 4+5

      reg.rb <- ifelse(for.rb12, 2, 4)

      rb4x_range <- rb4x.range(g, h04)


      if (print4)
      {
        catn('//// 4-onward debut')
        catn('rb4x_range', rb4x_range)
      }


      bg <- reg.predict(reg, reg.rb, h04)
      rb4 <- my.bs.init(rb4x_range, bg$pred, bg$margin, suggested.2nd.step=suggested.2nd.step, epsilon=g$epsilon)


      while (rb4$continue)
      {
        tmp <- parallelepiped(r=4, g, rb4$x, h04)
        s4 <- tmp$s
        h4 <- tmp$h[-1]
        if (h4[2] < 0 & h4[2] > -1e-8) h4[2] <- 0 # correction for numerical imprecision
        body4 <- list(h=c(h04, h4), s=s4, endrb.nSegments=c(0,0))


        # Riskband 5

        rb5 <- quick.close(G, h4[2])
        rb5$rightSide.closing.rb <- TRUE


        rough <- roughness(g, body4, rb5, incoming.slope=incoming.slope)
          out4 <- rough$score

          where2go.next <- ifelse(g$use.educated.guess, ask.the.wizard(rough$key), 0)

          tmp <- list(rb5=rb5, rb4=list(h=h4, s=s4))


          if (print4)
          {
            catn('---')
            catn('step4#', rb4$n.steps)
            catn('rb4$x', rb4$x)
            catn('h4', h4)
            h05 <- h4[2]
            catn('out4', out4)
            catn('h05', h05)
            catn('where2go.next', where2go.next)
            catn('key', unlist(rough$key))

            catn('body4-s', body4$s)
            catn('body4-h', body4$h)
            catn('rb5-s', rb5$s)
            catn('rb5-h', rb5$h)
          }


        rb4 <- bs.update(rb4, out4, tmp, where2go.next)


        if (print4 && rb4$last.score.was.the.best) catn('<---- best h04 so far ---------------//')
      } # end of while-rb4

      if (print4) catn('4-onward fin ////')


      on.the.edge <- any(rb4$visited[rb4$which.best] == rb4x_range)

      rb5 <- rb4$best.parms$rb5
      rb4 <- rb4$best.parms$rb4
        rb4$nSegments <- -1

                       
      rb4 <- augmented.body(rb5, rb4)
        rb4$h <- c(h04, rb4$h)
        rb4$rightSide.closing.rb = G$rightSide.closing.rb
        rb4$on.the.edge <- on.the.edge

      reg <- reg.update(reg, reg.rb, rb4)

      tmp <- corrected.onePiece.rb4(rb4, G)
      if (tmp$corrected) rb4 <- tmp$rb4

      list(rb4=rb4, reg=reg)
    } # end of riskband4.onward


    roughness <- function(g, pb, closing.rb=list(), incoming.slope=NULL, incoming.slope.sign=1, ssc.threshold=-2)
    {
      # g: set of geometrics -- can be the complete g list, or one of its subsets [g.left or g.right]
      # pb [partial body]: a list with dimensions(h, s, endrb.nSegments)
      # closing.rb: provided if one of the closing riskbands is re-estimated 
      #             (that is, AFTER the curve was estimated throughout, to see if we can find a 
      #              better density closing).
      #             > When provided, closing.rb is a list with dimensions $s, $h & $rightSide.closing.rb
      #
      # ssc.threshold: used when the value of sharpest.slope.change [ssc] is of no interest when smaller to that value
      #                --- indicating that the current distrn is already better than the reference distrn [or last considered best distrn])
      #                and that the 'beauty criteria' used after ssc should be considered rather than ssc --- this will be done by setting the ssc value
      #                to a fixed low value.


      if (!is.null(incoming.slope))  incoming.slope.sign <- sign(incoming.slope)
      
      s <- pb$s
      h <- pb$h
      
      endrb.nSegments <- pb$endrb.nSegments
        if (is.null(endrb.nSegments))  endrb.nSegments <- rep(0, 2)
      closing.rb.provided <- length(closing.rb) > 0


      if (closing.rb.provided)
      {
        rightSide.closing.rb <- closing.rb$rightSide.closing.rb  # will be used later on, so defining it here is not useless (!!)


        if (rightSide.closing.rb)
        {  
          s <- c(s, closing.rb$s)
          h <- c(h, closing.rb$h)
        }
        else
        {                              
          s <- c(closing.rb$s, s)
          h <- c(closing.rb$h, h)
        }
        
        
        j <- ifelse(rightSide.closing.rb, 2, 1)
        endrb.nSegments[j] <- ifelse(length(closing.rb$s)%%2 == 1, 1, 2)
      }
      else rightSide.closing.rb <- FALSE

      
      s.len <- length(s)
      h.len <- length(h)
      
      h.2nd.peak <- second.higher.peak.height(h, s, incoming.slope.sign)
             
  
      # Correct for (most likely) numerical imprecision

      w <- which(abs(s) < 1e-12)
      s[w] <- 0
        

      # Final slope in closing riskbands should not be below some value, as indicated by s.fade.out
      
      if      (                           h[h.len] == 0 && s[s.len] > g$s.fade.out[2])  end.too.brutal <- TRUE
      else if (endrb.nSegments[1] > 0  && h[1]     == 0 && s[1]     < g$s.fade.out[1])  end.too.brutal <- TRUE
      else end.too.brutal <- FALSE
      

      # See if density fct is going too deep (when compared to the piecewise uniform prior)

      if (!g$impossible)
      {        
        h.last.index  <- h.len - endrb.nSegments[2]
        h.first.index <- max(h.last.index - 6, 1)

        h.indices <- seq(from=h.first.index, to=h.last.index)
  
        ulow.indices <- seq(to=7, by=1, length=length(h.indices))
        ratios <- h[h.indices] / g$u.low[ulow.indices]
        too.deep <- any(ratios < g$max.deepness & (g$max.deepness - ratios) > 1e-6) 
      }
      else  too.deep <- -1
      
      
      # Number of sign changes -- number of Summits[up/down]

      augmented.s <- c(s, -1)
      if (endrb.nSegments[1] == 0) augmented.s <- c(incoming.slope.sign, augmented.s)
      else                         augmented.s <- c(1, augmented.s) # imaginary incoming positive slope
      
      number.of.sSign.changes <- n.sign.changes(augmented.s)   
      
  
      # Remove leading/trailing slopes from closing riskbands (if there are two slopes in them)
      
      s.complete <- s

      if (s.len > 2)  tails <- tail.angles(g, endrb.nSegments, s, h, incoming.slope)
      else            tails <- numeric(0)


      if (endrb.nSegments[1] > 0)
      {
        if (endrb.nSegments[1] > 1)  s <- s[-1]
      }
      else if (!is.null(incoming.slope))  s <- c(incoming.slope, s)


      if (endrb.nSegments[2] > 1)  s <- s[-length(s)]



      # Number of inflexion points
      #   (ignored when we have 2+ peaks+hollows)

      if (number.of.sSign.changes > 1 || (number.of.sSign.changes == 1 && too.deep == 0))
      {
        n.inflex <- -1
      }
      else
      {
        n.inflex <- n.inflexion.points(s.complete)
        if (n.inflex == 1)  n.inflex <- 0 # one inflexion point is as acceptable as no inflexion pts (but more than that is bad)
      }
      
      

      # Find the sharpest turn / slope change

      ssc <- sharpest.slope.change(s)
      ssc.value <- ssc$value
        if (length(ssc.value) > 0 && ssc.value < ssc.threshold)  ssc.value <- -2


      # Max height of h (excluding the two endpoint values)

      key.htrim <- h
      if (endrb.nSegments[1] > 1)  key.htrim <- key.htrim[-1]
      if (endrb.nSegments[2] > 1)  key.htrim <- key.htrim[-length(key.htrim)]
      
      h <- h[-length(h)]
      if (endrb.nSegments[1] > 0) h <- h[-1]
      h.max <- max(h)
      

      h.max       <- round(h.max,       digits=5)
      h.2nd.peak  <- round(h.2nd.peak,  digits=5)


      score <- c(number.of.sSign.changes,
                 too.deep,
                 ssc.value,
                 n.inflex,
                 tails,
                 h.max,
                 end.too.brutal,
                 h.2nd.peak)


      key <- list(w=ssc$w, sign=ssc$sign, sSum.sign=ssc$sSum.sign, null.s1=ssc$null.s1,
                  incoming.slope.included=!is.null(incoming.slope),
                  h=key.htrim,
                  u=g$u, ground=g$ground, is.local.min3=g$is.local.min[3],
                  middle.step=g$middle.step, endrb.nSegments=endrb.nSegments)

      return(list(score=score, key=key))
    } # end of roughness

    
    second.higher.peak.height <- function(h, s, incoming.slope.sign)
    {
      s <- c(incoming.slope.sign, s)
      is.local.max <- s[-length(s)] > 0 & s[-1] <= 0
      w.local.max <- which(is.local.max)
      
      if (length(w.local.max) < 2) return(-1)
      
      local.max.heights <- h[w.local.max]
      local.max.heights <- sort(local.max.heights, decreasing=TRUE)
      local.max.heights[2]
    } # end of second.higher.peak.height


    sharpest.slope.change <- function(s)
    {
      S <- length(s)
      s.module <- sqrt(1+s^2)

      s.j <- s[-S]
      s.module.j <- s.module[-S]
      # following slope
      s.k <- s[-1]
      s.module.k <- s.module[-1]

      scalar.product <- 1 + s.j*s.k
      tmp <- scalar.product / s.module.j / s.module.k

      w <- which.min(tmp)
      value <- -tmp[w]
      slopeChange.sign <- sign(s[w+1] - s[w])
      sSum.sign        <- sign(s[w+1] + s[w])
      null.s1 <- abs(s[w]) < 1e-10


      return(list(value=value, w=w, sign=slopeChange.sign, sSum.sign=sSum.sign, null.s1=null.s1))
    } # end of sharpest.slope.change


    Sharpest.slope.change <- function(d)
    {
      s <- d$s

      if (d$endrb.nSegments[1] > 1)  s <- s[-1]
      if (d$endrb.nSegments[2] > 1)  s <- s[-length(s)]

      return(sharpest.slope.change(s)$value)
    } # end of Sharpest.slope.change


    side.info <- function(G)
    {
      list(h0.SE=G$h0.SE, SE.connexion=G$SE.connexion, cl=G$cl)
    } # end of side.info


    slope2angle <- function(s, sp=0, ss=0)
    {
      # Return the angle f=phi below a line of slope s
      # as well as fp = f' = df/dt and fs = f'' = d2f/dt2
      # where
      # sp = ds/dt
      # ss = d2s/dt2

      f <- atan(s)

      # ___________________________________________________________
      #
      # d{atan(s(t))} / dt = 1 / (1 + s^2) * ds/dt  = s' / (1+s^2)
      #
      # where f = atan(s(t))
      # ___________________________________________________________


      v <- 1 + s^2

      fp <- sp / v
      fs <- ss/v - 2 * s * sp^2 / v^2

      return(list(f=f, fp=fp, fs=fs))
    } # end of slope2angle


    slopes2angle <- function(phi1, phi2)
    {
      # Return the angle between two slopes s1(t) and s2(t)
      # each described by phi_i = slope2angle(s_i, dsi_dt, d2si_dt2)

      f <- phi2$f - phi1$f
      m <- sign(f)

      return(list(f=abs(f), fp=(phi2$fp - phi1$fp)*m, fs=(phi2$fs-phi1$fs)*m))
    } # end of slopes2angle


    southEast.connexion <- function(G, h0, s)
    {
      h <- c(h0 + s * G$b, 0)
      s2 <- - h[1] / G$d
      s <- c(s, s2)

      if (G$side < 0)
      {
        h <- rev(h)
        s <- -rev(s)
      }

      return(list(h=h, s=s, x95b=G$x95b2, rightSide.closing.rb=G$rightSide.closing.rb, nSegments=2))
    } # end of southEast.connexion

    
    southEast.connexion.slope <- function(G, h0)
    {
      # G: closing riskband geometrics
      
      # Return the slope of the 1st segment of the line which, broken at t=b, goes south-east to (D,0)
      # and leads to the appropriate riskband prior probability
            
      G$SE.connexion$a - h0 * G$SE.connexion$b
    } # end of southEast.connexion.slope     


    tail.angles <- function(g, endrb.nSegments, s, h, incoming.slope)
    {
      # Complete s

      if (!is.null(incoming.slope))
      {
        if (endrb.nSegments[1] > 0)  s <- c(s, incoming.slope)
        else                         s <- c(incoming.slope, s)
      }


      tails <- numeric(0)

      # Riskband 5
      if (endrb.nSegments[2] > 0)  tails <- endrb.cost(endrb.nSegments[2], s, h[length(h)])

      # Riskband 1
      if (endrb.nSegments[1] > 0)
      {
        tmp <- endrb.cost(endrb.nSegments[1], s, h[1], FALSE)
        tails <- c(tails, tmp)
      }


      if (!g$forward)  tails <- rev(tails)

      return(tails)
    } # end of tail.angles


    t.equal_angles <- function(t.parms, s0, s0.is.outgoing.slope=FALSE, NR.precision=1e-6)
    {
      # s0: either incoming slope  [if s0.is.outgoing.slope=FALSE]
      #         or outgoing slope  [if s0.is.outgoing.slope=TRUE]
      #
      # Return the value of t for which the angle between incoming slope & s1 is equal to the angle between s1 & s2
      #  (where s1 and s2 meet at t and provide the expected riskband probability)
      #  -- solved by Newton-Raphson
      #
      # If the algorithm fails to converge => soln does not seem to exist: we then find the point where the angles are the closest
      # (that is, where f = angle2 - angle1 reaches a minimum value)

      # The above is true when incoming slope is specified [that is, when s0.is.outgoing.slope=FALSE];
      #   if it is rather outgoing slope that is specified [s0.is.outgoing.slope=TRUE], then we try to equalize the angle between s1 & s2
      #   and the angle between s2 & outgoing slope s0.


      better.t <- function(t.parms, s0, left.side.angle, epsilon)
      {
        # When this fct is called, one of the t.range interval limits must be different from 0 & b

        t.range <- t.parms$range
        b <- t.parms$b


        t <- setdiff(t.range, c(0, b))

        if (length(t) == 0)  return(list(t=NA))


        if (length(t) == 2)
        {
          f  <- rep(NA, 2)
          fp <- rep(NA, 2)

          for (i in 1:2)
          {
            tmp <- h(t[i], t.parms, s0, left.side.angle)
            f[i]  <- tmp$f
            fp[i] <- tmp$fp
          }

          zero.exists <- sign(f[1]) != sign(f[2])

          if (!zero.exists)
          {
            j <- which.min(abs(f))
            t <- t[j]

            s <- s.t(t, t.parms)
            return(list(t=t, s=s))
          }


          j <- which.max(abs(fp))
          t <- t[j]
          fp.sign <- sign(fp[j])
        }
        else
        {
          tmp <- h(t, t.parms, s0, left.side.angle)

          t.side <- which(t.range==t)
          m <- ifelse(t.side==1, -1, 1)

          fp.sign <- sign(tmp$fp)
          zero.exists <- sign(tmp$f) * m == fp.sign

          if (!zero.exists)  return(list(t=t, s=tmp$s))
        }


        continue <- TRUE


        while (continue)
        {
          tmp <- h(t, t.parms, s0, left.side.angle)
          change <- - fp.sign * tmp$f / abs(tmp$fp)

          j <- ifelse(change > 0, 1, 2)
          t.range[j] <- t

          t <- t + change
          if (t <= t.range[1] || t >= t.range[2])  t <- mean(t.range)

          converged <- abs(change) < epsilon
          squeezed  <- diff(t.range) < epsilon
          continue <- !converged && !squeezed
        }

        s <- s.t(t, t.parms)
        return(list(t=t, s=s))
      } # end of better.t


      h <- function(t, t.parms, s0, left.side.angle)
      {
        phi0 <- slope2angle(s0)

        k1 <- t.parms$k1
        k2 <- t.parms$k2
        k3 <- t.parms$k3
        alpha <- t.parms$alpha
        b <- t.parms$b


        s1 <- alpha + k1 / t
        s2 <- (k2 - alpha*t) / (b - t)

        s1p <- -k1 / t^2             # First derivatives
        s2p <-  k3 / (b - t)^2

        s1s <- 2 * k1 / t^3          # Second derivatives
        s2s <- 2 * k3 / (b-t)^3


        phi1 <- slope2angle(s1, s1p, s1s)
        phi2 <- slope2angle(s2, s2p, s2s)


        if (left.side.angle)  theta.ext <- slopes2angle(phi1, phi0)
        else                  theta.ext <- slopes2angle(phi2, phi0)


        theta.int <- slopes2angle(phi1, phi2)


        f  <- theta.ext$f  - theta.int$f    # This is what we are trying to minimize!
        fp <- theta.ext$fp - theta.int$fp   # <- its first two derivatives are:
        fs <- theta.ext$fs - theta.int$fs

        return(list(f=f, fp=fp, fs=fs, s=c(s1,s2)))
      } # end of h


      NR <- function(t.init, qA.sign, t.parms, s0, range, left.side.angle, epsilon, looking4left.side.soln=TRUE)
      {
        expected.fp.sign <- qA.sign * ifelse(looking4left.side.soln, -1, 1)
        t <- t.init

        converged <- FALSE
        squeezed  <- FALSE
        continue  <- TRUE


        while (continue)
        {
          tmp <- h(t, t.parms, s0, left.side.angle)

          change <- - tmp$f * expected.fp.sign / abs(tmp$fp)
            j <- ifelse(change > 0, 1, 2)
            range[j] <- t

          t <- t + change

          if (t < range[1] || t > range[2])  t <- mean(range)

          converged <- abs(change) < epsilon

          if (!converged)  squeezed <- diff(range) < epsilon
          continue <- !converged && !squeezed
        }


        tmp <- h(t, t.parms, s0, left.side.angle)

        return(list(t=t, f=tmp$f, converged=converged, s=tmp$s))
      } # end of NR


      s.t <- function(t, t.parms)
      {
        k1 <- t.parms$k1
        k2 <- t.parms$k2
        alpha <- t.parms$alpha
        b <- t.parms$b


        s1 <- alpha + k1 / t
        s2 <- (k2 - alpha*t) / (b - t)
        s <- c(s1, s2)

        return(s)
      } # end of s.t


      t.opt <- function(t.parms, s0, left.side.angle, epsilon)
      {
        # Find the nearby min/max

        range <- t.parms$range
        t     <- t.parms$soft

        continue <- TRUE


        while (continue && t >= range[1] && t <= range[2])
        {
          tmp <- h(t, t.parms, s0, left.side.angle)

          if (tmp$fp == 0 || is.nan(tmp$fp)) break

          change <- - tmp$fp / tmp$fs

          j <- ifelse(change > 0, 1, 2)
          range[j] <- t

          t <- t + change
          if (t < range[1] || t > range[2])  t <- mean(range)


          converged <- abs(change) < epsilon
          continue <- !converged

          if (continue)
          {
            squeezed <- diff(range) < epsilon
            continue <- !squeezed
          }
        }


        if (t < range[1] || t > range[2])  return(list(t=NA))


        tmp <- h(t, t.parms, s0, left.side.angle)
        tmp$t <- t


        # Estimate the parameters of the quadratic approximation of 'angle' around t_opt

        A <- tmp$fs / 2
        B <- -2 * A * t
        C <- tmp$f - A*t^2 - B*t

        tmp$has.zeros <- sign(tmp$f) != sign(A) && A != 0 && !is.nan(A)


        if (tmp$has.zeros)
        {
          m <- sqrt(B^2 - 4*A*C) / 2 / abs(A)
          tmp$zeros.approx <- t + c(-1, 1) * m
          tmp$qA.sign <- sign(A)

          # Make sure to use zeros.approx that are within range

          if (tmp$zeros.approx[1] <= range[1])  tmp$zeros.approx[1] <- (range[1] + t) / 2
          if (tmp$zeros.approx[2] >= range[2])  tmp$zeros.approx[2] <- (range[2] + t) / 2
        }


        return(tmp)
      } # end of t.opt


      # ------------------------------------------------------------------------
      # Start of t.equal_angles

      delta.h <- t.parms$delta.h
      b       <- t.parms$b


      left.side.angle <- !s0.is.outgoing.slope


      if (t.parms$one.piece)
      {
        s <- delta.h / b
        return(list(t=b/2, s=rep(s, 2), accepted=TRUE, external.slope=s0))
      }



      tmp <- t.opt(t.parms, s0, left.side.angle, NR.precision)

        if (is.na(tmp$t))
        {
          tmp <- better.t(t.parms, s0, left.side.angle, NR.precision)
          tmp$accepted <- !is.na(tmp$t)
          tmp$external.slope <- s0
          return(tmp)
        }


      t_opt <- tmp$t

      if (!tmp$has.zeros)
      {
        s <- s.t(t_opt, t.parms)
        return(list(t=t_opt, s=s, accepted=TRUE, external.slope=s0))
      }

      qA.sign <- tmp$qA.sign


      # We run a search to the left of t_opt and one to its right and keep the best solution

      t.start <- tmp$zeros.approx

      range <- t.parms$range
      if (t_opt > range[1] && t_opt < range[2])  range[2] <- t_opt
      soln1 <- NR(t.start[1], qA.sign, t.parms, s0, range, left.side.angle, NR.precision)

      range <- t.parms$range
      if (t_opt > range[1] && t_opt < range[2])  range[1] <- t_opt
      soln2 <- NR(t.start[2], qA.sign, t.parms, s0, range, left.side.angle, NR.precision, FALSE)


      # Now we have to compare these two solns

      s1 <- s.t(soln1$t, t.parms)
      s2 <- s.t(soln2$t, t.parms)


      if (left.side.angle)
      {
        s1 <- c(s0, s1)
        s2 <- c(s0, s2)
      }
      else
      {
        s1 <- c(s1, s0)
        s2 <- c(s2, s0)
      }


      a1 <- max(angles(s1))
      a2 <- max(angles(s2))

      if (a1 <= a2)  return(list(t=soln1$t, s=soln1$s, accepted=TRUE, external.slope=s0))

                     return(list(t=soln2$t, s=soln2$s, accepted=TRUE, external.slope=s0))
    } # end of t.equal_angles


    theta4.limits <- function(g, h04)
    {
      r <- g$R - 1  # riskband r = 4
      
      s.max <- parallelepiped.connexion(r, g, h04, g$u[g$R])
              
      # Plateau
      s.plateau <- 8/3/g$b[r] * (g$hbar[r] - h04)
      
      s <- c(s.plateau, s.max)
      
      atan(s)
    } # end of theta4.limits


    t.geometrix <- function(b, p, h.bar, delta.h)
    {
      # Calculate 't soft spot'
      #   (the value of t for which the middle angle is the smoothest, that is, directly
      #    perpendicular/above the h0-h1 straight line midpoint)

      u <- p / b
      convexity <- ifelse(h.bar < u, 1, -1)
      one.piece <- abs(h.bar-u) < 1e-6

      if (delta.h == 0)  return(list(soft=b/2, convexity=convexity, one.piece=one.piece))


      ym <- h.bar
      direct.slope <- delta.h / b
      gamma <- -1 / direct.slope
      xm <- b / 2
      h1 <- delta.h/2 + h.bar

      soft <- (ym - gamma*xm - 2*p/b + h1) / (direct.slope - gamma)


      return(list(soft=soft, convexity=convexity, one.piece=one.piece))
    } # end of t.geometrix


    t.info <- function(my.g, h0, h1, r=NULL, avoid.entry.bump=FALSE)
    {
      # my.g: either g or g.left$n2l or g.right$n2l


      t.flat <- function(h0, h1, b, p)
      {
        delta.h <- h1 - h0

        u <- p / b
        h.bar <- (h0 + h1) / 2

        convexity <- ifelse(h.bar < u, 1, -1)

        right.side <- (convexity > 0 && h1 > h0) || (convexity < 0 && h1 < h0)

        if (right.side)  return(2*(h1*b-p)/delta.h)
                         return(2*(h0*b-p)/delta.h + b)
      } # end of t.flat



      delta.h <- h1 - h0

      middle.step    <- my.g$middle.step
      ground         <- my.g$ground
      b              <- my.g$b
      p              <- my.g$p
      u              <- my.g$u
      going.up       <- my.g$going.up
      going.dn       <- my.g$going.dn
      going.west.up  <- my.g$going.west.up
      going.west.dn  <- my.g$going.west.dn


        if (!is.null(r))
        {
          middle.step <- middle.step[r]
          ground <- ground[r-1]
          b <- b[r]
          p <- p[r]
          u <- u[r]

          going.up       <- going.up[r]
          going.dn       <- going.dn[r]
          going.west.up  <- going.west.up[r]
          going.west.dn  <- going.west.dn[r]
        }


      # Make sure that t-range does not allow values that would need to go underground (for h)



      if ((!is.null(r) && r == 4) || my.g$side != 0)
      {
        t.min <- ifelse(h0 < u & h1 > u, ((h0 + h1)*b - 2*p) / (h1 - h0), 0)
        if (t.min < 0)  t.min <- 0
      }
      else t.min <- 0


      range <- c(t.min, b)
        untouchable.bounds <- rep(TRUE, 2)



      delta.h.sign <- sign(delta.h)
      h.bar <- (h0 + h1) / 2


      tmp <- t.geometrix(b, p, h.bar, delta.h)
        convexity  <- tmp$convexity
        soft       <- tmp$soft
        one.piece  <- tmp$one.piece


      if (one.piece)  straight <- list(t=b/2, s=rep(delta.h/b, 2), accepted=TRUE)
      else            straight <- NULL


      if (delta.h == 0 && h0 == u)  return(list(range=rep(b/2, 2), soft=b/2, untouchable.bounds=rep(FALSE,2), delta.h=0, convexity=0, h1=h1, b=b, p=p, one.piece=TRUE, straight=straight))



      if (avoid.entry.bump && convexity > 0)
      {
        t.critic <- t.flat(h0, h1, b, p)

        if (r == 2)
        {
          if (t.critic > range[1])
          {
            range[1] <- t.critic
            untouchable.bounds[1] <- FALSE
          }
        }
        else if (r == 3)
        {
          if (my.g$middle.step.up[3])
          {
            if (t.critic > range[1])
            {
              range[1] <- t.critic
              untouchable.bounds[1] <- FALSE
            }
          }
          else if (my.g$middle.step.dn[3])
          {
            if (t.critic < range[2])
            {
              range[2] <- t.critic
              untouchable.bounds[2] <- FALSE
            }
          }
        }
        else
        {
          # r == 4

          if (t.critic < range[2])
          {
            range[2] <- t.critic
            untouchable.bounds[2] <- FALSE
          }
        }
      }



      if (delta.h.sign != 0)
      {
        t.critic <- - (2*p - (h1 + ground)*b) / delta.h

        if (t.critic > range[1] & t.critic < range[2])
        {
          j <- ifelse(delta.h.sign < 0, 2, 1)
          range[j] <- t.critic
          untouchable.bounds[j] <- FALSE
        }
      }



      t.critic <- list(side=0)
        # side = 1 for new min value
        #      = 2 for new max value
        #      = 0 in the absence of new min/max value


      if (delta.h > 0)
      {
        if (convexity > 0)
        {
          if (going.up)          t.critic <- list(side=2, t=t.flat(h0, h1, b, p))  # Max
        }
        else if (going.west.dn)  t.critic <- list(side=1, t=t.flat(h0, h1, b, p))  # Min
      }
      else if (convexity < 0)
      {
        if (going.dn)            t.critic <- list(side=2, t=t.flat(h0, h1, b, p))  # Max
      }
      else if (going.west.up)    t.critic <- list(side=1, t=t.flat(h0, h1, b, p))  # Min



      if (t.critic$side > 0 && t.critic$t > range[1] && t.critic$t < range[2])
      {
        j <- t.critic$side
        range[j] <- t.critic$t
        untouchable.bounds[j] <- FALSE
      }



      if      (soft > range[2])  soft <- range[2]
      else if (soft < range[1])  soft <- range[1]


      alpha <- delta.h / b
      beta  <- 2*u - h1

      k1 <- beta - h0
      k2 <- 2 * (h1-u)
      k3 <- k2 - alpha*b



      return(list(b=b, convexity=convexity, delta.h=delta.h, h1=h1,
                  p=p, range=range, soft=soft,
                  one.piece=one.piece, straight=straight,
                  h0=h0, alpha=alpha, beta=beta, k1=k1, k2=k2, k3=k3,
                  untouchable.bounds=untouchable.bounds))
    } # end of t.info


    which1 <- function(cond)
    {
      w <- which(cond)

      if (length(w) == 0) return(0)
      if (length(w) > 1)  return(w[1])
      return(w)
    } # end of which1
    
    
    wrap <- function(best, fwd)
    {
      best$r <- rep(seq(5), c(length(best$R1b)-1, 2, 2, 2, length(best$R5b)-1))
      best$x95 <- c(best$R1b, best$inner.b[-1], best$R5b[-1])
      best$forward <- fwd
      
      return(best)
    } # end of wrap


    # --------------------------------------------------------------------------
    # Start of smooth

    g.left  <- geometry.closing.riskband(g)
    g.right <- geometry.closing.riskband(g, TRUE)

    g.bwd.left  <- geometry.closing.riskband(g.bwd)
    g.bwd.right <- geometry.closing.riskband(g.bwd, TRUE)

      # Store some side-specific info into g & g.bwd

      g$left  <- side.info(g.left)
      g$right <- side.info(g.right)

      g.bwd$left  <- side.info(g.bwd.left)
      g.bwd$right <- side.info(g.bwd.right)


    # --------------------------------------------------------------------------
    tmp <- first.draft(g, g.right, g.bwd, g.bwd.right, reg)
      best <- tmp$d
      reg  <- tmp$reg


    # Revisit the first/last two riskbands if necessary

    # First two riskbands

    if (enabled$revisit.first2rb)
    {
      tmp <- revisit.closing2riskbands(g.left, g, best)
      if (tmp$better)  best <- tmp$d
    }


    # Last two riskbands

    if (enabled$revisit.last2rb)
    {
      tmp <- revisit.closing2riskbands(g.right, g, best)
      if (tmp$better)  best <- tmp$d
    }


    # Revisit the inner riskbands breakpoints & their corresponding heights
          
    if (enabled$revisit.innerrb)
    {
      tmp <- corrected.onePiece.innerRiskbands(best, g)
      if (tmp$corrected)  best <- tmp$d

      tmp <- revisit.inner.riskbands(g, best) # score is calculated herein for the last time
      if (tmp$better)  best <- tmp$d

      tmp <- corrected.onePiece.innerRiskbands(best, g) # Yes, we do it again!
      if (tmp$corrected)  best <- tmp$d


      # Also revisit riskbands 1+2 & 4+5, with h values fixed (loop over t's only to make curves smoother)

      best <- revisit.closing2riskbands.fixedh(g.left,  g, best)  # Riskbands 1+2
      best <- revisit.closing2riskbands.fixedh(g.right, g, best)  # Riskbands 4+5
    }



    tmp <- corrected.onePiece.innerRiskbands(best, g) # Yes, we do it once again!
    if (tmp$corrected)  best <- tmp$d


    # Recalculate score (to include the above-corrected slopes)

    best$score <- roughness(g, best)$score  # score update necessary after possible corrections to slopes in one-piece
                                            # riskbands, for example
    

    # Rearrange best's dimensions in sensible order
    
    if (!g$forward)
    {
      best$h <-  rev(best$h)
      best$s <- -rev(best$s)
          
      tmp      <-  best$R5b
      best$R5b <- -rev(best$R1b)
      best$R1b <- -rev(tmp)

      best$endrb.nSegments <- rev(best$endrb.nSegments)
      
      best$inner.b <- -rev(best$inner.b)
    }
            
    
    best <- wrap(best, g$forward)

    return(list(best=best, reg=reg))
  } # end of smooth


  truncated.s <- function(d)
  {
    s <- d$s
    if (d$endrb.nSegments[2] > 1)  s <- s[-length(s)]
    if (d$endrb.nSegments[1] > 1)  s <- s[-1]

    return(s)
  } # end of truncated.s
    

  #                                                            (.)~(.)
  #                                                           (-------)
  #  --------------------------------------------------------ooO-----Ooo--------
  #
  #                                   Beginning of Function
  #                                               CPL.prior
  #
  #  ---------------------------------------------------------------------------
  #                                                           ( )   ( )
  #                                                           /|\   /|\


  USE.EDUCATED.GUESS <- TRUE
  prior.elicitation.enabling.level <- 4  # can be set to {0, 1, 2, 3, 4}
                                         # for debugging the different layers
                                         # of the algorithm, but the algorithm will
                                         # be fully working as intended only when set to '4'


  enabled <- list(revisit.first2rb = prior.elicitation.enabling.level >= 1,
                  revisit.last2rb  = prior.elicitation.enabling.level >= 2,
                  revisit.innerrb  = prior.elicitation.enabling.level >= 3,
                  backward         = prior.elicitation.enabling.level >= 4)
      

  reg <- reg.init(precision)

  
  # Prepare geometry
  
  g <- geometry(riskband.prior.prob, lim, z, A, outcome.is.logNormally.distributed)

    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    # Verify geometric assumptions:
    # i) Five regions must be defined
    # ii) The two limit risbands/regions must be trapezoid-shaped.

    if (length(A) != 4) stop("The fct CPL.prior was developped to find the prior for\nfive (5) regions [hence 'A' must be of length 4].")

    NW.SE.corner.values <- lim$mu + z*rev(lim$sigma)
    tmp <- range(g$my.A)
    cond <- c(NW.SE.corner.values[1] < tmp[1], NW.SE.corner.values[2] > tmp[2])

    if (!all(cond))
    {
      if (outcome.is.logNormally.distributed)  NW.SE.corner.values <- exp(NW.SE.corner.values)

      cat('This function was developped to deal with 5 riskbands where the two end/limit riskbands (painted in figure) must be trapezoid-shaped.\n')
      cat('To do so, you must either\n')
      cat('\ta) change the cut-off values, making sure that\n')

      if (!cond[1]) cat('\t\tmin(A) > ', NW.SE.corner.values[1], '\n')
      if (!cond[2]) cat('\t\tmax(A) < ', NW.SE.corner.values[2], '\n')

      cat('or, maybe more naturally,\n\tb) change the (mu, sigma) lower & upper limits, making sure that\n')

      if (!cond[1]) cat('\t\tmu_lower + z * sigma_upper <', tmp[1], '\n')
      if (!cond[2]) cat('\t\tmu_upper + z * sigma_lower >', tmp[2], '\n')
      cat('\t\twhere z =', z, '\n')


      riskband.plot.regions(g)
      stop('\n')
    }

    # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    g.bwd <- geometry.backward(g)
  
    g     <- geometry.completeDefn(g)    
    g.bwd <- geometry.completeDefn(g.bwd)

    g$epsilon     <- precision
    g.bwd$epsilon <- precision

    g$use.educated.guess <- USE.EDUCATED.GUESS
    g.bwd$use.educated.guess <- g$use.educated.guess
                
    
  # Actually estimate the best density curve
  
  tmp <- smooth(g, g.bwd, enabled, reg)
    best <- tmp$best

  
  if (enabled$backward)
  {
    reg <- tmp$reg
    reg$forward <- FALSE
    reg$n   <- reg$backward$n
    reg$sum <- reg$backward$sum
    reg$backward <- NULL # drop backward object (useless after this point)


    bwd <- smooth(g.bwd, g, enabled, reg)$best


    if (bwd$score[1] == best$score[1] && bwd$score[2] == best$score[2] && abs(bwd$score[3] - best$score[3]) < 1e-6)
    {
      s.fwd <- truncated.s(best)
      s.bwd <- truncated.s(bwd)

      a.fwd <- angles(s.fwd)
      a.bwd <- angles(s.bwd)

      bwd.is.better <- sum(a.bwd^2) < sum(a.fwd^2)
    }
    else bwd.is.better <- better.score(bwd$score, best$score)

    if (bwd.is.better) best <- bwd
  }  
    

  # Wrap-up

  tmp <- simplified.onePiece.innerRiskbands(best, g)
  if (tmp$simplified)  best <- tmp$d

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
    legend.size <- legend('bottom', legend.txt, lty=legend.lty, bty='n', plot=FALSE)$rect # do not plot -- just note legend size
  
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
