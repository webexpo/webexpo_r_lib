

riskband.area <- function(region.prior.prob, lim, z, A, outcome.is.logNormally.distributed=T)

  {
  
  max.deepness=1/3
  
  if (outcome.is.logNormally.distributed) A <- log(A)
  
  theta <- atan(z)
  cos.theta <- cos(theta)
  R <- length(A) + 1
  
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
  mu.top <- c(mu0, A - z*sigma1, mu1 - Bt)
  
  B <- diff(mu.top)
  b <- B*cos.theta
  
  area <- B*H
  w <- c(1,5)
  area[w] <- area[w] + Bt * H / 2
  
  u <- region.prior.prob / area
  p <- region.prior.prob / L
  hbar <- p / b  # meaningless for the left-end and right-end regions
  
  
  # regression coefficients that permit to calculate s1 from s2 in paralleliped regions
  
  s1.intercept <- 8/b * hbar
  s1.h0 <- -8/b
  
  
  # Calculate the floor level to respect max deepness
  # and the highest height (on either side of the riskband)
  # to be able to respect max deepness
  
  lo <- u * max.deepness    # max deepness = 1/3 of the piecewise-uniform level
  hi <- 4 * hbar - 3 * lo
  hi[c(1,5)] <- Inf # meaningless for the left-end and right-end regions
  
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
  
  # inner-region lower points 
  #   (for the calculation of density curves appropriateness, 
  #   adding a cost to curves going too low when compared to piecewise-uniform prior inner-regions)
  
  u.low <- pmin(u[-R], u[-1])
  u.low <- c(u.low[1], u[2], u.low[2], 
             u[3], u.low[3],
             u[4], u.low[4])
  
  # Region bounds
  
  x95b <- c(mu.lim[1] + z*sigma.lim[1], A, mu.lim[2] + z*sigma.lim[2])
  
  
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
  
  
  endregion.is.local.min <- rep(NA, 2)
  uniq.u <- u[!duplicated(u)]
  endregion.is.local.min[1] <- uniq.u[1] < uniq.u[2]
  
  uniq.u <- rev(u)
  uniq.u <- uniq.u[!duplicated(uniq.u)]
  endregion.is.local.min[2] <- uniq.u[1] < uniq.u[2]
  
  j <- c(1,5)
  s.fade.out <- u[j] * max.deepness / (d + b[j]) * c(1, -1)
  
  
  list(area=area)
  
} # end of geometry
