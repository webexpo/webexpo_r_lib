##############################################################
#
#   WEBEXPO official R scripts
#
#   SEG ANALYSIS
#  
#   BAYESIAN MODELS IN JAGS 
#
#   RISKBAND FUNCTION - JAGS MODEL STATEMENTS and support functions written by
#   the McGill TEAM
#   
#   V1.0 Sept 2018 
#
#
#
#
#
############################################################# 

jags.model.riskband <-function(# Distribution
                                is.lognormal=TRUE, 
                                
                                # measurement error structure
                                error.type="none",
                                me.lower=0.3,
                                me.upper=0.3)
  
{


###### LOGNORMAL MODELS
  
if (is.lognormal)  
  
{

  ### measurement error as CV  

  if (error.type=="CV") model.jags <-paste("
                          
                                     
            mode{
               
            ##### MCGILL CUSTOM SECTION
          
            # Sample from marginal prior f(sigma)

             sigma.U ~ dunif(0, 1)
             
             for (j in 1:J)
             {
             delta[j] <- pow(q[j,2], 2) - 4*q[j,1]*(q[j,3] - sigma.U)
             sigma.tmp.soln[j] <- (-q[j,2] + sqrt(abs(delta[j]))) / (2 * q[j,1])
             # accept tmp soln if within range, otherwise score as 0
             sigma.soln[j] <- step(delta[j]) * step(sigma.tmp.soln[j] - sigma0[j]) * step(sigma1[j] - sigma.tmp.soln[j]) * sigma.tmp.soln[j]
             }
             
             sigma <- sum(sigma.soln[]) # sum of accepted solutions (there is only one)
             
             
             # Sample from conditional prior f(mu|sigma)
             
             mu.U ~ dunif(0, 1)
             
             # T_j(sigma) interval endpoints
             
             mu0[1] <- mu.range[1]
             mu1[1] <- min(max(A[1] - z*sigma, mu.range[1]), mu.range[2])
             
             for (r in 2:(R-1))
             {
             mu0[r] <- min(max(A[r-1] - z*sigma, mu.range[1]), mu.range[2]) # left-hand  side limit
             mu1[r] <- min(max(A[r]   - z*sigma, mu.range[1]), mu.range[2]) # right-hand side limit
             }
             
             mu0[R] <- min(max(A[R-1] - z*sigma, mu.range[1]), mu.range[2])
             mu1[R] <- mu.range[2]
             
             # T_j(sigma) interval lengths & weights
             
             for (r in 1:R)
             {
             mu.len[r] <- mu1[r] - mu0[r]
             mu.riskband.wt[r] <- mu.len[r] * f[r]
             }
             
             mu.totwt <- sum(mu.riskband.wt[])
             
             # cdf at regions endpoints
             
             mu.fcum[1] <- 0
             
             for (r in 2:R)
             {
             mu.fcum[r] <- mu.fcum[r-1] + mu.riskband.wt[r-1]
             }
             
             # solve inverse cdf
             
             for (r in 1:R)
             {
             mu.tmp.soln[r] <- mu0[r] + (mu.U * mu.totwt - mu.fcum[r]) / f[r]
             # accept tmp soln if within range, otherwise score as 0
             mu.soln[r] <- step(mu1[r] - mu.tmp.soln[r]) * step(mu.tmp.soln[r] - mu0[r]) * mu.tmp.soln[r]
             }
             
             mu <- sum(mu.soln[]) # sum of accepted solutions (there is only one)
                                         
                                         
            ##### WEBEXPO SECTION   
              
              ## prior on CV
                                      
              CV ~ dunif(",me.lower,",", me.upper,")

              # Likelihood
              
              for (i in 1:N) 
              {
             
              
               true.value[i] ~ dlnorm(mu, precision)
                                      
               x[i] ~ dnorm (true.value[i], MeasErrorPrecision[i])
                                         
               MeasErrorPrecision[i] <- pow(CV * true.value[i], -2)

                CensorType[i] ~ dinterval( x[i] , censorLimitMat[i,] )
              }
              
              precision <- pow(sigma, -2)
            }
                                      
                                      ",sep='')

  ### measurement error as SD
  
  if (error.type=="SD") model.jags <-paste("
                          
                                     
            model {

     ##### MCGILL CUSTOM SECTION

    # Sample from marginal prior f(sigma)

     sigma.U ~ dunif(0, 1)
     
     for (j in 1:J)
     {
     delta[j] <- pow(q[j,2], 2) - 4*q[j,1]*(q[j,3] - sigma.U)
     sigma.tmp.soln[j] <- (-q[j,2] + sqrt(abs(delta[j]))) / (2 * q[j,1])
     # accept tmp soln if within range, otherwise score as 0
     sigma.soln[j] <- step(delta[j]) * step(sigma.tmp.soln[j] - sigma0[j]) * step(sigma1[j] - sigma.tmp.soln[j]) * sigma.tmp.soln[j]
     }
     
     sigma <- sum(sigma.soln[]) # sum of accepted solutions (there is only one)
     
     
     # Sample from conditional prior f(mu|sigma)
     
     mu.U ~ dunif(0, 1)
     
     # T_j(sigma) interval endpoints
     
     mu0[1] <- mu.range[1]
     mu1[1] <- min(max(A[1] - z*sigma, mu.range[1]), mu.range[2])
     
     for (r in 2:(R-1))
     {
     mu0[r] <- min(max(A[r-1] - z*sigma, mu.range[1]), mu.range[2]) # left-hand  side limit
     mu1[r] <- min(max(A[r]   - z*sigma, mu.range[1]), mu.range[2]) # right-hand side limit
     }
     
     mu0[R] <- min(max(A[R-1] - z*sigma, mu.range[1]), mu.range[2])
     mu1[R] <- mu.range[2]
     
     # T_j(sigma) interval lengths & weights
     
     for (r in 1:R)
     {
     mu.len[r] <- mu1[r] - mu0[r]
     mu.riskband.wt[r] <- mu.len[r] * f[r]
     }
     
     mu.totwt <- sum(mu.riskband.wt[])
     
     # cdf at regions endpoints
     
     mu.fcum[1] <- 0
     
     for (r in 2:R)
     {
     mu.fcum[r] <- mu.fcum[r-1] + mu.riskband.wt[r-1]
     }
     
     # solve inverse cdf
     
     for (r in 1:R)
     {
     mu.tmp.soln[r] <- mu0[r] + (mu.U * mu.totwt - mu.fcum[r]) / f[r]
     # accept tmp soln if within range, otherwise score as 0
     mu.soln[r] <- step(mu1[r] - mu.tmp.soln[r]) * step(mu.tmp.soln[r] - mu0[r]) * mu.tmp.soln[r]
     }
     
     mu <- sum(mu.soln[]) # sum of accepted solutions (there is only one)
     
     
                ##### WEBEXPO SECTION                          
                                         
                     ## prior on SD
                     
                     SD ~ dunif(",me.lower,",", me.upper,")
                     
                     # Likelihood
                     
                     for (i in 1:N) 
                     {
                     
                     
                     true.value[i] ~ dlnorm(mu, precision)
                     
                     x[i] ~ dnorm (true.value[i], MeasErrorPrecision[i])
                     
                     MeasErrorPrecision[i] <- pow( SD , -2)

                          CensorType[i] ~ dinterval( x[i] , censorLimitMat[i,] )
                     }
                     
                     precision <- pow(sigma, -2)
}                                    
                                         ",sep='')



  

  ### no measurement error
  
  if (error.type=="none") model.jags <-paste("
                          
                                           
      model {

           ##### MCGILL CUSTOM SECTION

  # Sample from marginal prior f(sigma)

   sigma.U ~ dunif(0, 1)
   
   for (j in 1:J)
   {
   delta[j] <- pow(q[j,2], 2) - 4*q[j,1]*(q[j,3] - sigma.U)
   sigma.tmp.soln[j] <- (-q[j,2] + sqrt(abs(delta[j]))) / (2 * q[j,1])
   # accept tmp soln if within range, otherwise score as 0
   sigma.soln[j] <- step(delta[j]) * step(sigma.tmp.soln[j] - sigma0[j]) * step(sigma1[j] - sigma.tmp.soln[j]) * sigma.tmp.soln[j]
   }
   
   sigma <- sum(sigma.soln[]) # sum of accepted solutions (there is only one)
   
   
   # Sample from conditional prior f(mu|sigma)
   
   mu.U ~ dunif(0, 1)
   
   # T_j(sigma) interval endpoints
   
   mu0[1] <- mu.range[1]
   mu1[1] <- min(max(A[1] - z*sigma, mu.range[1]), mu.range[2])
   
   for (r in 2:(R-1))
   {
   mu0[r] <- min(max(A[r-1] - z*sigma, mu.range[1]), mu.range[2]) # left-hand  side limit
   mu1[r] <- min(max(A[r]   - z*sigma, mu.range[1]), mu.range[2]) # right-hand side limit
   }
   
   mu0[R] <- min(max(A[R-1] - z*sigma, mu.range[1]), mu.range[2])
   mu1[R] <- mu.range[2]
   
   # T_j(sigma) interval lengths & weights
   
   for (r in 1:R)
   {
   mu.len[r] <- mu1[r] - mu0[r]
   mu.riskband.wt[r] <- mu.len[r] * f[r]
   }
   
   mu.totwt <- sum(mu.riskband.wt[])
   
   # cdf at regions endpoints
   
   mu.fcum[1] <- 0
   
   for (r in 2:R)
   {
   mu.fcum[r] <- mu.fcum[r-1] + mu.riskband.wt[r-1]
   }
   
   # solve inverse cdf
   
   for (r in 1:R)
   {
   mu.tmp.soln[r] <- mu0[r] + (mu.U * mu.totwt - mu.fcum[r]) / f[r]
   # accept tmp soln if within range, otherwise score as 0
   mu.soln[r] <- step(mu1[r] - mu.tmp.soln[r]) * step(mu.tmp.soln[r] - mu0[r]) * mu.tmp.soln[r]
   }
   
   mu <- sum(mu.soln[]) # sum of accepted solutions (there is only one)
   
   
            ##### WEBEXPO SECTION                                

                   # Likelihood
                   
                   for (i in 1:N) 
                   {
                   
                   
                   x[i] ~ dlnorm(mu, precision)

                      CensorType[i] ~ dinterval( x[i] , censorLimitMat[i,] )
                   
                                                              }
                   
                   precision <- pow(sigma, -2)

                   }                                    
                   ",sep='')
  
  
}  
  
  
###### NORMAL MODELS

if (!is.lognormal)  
  
{
  
  
  ### measurement error as CV
  
  if (error.type=="CV") model.jags <-paste("
                                           
                                     
           model {

          ##### MCGILL CUSTOM SECTION

  # Sample from marginal prior f(sigma)

   sigma.U ~ dunif(0, 1)
   
   for (j in 1:J)
   {
   delta[j] <- pow(q[j,2], 2) - 4*q[j,1]*(q[j,3] - sigma.U)
   sigma.tmp.soln[j] <- (-q[j,2] + sqrt(abs(delta[j]))) / (2 * q[j,1])
   # accept tmp soln if within range, otherwise score as 0
   sigma.soln[j] <- step(delta[j]) * step(sigma.tmp.soln[j] - sigma0[j]) * step(sigma1[j] - sigma.tmp.soln[j]) * sigma.tmp.soln[j]
   }
   
   sigma <- sum(sigma.soln[]) # sum of accepted solutions (there is only one)
   
   
   # Sample from conditional prior f(mu|sigma)
   
   mu.U ~ dunif(0, 1)
   
   # T_j(sigma) interval endpoints
   
   mu0[1] <- mu.range[1]
   mu1[1] <- min(max(A[1] - z*sigma, mu.range[1]), mu.range[2])
   
   for (r in 2:(R-1))
   {
   mu0[r] <- min(max(A[r-1] - z*sigma, mu.range[1]), mu.range[2]) # left-hand  side limit
   mu1[r] <- min(max(A[r]   - z*sigma, mu.range[1]), mu.range[2]) # right-hand side limit
   }
   
   mu0[R] <- min(max(A[R-1] - z*sigma, mu.range[1]), mu.range[2])
   mu1[R] <- mu.range[2]
   
   # T_j(sigma) interval lengths & weights
   
   for (r in 1:R)
   {
   mu.len[r] <- mu1[r] - mu0[r]
   mu.riskband.wt[r] <- mu.len[r] * f[r]
   }
   
   mu.totwt <- sum(mu.riskband.wt[])
   
   # cdf at regions endpoints
   
   mu.fcum[1] <- 0
   
   for (r in 2:R)
   {
   mu.fcum[r] <- mu.fcum[r-1] + mu.riskband.wt[r-1]
   }
   
   # solve inverse cdf
   
   for (r in 1:R)
   {
   mu.tmp.soln[r] <- mu0[r] + (mu.U * mu.totwt - mu.fcum[r]) / f[r]
   # accept tmp soln if within range, otherwise score as 0
   mu.soln[r] <- step(mu1[r] - mu.tmp.soln[r]) * step(mu.tmp.soln[r] - mu0[r]) * mu.tmp.soln[r]
   }
   
   mu <- sum(mu.soln[]) # sum of accepted solutions (there is only one)
   

             
             ##### WEBEXPO SECTION

             ## prior on CV
             
             CV ~ dunif(",me.lower,",", me.upper,")
             
             # Likelihood
             
             for (i in 1:N) 
             {
             
             
             true.value[i] ~ dnorm(mu, precision)T(0,)
             
             x[i] ~ dnorm (true.value[i], MeasErrorPrecision[i])
             
             MeasErrorPrecision[i] <- pow(CV * true.value[i], -2)

                CensorType[i] ~ dinterval( x[i] , censorLimitMat[i,] )
             }
             
             precision <- pow(sigma, -2)
}                                           
                                           ",sep='')
  
  ### measurement error as SD
  
  if (error.type=="SD") model.jags <-paste("
                                           
               model {

                  ##### MCGILL CUSTOM SECTION

  # Sample from marginal prior f(sigma)

   sigma.U ~ dunif(0, 1)
   
   for (j in 1:J)
   {
   delta[j] <- pow(q[j,2], 2) - 4*q[j,1]*(q[j,3] - sigma.U)
   sigma.tmp.soln[j] <- (-q[j,2] + sqrt(abs(delta[j]))) / (2 * q[j,1])
   # accept tmp soln if within range, otherwise score as 0
   sigma.soln[j] <- step(delta[j]) * step(sigma.tmp.soln[j] - sigma0[j]) * step(sigma1[j] - sigma.tmp.soln[j]) * sigma.tmp.soln[j]
   }
   
   sigma <- sum(sigma.soln[]) # sum of accepted solutions (there is only one)
   
   
   # Sample from conditional prior f(mu|sigma)
   
   mu.U ~ dunif(0, 1)
   
   # T_j(sigma) interval endpoints
   
   mu0[1] <- mu.range[1]
   mu1[1] <- min(max(A[1] - z*sigma, mu.range[1]), mu.range[2])
   
   for (r in 2:(R-1))
   {
   mu0[r] <- min(max(A[r-1] - z*sigma, mu.range[1]), mu.range[2]) # left-hand  side limit
   mu1[r] <- min(max(A[r]   - z*sigma, mu.range[1]), mu.range[2]) # right-hand side limit
   }
   
   mu0[R] <- min(max(A[R-1] - z*sigma, mu.range[1]), mu.range[2])
   mu1[R] <- mu.range[2]
   
   # T_j(sigma) interval lengths & weights
   
   for (r in 1:R)
   {
   mu.len[r] <- mu1[r] - mu0[r]
   mu.riskband.wt[r] <- mu.len[r] * f[r]
   }
   
   mu.totwt <- sum(mu.riskband.wt[])
   
   # cdf at regions endpoints
   
   mu.fcum[1] <- 0
   
   for (r in 2:R)
   {
   mu.fcum[r] <- mu.fcum[r-1] + mu.riskband.wt[r-1]
   }
   
   # solve inverse cdf
   
   for (r in 1:R)
   {
   mu.tmp.soln[r] <- mu0[r] + (mu.U * mu.totwt - mu.fcum[r]) / f[r]
   # accept tmp soln if within range, otherwise score as 0
   mu.soln[r] <- step(mu1[r] - mu.tmp.soln[r]) * step(mu.tmp.soln[r] - mu0[r]) * mu.tmp.soln[r]
   }
   
   mu <- sum(mu.soln[]) # sum of accepted solutions (there is only one)
   

                  ##### WEBEXPO SECTION
                   
                   ## prior on CV
                   
                   SD ~ dunif(",me.lower,",", me.upper,")
                   
                   # Likelihood
                   
                   for (i in 1:N) 
                   {
                   
                   
                   true.value[i] ~ dnorm(mu, precision)
                   
                   x[i] ~ dnorm (true.value[i], MeasErrorPrecision[i])
                   
                   MeasErrorPrecision[i] <- pow( SD , -2)

                   CensorType[i] ~ dinterval( x[i] , censorLimitMat[i,] )
                   }
                   
                   precision <- pow(sigma, -2)
}                                       
                   ",sep='')
  
  ### no measurement error
  
  if (error.type=="none") model.jags <-paste("
                                             
                model  {

                  ##### MCGILL CUSTOM SECTION

              # Sample from marginal prior f(sigma)

                 sigma.U ~ dunif(0, 1)
                 
                 for (j in 1:J)
                 {
                 delta[j] <- pow(q[j,2], 2) - 4*q[j,1]*(q[j,3] - sigma.U)
                 sigma.tmp.soln[j] <- (-q[j,2] + sqrt(abs(delta[j]))) / (2 * q[j,1])
                 # accept tmp soln if within range, otherwise score as 0
                 sigma.soln[j] <- step(delta[j]) * step(sigma.tmp.soln[j] - sigma0[j]) * step(sigma1[j] - sigma.tmp.soln[j]) * sigma.tmp.soln[j]
                 }
                 
                 sigma <- sum(sigma.soln[]) # sum of accepted solutions (there is only one)
                 
                 
                 # Sample from conditional prior f(mu|sigma)
                 
                 mu.U ~ dunif(0, 1)
                 
                 # T_j(sigma) interval endpoints
                 
                 mu0[1] <- mu.range[1]
                 mu1[1] <- min(max(A[1] - z*sigma, mu.range[1]), mu.range[2])
                 
                 for (r in 2:(R-1))
                 {
                 mu0[r] <- min(max(A[r-1] - z*sigma, mu.range[1]), mu.range[2]) # left-hand  side limit
                 mu1[r] <- min(max(A[r]   - z*sigma, mu.range[1]), mu.range[2]) # right-hand side limit
                 }
                 
                 mu0[R] <- min(max(A[R-1] - z*sigma, mu.range[1]), mu.range[2])
                 mu1[R] <- mu.range[2]
                 
                 # T_j(sigma) interval lengths & weights
                 
                 for (r in 1:R)
                 {
                 mu.len[r] <- mu1[r] - mu0[r]
                 mu.riskband.wt[r] <- mu.len[r] * f[r]
                 }
                 
                 mu.totwt <- sum(mu.riskband.wt[])
                 
                 # cdf at regions endpoints
                 
                 mu.fcum[1] <- 0
                 
                 for (r in 2:R)
                 {
                 mu.fcum[r] <- mu.fcum[r-1] + mu.riskband.wt[r-1]
                 }
                 
                 # solve inverse cdf
                 
                 for (r in 1:R)
                 {
                 mu.tmp.soln[r] <- mu0[r] + (mu.U * mu.totwt - mu.fcum[r]) / f[r]
                 # accept tmp soln if within range, otherwise score as 0
                 mu.soln[r] <- step(mu1[r] - mu.tmp.soln[r]) * step(mu.tmp.soln[r] - mu0[r]) * mu.tmp.soln[r]
                 }
                 
                 mu <- sum(mu.soln[]) # sum of accepted solutions (there is only one)
                 

                  ##### WEBEXPO SECTION

                   # Likelihood
                   
                   for (i in 1:N) 
                   {
                   
                   
                   x[i] ~ dnorm(mu, precision)

                    CensorType[i] ~ dinterval( x[i] , censorLimitMat[i,] )
                   
                   }
                   
                   precision <- pow(sigma, -2)
}                                             
                   ",sep='')
  
  
}


return(model.jags)
     
}




######################## support function for the riskband JAGS models

###################   Written by Patrick Belisle July 2018
riskband.rjags.parms <- function(m.range, s.range, # See note below
                                 A,
                                 region.prior.prob=rep(1/(length(A)+1), length(A)+1), 
                                 level=0.95, outcome.is.logNormally.distributed=T)
{
  # Version 0.3
  # (5 september 2018)
  
  
  # IMPORTANT:
  # ----------
  #
  # - if outcome.is.logNormally.distributed = TRUE,
  #     then arguments m.range and s.range MUST represent gm.range and gsd.range
  #
  # - if outcome.is.logNormally.distributed = FALSE,
  #     then arguments m.range and s.range MUST represent mu.range and sd.range
  
  
  z <- qnorm(level)
  
  if (outcome.is.logNormally.distributed)
  { 
    A <- log(A)
    
    gm.range <- m.range
    gsd.range <- s.range
    
    mu.range <- log(gm.range)
    sigma.range <- log(gsd.range)
  }
  else
  {
    mu.range <- m.range
    sigma.range <- s.range
  }
  
  
  pr <- region.prior.prob / sum(region.prior.prob) 
  
  
  # --- Define useful functions ------------------------------------------------
  
  in.range <- function(x, r)
  {
    out.of.range <- x < min(r) | x > max(r)
    x[out.of.range] <- NA
    x
  } # end of in.range
  
  
  lengths <- function(l, mu.range, diff.A, R, r1)
  {
    one <- matrix(1, nrow=1, ncol=R)
    mu.width <- diff(mu.range)
    
    len <- matrix(NA, nrow=R, ncol=2)
    len[r1,] <- l
    
    
    len[is.na(len)] <- 0
    
    seq.from <- r1 + 1
    seq.to <- R-1
    
    if (seq.to >= seq.from)
    {
      
      r.seq <- seq(from=seq.from, to=seq.to)
      
      
      for (r in r.seq)
      {
        len[r,] <- rep(diff(A)[r-1], 2)
        col.sums <- as.vector(one %*% len)
        
        W <- which(col.sums > mu.width)
        for (w in W)
        {
          col.sum <- sum(len[seq(r-1), w])
          len[r,w] <- mu.width - col.sum
        }
      }
    }
    
    col.sums <- as.vector(one %*% len)
    w <- which(col.sums < mu.width)
    len[R,w] <- mu.width - col.sums[w]
    
    len
  } # end of lengths
  
  
  sigma.regions <- function(mu.range, sigma.range, y.left, y.right, x.top, z)
  {
    sigma.limits <- sort(c(y.left[!is.na(y.left)], sigma.range))
    y.right <- y.right[!is.na(y.right)]
    
    df <- data.frame(from=sigma.limits[-length(sigma.limits)])
    df$to <- sigma.limits[-1]
    df$len <- df$to - df$from
    
    df$mu.len.min <- 0
    J <- nrow(df)
    if (J > 1)
    { 
      if (any(!is.na(x.top))) df$mu.len.min[J] <- min(x.top, na.rm=T) - min(mu.range)
      else df$mu.len.min[J] <- diff(mu.range)
    }
    
    df$mu.len.max <- df$mu.len.min + z*df$len
    df$mu.len.max <- pmin(df$mu.len.max, diff(mu.range))
    df$r1 <- seq(J)
    
    # Break lines where y.right values are found in (from, to)
    if (length(y.right) > 0)
    {
      for (yr in y.right)
      {
        w <- which(yr > df$from & yr < df$to)
        
        if (length(w) > 0)
        {
          tmp <- df[w,]
          df <- df[-w,]
          tmp.new <- data.frame(from=c(tmp$from, yr), to=c(yr, tmp$to))
          tmp.new$len <- tmp.new$to - tmp.new$from
          new.mu <- tmp$mu.len.max - z*(yr-tmp$from)
          tmp.new$mu.len.min <- c(new.mu, tmp$mu.len.min)
          tmp.new$mu.len.max <- c(tmp$mu.len.max, new.mu)
          tmp.new$r1 <- tmp$r1
          
          df <- rbind.data.frame(df, tmp.new)
        }
      }
    }
    
    o <- order(df$from)
    df <- df[o,]
    
    df
  } # end of sigma.regions
  
  
  wts <- function(mu.range, sigma.range, A, pr, z)
  {
    leftZone.area <- function(xbottom, yleft, mu.range, sigma.range, A.r, z, left=T)
    {
      if (is.na(xbottom))
      {
        if (is.na(yleft))
        {        
          if (is.na(A.r))
          {
            S <- diff(mu.range) * diff(sigma.range)
          }
          else
          {
            xtop <- A.r - z*max(sigma.range)
            yright <- (A.r - max(mu.range)) / z
            S <- diff(mu.range) * diff(sigma.range) - (max(mu.range) - xtop) * (max(sigma.range) - yright) / 2
          }
        }
        else
        {
          yright <- (A.r - max(mu.range)) / z
          S <- ((yleft + yright) / 2 - min(sigma.range)) * diff(mu.range)
        }
      }
      else
      {
        if (is.na(yleft))
        {
          xtop <- A.r - z*max(sigma.range)
          S <- ((xtop + xbottom) / 2 - min(mu.range)) * diff(sigma.range)
        }
        else
        {
          S <- (xbottom - min(mu.range)) * (yleft-min(sigma.range)) / 2
        } 
      } 
      
      if (!left) S <- diff(mu.range) * diff(sigma.range) - S
      
      S
    } # end of leftZone.area
    
    
    # Begining of fct wts
    
    # Where the lines intersect with rectangle's bottom & left borders
    x.bottom <- A - z*min(sigma.range)
    y.left <- (A - min(mu.range)) / z
    
    within.range <- x.bottom > min(mu.range) & x.bottom < max(mu.range)
    x.bottom[!within.range] <- NA
    
    within.range <- y.left > min(sigma.range) & y.left < max(sigma.range)
    y.left[!within.range] <- NA
    
    mu.len <- diff(mu.range)
    sigma.len <- diff(sigma.range)
    
    R <- length(A) + 1 # Number of regions
    S <- rep(NA, R)
    
    S[1] <- leftZone.area(x.bottom[1], y.left[1], mu.range, sigma.range, A[1], z)
    S[R] <- leftZone.area(x.bottom[R-1], y.left[R-1], mu.range, sigma.range, A[R-1], z, left=F)  
    
    if (R > 2)
    {
      left.area <- rep(NA, R-1)
      left.area[1] <- S[1] # already computed above
      
      for (r in seq(from=2, to=R-1))
      {
        left.area[r] <- leftZone.area(x.bottom[r], y.left[r], mu.range, sigma.range, A[r], z) 
      }      
      
      S[2:(R-1)] <- diff(left.area)
    }
    
    f <- pr/S * S[1]/pr[1]
    f <- f/sum(f*S)
    
    list(S=S, f=f)
  } # end of wts
  
  
  # --- Beginning of function --------------------------------------------------
  
  
  w <- wts(mu.range, sigma.range, A, pr, z)
  R <- length(w$f)
  
  # Compute f(sigma)  [marginal prior]
  
  y.left <- (A - min(mu.range)) / z
  y.left <- in.range(y.left, sigma.range)
  
  y.right <- (A - max(mu.range)) / z
  y.right <- in.range(y.right, sigma.range)
  
  x.bottom <- A - z*min(sigma.range)
  x.top <- A - z*max(sigma.range)
  x.bottom <- in.range(x.bottom, mu.range)
  x.top <- in.range(x.top, mu.range)
  
  
  # --- Safety check ---
  #
  # Verify that values in A all lead to lines that DO cross the area
  
  boundaries.cutoffs <- matrix(c(y.left, y.right, x.bottom, x.top), ncol=4)
  any.cutoff <- apply(!is.na(boundaries.cutoffs), 1, any)
  
  if (any(!any.cutoff))
  {
    w <- which(!any.cutoff)
    msg <- paste("Boundary defined by A[", w, "] does not cross the domain m.range X s.range. Sorry.", sep='')
    stop(msg)
  }
  
  
  
  df <- sigma.regions(mu.range, sigma.range, y.left, y.right, x.top, z)
  
  df$f0 <- NA
  df$f1 <- NA
  df$fp <- NA
  
  for (r in seq(nrow(df)))
  {
    tmp <- df[r,]
    s <- c(tmp$from, tmp$to)
    l1 <- c(tmp$mu.len.max, tmp$mu.len.min)
    len <- lengths(l1, mu.range, diff(A), R, tmp$r1)
    f <- as.vector(matrix(w$f, nrow=1) %*% len)
    df$f0[r] <- f[1]
    df$f1[r] <- f[2]
    df$fp[r] <- diff(f)/diff(s)
  }
  
  df$ftot <- (df$f0 + df$f1) / 2 * df$len
  fcum <- c(0, cumsum(df$ftot))
  fcum <- fcum[-length(fcum)]
  
  
  df.q <- data.frame(A=df$fp/2, B=df$f0-df$fp*df$from, C=df$fp*df$from^2/2-df$f0*df$from)
  df.q$C <- df.q$C + fcum
  
  
  list(q=df.q, mu.range=mu.range, sigma0=df$from, sigma1=df$to,
       A=A, z=z, f=w$f, R=R, J=nrow(df))
} # end of riskband.rjags.parms