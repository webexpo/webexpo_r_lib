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
#   V1.0 Sept2018
#   V1.1 Jan 2020 # implementation by JL of the latest McGill modifications 
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
            
              for (j in Quadratic.sigma.cdf.segments)
              {
                delta[j] <- pow(q[j,2], 2) - 4*q[j,1]*(q[j,3] - sigma.U)
                sigma.tmp.soln[j] <- (-q[j,2] + sqrt(abs(delta[j]))) / (2 * q[j,1])
                # accept tmp soln if within range, otherwise score as 0
                sigma.soln[j] <- step(delta[j]) * step(sigma.tmp.soln[j] - sigma0[j]) * step(sigma1[j] - sigma.tmp.soln[j]) * sigma.tmp.soln[j]
              }
              
              for (j in Linear.sigma.cdf.segments)
              {
                sigma.tmp.soln[j] <- (sigma.U - q[j,3]) / q[j,2]
                sigma.soln[j] <- step(sigma.tmp.soln[j] - sigma0[j]) * step(sigma1[j] - sigma.tmp.soln[j]) * sigma.tmp.soln[j]
              }
              
              
              sigma <- sum(sigma.soln[]) # sum of accepted solutions (there is only one)
              
              # Sample from conditional prior f(mu|sigma)
              
              mu.U ~ dunif(0, 1)
              
                # T_j(sigma) interval endpoints
              
                mu0[1] <- mu.range[1]
                mu1[1] <- min(max(B[1] - z*sigma, mu.range[1]), mu.range[2])
              
                for (r in 2:(R-1))
                {
                  mu0[r] <- min(max(B[r-1] - z*sigma, mu.range[1]), mu.range[2]) # left-hand  side limit
                  mu1[r] <- min(max(B[r]   - z*sigma, mu.range[1]), mu.range[2]) # right-hand side limit
                }
            
                mu0[R] <- min(max(B[R-1] - z*sigma, mu.range[1]), mu.range[2])
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
    
      for (j in Quadratic.sigma.cdf.segments)
      {
        delta[j] <- pow(q[j,2], 2) - 4*q[j,1]*(q[j,3] - sigma.U)
        sigma.tmp.soln[j] <- (-q[j,2] + sqrt(abs(delta[j]))) / (2 * q[j,1])
        # accept tmp soln if within range, otherwise score as 0
        sigma.soln[j] <- step(delta[j]) * step(sigma.tmp.soln[j] - sigma0[j]) * step(sigma1[j] - sigma.tmp.soln[j]) * sigma.tmp.soln[j]
      }
      
      for (j in Linear.sigma.cdf.segments)
      {
        sigma.tmp.soln[j] <- (sigma.U - q[j,3]) / q[j,2]
        sigma.soln[j] <- step(sigma.tmp.soln[j] - sigma0[j]) * step(sigma1[j] - sigma.tmp.soln[j]) * sigma.tmp.soln[j]
      }
      
      
      sigma <- sum(sigma.soln[]) # sum of accepted solutions (there is only one)
      
      # Sample from conditional prior f(mu|sigma)
      
      mu.U ~ dunif(0, 1)
      
        # T_j(sigma) interval endpoints
      
        mu0[1] <- mu.range[1]
        mu1[1] <- min(max(B[1] - z*sigma, mu.range[1]), mu.range[2])
      
        for (r in 2:(R-1))
        {
          mu0[r] <- min(max(B[r-1] - z*sigma, mu.range[1]), mu.range[2]) # left-hand  side limit
          mu1[r] <- min(max(B[r]   - z*sigma, mu.range[1]), mu.range[2]) # right-hand side limit
        }
    
        mu0[R] <- min(max(B[R-1] - z*sigma, mu.range[1]), mu.range[2])
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
        
          for (j in Quadratic.sigma.cdf.segments)
          {
            delta[j] <- pow(q[j,2], 2) - 4*q[j,1]*(q[j,3] - sigma.U)
            sigma.tmp.soln[j] <- (-q[j,2] + sqrt(abs(delta[j]))) / (2 * q[j,1])
            # accept tmp soln if within range, otherwise score as 0
            sigma.soln[j] <- step(delta[j]) * step(sigma.tmp.soln[j] - sigma0[j]) * step(sigma1[j] - sigma.tmp.soln[j]) * sigma.tmp.soln[j]
          }
          
          for (j in Linear.sigma.cdf.segments)
          {
            sigma.tmp.soln[j] <- (sigma.U - q[j,3]) / q[j,2]
            sigma.soln[j] <- step(sigma.tmp.soln[j] - sigma0[j]) * step(sigma1[j] - sigma.tmp.soln[j]) * sigma.tmp.soln[j]
          }
          
          
          sigma <- sum(sigma.soln[]) # sum of accepted solutions (there is only one)
          
          # Sample from conditional prior f(mu|sigma)
          
          mu.U ~ dunif(0, 1)
          
            # T_j(sigma) interval endpoints
          
            mu0[1] <- mu.range[1]
            mu1[1] <- min(max(B[1] - z*sigma, mu.range[1]), mu.range[2])
          
            for (r in 2:(R-1))
            {
              mu0[r] <- min(max(B[r-1] - z*sigma, mu.range[1]), mu.range[2]) # left-hand  side limit
              mu1[r] <- min(max(B[r]   - z*sigma, mu.range[1]), mu.range[2]) # right-hand side limit
            }
        
            mu0[R] <- min(max(B[R-1] - z*sigma, mu.range[1]), mu.range[2])
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
        
          for (j in Quadratic.sigma.cdf.segments)
          {
            delta[j] <- pow(q[j,2], 2) - 4*q[j,1]*(q[j,3] - sigma.U)
            sigma.tmp.soln[j] <- (-q[j,2] + sqrt(abs(delta[j]))) / (2 * q[j,1])
            # accept tmp soln if within range, otherwise score as 0
            sigma.soln[j] <- step(delta[j]) * step(sigma.tmp.soln[j] - sigma0[j]) * step(sigma1[j] - sigma.tmp.soln[j]) * sigma.tmp.soln[j]
          }
          
          for (j in Linear.sigma.cdf.segments)
          {
            sigma.tmp.soln[j] <- (sigma.U - q[j,3]) / q[j,2]
            sigma.soln[j] <- step(sigma.tmp.soln[j] - sigma0[j]) * step(sigma1[j] - sigma.tmp.soln[j]) * sigma.tmp.soln[j]
          }
          
          
          sigma <- sum(sigma.soln[]) # sum of accepted solutions (there is only one)
          
          # Sample from conditional prior f(mu|sigma)
          
          mu.U ~ dunif(0, 1)
          
            # T_j(sigma) interval endpoints
          
            mu0[1] <- mu.range[1]
            mu1[1] <- min(max(B[1] - z*sigma, mu.range[1]), mu.range[2])
          
            for (r in 2:(R-1))
            {
              mu0[r] <- min(max(B[r-1] - z*sigma, mu.range[1]), mu.range[2]) # left-hand  side limit
              mu1[r] <- min(max(B[r]   - z*sigma, mu.range[1]), mu.range[2]) # right-hand side limit
            }
        
            mu0[R] <- min(max(B[R-1] - z*sigma, mu.range[1]), mu.range[2])
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
              
                for (j in Quadratic.sigma.cdf.segments)
                {
                  delta[j] <- pow(q[j,2], 2) - 4*q[j,1]*(q[j,3] - sigma.U)
                  sigma.tmp.soln[j] <- (-q[j,2] + sqrt(abs(delta[j]))) / (2 * q[j,1])
                  # accept tmp soln if within range, otherwise score as 0
                  sigma.soln[j] <- step(delta[j]) * step(sigma.tmp.soln[j] - sigma0[j]) * step(sigma1[j] - sigma.tmp.soln[j]) * sigma.tmp.soln[j]
                }
                
                for (j in Linear.sigma.cdf.segments)
                {
                  sigma.tmp.soln[j] <- (sigma.U - q[j,3]) / q[j,2]
                  sigma.soln[j] <- step(sigma.tmp.soln[j] - sigma0[j]) * step(sigma1[j] - sigma.tmp.soln[j]) * sigma.tmp.soln[j]
                }
                
                
                sigma <- sum(sigma.soln[]) # sum of accepted solutions (there is only one)
                
                # Sample from conditional prior f(mu|sigma)
                
                mu.U ~ dunif(0, 1)
                
                  # T_j(sigma) interval endpoints
                
                  mu0[1] <- mu.range[1]
                  mu1[1] <- min(max(B[1] - z*sigma, mu.range[1]), mu.range[2])
                
                  for (r in 2:(R-1))
                  {
                    mu0[r] <- min(max(B[r-1] - z*sigma, mu.range[1]), mu.range[2]) # left-hand  side limit
                    mu1[r] <- min(max(B[r]   - z*sigma, mu.range[1]), mu.range[2]) # right-hand side limit
                  }
              
                  mu0[R] <- min(max(B[R-1] - z*sigma, mu.range[1]), mu.range[2])
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
              
                for (j in Quadratic.sigma.cdf.segments)
                {
                  delta[j] <- pow(q[j,2], 2) - 4*q[j,1]*(q[j,3] - sigma.U)
                  sigma.tmp.soln[j] <- (-q[j,2] + sqrt(abs(delta[j]))) / (2 * q[j,1])
                  # accept tmp soln if within range, otherwise score as 0
                  sigma.soln[j] <- step(delta[j]) * step(sigma.tmp.soln[j] - sigma0[j]) * step(sigma1[j] - sigma.tmp.soln[j]) * sigma.tmp.soln[j]
                }
                
                for (j in Linear.sigma.cdf.segments)
                {
                  sigma.tmp.soln[j] <- (sigma.U - q[j,3]) / q[j,2]
                  sigma.soln[j] <- step(sigma.tmp.soln[j] - sigma0[j]) * step(sigma1[j] - sigma.tmp.soln[j]) * sigma.tmp.soln[j]
                }
                
                
                sigma <- sum(sigma.soln[]) # sum of accepted solutions (there is only one)
                
                # Sample from conditional prior f(mu|sigma)
                
                mu.U ~ dunif(0, 1)
                
                  # T_j(sigma) interval endpoints
                
                  mu0[1] <- mu.range[1]
                  mu1[1] <- min(max(B[1] - z*sigma, mu.range[1]), mu.range[2])
                
                  for (r in 2:(R-1))
                  {
                    mu0[r] <- min(max(B[r-1] - z*sigma, mu.range[1]), mu.range[2]) # left-hand  side limit
                    mu1[r] <- min(max(B[r]   - z*sigma, mu.range[1]), mu.range[2]) # right-hand side limit
                  }
              
                  mu0[R] <- min(max(B[R-1] - z*sigma, mu.range[1]), mu.range[2])
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
  # Version 0.6
  # (30 september 2018)
  
  
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
    B <- log(A)
    
    gm.range <- m.range
    gsd.range <- s.range
    
    mu.range <- log(gm.range)
    sigma.range <- log(gsd.range)
  }
  else
  {
    B <- A
    
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
  
  
  sigma.regions <- function(mu.range, sigma.range, B, f, z)
  {
    f.Sigma <- function(sigma, mu.range, B, r0, r1, f, z)
    {
      if (r0 == r1)
      {
        f <- f[r0] * diff(mu.range)
        linear <- T
      }
      else
      {
        r <- seq(from=r0, to=r1)
        rB <- seq(from=r0, to=r1-1)
        mu.breaks <- B[rB] - z*sigma
        f <- sum(f[r] * diff(c(mu.range[1], mu.breaks, mu.range[2])))
        linear <- F
      }
      
      list(f=f, linear=linear)
    } # end of f.Sigma
    
    
    tmp <- data.frame(j=seq(along=B), left=(B - min(mu.range)) / z, right=(B - max(mu.range)) / z)
    
    r0 <- 1
    r1 <- sum(tmp$right < sigma.range[1]) + 1
    
    # Lower sigma value
    
    sigma.breaks <- sigma.range[1]
    g <- f.Sigma(sigma.breaks, mu.range, B, r0, r1, f, z)
    f.sigma <- g$f
    linear <- g$linear
    
    # Intermediate sigma values
    
    tmp0 <- data.frame(j=tmp$j, sigma=tmp$left, left.side=T)
    tmp1 <- data.frame(j=tmp$j, sigma=tmp$right, left.side=F)
    tmp <- rbind.data.frame(tmp0, tmp1)
    tmp <- subset(tmp, sigma > sigma.range[1] & sigma < sigma.range[2])
    o <- order(tmp$sigma)
    tmp <- tmp[o,]
    
    while (nrow(tmp) > 0)
    {
      sigma <- tmp$sigma[1]
      sigma.breaks <- c(sigma.breaks, sigma)
      w <- which(tmp$sigma == sigma)
      
      if (any(tmp$left.side[w])) r0 <- r0 + 1
      
      g <- f.Sigma(sigma, mu.range, B, r0, r1, f, z)
      f.sigma <- c(f.sigma, g$f)
      linear <- c(linear, g$linear)
      
      if (any(!tmp$left.side[w])) r1 <- r1 + 1
      tmp <- tmp[-w,]
    }
    
    # Higher sigma value
    
    g <- f.Sigma(sigma.range[2], mu.range, B, r0, r1, f, z)
    sigma.breaks <- c(sigma.breaks, sigma.range[2])
    f.sigma <- c(f.sigma, g$f)
    
    
    fp <- diff(f.sigma) / diff(sigma.breaks)
    ftot <- (f.sigma[-length(f.sigma)] + f.sigma[-1]) / 2 * diff(sigma.breaks)
    fcum <- cumsum(ftot)
    fcum.max <- max(fcum)
    fp <- fp / fcum.max
    fcum <- fcum / fcum.max
    fcum <- c(0, fcum)
    
    
    qA <- fp/2
    qA[which(linear)] <- 0
    qB <- rep(NA, length(qA))
    qC <- qB
    
    
    for (i in seq(length(qA)))
    {
      u <- matrix(c(sigma.breaks[i+c(0,1)], 1, 1), ncol=2)
      v <- matrix(fcum[i+c(0,1)] - qA[i]*(sigma.breaks[i+c(0,1)]^2), ncol=1)
      x <- as.vector(solve(u) %*% v)
      qB[i] <- x[1]
      qC[i] <- x[2]
    }
    
    
    out <- list(q=matrix(c(qA, qB, qC), ncol=3),
                from=sigma.breaks[-length(sigma.breaks)],
                to=sigma.breaks[-1])
    
    
    if (any(linear))
    {
      Linear.sigma.cdf.segments <- which(linear)
      Quadratic.sigma.cdf.segments <- which(!linear)
      out$Quadratic.sigma.cdf.segments <- Quadratic.sigma.cdf.segments
      out$Linear.sigma.cdf.segments <- Linear.sigma.cdf.segments
    }
    else
    {
      out$J <- length(qA)
    }
    
    out
  } # end of sigma.regions
  
  
  wts <- function(mu.range, sigma.range, B, pr, z)
  {
    leftZone.area <- function(xbottom, yleft, mu.range, sigma.range, B.r, z, left=T)
    {
      if (is.na(xbottom))
      {
        if (is.na(yleft))
        {        
          if (is.na(B.r))
          {
            S <- diff(mu.range) * diff(sigma.range)
          }
          else
          {
            xtop <- B.r - z*max(sigma.range)
            yright <- (B.r - max(mu.range)) / z
            S <- diff(mu.range) * diff(sigma.range) - (max(mu.range) - xtop) * (max(sigma.range) - yright) / 2
          }
        }
        else
        {
          yright <- (B.r - max(mu.range)) / z
          S <- ((yleft + yright) / 2 - min(sigma.range)) * diff(mu.range)
        }
      }
      else
      {
        if (is.na(yleft))
        {
          xtop <- B.r - z*max(sigma.range)
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
    x.bottom <- B - z*min(sigma.range)
    y.left <- (B - min(mu.range)) / z
    
    within.range <- x.bottom > min(mu.range) & x.bottom < max(mu.range)
    x.bottom[!within.range] <- NA
    
    within.range <- y.left > min(sigma.range) & y.left < max(sigma.range)
    y.left[!within.range] <- NA
    
    mu.len <- diff(mu.range)
    sigma.len <- diff(sigma.range)
    
    R <- length(B) + 1 # Number of regions
    S <- rep(NA, R)
    
    S[1] <- leftZone.area(x.bottom[1], y.left[1], mu.range, sigma.range, B[1], z)
    S[R] <- leftZone.area(x.bottom[R-1], y.left[R-1], mu.range, sigma.range, B[R-1], z, left=F)  
    
    if (R > 2)
    {
      left.area <- rep(NA, R-1)
      left.area[1] <- S[1] # already computed above
      
      for (r in seq(from=2, to=R-1))
      {
        left.area[r] <- leftZone.area(x.bottom[r], y.left[r], mu.range, sigma.range, B[r], z) 
      }      
      
      S[2:(R-1)] <- diff(left.area)
    }
    
    f <- pr/S * S[1]/pr[1]
    f <- f/sum(f*S)
    
    list(S=S, f=f)
  } # end of wts
  
  
  # --- Beginning of function --------------------------------------------------
  
  
  w <- wts(mu.range, sigma.range, B, pr, z)
  R <- length(w$f)
  
  # Compute f(sigma)  [marginal prior]
  
  y.left <- (B - min(mu.range)) / z
  y.left <- in.range(y.left, sigma.range)
  
  y.right <- (B - max(mu.range)) / z
  y.right <- in.range(y.right, sigma.range)
  
  x.bottom <- B - z*min(sigma.range)
  x.top <- B - z*max(sigma.range)
  x.bottom <- in.range(x.bottom, mu.range)
  x.top <- in.range(x.top, mu.range)
  
  
  # --- Safety check ---
  #
  # Verify that values in B all lead to lines that DO cross the area
  
  boundaries.cutoffs <- matrix(c(y.left, y.right, x.bottom, x.top), ncol=4)
  any.cutoff <- apply(!is.na(boundaries.cutoffs), 1, any)
  
  if (any(!any.cutoff))
  {
    w <- which(!any.cutoff)
    msg <- paste("Boundary defined by A[", w, "] does not cross the domain m.range X s.range. Sorry.", sep='')
    stop(msg)
  }
  
  
  df <- sigma.regions(mu.range, sigma.range, B, w$f, z)
  out <- list(q=df$q, mu.range=mu.range, sigma0=df$from, sigma1=df$to, B=B, z=z, f=w$f, R=R)
  
  df.names2keep <- c('Quadratic.sigma.cdf.segments', 'Linear.sigma.cdf.segments', 'J')
  df.names <- intersect(names(df), df.names2keep)
  out[df.names] <- df[df.names] 
  
  out
} # end of riskband.rjags.parms



#*****************************************

riskband.rjags.model <- function(file, parms)
{
  # Version 0.1 (29 september 2018)
  
  # see if Linear.sigma.cdf.segments is required
  
  required.Linear <- !is.na(match('Linear.sigma.cdf.segments', names(parms)))
  
  # Read model (and edit if necessary -- that is, if Linear.sigma.cdf.segments is not required)
  
  model <- readLines(file)
  
  if (!required.Linear)
  {
    Linear <- grep('Linear', model)
    closing.brackets <- grep('}', model)
    closing.bracket <- min(closing.brackets[closing.brackets > Linear])
    model <- model[-seq(from=Linear, to=closing.bracket)]
    model <- gsub('Quadratic.sigma.cdf.segments', '1:J', model)
  }
  
  model
} # end of riskband.rjags.model
