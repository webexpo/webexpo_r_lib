##############################################################
#
#   WEBEXPO official R scripts
#
#   SEG ANALYSIS
#  
#   BAYESIAN MODELS IN JAGS 
#
#   INFORMEDVAR FUNCTION - JAGS MODEL STATEMENTS
#
#   V1.0 Sept 2018 
#   
#
#
#
############################################################# 


jags.model.informedvar <-function(# Distribution
                                  is.lognormal=TRUE, 
                                  
                                  # measurement error structure
                                  error.type="none",
                                  me.lower=0.3,
                                  me.upper=0.3,
                                  
                                  #Past data                                
                                  past.data=numeric(0),
                                  
                                  #prior informedvar (default values valid for the lognormal distribution)
                                  mu.lower = -20 ,
                                  mu.upper = 20  ,
                                  log.sigma.mu = -0.1744 ,
                                  log.sigma.prec = 2.5523
                                  )
  
{

  
###### past.data
  
if (length(past.data)==3) past.data.chunk <- paste("mean.past ~ dnorm(mu, mean.past.precision)
                                                    mean.past.precision <- n.past / pow(sigma, 2)
                                                    ns2 ~ dgamma(ns2.a, b)
                                                    b <- pow(sigma, -2) / 2",sep='')



if (length(past.data)!=3) past.data.chunk <- ''
  
###### LOGNORMAL MODELS
  
if (is.lognormal)  
  
{

  ### measurement error as CV
  
  if (error.type=="CV") model.jags <-paste("
                          
                                      model {
                                      
                                      
                                      ##prior on mu
                                      
                                      mu ~dunif(",mu.lower,",",mu.upper,")
                                      
                                      ##prior on sigma
                                      
                                      precision <-1/(sigma*sigma)
                                      
                                      sigma <-exp(log.sigma)
                                      
                                      log.sigma ~dnorm(",log.sigma.mu,",",log.sigma.prec,")
                                      
                                      ## prior on CV
                                      
                                      CV ~ dunif(",me.lower,",", me.upper,")
                                      
                                      ###likelihood
                                      
                                      for (i in 1:n.obs) {
                                      
                                      true.value[i] ~ dlnorm(mu, precision)
                                      
                                      y[i] ~ dnorm (true.value[i], MeasErrorPrecision[i])
                                      
                                      MeasErrorPrecision[i] <- pow(CV * true.value[i], -2)
                                      
                                      CensorType[i] ~ dinterval( y[i] , censorLimitMat[i,] )
                                      }
                                      ",past.data.chunk,"
                                      }
                                      
                                      ",sep='')
  ### measurement error as SD

  if (error.type=="SD") model.jags <-paste("
                          
                                         model {
                                         
                                         
                                         ##prior on mu
                                         
                                         mu ~dunif(",mu.lower,",",mu.upper,")
                                         
                                         ##prior on sigma
                                         
                                         precision <-1/(sigma*sigma)
                                         
                                         sigma <-exp(log.sigma)
                                         
                                         log.sigma ~dnorm(",log.sigma.mu,",",log.sigma.prec,")
                                         
                                         ## prior on SD

                                          MeasErrorPrecision <-1/(SD*SD)
                                         
                                         SD ~ dunif(",me.lower,",", me.upper,")
                                         
                                         ###likelihood
                                         
                                         for (i in 1:n.obs) {
                                         
                                         true.value[i] ~ dlnorm(mu, precision)
                                         
                                         y[i] ~ dnorm (true.value[i], MeasErrorPrecision)
                                         
                                         CensorType[i] ~ dinterval( y[i] , censorLimitMat[i,] )

                                         }
                                         ",past.data.chunk,"
                                         }
                                         
                                         ",sep='')
  ### no measurement error
  
  if (error.type=="none") model.jags <-paste("
                          
                                         model {
                                         
                                         
                                         ##prior on mu
                                         
                                         mu ~dunif(",mu.lower,",",mu.upper,")
                                         
                                         ##prior on sigma
                                         
                                         precision <-1/(sigma*sigma)
                                         
                                         sigma <-exp(log.sigma)
                                         
                                         log.sigma ~dnorm(",log.sigma.mu,",",log.sigma.prec,")
                                         
                                        
                                         ###likelihood
                                         
                                         for (i in 1:n.obs) {
                                         
                                         y[i] ~ dlnorm(mu, precision)
                                         
                                         CensorType[i] ~ dinterval( y[i] , censorLimitMat[i,] )
                                         
                                         }
                                         ",past.data.chunk,"
                                         }
                                         
                                         ",sep='')


}
  

###### NORMAL MODELS

if (!is.lognormal)  
  
{
  
  
  ### measurement error as CV
  
  if (error.type=="CV") model.jags <-paste("
                                           
                                           model {
                                           
                                           
                                           ##prior on mu
                                           
                                           mu ~dunif(",mu.lower,",",mu.upper,")
                                           
                                           ##prior on sigma
                                           
                                           precision <-1/(sigma*sigma)
                                           
                                           sigma <-exp(log.sigma)
                                           
                                           log.sigma ~dnorm(",log.sigma.mu,",",log.sigma.prec,")
                                           
                                           ## prior on CV
                                           
                                           CV ~ dunif(",me.lower,",", me.upper,")
                                           
                                           ###likelihood
                                           
                                           for (i in 1:n.obs) {
                                           
                                           true.value[i] ~ dnorm(mu, precision)T(0,)
                                           
                                           y[i] ~ dnorm (true.value[i], MeasErrorPrecision[i])
                                           
                                           MeasErrorPrecision[i] <- pow(CV * true.value[i], -2)
                                           
                                           CensorType[i] ~ dinterval( y[i] , censorLimitMat[i,] )
                                           }
                                           ",past.data.chunk,"
                                           }
                                           
                                           ",sep='')
  
  ### measurement error as SD
  
  if (error.type=="SD") model.jags <-paste("
                                           
                                           model {
                                           
                                           
                                           ##prior on mu
                                           
                                           mu ~dunif(",mu.lower,",",mu.upper,")
                                           
                                           ##prior on sigma
                                           
                                           precision <-1/(sigma*sigma)
                                           
                                           sigma <-exp(log.sigma)
                                           
                                           log.sigma ~dnorm(",log.sigma.mu,",",log.sigma.prec,")
                                           
                                           ## prior on SD
                                           
                                           MeasErrorPrecision <-1/(SD*SD)
                                           
                                           SD ~ dunif(",me.lower,",", me.upper,")
                                           
                                           ###likelihood
                                           
                                           for (i in 1:n.obs) {
                                           
                                           true.value[i] ~ dnorm(mu, precision)
                                           
                                           y[i] ~ dnorm (true.value[i], MeasErrorPrecision)
                                           
                                           CensorType[i] ~ dinterval( y[i] , censorLimitMat[i,] )
                                           
                                           }
                                           ",past.data.chunk,"
                                           }
                                           
                                           ",sep='')
  
  ### no measurement error
  
  if (error.type=="none") model.jags <-paste("
                                             
                                             model {
                                             
                                             
                                             ##prior on mu
                                             
                                             mu ~dunif(",mu.lower,",",mu.upper,")
                                             
                                             ##prior on sigma
                                             
                                             precision <-1/(sigma*sigma)
                                             
                                             sigma <-exp(log.sigma)
                                             
                                             log.sigma ~dnorm(",log.sigma.mu,",",log.sigma.prec,")
                                             
                                             
                                             ###likelihood
                                             
                                             for (i in 1:n.obs) {
                                             
                                             y[i] ~ dnorm(mu, precision)
                                             
                                             CensorType[i] ~ dinterval( y[i] , censorLimitMat[i,] )
                                             
                                             }
                                             ",past.data.chunk,"
                                             }
                                             
                                             ",sep='')
  
  
}


return(model.jags)
     
}