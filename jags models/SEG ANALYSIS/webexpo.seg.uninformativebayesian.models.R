##############################################################
#
#   WEBEXPO official R scripts
#
#   SEG ANALYSIS
#  
#   BAYESIAN MODELS IN JAGS 
#
#   UNINFORMATIVE FUNCTION - JAGS MODEL STATEMENTS
#   
#   V1.0 Sept 2018  
#
############################################################# 

jags.model.uninformative <-function(
  
                                  # Distribution
                                  is.lognormal=TRUE, 
                                  
                                  # measurement error structure
                                  error.type="none",
                                  me.lower=0.3,
                                  me.upper=0.3,
                                  
                                  # prior definition  (default values valid for the lognormal distribution)
                                  mu.lower.uninf = -20, 
                                  mu.upper.uninf = 20 ,
                                  sd.range=c(0, 2.3) 
                                 )
  
{

  
##formatting  sd.range
  
  sd.lower <-sd.range[1]
  sd.upper <-sd.range[2]
  
###### LOGNORMAL MODELS
  
if (is.lognormal)  
  
{

### measurement error as CV
  
if (error.type=="CV") model.jags <-paste("
                          
                                      model {
                                      
                                      
                                      ##prior on mu
                                      
                                      mu ~ dunif(",mu.lower.uninf,",",mu.upper.uninf,")
                                      
                                      ##prior on sigma
                                      
                                      sigma ~ dunif(",sd.lower,",",sd.upper,")
                                      
                                      precision <- pow( sigma , -2)

                                      ## prior on CV
                                      
                                      CV ~ dunif(",me.lower,",", me.upper,")

                                      ###likelihood
                                      
                                      for (i in 1:n.obs) {
                                      
                                      true.value[i] ~ dlnorm(mu, precision)
                                      
                                      y[i] ~ dnorm (true.value[i], MeasErrorPrecision[i])
                                      
                                      MeasErrorPrecision[i] <- pow(CV * true.value[i], -2)
                                      
                                      CensorType[i] ~ dinterval( y[i] , censorLimitMat[i,] )
                                      }
                                     
                                      }
                                      
                                      ",sep='')

### measurement error as SD

if (error.type=="SD") model.jags <-paste("
                          
                                         model {
                                         
                                         
                                         ##prior on mu
                                         
                                         mu ~dunif(",mu.lower.uninf,",",mu.upper.uninf,")
                                         
                                         ##prior on sigma
                                         
                                         sigma ~ dunif(",sd.lower,",",sd.upper,")
                                      
                                         precision <- pow( sigma , -2)
                                         
                                         ## prior on SD

                                          MeasErrorPrecision <-1/(SD*SD)
                                         
                                         SD ~ dunif(",me.lower,",", me.upper,")
                                         
                                         ###likelihood
                                         
                                         for (i in 1:n.obs) {
                                         
                                         true.value[i] ~ dlnorm(mu, precision)
                                         
                                         y[i] ~ dnorm (true.value[i], MeasErrorPrecision)
                                         
                                         CensorType[i] ~ dinterval( y[i] , censorLimitMat[i,] )

                                         }
                   
                                         }
                                         
                                         ",sep='')

### no measurement error

if (error.type=="none") model.jags <-paste("
                          
                                         model {
                                         
                                         
                                         ##prior on mu
                                         
                                         mu ~dunif(",mu.lower.uninf,",",mu.upper.uninf,")
                                         
                                         ##prior on sigma
                                         
                                         sigma ~ dunif(",sd.lower,",",sd.upper,")
                                      
                                         precision <- pow( sigma , -2)
                                         
                                        
                                         ###likelihood
                                         
                                         for (i in 1:n.obs) {
                                         
                                         y[i] ~ dlnorm(mu, precision)
                                         
                                         CensorType[i] ~ dinterval( y[i] , censorLimitMat[i,] )
                                         
                                         }
                         
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
                                           
                                           mu ~dunif(",mu.lower.uninf,",",mu.upper.uninf,")
                                           
                                           ##prior on sigma
                                           
                                           sigma ~ dunif(",sd.lower,",",sd.upper,")
                                      
                                           precision <- pow( sigma , -2)
                                           
                                           ## prior on CV
                                           
                                           CV ~ dunif(",me.lower,",", me.upper,")
                                           
                                           ###likelihood
                                           
                                           for (i in 1:n.obs) {
                                           
                                           true.value[i] ~ dnorm(mu, precision)T(0,)
                                           
                                           y[i] ~ dnorm (true.value[i], MeasErrorPrecision[i])
                                           
                                           MeasErrorPrecision[i] <- pow(CV * true.value[i], -2)
                                           
                                           CensorType[i] ~ dinterval( y[i] , censorLimitMat[i,] )
                                           }
                                   
                                           }
                                           
                                           ",sep='')
  
  ### measurement error as SD
  
  if (error.type=="SD") model.jags <-paste("
                                           
                                           model {
                                           
                                           
                                           ##prior on mu
                                           
                                           mu ~dunif(",mu.lower.uninf,",",mu.upper.uninf,")
                                           
                                           ##prior on sigma
                                           
                                           sigma ~ dunif(",sd.lower,",",sd.upper,")
                                      
                                           precision <- pow( sigma , -2)
                                           
                                           ## prior on SD
                                           
                                           MeasErrorPrecision <-1/(SD*SD)
                                           
                                           SD ~ dunif(",me.lower,",", me.upper,")
                                           
                                           ###likelihood
                                           
                                           for (i in 1:n.obs) {
                                           
                                           true.value[i] ~ dnorm(mu, precision)
                                           
                                           y[i] ~ dnorm (true.value[i], MeasErrorPrecision)
                                           
                                           CensorType[i] ~ dinterval( y[i] , censorLimitMat[i,] )
                                           
                                           }
                            
                                           }
                                           
                                           ",sep='')
  
  ### no measurement error
  
  if (error.type=="none") model.jags <-paste("
                                             
                                             model {
                                             
                                             
                                             ##prior on mu
                                             
                                             mu ~dunif(",mu.lower.uninf,",",mu.upper.uninf,")
                                             
                                             ##prior on sigma
                                             
                                             sigma ~ dunif(",sd.lower,",",sd.upper,")
                                      
                                             precision <- pow( sigma , -2)
                                             
                                             
                                             ###likelihood
                                             
                                             for (i in 1:n.obs) {
                                             
                                             y[i] ~ dnorm(mu, precision)
                                             
                                             CensorType[i] ~ dinterval( y[i] , censorLimitMat[i,] )
                                             
                                             }
                                          
                                             }
                                             
                                             ",sep='')
  
  
}


return(model.jags)
     
}