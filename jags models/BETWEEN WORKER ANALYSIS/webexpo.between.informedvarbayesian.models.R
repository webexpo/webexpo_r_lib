##############################################################
#
#   WEBEXPO official R scripts
#
#   BETWEEN WORKER ANALYSIS
#  
#   BAYESIAN MODELS IN JAGS 
#
#   INFORMEDVAR FUNCTION - JAGS MODEL STATEMENTS
#   
#
#   
#   V1.0 Sept 2018  
#
############################################################# 

jags.model.informedvar.between <-function(
                                  #Distribution
                                  is.lognormal=TRUE,
                                  #Measurement error structure
                                  error.type="none",
                                  me.lower = 0.3,
                                  me.upper = 0.31,
                                  #priors (default values valid for lognormal distribution)
                                  mu.overall.lower = -20, 
                                  mu.overall.upper = 20,
                                  log.sigma.between.mu=-0.8786,
                                  log.sigma.between.prec=1.634, 
                                  log.sigma.within.mu=-0.4106,
                                  log.sigma.within.prec=1.9002)
  
{

  
###### LOGNORMAL MODELS
  
if (is.lognormal)  
  
{


if (error.type=="CV") model.jags <-paste("
                          
                                      model {
                                      
                                      
                                      ##prior on fixed effects
                                      
                                      mu.overall ~ dunif(",mu.overall.lower,",",mu.overall.upper,")
                                      
                                      precision.within <- pow( sigma.within , -2)
                                      
                                      sigma.within <-exp(log.sigma.within)
                                      
                                      log.sigma.within ~dnorm(",log.sigma.within.mu,",",log.sigma.within.prec,")


                                      #random effect / hierarchy

                                      for (i in 1:n.group)
                                         {
                                         
                                        worker.effect[i] ~ dnorm(0, precision.between) 
                                         }
                                         
                                         precision.between <- pow( sigma.between , -2)
                                      
                                         sigma.between <-exp(log.sigma.between)
                                         
                                         log.sigma.between ~dnorm(",log.sigma.between.mu,",",log.sigma.between.prec,")
                                         



                                      ## prior on CV
                                      
                                      CV ~ dunif(",me.lower,",", me.upper,")

                                     
                                       ###likelihood
                                      
                                         for (i in 1:n.obs) 

                                        {

                                         true.value[i] ~ dlnorm(mu[i], precision.within)

                                         mu[i] <- mu.overall + worker.effect[groups[i]]
                                         
                                         y[i] ~ dnorm(true.value[i], MeasErrorPrecision[i])

                                         MeasErrorPrecision[i] <- pow(CV * true.value[i], -2)
                                         
                                         CensorType[i] ~ dinterval(y[i], censorLimitMat[i,])
                                      }
                                     
                                      }
                                      
                                      ",sep='')


if (error.type=="SD") model.jags <-paste("
                          
                                      model {
                                      
                                         
                                      ##prior on fixed effects
                                      
                                         mu.overall ~ dunif(",mu.overall.lower,",",mu.overall.upper,")
                                         
                                         precision.within <- pow( sigma.within , -2)
                                         
                                         sigma.within <-exp(log.sigma.within)
                                         
                                         log.sigma.within ~dnorm(",log.sigma.within.mu,",",log.sigma.within.prec,")
                                         
                                         
                                         #random effect / hierarchy
                                         
                                         for (i in 1:n.group)
                                         {
                                         
                                         worker.effect[i] ~ dnorm(0, precision.between) 
                                         }
                                         
                                         precision.between <- pow( sigma.between , -2)
                                         
                                         sigma.between <-exp(log.sigma.between)
                                         
                                         log.sigma.between ~dnorm(",log.sigma.between.mu,",",log.sigma.between.prec,")
                                         
                                         
                                         
                                         ## prior on SD
                                         
                                         SD ~ dunif(",me.lower,",", me.upper,")
                                         
                                         
                                         ###likelihood
                                         
                                         for (i in 1:n.obs) 
                                         
                                         {
                                         
                                         true.value[i] ~ dlnorm(mu[i], precision.within)
                                         
                                         mu[i] <- mu.overall + worker.effect[groups[i]]
                                         
                                         y[i] ~ dnorm(true.value[i], MeasErrorPrecision[i])
                                         
                                         MeasErrorPrecision[i] <- pow( SD , -2)
                                         
                                         CensorType[i] ~ dinterval(y[i], censorLimitMat[i,])
                                         }
                                         
}                                         
                                         ",sep='')

if (error.type=="none") model.jags <-paste("
                          
                                      model {
                                      
                                           
                                      ##prior on fixed effects
                                      
                                      mu.overall ~ dunif(",mu.overall.lower,",",mu.overall.upper,")
                                      
                                      precision.within <- pow( sigma.within , -2)
                                      
                                      sigma.within <-exp(log.sigma.within)
                                      
                                      log.sigma.within ~dnorm(",log.sigma.within.mu,",",log.sigma.within.prec,")


                                      #random effect / hierarchy

                                      for (i in 1:n.group)
                                         {
                                         
                                        worker.effect[i] ~ dnorm(0, precision.between) 
                                         }
                                         
                                         precision.between <- pow( sigma.between , -2)
                                      
                                         sigma.between <-exp(log.sigma.between)
                                         
                                         log.sigma.between ~dnorm(",log.sigma.between.mu,",",log.sigma.between.prec,")
                                         

                                           
                                           ###likelihood
                                           
                                           for (i in 1:n.obs) 
                                           
                                           {
                                           
                                           y[i] ~ dlnorm(mu[i], precision.within)
                                           
                                           mu[i] <- mu.overall + worker.effect[groups[i]]
                                           
                                           CensorType[i] ~ dinterval(y[i], censorLimitMat[i,])
                                           }
                                           
}                                         
                                         ",sep='')


}
  

###### NORMAL MODELS

if (!is.lognormal)  
{
  
  
  if (error.type=="CV") model.jags <-paste("
                                           
                                           model {
                                           
                                           
                                      ##prior on fixed effects
                                      
                                           mu.overall ~ dunif(",mu.overall.lower,",",mu.overall.upper,")
                                           
                                           precision.within <- pow( sigma.within , -2)
                                           
                                           sigma.within <-exp(log.sigma.within)
                                           
                                           log.sigma.within ~dnorm(",log.sigma.within.mu,",",log.sigma.within.prec,")
                                           
                                           
                                           #random effect / hierarchy
                                           
                                           for (i in 1:n.group)
                                           {
                                           
                                           worker.effect[i] ~ dnorm(0, precision.between) 
                                           }
                                           
                                           precision.between <- pow( sigma.between , -2)
                                           
                                           sigma.between <-exp(log.sigma.between)
                                           
                                           log.sigma.between ~dnorm(",log.sigma.between.mu,",",log.sigma.between.prec,")
                                           
                                           
                                           
                                           
                                           ## prior on CV
                                           
                                           CV ~ dunif(",me.lower,",", me.upper,")
                                           
                                           
                                           ###likelihood
                                           
                                           for (i in 1:n.obs) 
                                           
                                           {
                                           
                                           true.value[i] ~ dnorm(mu[i], precision.within)T(0,)
                                           
                                           mu[i] <- mu.overall + worker.effect[groups[i]]
                                           
                                           y[i] ~ dnorm(true.value[i], MeasErrorPrecision[i])
                                           
                                           MeasErrorPrecision[i] <- pow(CV * true.value[i], -2)
                                           
                                           CensorType[i] ~ dinterval(y[i], censorLimitMat[i,])
                                           }
                                           
                                           }
                                           
                                           ",sep='')
  
  
  if (error.type=="SD") model.jags <-paste("
                                           
                                           model {
                                           
                                           
                                      ##prior on fixed effects
                                      
                                           mu.overall ~ dunif(",mu.overall.lower,",",mu.overall.upper,")
                                           
                                           precision.within <- pow( sigma.within , -2)
                                           
                                           sigma.within <-exp(log.sigma.within)
                                           
                                           log.sigma.within ~dnorm(",log.sigma.within.mu,",",log.sigma.within.prec,")
                                           
                                           
                                           #random effect / hierarchy
                                           
                                           for (i in 1:n.group)
                                           {
                                           
                                           worker.effect[i] ~ dnorm(0, precision.between) 
                                           }
                                           
                                           precision.between <- pow( sigma.between , -2)
                                           
                                           sigma.between <-exp(log.sigma.between)
                                           
                                           log.sigma.between ~dnorm(",log.sigma.between.mu,",",log.sigma.between.prec,")
                                           
                                           
                                           
                                           
                                           ## prior on SD
                                           
                                           SD ~ dunif(",me.lower,",", me.upper,")
                                           
                                           
                                           ###likelihood
                                           
                                           for (i in 1:n.obs) 
                                           
                                           {
                                           
                                           true.value[i] ~ dnorm(mu[i], precision.within)
                                           
                                           mu[i] <- mu.overall + worker.effect[groups[i]]
                                           
                                           y[i] ~ dnorm(true.value[i], MeasErrorPrecision[i])
                                           
                                           MeasErrorPrecision[i] <- pow( SD , -2)
                                           
                                           CensorType[i] ~ dinterval(y[i], censorLimitMat[i,])
                                           }
                                           
                                           }                                         
                                           ",sep='')
  
  if (error.type=="none") model.jags <-paste("
                                             
                                             model {
                                             
                                             
                                      ##prior on fixed effects
                                      
                                      mu.overall ~ dunif(",mu.overall.lower,",",mu.overall.upper,")
                                      
                                      precision.within <- pow( sigma.within , -2)
                                      
                                      sigma.within <-exp(log.sigma.within)
                                      
                                      log.sigma.within ~dnorm(",log.sigma.within.mu,",",log.sigma.within.prec,")


                                      #random effect / hierarchy

                                      for (i in 1:n.group)
                                         {
                                         
                                        worker.effect[i] ~ dnorm(0, precision.between) 
                                         }
                                         
                                         precision.between <- pow( sigma.between , -2)
                                      
                                         sigma.between <-exp(log.sigma.between)
                                         
                                         log.sigma.between ~dnorm(",log.sigma.between.mu,",",log.sigma.between.prec,")
                                         

                                             
                                             ###likelihood
                                             
                                             for (i in 1:n.obs) 
                                             
                                             {
                                             
                                             y[i] ~ dnorm(mu[i], precision.within)
                                             
                                             mu[i] <- mu.overall + worker.effect[groups[i]]
                                             
                                             CensorType[i] ~ dinterval(y[i], censorLimitMat[i,])
                                             }
                                             
                                             }                                         
                                             ",sep='')
  
  
}

return(model.jags)
     
}