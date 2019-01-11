##############################################################
#
#   WEBEXPO official R scripts
#
#   BETWEEN WORKER ANALYSIS
#  
#   BAYESIAN MODELS IN JAGS 
#
#   UNINFORMATIVE FUNCTION - JAGS MODEL STATEMENTS
#   
#
#   
#   V1.0 Sept 2018  
#
############################################################# 

jags.model.uninformative.between <-function(
                                  # Distribution
                                  is.lognormal=TRUE,
                                  # Measurement error
                                  error.type="none",
                                  me.lower = 0.3,
                                  me.upper = 0.31,
                                  # Priors (default values valid for the lognormal model)
                                  mu.overall.lower = -20, 
                                  mu.overall.upper = 20,
                                  sigma.between.range=c(0.095,2.3),
                                  sigma.within.range=c(0.095,2.3))
  
{

  
##formatting  sd.range
  
  sd.between.lower <-sigma.between.range[1]
  sd.between.upper <-sigma.between.range[2]
  
  sd.within.lower <-sigma.within.range[1]
  sd.within.upper <-sigma.within.range[2]
  
###### LOGNORMAL MODELS
  
if (is.lognormal)  
  
{


if (error.type=="CV") model.jags <-paste("
                          
                                      model {
                                      
                                      
                                      ##prior on fixed effects
                                      
                                      mu.overall ~ dunif(",mu.overall.lower,",",mu.overall.upper,")
                                      
                                      sigma.within ~ dunif(",sd.within.lower,",",sd.within.upper,")

                                      precision.within <- pow( sigma.within , -2)

                                      #random effect / hierarchy

                                      for (i in 1:n.group)
                                         {
                                         
                                        worker.effect[i] ~ dnorm(0, precision.between) 
                                         }
                                         
                                      sigma.between ~ dunif(",sd.between.lower,",",sd.between.upper,")

                                      precision.between <- pow( sigma.between , -2)



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
                                         
                                         sigma.within ~ dunif(",sd.within.lower,",",sd.within.upper,")
                                         
                                         precision.within <- pow( sigma.within , -2)
                                         
                                         #random effect / hierarchy
                                         
                                         for (i in 1:n.group)
                                         {
                                         
                                         worker.effect[i] ~ dnorm(0, precision.between) 
                                         }
                                         
                                         sigma.between ~ dunif(",sd.between.lower,",",sd.between.upper,")
                                         
                                         precision.between <- pow( sigma.between , -2)
                                         
                                         
                                         
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
                                           
                                           sigma.within ~ dunif(",sd.within.lower,",",sd.within.upper,")
                                           
                                           precision.within <- pow( sigma.within , -2)
                                           
                                           #random effect / hierarchy
                                           
                                           for (i in 1:n.group)
                                           {
                                           
                                           worker.effect[i] ~ dnorm(0, precision.between) 
                                           }
                                           
                                           sigma.between ~ dunif(",sd.between.lower,",",sd.between.upper,")
                                           
                                           precision.between <- pow( sigma.between , -2)

                                           
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
                                           
                                           sigma.within ~ dunif(",sd.within.lower,",",sd.within.upper,")
                                           
                                           precision.within <- pow( sigma.within , -2)
                                           
                                           #random effect / hierarchy
                                           
                                           for (i in 1:n.group)
                                           {
                                           
                                           worker.effect[i] ~ dnorm(0, precision.between) 
                                           }
                                           
                                           sigma.between ~ dunif(",sd.between.lower,",",sd.between.upper,")
                                           
                                           precision.between <- pow( sigma.between , -2)
                                           
                                           
                                           
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
                                           
                                           sigma.within ~ dunif(",sd.within.lower,",",sd.within.upper,")
                                           
                                           precision.within <- pow( sigma.within , -2)
                                           
                                           #random effect / hierarchy
                                           
                                           for (i in 1:n.group)
                                           {
                                           
                                           worker.effect[i] ~ dnorm(0, precision.between) 
                                           }
                                           
                                           sigma.between ~ dunif(",sd.between.lower,",",sd.between.upper,")
                                           
                                           precision.between <- pow( sigma.between , -2)
                                           
                                           
                                           
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
                                             
                                             sigma.within ~ dunif(",sd.within.lower,",",sd.within.upper,")
                                             
                                             precision.within <- pow( sigma.within , -2)
                                             
                                             #random effect / hierarchy
                                             
                                             for (i in 1:n.group)
                                             {
                                             
                                             worker.effect[i] ~ dnorm(0, precision.between) 
                                             }
                                             
                                             sigma.between ~ dunif(",sd.between.lower,",",sd.between.upper,")
                                             
                                             precision.between <- pow( sigma.between , -2)
                                             
                                             
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