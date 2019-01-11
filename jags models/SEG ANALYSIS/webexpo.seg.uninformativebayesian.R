##############################################################
#
#   WEBEXPO official R scripts
#
#   SEG ANALYSIS
#  
#   BAYESIAN MODELS IN JAGS 
#
#   UNINFORMATIVE FUNCTION
#
#   V1.0 Sept 2018  
#   
#
#############################################################

fun.jags.uninformative <-function(#observations
                                  y , CensorType , censorLimitMat , 
                                  
                                  #distribution
                                  is.lognormal=TRUE,
                                  
                                  #measurement error
                                  error.type = "none",
                                  me.range = c(0.3,0.3),
                                
                                  # Prior (default values valid for the lognormal distribution)
                                  mu.lower.uninf = -20, 
                                  mu.upper.uninf = 20 ,
                                  sd.range=c(0, 2.3),
                                  
                                  #MCMC
                                  n.iter=25000, 
                                  n.burnin=5000, 
                                  init.mu = log(0.3), #(default values valid for the lognormal distribution)
                                  init.sigma = log(2.5))

{
  
##### preparation of measurement error parameters

me.lower <-me.range[1]
me.upper <-me.range[2]

if (me.lower==me.upper) me.upper <- me.upper+0.0001 # TRICK TO AVOID ERRORS IN THE JAGS MODEL

###JAGS MODELS

jags.model <- jags.model.uninformative(is.lognormal=is.lognormal,
                                     error.type=error.type,
                                     mu.lower.uninf = mu.lower.uninf, 
                                     mu.upper.uninf = mu.upper.uninf,
                                     sd.range = sd.range,
                                     me.lower=me.lower,
                                     me.upper=me.upper)


### DATA LIST FOR JAGS

dataList = list( 
  
  y = y ,
  
  n.obs = length(y) ,
  
  CensorType = CensorType ,
  
  censorLimitMat = censorLimitMat
)


#### initialising chains for random variables

initsList = list( mu = init.mu , sigma = init.sigma )

### Running the model

jagsModel = jags.model( file=textConnection(jags.model), data=dataList , inits=initsList ,
                        n.chains=1 , n.adapt=500 )

update( jagsModel , n.iter=n.burnin ) #### Burn-in

codaSamples.M1 = coda.samples( jagsModel , variable.names=c("mu", "sigma") ,
                               n.iter=n.iter, thin=1 )


#####extracting the chains

mu.chain <- codaSamples.M1[[1]][,"mu"]
sigma.chain <- codaSamples.M1[[1]][,"sigma"]


return(list(mu.chain=mu.chain , sigma.chain=sigma.chain))

}


