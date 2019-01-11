##############################################################
#
#   WEBEXPO official R scripts
#
#   SEG ANALYSIS
#  
#   BAYESIAN MODELS IN JAGS 
#
#   INFORMEDVAR FUNCTION
#
#   V1.0 Sept 2018  
#   
#
#############################################################

fun.jags.informedvar <-function(#observations
                                y , CensorType , censorLimitMat , 
                                
                                #past.data
                                past.data=numeric(0) ,
                                
                                #distribution
                                is.lognormal=TRUE,
                                
                                #measurement error
                                error.type = "none",
                                me.range = c(0.3,0.3),
                                
                                #Prior (default values valid for the lognormal distribution)
                                mu.lower = -20, 
                                mu.upper = 20, 
                                log.sigma.mu = -0.1744, 
                                log.sigma.prec = 2.5523, 
                                
                                #MCMC parameters
                                n.iter=25000, 
                                n.burnin=5000, 
                                init.mu = log(0.3),   #(default values valid for the lognormal distribution)
                                init.sigma = log(2.5))

{
  
##### preparation of measurement error parameters

me.lower <-me.range[1]
me.upper <-me.range[2]

if (me.lower==me.upper) me.upper <- me.upper+0.0001   # TRICK TO AVOID ERRORS IN THE JAGS MODEL

###JAGS MODELS

jags.model <- jags.model.informedvar(is.lognormal=is.lognormal,
                                     error.type=error.type,
                                     past.data=past.data,
                                     mu.lower = mu.lower, 
                                     mu.upper = mu.upper,
                                     log.sigma.mu = log.sigma.mu, 
                                     log.sigma.prec = log.sigma.prec,
                                     me.lower=me.lower,
                                     me.upper=me.upper)


### DATA LIST FOR JAGS

dataList = list( 
  
  y = y ,
  
  n.obs = length(y) ,
  
  CensorType = CensorType ,
  
  censorLimitMat = censorLimitMat
)


# IF past data is populated

if (length(past.data)==3) {
  
  dataList$mean.past <- past.data$mean
  dataList$sd.past <- past.data$sd
  dataList$n.past <- past.data$n
  dataList$ns2 <- (dataList$n.past - 1) * (dataList$sd.past^2)
  dataList$ns2.a <- (dataList$n.past - 1) / 2
}




#### initialising chains for random variables

initsList = list( mu = init.mu , log.sigma = log(init.sigma) )

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


