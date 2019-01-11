##############################################################
#
#   WEBEXPO official R scripts
#
#   SEG ANALYSIS
#  
#   BAYESIAN MODELS IN JAGS 
#
#   RISKBAND FUNCTION
#   
#   V1.0 Sept 2018 
#
#
#
#############################################################

fun.jags.riskband <-function(#observations
                              y , CensorType , censorLimitMat , 
                              
                              #past.data
                              past.data=numeric(0) ,
                              
                              #distribution
                              is.lognormal=TRUE,
                              
                              #measurement error
                              error.type = "none",
                              me.range = c(0.3,0.3),
                              
                              #MCMC
                              n.iter = 25000 , 
                              n.burnin = 5000 ,
                               
                              #PRIOR (default values valid for the lognormal distribution)
                               A = c(0.01,0.1,0.5,1) ,
                              
                              # mu and sigma are the log transformed gm and gsds if is.lognormal=TRUE
                              region.prior.prob = rep(0.2, 5) ,
                              mu.lower.riskb = log(0.001) ,
                              mu.upper.riskb = log(1000) ,
                              sigma.lower = log(1.05) ,
                              sigma.upper = log(10) ,
                              target_perc = 95)

{
  
##### preparation of measurement error parameters

me.lower <-me.range[1]
me.upper <-me.range[2]

if (me.lower==me.upper) me.upper <- me.upper+0.0001 # TRICK TO AVOID ERRORS IN THE JAGS MODEL


##### applying the PARMS function


if(is.lognormal) parms <- riskband.rjags.parms(m.range = exp(c(  mu.lower.riskb , mu.upper.riskb )), 
                                               s.range = exp(c( sigma.lower , sigma.upper )),
                                               A = A,
                                               region.prior.prob = region.prior.prob, 
                                               level = target_perc/100,
                                               outcome.is.logNormally.distributed = is.lognormal )

if(!is.lognormal) parms <- riskband.rjags.parms(m.range = c( mu.lower.riskb , mu.upper.riskb  ), 
                                                s.range = c( sigma.lower , sigma.upper ),
                                                A = A,
                                                region.prior.prob = region.prior.prob, 
                                                level = target_perc/100,
                                                outcome.is.logNormally.distributed = is.lognormal )


###JAGS MODELS

jags.model <- jags.model.riskband(is.lognormal=is.lognormal,
                                     error.type=error.type,
                                     me.lower=me.lower,
                                     me.upper=me.upper)


## Applyin the function testing the bands and potentially modifying the model


jags.model <- riskband.rjags.model(jags.model, parms)


### DATA LIST FOR JAGS

dataList = list( 
  
  x = y ,
  
  N = length(y) ,
  
  CensorType = CensorType ,
  
  censorLimitMat = censorLimitMat
)

dataList <- c( dataList, parms )



#### initialising chains for random variables

initsList = list(sigma.U=0.5, mu.U=0.4)

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


###########################################################


#Function to verify pertinence of A given the central and variance parameter ranges


webexpo.riksbandcheck <- function(is.lognormal,
                                  A = c(0.01,0.1,0.5,1) ,
                                  mu.lower.riskb = log(0.001) ,
                                  mu.upper.riskb = log(1000) ,
                                  sigma.lower = log(1.05) ,
                                  sigma.upper = log(10) ,
                                  target_perc = 95) {
  
  # quantile defining the prior
  z <- qnorm( target_perc / 100 )
  
  if (is.lognormal) {
    
    mu.range <- c( mu.lower.riskb , mu.upper.riskb )
    sigma.range <-  c( sigma.lower , sigma.upper )
    
    #rectangle
    plot( mu.range , sigma.range, col='white')
    rect( min( mu.range ) , min( sigma.range ) , max( mu.range ) , max( sigma.range ) )
    
    #zones defined by A
    for (a in A) abline( log(a)/z , -1/z , col = 'gray')
  }
  
  if (!is.lognormal) {
  
  mu.range <- c( mu.lower.riskb , mu.upper.riskb )
  sigma.range <- c( sigma.lower , sigma.upper )
  
  #rectangle
  plot( mu.range , sigma.range, col='white')
  rect( min( mu.range ) , min( sigma.range ) , max( mu.range ) , max( sigma.range ) )
  
  #zones defined by A
  for (a in A) abline( a/z , -1/z , col='gray')
  
  
  }
  
  
  
}
