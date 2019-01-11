##############################################################
#
#   WEBEXPO official R scripts
#
#   BETWEEN WORKER ANALYSIS
#  
#   BAYESIAN MODELS IN JAGS 
#
#   INFORMEDVAR FUNCTION
#   
#
#   
#   V1.0 Sept 2018  
#
#############################################################

fun.jags.informedvar.between <-function(
                                # Observations (transformed for the model format)
                                y , workers , CensorType , censorLimitMat , 
                                # Distribution
                                is.lognormal=TRUE,
                                # Measurement error
                                error.type = "CV",
                                me.range = c(0.2,0.2),
                                # Priors (default values valid for the lognormal model)
                                mu.overall.lower = -20, 
                                mu.overall.upper = 20,
                                log.sigma.between.mu=-0.8786,
                                log.sigma.between.prec=1.634, 
                                log.sigma.within.mu=-0.4106,
                                log.sigma.within.prec=1.9002,
                                # MCMC
                                n.iter=50000, 
                                n.burnin=5000, 
                                # Default values valid for the lognormal model
                                init.mu.overall = log(0.3), 
                                init.log.sigma.within = log(log(2.5)),
                                init.log.sigma.between = log(log(2.5)))

{
  
##### preparation of measurement error parameters

me.lower <-me.range[1]
me.upper <-me.range[2]

if (me.lower==me.upper) me.upper <- me.upper+0.0001

###JAGS MODELS

jags.model <- jags.model.informedvar.between(is.lognormal=is.lognormal,
                                     error.type=error.type,
                                     mu.overall.lower = mu.overall.lower, 
                                     mu.overall.upper =  mu.overall.upper,
                                     log.sigma.between.mu = log.sigma.between.mu,
                                     log.sigma.between.prec = log.sigma.between.prec, 
                                     log.sigma.within.mu = log.sigma.within.mu,
                                     log.sigma.within.prec = log.sigma.within.prec, 
                                     me.lower = me.lower,
                                     me.upper = me.upper)


#### initialising chains for random variables

#some inputs

#some inputs

n.obs<-length(y)

n.group <-  length(unique(workers))

groups <-as.integer(factor(workers))


#initial list

InitList <- list(log.sigma.within= init.log.sigma.within,
                 log.sigma.between=init.log.sigma.between,
                  mu.overall=init.mu.overall)


# data list

DataList <-list(y = y,
     n.obs = n.obs,
     n.group = n.group,
     groups = groups,
     CensorType = CensorType,
     censorLimitMat = censorLimitMat)

### Running the model

jagsModel = jags.model( file=textConnection(jags.model), data=DataList , inits=InitList ,
                        n.chains=1 , n.adapt=500 )

update( jagsModel , n.iter=n.burnin ) #### Burn-in

j.out <-coda.samples(model=jagsModel,
                     variable.names=c('mu.overall','sigma.between','sigma.within','worker.effect'),
                     n.iter=n.iter,
                     thin=1)


#####extracting the chains

mu.overall.chain <- j.out[[1]][,"mu.overall"]
sigma.between.chain <- j.out[[1]][,"sigma.between"]
sigma.within.chain <- j.out[[1]][,"sigma.within"]

mu.workers <-matrix(nrow=n.group,ncol=n.iter)

for (i in 1:n.group) mu.workers[i,] <-  j.out[[1]][,i+3]+j.out[[1]][,"mu.overall"]

worker.index <-data.frame(worker=names(table(factor(workers))),stringsAsFactors=F)

worker.index$numberinoutput <-1:length(worker.index[,1])

return(list(mu.overall.chain = mu.overall.chain , 
            sigma.between.chain = sigma.between.chain,
            sigma.within.chain = sigma.within.chain,
            mu.workers = mu.workers,
            worker.index = worker.index))

}


