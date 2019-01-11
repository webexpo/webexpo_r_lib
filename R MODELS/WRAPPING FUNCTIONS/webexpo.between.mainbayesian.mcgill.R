##############################################################
#
#   WEBEXPO official R scripts
#
#   BETWEEN WORKER ANALYSIS
#  
#   BAYESIAN MODELS WRITTEN BY L. JOSEPH and P. BELISLE 
#
#   Global function
#
#   V1.0 Sept 2018 
#
#############################################################




##
#
#  INPUT : 
#
# data.frame of observations (exposures and workers) 
# OEL
# Distribution
# Measurement error specification
# MCMC parameters
# prior specification
# 
##

#required scripts
# The following scripts must be sourced prior to sourcing this script

    #sourcing data preparation function 
    
      # source("DATA PREPARATION/BETWEEN WORKER ANALYSIS/webexpo.between.dataprep.R")
    
    #sourcing bayesian functions  
    
      # source("R MODELS/McGILL FUNCTIONS/data-summary.R")
      
      # source("R MODELS/McGILL FUNCTIONS/fcts.R")
      
      # source("R MODELS/McGILL FUNCTIONS/model-Between-worker.R")




Webexpo.between.globalbayesian.mcgill <-function( 
                            

  data.sample = sample.A,
  
  #distribution
  is.lognormal = TRUE,  
  
  #model options
  error.type = "none",  
  me.range = c(0.3,0.3),  
  
  #input parameters
  oel = 100, 

  #prior information
  prior.model = "informedvar",    ## options : ("informedvar"  / "uninformative" )
  
  #(default values valid for the lognormal distribution)
  mu.overall.lower = -20,
  mu.overall.upper = 20,
  
  sigma.between.range=c(0.095,2.3),
  sigma.within.range=c(0.095,2.3),
  
  log.sigma.between.mu=-0.8786,
  log.sigma.between.prec=1.634,
  log.sigma.within.mu=-0.4106,
  log.sigma.within.prec=1.9002,  
  
  
  #MCMC parameters  
  n.iter = 50000, 
  n.burnin = 5000, 
  #(default values valid for the lognormal distribution)
  init.mu.overall = log(0.3),
  init.sigma.within = log(2.5))
  
  
 {
  

#input preparation
  
        formatted.data <-webexpo.between.datapreparation(data.in = data.sample)
        
        x <- formatted.data$data$x
        
  
#censored data treatment preparation
  
        ##lognormal
        
        if (is.lognormal) {
          
          
          #observed values
          
          observed_values <-as.numeric(x[formatted.data$notcensored])/oel
          
          observed_workers <-formatted.data$data$worker[formatted.data$notcensored]
          
          #left censored values
          
          leftcensored_values <-as.numeric(substring(x[formatted.data$leftcensored],2))/oel
          
          leftcensored_workers <-formatted.data$data$worker[formatted.data$leftcensored]
          
          #right censored values
          
          rightcensored_values <-as.numeric(substring(x[formatted.data$rightcensored],2))/oel
          
          rightcensored_workers <-formatted.data$data$worker[formatted.data$rightcensored]
          
          #interval censored values
          
          intcensored_left_values <-as.numeric(substring(x[formatted.data$intcensored],2,regexpr("-",x[formatted.data$intcensored],fixed=TRUE)-1))/oel
          
          intcensored_right_values <-as.numeric(substring(x[formatted.data$intcensored],regexpr("-",x[formatted.data$intcensored],fixed=TRUE)+1,nchar(x[formatted.data$intcensored])-1))/oel
          
          intcensored_workers <-formatted.data$data$worker[formatted.data$intcensored] 
        }
        
        ##normal
        
        if (!is.lognormal) {
          
          #observed values
          
          observed_values <-as.numeric(x[formatted.data$notcensored])
          
          observed_workers <-formatted.data$data$worker[formatted.data$notcensored]
          
          #left censored values
          
          leftcensored_values <-as.numeric(substring(x[formatted.data$leftcensored],2))
          
          leftcensored_workers <-formatted.data$data$worker[formatted.data$leftcensored]
          
          #right censored values
          
          rightcensored_values <-as.numeric(substring(x[formatted.data$rightcensored],2))
          
          rightcensored_workers <-formatted.data$data$worker[formatted.data$rightcensored]
          
          #interval censored values
          
          intcensored_left_values <-as.numeric(substring(x[formatted.data$intcensored],2,regexpr("-",x[formatted.data$intcensored],fixed=TRUE)-1))
          
          intcensored_right_values <-as.numeric(substring(x[formatted.data$intcensored],regexpr("-",x[formatted.data$intcensored],fixed=TRUE)+1,nchar(x[formatted.data$intcensored])-1))
          
          intcensored_workers <-formatted.data$data$worker[formatted.data$intcensored] 
          
        }
        
## preparation of measurement error 
        
        if (error.type=="CV") {
          
          me.sd.range <- numeric(0)
          cv.range <- me.range
          
        }
        
        if (error.type=="SD") {
          
          me.sd.range <- me.range
          cv.range <- numeric(0)
          
          if (is.lognormal) me.sd.range <-me.sd.range/oel
          
        }
        
        if (error.type=="none") {
          
          me.sd.range <- numeric(0)
          cv.range <- numeric(0)
          
        }        
        
        
                
## bayesian calculations
  
  if (prior.model == "informedvar") use.uniform.prior.on.sds <- FALSE
        
  if (prior.model == "uninformative") use.uniform.prior.on.sds <- TRUE
    
    res <-Between.worker(y = observed_values , 
                         worker = observed_workers ,
                         lt = leftcensored_values ,
                         worker.lt = leftcensored_workers ,
                         gt = rightcensored_values ,
                         worker.gt = rightcensored_workers ,
                         interval.lower = intcensored_left_values ,
                         interval.upper = intcensored_right_values ,
                         worker.interval = intcensored_workers,
                         n.chains = 1 ,
                         n.iter = n.iter , 
                         n.burnin = n.burnin ,
                         n.thin = 1 ,
                         monitor.burnin = FALSE ,
                         log.sigma.between.mu = log.sigma.between.mu ,
                         log.sigma.between.prec = log.sigma.between.prec ,
                         log.sigma.within.mu = log.sigma.within.mu ,
                         log.sigma.within.prec = log.sigma.within.prec ,
                         mu.overall.lower = mu.overall.lower ,
                         mu.overall.upper = mu.overall.upper ,
                         outcome.is.logNormally.distributed = is.lognormal ,
                         use.uniform.prior.on.sds = use.uniform.prior.on.sds ,
                         sigma.between.range = sigma.between.range ,
                         sigma.within.range = sigma.within.range,
                         me.sd.range = me.sd.range ,
                         cv.range = cv.range ,
                         init.mu.overall = init.mu.overall , 
                         init.sigma.within = init.sigma.within ,
                         save.RData=FALSE)
  
    
    sw.chain <-res$sample$sigma.within
    sb.chain <-res$sample$sigma.between
    mu.chain <-res$sample$mu.overall
    
    worker.index <-data.frame(worker=names(table(factor(formatted.data$data$worker))),stringsAsFactors=F)
    
    worker.index$numberinoutput <-1:length(worker.index[,1])
    
    mu.workers <-matrix(res$sample$mu.worker , nrow=length(worker.index[,1]),ncol=n.iter)
    
    results <-list(mu.overall.chain = mu.chain,
                   sigma.within.chain = sw.chain,
                   sigma.between.chain = sb.chain,
                   mu.workers = mu.workers,
                   worker.index = worker.index)  
    
      
   
  if (is.lognormal) {
    
    results$mu.overall.chain <- results$mu.overall.chain + log(oel)
    results$mu.workers <- results$mu.workers + log(oel)
  
  }
  
  return(results)
  
  
}    


  



  
  
  
  