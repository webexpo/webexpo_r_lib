##############################################################
#
#   WEBEXPO official STAN scripts
#
#   SEG ANALYSIS
#  
#   STAN MODELS WRITTEN BY P. BELISLE 
#
#   Global function
#   
#   V0.5 Sept 2024 
#
#############################################################

#required scripts
# The following scripts must be sourced prior to sourcing this script

        
        #data preparation function 
        
        #        source("DATA PREPARATION/SEG ANALYSIS/webexpo.seg.dataprep.R")
        
        #bayesian stan functions  
        

        #        source("STAN MODELS/McGILL FUNCTIONS/stan-fcts.R")
        
        #        source("STAN MODELS/McGILL FUNCTIONS/SEG-informedvar-stan.R")
        
        #        source("STAN MODELS/McGILL FUNCTIONS/SEG-uninformative-stan.R")
        
        #        source("STAN MODELS/McGILL FUNCTIONS/SEG-informedMean-stan.R")


##
#
#  INPUT : 
#
# vector of  observations 
# OEL
# Distribution
# Measurement error specification
# MCMC parameters
# prior specification
# model list (see compile-stan-models.R and Webexpo_examples-Stan.R)
# WARNING : Prior to running the function, stan model object(s) must be created and stored in a list.



Webexpo.seg.globalbayesian.stan <-function( 
                            

                  data.sample = c("27.8","41","29","<36.3","<41.7","<10.8","54.4","30","[40.1-90.2]","200"), 
                  
                  #distribution
                  is.lognormal = TRUE , 
                  
                  #measurement error options
                  error.type = "none" , 
                  me.range = c(0.3,0.3) , 
                  
                  #input parameters
                  oel = 100 ,
                
                  #prior information
                  prior.model = "informedvar" ,   ## options : ("informedvar"  / "uninformative" / "riskband")
                  
                  #informed var prior (default values valid for the lognormal distribution)
                  mu.lower = -20 ,
                  mu.upper = 20  ,
                  
                  log.sigma.mu = -0.1744 ,
                  log.sigma.prec = 2.5523 , 
                  
                  #uninformative prior
                  mu.lower.uninf = -20 ,
                  mu.upper.uninf = 20  ,
                  sd.range = c(0, 2.3)  ,
                  
                  past.data = NULL ,
                  
                  #riskband prior  (default values valid for the lognormal distribution)
                  
                  # mu and sigma are the log transformed gm and gsds if is.lognormal=TRUE
                  
                  A = c(0.01,0.1,0.5,1) ,
                  region.prior.prob = rep(0.2, 5) ,
                  mu.lower.riskb = log(0.001) ,
                  mu.upper.riskb = log(1000) ,
                  sigma.lower = log(1.05) ,
                  sigma.upper = log(10) ,
                  target_perc = 95 , 
                  
                  
                  #MCMC parameters  
                  n.iter = 25000 ,
                  n.burnin = 5000 , 
                  
                  init.mu = log(0.3) ,  #(default values valid for the lognormal distribution)
                  init.sigma = log(2.5),
                  
                  models.list=models.list
                  
) {
  

  
  
#input preparation
  
        formatted.data <-webexpo.seg.datapreparation(data.in = data.sample)
        
        x <- formatted.data$data
        
  
#censored data treatment preparation
  
        ##lognormal
        
        if (is.lognormal) {
          
          observed_values <-as.numeric(x[formatted.data$notcensored])/oel
          
          #left censored values
          
          leftcensored_values <-as.numeric(substring(x[formatted.data$leftcensored],2))/oel
          
          #right censored values
          
          rightcensored_values <-as.numeric(substring(x[formatted.data$rightcensored],2))/oel
          
          #interval censored values
          
          intcensored_left_values <-as.numeric(substring(x[formatted.data$intcensored],2,regexpr("-",x[formatted.data$intcensored],fixed=TRUE)-1))/oel
          
          intcensored_right_values <-as.numeric(substring(x[formatted.data$intcensored],regexpr("-",x[formatted.data$intcensored],fixed=TRUE)+1,nchar(x[formatted.data$intcensored])-1))/oel
          
        }
        
        ##normal
        
        if (!is.lognormal) {
          
          observed_values <-as.numeric(x[formatted.data$notcensored])
          
          #left censored values
          
          leftcensored_values <-as.numeric(substring(x[formatted.data$leftcensored],2))
          
          #right censored values
          
          rightcensored_values <-as.numeric(substring(x[formatted.data$rightcensored],2))
          
          #interval censored values
          
          intcensored_left_values <-as.numeric(substring(x[formatted.data$intcensored],2,regexpr("-",x[formatted.data$intcensored],fixed=TRUE)-1))
          
          intcensored_right_values <-as.numeric(substring(x[formatted.data$intcensored],regexpr("-",x[formatted.data$intcensored],fixed=TRUE)+1,nchar(x[formatted.data$intcensored])-1))
          
        }
        

## preparation of past.data in case of lognormal model
        
        if (length(past.data)==3 & is.lognormal) past.data$mean <- past.data$mean-log(oel)
        
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
  
  if (prior.model == "informedvar") {
    
    
    res <-SEG.informedvar.stan(y = observed_values , 
                                lt = leftcensored_values ,
                                gt = rightcensored_values ,
                                interval.lower = intcensored_left_values ,
                                interval.upper = intcensored_right_values ,
                                n.iter = n.iter , 
                                n.burnin = n.burnin ,
                                mu.lower = mu.lower ,
                                mu.upper = mu.upper ,
                                log.sigma.mu = log.sigma.mu ,
                                log.sigma.prec = log.sigma.prec ,
                                init.mu = init.mu ,
                                init.sd = init.sigma ,
                                outcome.is.logNormally.distributed = is.lognormal ,
                                past.data = past.data ,
                                me.sd.range = me.sd.range , 
                                cv.range = cv.range ,
                                models.list=models.list) 
    
    
  }
  
  
  if (prior.model == "uninformative") {
    
    
    res <-SEG.uninformative(y = observed_values , 
                            lt = leftcensored_values ,
                            gt = rightcensored_values ,
                            interval.lower = intcensored_left_values ,
                            interval.upper = intcensored_right_values ,
                            n.iter = n.iter , 
                            n.burnin = n.burnin ,
                            mu.lower = mu.lower.uninf ,
                            mu.upper = mu.upper.uninf ,
                            init.mu = init.mu ,
                            init.sd = init.sigma ,
                            sd.range = sd.range ,
                            outcome.is.logNormally.distributed = is.lognormal ,
                            me.sd.range = me.sd.range , 
                            cv.range = cv.range ,
                            models.list=models.list) 
    
    
  }
  


  if (is.lognormal) results <- list(mu.chain=res$mu+log(oel),sigma.chain=res$sigma)
  
  if (!is.lognormal) results <- list(mu.chain=res$mu,sigma.chain=res$sigma)
  
  
  return(results)
  
  
}    


  



  
  
  
  
