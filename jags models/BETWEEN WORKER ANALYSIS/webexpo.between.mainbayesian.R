##############################################################
#
#   WEBEXPO official R scripts
#
#   BETWEEN WORKER ANALYSIS
#  
#   BAYESIAN MODELS IN JAGS 
#
#   Global function
#   
#   V1.0 Sept 2018  
# 
#############################################################

## required libraries

library(rjags)



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

#    source("DATA PREPARATION/BETWEEN WORKER ANALYSIS/webexpo.between.dataprep.R")

          #sourcing bayesian functions
          
          #    source("JAGS MODELS/BETWEEN WORKER ANALYSIS/webexpo.between.informedvarbayesian.R")
          
          #    source("JAGS MODELS/BETWEEN WORKER ANALYSIS/webexpo.between.informedvarbayesian.models.R")
          
          #    source("JAGS MODELS/BETWEEN WORKER ANALYSIS/webexpo.between.uninformativebayesian.R")
          
          #   source("JAGS MODELS/BETWEEN WORKER ANALYSIS/webexpo.between.uninformativebayesian.models.R")


# FUNCTION

Webexpo.between.globalbayesian.jags <- function(
  
                      data.sample = sample.A, 

                      #distribution
                      is.lognormal = TRUE, 
                      
                      #model options
                      error.type = "none", 
                      me.range = c(0.3,0.3), 
                      
                      #input parameters
                      oel = 100,
                      
                      #prior information
                      prior.model = "informedvar",   ## options : ("informedvar"  / "uninformative" )
                      
                      #prior for group mean, valid for informedvar and uninformative
                      #(default values valid for the lognormal distribution)
                      mu.overall.lower = -20,
                      mu.overall.upper = 20,
                      
                      #prior for variance for uninformative
                      #(default values valid for the lognormal distribution)
                      sigma.between.range=c(0.095,2.3),
                      sigma.within.range=c(0.095,2.3),
                      
                      #prior for variance for informedvar
                      #(default values valid for the lognormal distribution)
                      log.sigma.between.mu=-0.8786,
                      log.sigma.between.prec=1.634,
                      log.sigma.within.mu=-0.4106,
                      log.sigma.within.prec=1.9002, 
                      
                      #MCMC parameters  
                      n.iter = 50000, 
                      n.burnin = 5000, 
                      
                      #(default values valid for the lognormal distribution)
                      init.mu.overall = log(0.3),
                      init.sigma.within = log(2.5),
                      init.sigma.between = log(2.5),
                      init.log.sigma.within = log(log(2.5)),
                      init.log.sigma.between = log(log(2.5)))
                                               
{


#input preparation (dealing with censorship)

      formatted.data <-webexpo.between.datapreparation(data.in = data.sample)
      
      x <- formatted.data$data$x
        
        ##lognormal
        
        if (is.lognormal) {
          
          
          ###creating the input for the censoring functions in JAGS
          
          #Explanations about the use of dinterval is detailed in the JAGS user manuel (Plummer and Northcott, 2013) and in Kruschke (Kruschke, 2014).
          
          #CensorType is a vector of the length of y, which indicates the censoring status of each observation with regard to the limits indicated by the matrix of censoring points (censorLimitMat)
          
          
          ##must create a matrix with n rows and 2 columnns
          ##column 1 is for the left censoring points (<)
          ##column 2 is for the right censoring points (>)
          
          ## values of the dinterval function
          # 0  : < to the left censoring point
          # 1  : between the 2 censoring points
          # 2 : > to the right censoring points
          
          ##coding the dinterval values
          
          #In practice :
          
          #for left censored points :
          #Censortype is 0
          #left limit of censorLimitMat is the censoring value
          #right limit of censorLimitMat is 50 (value we are sure any observation is smaller than)
          
          #for right censored points :
          #Censortype is 2
          #left limit of censorLimitMat is -50 (value we are sure any observation is greater than)
          #right limit of censorLimitMat is the censoring value
          
          #for interval censored points :
          #Censortype is 1
          #left limit of censorLimitMat is left censoring point
          #rightlimit of censorLimitMat is right censoring value
          
          #for observed points :
          #Censortype is 1
          #left limit of censorLimitMat is -50 / 1 (value we are sure any observation is greater than) lognormal/normal
          #right limit of censorLimitMat is 50 / 10000 (value we are sure any observation is smaller than) lognormal/normal                
          
          #censorLimitMat is the matrix of censoring points
          
          #x/Y should be set to NA when censored
          
         
          y <- as.numeric(x)/oel  ##Nas automatically generated
          
          ##creating the dinterval vector : CensorType
          
          CensorType <- rep(1,length(x)) 
          
          #left censored
          CensorType[formatted.data$leftcensored] <-0  
          
          #right censored
          CensorType[formatted.data$rightcensored] <-2  
          
          ## creating the matrix of limits
          censorLimitMat <- matrix(nrow = length(y) , ncol = 2)
          censorLimitMat[,1] <-rep(-50,length(y))
          censorLimitMat[,2] <-rep(50,length(y))
          
          
          #left censored
          censorLimitMat[formatted.data$leftcensored,1] <-as.numeric(substring(x[formatted.data$leftcensored],2))/oel
          #right censored
          censorLimitMat[formatted.data$rightcensored,2] <-as.numeric(substring(x[formatted.data$rightcensored],2))/oel
          #interval censored
          censorLimitMat[formatted.data$intcensored,1] <-as.numeric(substring(x[formatted.data$intcensored],2,regexpr("-",x[formatted.data$intcensored],fixed=TRUE)-1))/oel
          censorLimitMat[formatted.data$intcensored,2] <-as.numeric(substring(x[formatted.data$intcensored],regexpr("-",x[formatted.data$intcensored],fixed=TRUE)+1,nchar(x[formatted.data$intcensored])-1))/oel
          
          
        }
      

      
      ##normal
      
      if (!is.lognormal) {
        
        y <- as.numeric(x)  ##Nas automatically generated
        
        ##creating the dinterval vector : CensorType
        
        CensorType <- rep(1,length(x)) 
        
        #left censored
        CensorType[formatted.data$leftcensored] <-0  
        
        #right censored
        CensorType[formatted.data$rightcensored] <-2  
        
        ## creating the matrix of limits
        censorLimitMat <- matrix(nrow = length(y) , ncol = 2)
        censorLimitMat[,1] <-rep(1,length(y))
        censorLimitMat[,2] <-rep(10000,length(y))
        
        
        #left censored
        censorLimitMat[formatted.data$leftcensored,1] <-as.numeric(substring(x[formatted.data$leftcensored],2))
        #right censored
        censorLimitMat[formatted.data$rightcensored,2] <-as.numeric(substring(x[formatted.data$rightcensored],2))
        #interval censored
        censorLimitMat[formatted.data$intcensored,1] <-as.numeric(substring(x[formatted.data$intcensored],2,regexpr("-",x[formatted.data$intcensored],fixed=TRUE)-1))
        censorLimitMat[formatted.data$intcensored,2] <-as.numeric(substring(x[formatted.data$intcensored],regexpr("-",x[formatted.data$intcensored],fixed=TRUE)+1,nchar(x[formatted.data$intcensored])-1))
        
      }




## preparation of measurement error when lognormal model and SD measurement error

      if (error.type=="SD" & is.lognormal) me.range <-me.range/oel
  
## bayesian calculations

      if (prior.model == "informedvar") {
        
        
        res <-fun.jags.informedvar.between( y=y , workers=formatted.data$data$worker,
                                CensorType = CensorType ,
                                censorLimitMat = censorLimitMat ,
                                is.lognormal=is.lognormal,
                                error.type = error.type,
                                me.range = me.range,
                                n.iter = n.iter , 
                                n.burnin = n.burnin ,
                                mu.overall.lower = mu.overall.lower ,
                                mu.overall.upper = mu.overall.upper ,
                                log.sigma.between.mu = log.sigma.between.mu,
                                log.sigma.between.prec = log.sigma.between.prec, 
                                log.sigma.within.mu = log.sigma.within.mu,
                                log.sigma.within.prec = log.sigma.within.prec, 
                                init.log.sigma.within =  init.log.sigma.within,
                                init.log.sigma.between = init.log.sigma.between,
                                init.mu.overall=init.mu.overall) 
        
        
      }
      
      if (prior.model == "uninformative") {
        
        
        res <-fun.jags.uninformative.between( y=y , workers=formatted.data$data$worker,
                                    CensorType = CensorType ,
                                    censorLimitMat = censorLimitMat ,
                                    is.lognormal=is.lognormal,
                                    error.type = error.type,
                                    me.range = me.range,
                                    n.iter = n.iter , 
                                    n.burnin = n.burnin ,
                                    mu.overall.lower = mu.overall.lower, 
                                    mu.overall.upper =  mu.overall.upper,
                                    sigma.between.range = sigma.between.range,
                                    sigma.within.range = sigma.within.range,
                                    init.mu.overall = init.mu.overall,
                                    init.sigma.within = init.sigma.within,
                                    init.sigma.between = init.sigma.between) 
        
        
      }


## final chain preparation

      if (is.lognormal) results <-list(mu.overall.chain = res$mu.overall.chain + log(oel) , 
                                       sigma.between.chain = res$sigma.between.chain,
                                       sigma.within.chain = res$sigma.within.chain,
                                       mu.workers = res$mu.workers + log(oel),
                                       worker.index = res$worker.index)
      
      if (!is.lognormal) results <-list(mu.overall.chain = res$mu.overall.chain , 
                                        sigma.between.chain = res$sigma.between.chain,
                                        sigma.within.chain = res$sigma.within.chain,
                                        mu.workers = res$mu.workers,
                                        worker.index = res$worker.index)


return(results)


}    


  
  



  
  
  
  