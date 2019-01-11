##############################################################
#
#   WEBEXPO official R scripts
#
#   SEG ANALYSIS
#  
#   BAYESIAN MODELS IN JAGS 
#
#   Global function
#
#   V1.0 Sept 2018  
#   
#
#############################################################


## required libraries

library(rjags)


#required scripts
# The following scripts must be sourced prior to sourcing this script

          #data preparation  
          
            #   source("DATA PREPARATION/SEG ANALYSIS/webexpo.seg.dataprep.R")
          
          #bayesian functions
          
            #    source("JAGS MODELS/SEG ANALYSIS/webexpo.seg.informedvarbayesian.R")
          
            #    source("JAGS MODELS/SEG ANALYSIS/webexpo.seg.informedvarbayesian.models.R")
          
            #   source("JAGS MODELS/SEG ANALYSIS/webexpo.seg.uninformativebayesian.R")
          
            #   source("JAGS MODELS/SEG ANALYSIS/webexpo.seg.uninformativebayesian.models.R")
          
            #  source("JAGS MODELS/SEG ANALYSIS/webexpo.seg.riskbandbayesian.R")
          
            #  source("JAGS MODELS/SEG ANALYSIS/webexpo.seg.riskbandbayesian.models.R")
          


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
# 
##



Webexpo.seg.globalbayesian.jags <- function(
  
                      data.sample = c("27.8","41","29","<36.3","<41.7","<10.8","54.4","30","[40.1-90.2]","200") ,

                      #distribution
                      is.lognormal = TRUE , 
                      
                      #measurement error model options
                      error.type = "none" ,
                      me.range = c(0.3,0.3) , 
                      
                      #input parameters
                      oel = 100 ,
                      
                      #prior information
                      prior.model = "informedvar" ,  ## options : ("informedvar"  / "uninformative" / "riskband" )
                      
                      #prior informedvar (default values valid for the lognormal distribution)
                      mu.lower = -20 ,
                      mu.upper = 20  ,
                      log.sigma.mu = -0.1744 ,
                      log.sigma.prec = 2.5523 , 
                      
                      #prior uninformative (default values valid for the lognormal distribution)
                      mu.lower.uninf = -20 ,
                      mu.upper.uninf = 20 ,
                      sd.range = c(0, 2.3) ,
                      
                      #past.data
                      past.data = numeric(0) ,
                      
                      #prior riskband 
                      # mu and sigma are the log transformed gm and gsds if is.lognormal=TRUE
                      mu.lower.riskb = log(0.001) ,
                      mu.upper.riskb = log(1000) ,
                      sigma.lower = log(1.05) ,
                      sigma.upper = log(10) ,
                      
                      #params riskband
                      target_perc = 95 ,
                      A = c(0.01,0.1,0.5,1) ,
                      region.prior.prob = rep(0.2, 5) , 
                        
                      #MCMC parameters  
                      n.iter = 25000 ,
                      n.burnin = 5000 , 
                      
                      init.mu = log(0.3) ,  #(default values valid for the lognormal distribution)
                      init.sigma = log(2.5)
                                              ) {



  
#input preparation (dealing with censorship)

      formatted.data <-webexpo.seg.datapreparation(data.in = data.sample)
      
      x <- formatted.data$data
        
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
          #left limit of censorLimitMat is -50 (value we are sure any observation is greater than)
          #right limit of censorLimitMat is 50 (value we are sure any observation is smaller than)                
          
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


## preparation of past.data in case of lognormal model

      if (length(past.data)==3 & is.lognormal) past.data$mean <- past.data$mean-log(oel)

## preparation of measurement error when lognormal model and SD measurement error

      if (error.type=="SD" & is.lognormal) me.range <-me.range/oel
  
## bayesian calculations

      if (prior.model == "informedvar") {
        
        
        res <-fun.jags.informedvar( y=y ,
                                CensorType = CensorType ,
                                censorLimitMat = censorLimitMat ,
                                past.data = past.data,
                                is.lognormal=is.lognormal,
                                error.type = error.type,
                                me.range = me.range,
                                n.iter = n.iter , 
                                n.burnin = n.burnin ,
                                mu.lower = mu.lower ,
                                mu.upper = mu.upper ,
                                log.sigma.mu = log.sigma.mu ,
                                log.sigma.prec = log.sigma.prec ,
                                init.mu = init.mu ,
                                init.sigma = init.sigma ) 
        
        
      }
      
      if (prior.model == "uninformative") {
        
        
        res <-fun.jags.uninformative( y=y ,
                                    CensorType = CensorType ,
                                    censorLimitMat = censorLimitMat ,
                                    is.lognormal=is.lognormal,
                                    error.type = error.type,
                                    me.range = me.range,
                                    n.iter = n.iter , 
                                    n.burnin = n.burnin ,
                                    mu.lower.uninf = mu.lower.uninf, 
                                    mu.upper.uninf = mu.upper.uninf,
                                    sd.range = sd.range ,
                                    init.mu = init.mu ,
                                    init.sigma = init.sigma ) 
        
        
      }


      if (prior.model == "riskband") {
        
        
        res <-fun.jags.riskband( y=y ,
                                    CensorType = CensorType ,
                                    censorLimitMat = censorLimitMat ,
                                    is.lognormal=is.lognormal,
                                    error.type = error.type,
                                    me.range = me.range,
                                    n.iter = n.iter , 
                                    n.burnin = n.burnin ,
                                    mu.lower.riskb = mu.lower.riskb ,
                                    mu.upper.riskb = mu.upper.riskb ,
                                    sigma.lower = sigma.lower,
                                    sigma.upper = sigma.upper,
                                    target_perc = target_perc ,
                                    A = A ,
                                    region.prior.prob = region.prior.prob  ) 
        
        
      }
      
      
## final chain preparation

      if (is.lognormal) results <-list(mu.chain=res$mu.chain+log(oel),sigma.chain=res$sigma.chain)
      
      if (!is.lognormal) results <-list(mu.chain=res$mu.chain,sigma.chain=res$sigma.chain)


return(results)


}    


  
  



  
  
  
  