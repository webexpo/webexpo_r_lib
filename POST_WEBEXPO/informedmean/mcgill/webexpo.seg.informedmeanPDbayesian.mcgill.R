################## WEBEXPO bayesian function for the linmod prior #########################
#
##################  LINMOD PRIOR ###################################
#
#



fun.mcgill.informedmeanPD <-function( 
  
  ## data (observations)
  
  data.sample = NULL, 
  
  #distribution
  is.lognormal = TRUE , 
  
  
  ##prior parameters
  
  mu.mean = 0,
  mu.sd = 1,
  log.sigma.mu = -0.1744 ,
  log.sigma.prec = 2.5523 ,
  mu.lower = -20,
  mu.upper = 20 ,
  
  #nits
  init.mu = log(0.3), 
  init.sigma = log(2.5),
  
  
  #input parameters
  
  oel = 100 ,
  
  #MCMC parameters   (default values recommended for most analyses)
  n.iter = 25000 ,
  n.burnin = 5000 ,
  
  past.data = NULL 
  
  
) {
  

  ###### test
 
 # data.sample = samples[[1]]$data 
  
  #distribution
  #is.lognormal = TRUE  
  
  
  ##prior parameters
  
  #mu.mean = 1
  #mu.sd = 0.2
  #log.sigma.mu = -0.1744 
  #log.sigma.prec = 2.5523 
  
  
  #input parameters
  
  #oel = 100 
  
  #MCMC parameters   (default values recommended for most analyses)
  #n.iter = 50000 
  #n.burnin = 5000
  
  #init.mu = log(0.3)
  #init.sigma = log(2.5)
  
    
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
  
  
          
        #bayesian model
                  
                  
        res <-SEG.informedvar(y = observed_values , 
                                        lt = leftcensored_values ,
                                        gt = rightcensored_values ,
                                        interval.lower = intcensored_left_values ,
                                        interval.upper = intcensored_right_values ,
                                        n.chains = 1 ,
                                        n.iter = n.iter , 
                                        n.burnin = n.burnin ,
                                        n.thin = 1 ,
                                        monitor.burnin = FALSE ,
                                        log.sigma.mu = log.sigma.mu ,
                                        log.sigma.prec = log.sigma.prec , 
                                        init.mu = init.mu, 
                                        init.sigma = init.sigma,
                                        outcome.is.logNormally.distributed = is.lognormal ,
                                        past.data = past.data ,
                                        mu.lower=mu.lower, 
                                        mu.upper=mu.upper,
                                        mu.mean = (mu.mean-log(oel)),
                                        mu.sd = mu.sd,
                                        me.sd.range = numeric(0) , 
                                        cv.range = numeric(0))          
  
  #### final results presentation  
      
  if (is.lognormal==TRUE) results <- list(mu.chain=res$sample$mu+log(oel),sigma.chain=res$sample$sd)
  
  if (is.lognormal==FALSE) results <- list(mu.chain=res$sample$mu,sigma.chain=res$sample$sd)
  

  return(results)
  
  
}    
