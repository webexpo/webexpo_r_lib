################## MCGILL bayesian function for INFORMED VAR PRIOR

#############   LOGNORMAL DISTRIBUTION

library(rjags)


fun.mcgill.informedmean <-function( 
  
  input.option = "observations" , ## OR "summary"
  
  
  ## data (observations)
  
  data.sample , 
  
  ## data (summary)

  log.gm ,
  
  log.gsd , 
  
  n.obs ,
  
  #distribution
  is.lognormal = TRUE , 
  
  
  ##prior parameters
  
  mu.mean = 0,
  mu.sd = 1,
  log.sigma.mu = -0.1744 ,
  log.sigma.prec = 2.5523 ,
  
  
  #input parameters
  
  oel = 100 ,
  
  #MCMC parameters   (default values recommended for most analyses)
  n.iter = 25000 ,
  n.burnin = 5000
  
  
) {
  
  if (input.option=="observations") {
    
    
    #input preparation
    
    formatted.data <-webexpo.seg.datapreparation(data.in = data.sample)
    
    x <- formatted.data$data
    
    
    #input preparation (dealing with censorship)
    
    formatted.data <-webexpo.seg.datapreparation(data.in = data.sample)
    
    x <- formatted.data$data
    
    
    observed_values <-as.numeric(x[formatted.data$notcensored])/oel
    
    #left censored values
    
    leftcensored_values <-as.numeric(substring(x[formatted.data$leftcensored],2))/oel
    
    #right censored values
    
    rightcensored_values <-as.numeric(substring(x[formatted.data$rightcensored],2))/oel
    
    #interval censored values
    
    intcensored_left_values <-as.numeric(substring(x[formatted.data$intcensored],2,regexpr("-",x[formatted.data$intcensored],fixed=TRUE)-1))/oel
    
    intcensored_right_values <-as.numeric(substring(x[formatted.data$intcensored],regexpr("-",x[formatted.data$intcensored],fixed=TRUE)+1,nchar(x[formatted.data$intcensored])-1))/oel
    
    ########bayesian model
    
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
                          init.mu = 0, 
                          init.sigma = log(2.5),
                          outcome.is.logNormally.distributed = is.lognormal ,
                          past.data = NULL ,
                          mu.lower= -20, 
                          mu.upper= 20,
                          mu.mean = mu.mean-log(oel),
                          mu.sd = mu.sd,
                          me.sd.range = numeric(0) , 
                          cv.range = numeric(0))          
    
    #### final results presentation  
    
    mu.chain <- res$sample$mu+log(oel)
    sigma.chain <- res$sample$sd
    
    
    
  }
  
  
  if (input.option=="summary") {
    
    
    ########bayesian model
    
    res <-SEG.informedvar(y = NULL , 
                          lt = NULL ,
                          gt = NULL ,
                          interval.lower = NULL ,
                          interval.upper = NULL ,
                          n.chains = 1 ,
                          n.iter = n.iter , 
                          n.burnin = n.burnin ,
                          n.thin = 1 ,
                          monitor.burnin = FALSE ,
                          log.sigma.mu = log.sigma.mu ,
                          log.sigma.prec = log.sigma.prec , 
                          init.mu = 0, 
                          init.sigma = log(2.5),
                          outcome.is.logNormally.distributed = is.lognormal ,
                          past.data = list(mean = log.gm -log(oel), sd = log.gsd, n = n.obs ) ,
                          mu.lower= -20, 
                          mu.upper= 20,
                          mu.mean = mu.mean-log(oel),
                          mu.sd = mu.sd,
                          me.sd.range = numeric(0) , 
                          cv.range = numeric(0))          
    
    #### final results presentation  
    
    mu.chain <- res$sample$mu+log(oel)
    sigma.chain <- res$sample$sd
    
    
    
  }
  
  
  
  results <- list(mu.chain=as.numeric(mu.chain),sigma.chain=as.numeric(sigma.chain))
  
  
  return(results)
  
  
} 