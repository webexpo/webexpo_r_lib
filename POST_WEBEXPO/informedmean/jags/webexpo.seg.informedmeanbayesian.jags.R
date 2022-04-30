################## RJAGS bayesian function for INFORMED VAR PRIOR

#############   LOGNORMAL DISTRIBUTION



library(rjags)


fun.jags.informedmean <-function( 
  
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
    
    
    
    
    ## bayesian calculations
    
    ### DATA LIST FOR JAGS
    
    dataList = list( 
      
      y = y ,
      
      n.obs = length(y) ,
      
      CensorType = CensorType ,
      
      censorLimitMat = censorLimitMat
    )
    
    
    
    #### initialising chains for random variables
    
    initsList = list( mu = 0 , log.sigma = 1 )
    
    
    ####JAGS MODEL
    
    mu.mean <-mu.mean-log(oel)
    
    mu.prec <- 1/(mu.sd*mu.sd)
    
    jags.model <-paste("
                          
                     model {
                     
                     
                     ##prior on mu
                     
                     mu ~dnorm(",mu.mean,",",mu.prec,")
                     
                     ##prior on sigma
                     
                     precision <-1/(sigma*sigma)
                     
                     sigma <-exp(log.sigma)
                     
                     log.sigma ~dnorm(",log.sigma.mu,",",log.sigma.prec,")
                     
                     
                     ###likelihood
                     
                     for (i in 1:n.obs) {
                     
                     y[i] ~ dlnorm(mu, precision)
                     
                     CensorType[i] ~ dinterval( y[i] , censorLimitMat[i,] )
                     
                     }
                 
                     }
                     
                     ",sep='')
    
    ### Running the model
    
    jagsModel = jags.model( file=textConnection(jags.model), data=dataList , inits=initsList ,
                            n.chains=1 , n.adapt=500 , quiet=TRUE)
    
    update( jagsModel , n.iter=5000 ) #### Burn-in
    
    codaSamples.M1 = coda.samples( jagsModel , variable.names=c("mu", "sigma") ,
                                   n.iter=n.iter, thin=1 )
    
    
    #####extracting the chains
    
    mu.chain <- codaSamples.M1[[1]][,"mu"]
    sigma.chain <- codaSamples.M1[[1]][,"sigma"] 
    
    
  }
  
  
  if (input.option=="summary") {
    
    
    
    #MODEL
    
    mu.mean <-mu.mean-log(oel)
    
    mu.prec <- 1/(mu.sd*mu.sd)
    
    jags.model <-paste("
                       
                       model {
                       
                       
                     ##prior on mu
                     
                       mu ~dnorm(",mu.mean,",",mu.prec,")
                       
                       ##prior on sigma
                       
                       precision <-1/(sigma*sigma)
                       
                       sigma <-exp(log.sigma)
                       
                       log.sigma ~dnorm(",log.sigma.mu,",",log.sigma.prec,")
                       
                       
                       ###likelihood
                       
                       mean.past ~ dnorm(mu, mean.past.precision)
                       
                       mean.past.precision <- n.past / pow(sigma, 2)
                       
                       ns2 ~ dgamma(ns2.a, b)
                       
                       b <- pow(sigma, -2) / 2
                       
                       }
                       
                       ",sep='')
    
    
    ### DATA LIST FOR JAGS
    
    dataList = list() 
    
    dataList$mean.past <- log.gm-log(oel)
    sd.past <- log.gsd
    dataList$n.past <- n.obs
    dataList$ns2 <- (dataList$n.past - 1) * sd.past^2
    dataList$ns2.a <- (dataList$n.past - 1) / 2
    
    
    
    #### initialising chains for random variables
    
    initsList = list( mu = 0 , log.sigma = 1 )
    
    ### Running the model
    
    jagsModel = jags.model( file=textConnection(jags.model), data=dataList , inits=initsList ,
                            n.chains=1 , n.adapt=500 , quiet=TRUE )
    
    update( jagsModel , n.iter=5000 ) #### Burn-in
    
    codaSamples.M1 = coda.samples( jagsModel , variable.names=c("mu", "sigma") ,
                                   n.iter=n.iter, thin=1 )
    
    
    #####extracting the chains
    
    mu.chain <- codaSamples.M1[[1]][,"mu"]
    sigma.chain <- codaSamples.M1[[1]][,"sigma"]
    
    
  }
  
  
  
  results <- list(mu.chain=as.numeric(mu.chain+log(oel)),sigma.chain=as.numeric(sigma.chain))
  
  
  return(results)
  
  
} 