###################################################################################
#
#
#   Using the INFORMEDMEAN model - EXAMPLE FOR THE LOGNORMAL MODEL
#
#
#
###################################################################################


##### path

setwd( "C:/jerome/Dropbox/bureau/RStudio/Webexpo")


##### sourcing

      ## data preparation 
      
      source("scripts/WEBEXPO OFFICIAL R SCRIPTS/DATA PREPARATION/SEG ANALYSIS/webexpo.seg.dataprep.R")
      
      # MCGILL models - Webexpo folder
      
      source("scripts/WEBEXPO OFFICIAL R SCRIPTS/R MODELS/McGILL FUNCTIONS/data-summary.R")
      
      # MCGILL models - most recent for informedmean - Validation folder
      
      
      source("scripts/WEBEXPO OFFICIAL R SCRIPTS/PostWebexpo additions from McGill/informedmean/mcgill/fcts.R")
      source("scripts/WEBEXPO OFFICIAL R SCRIPTS/PostWebexpo additions from McGill/informedmean/mcgill/model-SEG-informedvar.R")


      ## JAGS wrapp function

      source("scripts/WEBEXPO OFFICIAL R SCRIPTS/PostWebexpo additions from McGill/informedmean/jags/webexpo.seg.informedmeanbayesian.jags.R")

      ## MCGILL wrapp function
      
      source("scripts/WEBEXPO OFFICIAL R SCRIPTS/PostWebexpo additions from McGill/informedmean/mcgill/webexpo.seg.informedmeanbayesian.mcgill.R")

###### EXEMPLE data

      # prior and OEL
      
      mu.mean = -2.64
      mu.sd = 0.28
      oel = 0.1
      
      # Data as observations
      
      x <- c( "0.0923" , "0.0595" , "0.222" , "0.1380" , "0.149" )
      
      # Data as summary
      
      obs.gm <-0.12
      obs.gsd <- 1.8
      obs.n <- 5

      # prior on GSD
      
      gsd <- 2.5
      
      prec <- 10
      
      log.sigma.mu = log( log( 2.5 ) ) 
      log.sigma.prec = 10 

      
######## JAGS function - observations      


res.jags.obs <- fun.jags.informedmean( 
    
    input.option = "observations",
    data.sample = x,
    log.gm = NULL ,
    log.gsd = NULL , 
    n.obs = NULL,
    mu.mean = mu.mean,
    mu.sd = mu.sd,
    log.sigma.mu = log.sigma.mu ,
    log.sigma.prec = log.sigma.prec ,
    oel = oel ,
    n.iter = 25000 ,
    n.burnin = 5000)

exp(median( res.jags.obs$mu.chain))  #0.089    
exp(median( res.jags.obs$sigma.chain)) # 2.14


############# JAGS function - summmary

######## JAGS function - observations      


res.jags.sum <- fun.jags.informedmean( 
  
  input.option = "summary",
  data.sample = NULL,
  log.gm = log( obs.gm ) ,
  log.gsd = log( obs.gsd ) , 
  n.obs = 5,
  mu.mean = mu.mean,
  mu.sd = mu.sd,
  log.sigma.mu = log.sigma.mu ,
  log.sigma.prec = log.sigma.prec ,
  oel = oel ,
  n.iter = 25000 ,
  n.burnin = 5000)

exp(median( res.jags.sum$mu.chain))  #0.089    
exp(median( res.jags.sum$sigma.chain)) # 2.14


######## MCGILL function - observations      


res.mcgill.obs <- fun.mcgill.informedmean( 
  
  input.option = "observations",
  data.sample = x,
  log.gm = NULL ,
  log.gsd = NULL , 
  n.obs = NULL,
  mu.mean = mu.mean,
  mu.sd = mu.sd,
  log.sigma.mu = log.sigma.mu ,
  log.sigma.prec = log.sigma.prec ,
  oel = oel ,
  n.iter = 25000 ,
  n.burnin = 5000)

exp(median( res.mcgill.obs$mu.chain))  #0.089    
exp(median( res.mcgill.obs$sigma.chain)) # 2.14



