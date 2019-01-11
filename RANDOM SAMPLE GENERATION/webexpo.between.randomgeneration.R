##############################################################
#
#   WEBEXPO official R scripts
#
#   BETWEEN WORKER ANALYSIS
#  
#   GENERATION OF RANDOM SAMPLES 
#
#   V1.0 Sept 2018  
#
#############################################################

##
#
## INPUT FOR CREATING RANDOM SAMPLES
#
# model : lognormal or normal
#
# sample size
# proportion of left censored
# proportion of right censored
# proportion of interval censored
# true distributional parameters
# measurement error structure
# between worker correlation

#NOTE

# This algorithms create censored samples with the exact % of data censored as selected.
# to achieve this, some values generated are arbitrarily censored, and therefore the 
# generated sample does not really come from the parent distribution anymore : these algorithms
# are not made to permit evaluation of the performance of censored data treatment options.
# they are used for illustrative purposes.


#######  A. LOGNORMAL


webexpo.between.gener.LN <-function(
                                               # Number of workers     
                                               n.worker = 20,
                                               # Number of days
                                               n.days = rep(20,20), #vector of length n.worker (allows unbalanced data)
                                               # Presence of censoring
                                               no.censoring = TRUE,
                                               # Percentage of left censoring
                                               perc.lowerthan = 20,
                                               # Percentage of right censoring
                                               perc.greaterthan = 10,
                                               # Percentage of interval cenoring
                                               perc.between = 20,
                                               # Ditributional parameters
                                               gm = 0.3,
                                               gsd = 2.5,
                                               rho = 0.5,   ## within worker correlation
                                               
                                               # Censoring parameters
                                               left_factor = 1.5,
                                               right_factor = 1/1.5,
                                               int_factor = 1.5,
                                               
                                               # Measurement error structure
                                               error = "none",   #(or "cv" or "sd")
                                               me.cv = 0.20 ,
                                               me.sd = 0.05) {
  
          #######generating the raw sample
          
                  raw.sample <- gener.sample.V2(mu.y=log(gm),gsd=gsd,rho=rho,n.t=n.worker,n.w=n.days)
          
          ###### adding measurement error
          
                  #### WARNING : error should be defined to minimize the probability of negative values. Sampling will be repeated until no such occurence is found.
                  
                  if (error=="cv")
                    
                  {
                    
                    m.error <- rnorm( sum(n.days) , 0 , gm * me.cv)
                    while (min(m.error+raw.sample$x)<=0) m.error <- rnorm( sum(n.days) , 0 , gm * me.cv)
                    raw.sample$x <- raw.sample$x + m.error
                    
                  }
                  
                  if (error=="sd")
                    
                  {
                    m.error <- rnorm( sum(n.days) , 0 , me.sd)
                    while (min(m.error+raw.sample$x)<=0) m.error <- rnorm( sum(n.days) , 0 , me.sd)
                    raw.sample$x <- raw.sample$x + m.error
                  }
                  
          
        ###### censoring
        #
        # left censoring : x becomes <<x*left_factor 
        # right censoring : x becomes x>x*right_factor
        # interval censoring : x becomes [x/int_factor - x*int_factor]
                  
        ####number of censored = perc*n rounded 
          
          #####left censored chosen on whole data ordered
          #####right censored chosen on whole data ordered
          #####intervall censored chosen on whole data randomly
          
          
          raw.sample$final.x <-as.character(signif(raw.sample$x,3))
          
          if (!no.censoring) {
            
            
            #####ordering data
            
            raw.sample <-raw.sample[order(raw.sample$x,decreasing = FALSE),]
            
            #####left censoring
            
            n.Lcensored <-ceiling(perc.lowerthan*sum(n.days)/100)
            
            if(n.Lcensored!=0) {
              
              raw.sample$final.x[1:n.Lcensored] <- paste("<",signif(raw.sample$x[1:n.Lcensored]*left_factor,3),sep="")
              
            }
            
            #####right censoring
            
            n.Rcensored <-ceiling(perc.greaterthan*sum(n.days)/100)
            
            if(n.Rcensored!=0) {
              
              
              
              raw.sample$final.x[(sum(n.days)-n.Rcensored+1):sum(n.days)] <- paste(">",signif(raw.sample$x[(sum(n.days)-n.Rcensored+1):sum(n.days)]*right_factor,3),sep="")
              
            }
            
            #####intervall censoring
            
            n.Icensored <-ceiling(perc.between*sum(n.days)/100)
            
            if(n.Icensored!=0) {
              
              
              index <-sample((n.Lcensored+1):(sum(n.days)-n.Rcensored), size=n.Icensored,replace = FALSE)
              
              raw.sample$final.x[index] <- paste("[",signif(raw.sample$x[index]/int_factor,3),
                                                 "-",signif(raw.sample$x[index]*int_factor,3),"]",sep="")
              
            }
            
            ########reordering by worker + random
            
            raw.sample <-raw.sample[order(raw.sample$worker),]
            
            new.index <-sample(1:n.days[1],size=n.days[1],replace=FALSE)
            
            for (i in 2:n.worker) new.index <-c(new.index,
                                                sample(sum(n.days[1:(i-1)]):sum(n.days[1:(i)]),size=n.days[i],replace=FALSE))
            
            raw.sample <-raw.sample[new.index,]
            
          }
          
          #######results
          
          results <-list()
          
          ###as a data frame
          
          results <-raw.sample[,c(2,3)]
          
          names(results) <-c("worker","x")
          
          
          return(results)
          
        }


############# B. NORMAL


webexpo.between.gener.N <-function(          
                                              # Number of workers
                                              n.worker = 20,
                                              # Number of days
                                              n.days = rep(20,20), #vector of length n.worker (allows unbalanced data)
                                              # Presence of censoring
                                              no.censoring = TRUE,
                                              # Percentage of left censoring
                                              perc.lowerthan = 20,
                                              # Percentage of right censoring
                                              perc.greaterthan = 10,
                                              # Percentage of interval censoring
                                              perc.between = 20,
                                              # Distributional parameters
                                              mu = 85,
                                              sd = 5,
                                              rho = 0.5, # within worker correlation
                                              
                                              # Censoring parameters
                                              left_factor = 5,
                                              right_factor = 5,
                                              int_factor = 5,
                                              # Measurement error structure
                                              error = "none",   #(or "cv" or "sd")
                                              me.cv = 0.20 ,
                                              me.sd = 0.05) {
  
        #######generating the raw sample
        
        
        raw.sample <- gener.sample.V3(mu=mu,sd=sd,rho=rho,n.t=n.worker,n.w=n.days)
        
        ###### adding measurement error
        
        #### WARNING : error should be defined to minimize the probability of negative values. Sampling will be repeated until no such occurence is found.
        
        if (error=="cv")
          
        {
          
          m.error <- rnorm( sum(n.days) , 0 , mu * me.cv)
          while (min(m.error+raw.sample$x)<=0) m.error <- rnorm( sum(n.days) , 0 , mu * me.cv)
          raw.sample$x <- raw.sample$x + m.error
          
        }
        
        if (error=="sd")
          
        {
          m.error <- rnorm( sum(n.days) , 0 , me.sd)
          while (min(m.error+raw.sample$x)<=0) m.error <- rnorm( sum(n.days) , 0 , me.sd)
          raw.sample$x <- raw.sample$x + m.error
        }
        
        
        ###### censoring
        #
        # left censoring : x becomes <<x*left_factor 
        # right censoring : x becomes x>x*right_factor
        # interval censoring : x becomes [x/int_factor - x*int_factor]
        
        ####number of censored = perc*n rounded 
        
        #####left censored chosen on whole data ordered
        #####right censored chosen on whole data ordered
        #####intervall censored chosen on whole data randomly
        
        
        raw.sample$final.x <-as.character(signif(raw.sample$x,3))
        
        if (!no.censoring) {
          
          
          #####ordering data
          
          raw.sample <-raw.sample[order(raw.sample$x,decreasing = FALSE),]
          
          #####left censoring
          
          n.Lcensored <-ceiling(perc.lowerthan*sum(n.days)/100)
          
          if(n.Lcensored!=0) {
            
            raw.sample$final.x[1:n.Lcensored] <- paste("<",signif(raw.sample$x[1:n.Lcensored]+left_factor,3),sep="")
            
          }
          
          #####right censoring
          
          n.Rcensored <-ceiling(perc.greaterthan*sum(n.days)/100)
          
          if(n.Rcensored!=0) {
            
            
            
            raw.sample$final.x[(sum(n.days)-n.Rcensored+1):sum(n.days)] <- paste(">",signif(raw.sample$x[(sum(n.days)-n.Rcensored+1):sum(n.days)]-right_factor,3),sep="")
            
          }
          
          #####intervall censoring
          
          n.Icensored <-ceiling(perc.between*sum(n.days)/100)
          
          if(n.Icensored!=0) {
            
            
            index <-sample((n.Lcensored+1):(sum(n.days)-n.Rcensored), size=n.Icensored,replace = FALSE)
            
            raw.sample$final.x[index] <- paste("[",signif(raw.sample$x[index]-int_factor,3),
                                               "-",signif(raw.sample$x[index]+int_factor,3),"]",sep="")
            
          }
          
          ########reordering by worker + random
          
          raw.sample <-raw.sample[order(raw.sample$worker),]
          
          new.index <-sample(1:n.days[1],size=n.days[1],replace=FALSE)
          
          for (i in 2:n.worker) new.index <-c(new.index,
                                              sample(sum(n.days[1:(i-1)]):sum(n.days[1:(i)]),size=n.days[i],replace=FALSE))
          
          raw.sample <-raw.sample[new.index,]
          
        }
        
        #######results
        
        results <-list()
        
        ###as a data frame
        
        results <-raw.sample[,c(2,3)]
        
        names(results) <-c("worker","x")
        
        return(results)
        
      }
      




############  C. SUPPORT FUNCTIONS


##### generating identifiers for workers

index.fun.V2 <-function(n.t,n.w=rep(3,n.t)) {   
  
  #n.t number of workers
  #n.w number of days per worker  (n.w is a vector of size n.t)
  
  res <-integer(0)
  
  for (i in 1:n.t) {
    
    res <-c(res,rep(i,n.w[i]))
    
  }
  return(res)
  
}


################# generating data from a hierarchical model - LOGNORMAL


gener.sample.V2 <-function(mu.y,gsd,rho,n.t,n.w=rep(3,n.t),worst=1)   
  
{
  
  #mu.y : mean of logtransformed exposures
  #gsd : geometric standard deviation
  #rho : within worker correlation
  #n.t : number of workers
  ##n.w : number of days per workers (vector of length n.t)
  #worst case: 1-no worst case (random worker within group), 
  #2-50% worst, worker in most exposed half of population
  #3 75%worst,	 worker in most exposed quarter of population 
  
  
  
  worker.effect <-numeric(n.t)
  
  #identifier for each worker
  
  worker.name <-paste(rep('worker',n.t),as.character(1:n.t),sep='-')
  
  #estimation of within and between worker SDs from rho
  
  sigB <-sqrt(rho*(log(gsd))^2)
  sigW <-sqrt((1-rho)*(log(gsd))^2)
  
  #worker effect 
  a<-numeric(1)
  
  ###if sigB!=0
  
  if (sigB!=0) {
    
    #sampling a random worker (with rejection for worst case options)
    for (i in 1:n.t) {	if (worst==1) {a<-rnorm(1,0,sigB)} 
      else if (worst==2) {		a<-rnorm(1,0,sigB)
      while (a<qnorm(0.5,0,sigB)) {a<-rnorm(1,0,sigB)}
      }
      else if (worst==3) {		a<-rnorm(1,0,sigB)
      while (a<qnorm(0.75,0,sigB)) {a<-rnorm(1,0,sigB)}
      }
      worker.effect[i] <-a
    }
  }
  ###if sigB==0				
  else {worker.effect <- c(rep(0,n.t)) }
  
  #repeating each worker effect n..w times						
  worker.effect <-worker.effect[index.fun.V2(n.t,n.w)]
  
  #adding the inter day effects
  sim.sample <-rnorm(sum(n.w),mu.y,sigW)+worker.effect
  
  #results
  
  result <-data.frame(x=exp(sim.sample),worker=worker.name[index.fun.V2(n.t,n.w)],stringsAsFactors=F)
  return(result)
}	



################# generating data from a hierarchical model - NORMAL


gener.sample.V3 <-function(mu,sd,rho,n.t,n.w=rep(3,n.t),worst=1)   
  
{
  
  #mu : mean exposure
  #gd : standard deviation
  #rho : within worker correlation
  #n.t : number of workers
  ##n.w : number of days per workers (vector of length n.t)
  #worst case: 1-no worst case (random worker within group), 
  #2-50% worst, worker in most exposed half of population
  #3 75%worst,	 worker in most exposed quarter of population 
  
  
  
  
  
  worker.effect <-numeric(n.t)
  
  #identifier for each worker
  
  worker.name <-paste(rep('worker',n.t),as.character(1:n.t),sep='-')
  
  #estimation of within and between worker SDs from rho
  
  sigB <-sqrt(rho*sd^2)
  sigW <-sqrt((1-rho)*sd^2)
  
  #worker effect 
  a<-numeric(1)
  
  ###if sigB!=0
  
  if (sigB!=0) {
    
    #sampling a random worker (with rejection for worst case options)
    for (i in 1:n.t) {	if (worst==1) {a<-rnorm(1,0,sigB)} 
      else if (worst==2) {		a<-rnorm(1,0,sigB)
      while (a<qnorm(0.5,0,sigB)) {a<-rnorm(1,0,sigB)}
      }
      else if (worst==3) {		a<-rnorm(1,0,sigB)
      while (a<qnorm(0.75,0,sigB)) {a<-rnorm(1,0,sigB)}
      }
      worker.effect[i] <-a
    }
  }
  ###if sigB==0				
  else {worker.effect <- c(rep(0,n.t)) }
  
  #repeating each worker effect n..w times						
  worker.effect <-worker.effect[index.fun.V2(n.t,n.w)]
  
  #adding the inter day effects
  sim.sample <-rnorm(sum(n.w),mu,sigW)+worker.effect
  
  #results	
  
  result <-data.frame(x=sim.sample,worker=worker.name[index.fun.V2(n.t,n.w)],stringsAsFactors=F)
  return(result)
}	

