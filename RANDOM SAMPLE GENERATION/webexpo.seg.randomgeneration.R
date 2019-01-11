##############################################################
#
#   WEBEXPO official R scripts
#
#   SEG ANALYSIS
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

#NOTE

# This algorithms create censored samples with the exact % of data censored as selected.
# to achieve this, some values generated are arbitrarily censored, and therefore the 
# generated sample does not really come from the parent distribution anymore : these algorithms
# are not made to permit evaluation of the performance of censored data treatment options.
# they are used for illustrative purposes.


#######  A. LOGNORMAL

      webexpo.seg.gener.LN <-function(
                                    # Sample size
                                    n = 10,
                                    # indicator for the presence of censoring
                                    no.censoring = FALSE,
                                    # Percentage Left censored 
                                    perc.lowerthan = 20,
                                    # Percentage right censored
                                    perc.greaterthan = 10,
                                    # Percentage interval censored
                                    perc.between = 20,
                                    # Distributional parameters
                                    gm = 0.3,
                                    gsd = 2.5,
                                    # Censoring parameters
                                    left_factor = 1.5,
                                    right_factor = 1/1.5,
                                    int_factor = 1.5,
                                    # Measurement error structure
                                    error = "none",   #(or "cv" or "sd")
                                    me.cv = 0.20 ,
                                    me.sd = 0.05) {
        
        #######generating the raw sample
        
        raw.sample <- exp(rnorm(n,log(gm),log(gsd)))
        
        
        ###### adding measurement error
        
        #### WARNING : error should be defined to minimize the probability of negative values. Sampling will be repeated until no such occurence is found.
        
        if (error=="cv")
          
        {
              
              for (i in 1:n) {
          
              m.error <- rnorm( 1 , 0 , raw.sample[i] * me.cv)
              
              while (min(m.error+raw.sample[i])<=0) m.error <- rnorm( 1 , 0 , raw.sample[i] * me.cv)
              
              raw.sample[i] <- raw.sample[i] + m.error
              }
        }
        
        if (error=="sd")
          
        {
          for (i in 1:n) {
          
          m.error <- rnorm( 1 , 0 , me.sd )
          
          while (min(m.error+raw.sample[i])<=0) m.error <- rnorm( 1 , 0 , me.sd)
          
          raw.sample[i] <- raw.sample[i] + m.error }
        }
        
        ###### censoring
        #
        # left censoring : x becomes <<x*left_factor 
        # right censoring : x becomes x>x*right_factor
        # interval censoring : x becomes [x/int_factor - x*int_factor]
        
        ####number of censored = perc*n rounded up
        
        #####left censored chosen on whole data ordered
        #####right censored chosen on whole data ordered
        #####intervall censored chosen on whole data randomly
        
        
        samplev1 <-as.character(signif(raw.sample,3))
        
        if (!no.censoring) {
          
          
          #####ordering data
          
          raw.sample <-raw.sample[order(raw.sample,decreasing = FALSE)]
          
          
          ######insuring a minimum of 3 uncensored
          
          if (floor((100-perc.lowerthan-perc.greaterthan-perc.between)*n/100)<3) {
            
            
            perc.greaterthan <-0
            perc.between <-0
            perc.lowerthan <- 100*(n-3)/n
          }
          
          
          
          #####left censoring
          
          n.Lcensored <-round(perc.lowerthan*sum(n)/100)
          
          if(n.Lcensored!=0) {
            
            samplev1[1:n.Lcensored] <- paste("<",signif(raw.sample[1:n.Lcensored]*left_factor,3),sep="")
            
          }
          
          #####right censoring
          
          n.Rcensored <-round(perc.greaterthan*sum(n)/100)
          
          if(n.Rcensored!=0) {
            
            
            
            samplev1[(sum(n)-n.Rcensored+1):sum(n)] <- paste(">",signif(raw.sample[(sum(n)-n.Rcensored+1):sum(n)]*right_factor,3),sep="")
            
          }
          
          #####intervall censoring
          
          n.Icensored <-round(perc.between*sum(n)/100)
          
          if(n.Icensored!=0) {
            
            
            index <-sample((n.Lcensored+1):(sum(n)-n.Rcensored), size=n.Icensored,replace = FALSE)
            
            samplev1[index] <- paste("[",signif(raw.sample[index]/int_factor,3),
                                     "-",signif(raw.sample[index]*int_factor,3),"]",sep="")
            
          }
          
          ########reordering  randomly
          
          new.index <-sample(1:n,size=n,replace=FALSE)
          
          samplev1 <-samplev1[new.index]
          
        }
        
        #######results
        
        results <- samplev1
        
        return(results)
        
        
      }

      
############# B. NORMAL

      
      webexpo.seg.gener.N <-function(
                                      # Sample size    
                                      n = 10,
                                      # Presense of censoring
                                      no.censoring = FALSE,
                                      # Percantage left censored
                                      perc.lowerthan = 20,
                                      # Percentage right censored
                                      perc.greaterthan = 10,
                                      # Percentage interval censored
                                      perc.between = 20,
                                      # distributional parameters
                                      mu = 85,
                                      sigma = 5,
                                      # Censoring parameters
                                      left_factor = 5,
                                      right_factor = 5,
                                      int_factor = 5,
                                      # Measurement error structure
                                      error = "none",   #(or "cv" or "sd")
                                      me.cv = 0.20 ,
                                      me.sd = 0.05) {
        
        #######generating the raw sample
        
        raw.sample <- rnorm(n,mu,sigma)
        
        ###### adding measurement error
        
        if (error=="cv")
          
        {
          
          for (i in 1:n) {
            
            m.error <- rnorm( 1 , 0 , raw.sample[i] * me.cv)
            
            while (min(m.error+raw.sample[i])<=0) m.error <- rnorm( 1 , 0 , raw.sample[i] * me.cv)
            
            raw.sample[i] <- raw.sample[i] + m.error
          }
        }
        
        if (error=="sd")
          
        {
          for (i in 1:n) {
            
            m.error <- rnorm( 1 , 0 , me.sd )
            
            while (min(m.error+raw.sample[i])<=0) m.error <- rnorm( 1 , 0 , me.sd)
            
            raw.sample[i] <- raw.sample[i] + m.error }
        }
        ###### Censoring
        #
        # left censoring : x becomes <x+left_factor 
        # right censoring : x becomes x>x+right_factor
        # interval censoring : x becomes [x-int_factor - x+int_factor]
        
        ####number of censored = perc*n rounded up
        
        #####left censored chosen on whole data ordered
        #####right censored chosen on whole data ordered
        #####intervall censored chosen on whole data randomly
        
        
        samplev1 <-as.character(signif(raw.sample,3))
        
        if (!no.censoring) {
          
          
          #####ordering data
          
          raw.sample <-raw.sample[order(raw.sample,decreasing = FALSE)]
          
          ######insuring a minimum of 3 uncensored
          
          
          if (floor((100-perc.lowerthan-perc.greaterthan-perc.between)*n/100)<3) {
            
            
            perc.greaterthan <-0
            perc.between <-0
            perc.lowerthan <- 100*(n-3)/n
          }
          
          #####left censoring
          
          n.Lcensored <-round(perc.lowerthan*sum(n)/100)
          
          if(n.Lcensored!=0) {
            
            samplev1[1:n.Lcensored] <- paste("<",signif(raw.sample[1:n.Lcensored]+left_factor,3),sep="")
            
          }
          
          #####right censoring
          
          n.Rcensored <-round(perc.greaterthan*sum(n)/100)
          
          if(n.Rcensored!=0) {
            
            
            
            samplev1[(sum(n)-n.Rcensored+1):sum(n)] <- paste(">",signif(raw.sample[(sum(n)-n.Rcensored+1):sum(n)]-right_factor,3),sep="")
            
          }
          
          #####intervall censoring
          
          n.Icensored <-round(perc.between*sum(n)/100)
          
          if(n.Icensored!=0) {
            
            
            index <-sample((n.Lcensored+1):(sum(n)-n.Rcensored), size=n.Icensored,replace = FALSE)
            
            samplev1[index] <- paste("[",signif(raw.sample[index]-int_factor,3),
                                     "-",signif(raw.sample[index]+int_factor,3),"]",sep="")
            
          }
          
          ########reordering randomly
          
          new.index <-sample(1:n,size=n,replace=FALSE)
          
          samplev1 <-samplev1[new.index]
          
        }
        
        #######results
        
        results <- samplev1
        
        return(results)
        
        
      }
      
      
      
      
      
