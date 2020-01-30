##############################################################
#
#   WEBEXPO official R scripts
#
#   BETWEEN WORKER ANALYSIS
#  
#   OFFICIAL REPORT EXAMPLES
#
#   V1.0 Sept 2018    
#
#############################################################


###### library

library(rjags)

library(here)

#####  sourcing files

      ## setting path
      
      
      ## data preparation / generation
      
      source(here("RANDOM SAMPLE GENERATION", "webexpo.between.randomgeneration.R"))
      
      source(here("DATA PREPARATION", "BETWEEN WORKER ANALYSIS", "webexpo.between.dataprep.R"))
      
      # JAGS models
      
      source(here("JAGS MODELS", "BETWEEN WORKER ANALYSIS", "webexpo.between.mainbayesian.R"))
      
      source(here("JAGS MODELS", "BETWEEN WORKER ANALYSIS", "webexpo.between.informedvarbayesian.R"))
      
      source(here("JAGS MODELS", "BETWEEN WORKER ANALYSIS", "webexpo.between.informedvarbayesian.models.R"))
      
      # McGill models
      
      source(here("R MODELS", "WRAPPING FUNCTIONS", "webexpo.between.mainbayesian.mcgill.R"))
      
      source(here("R MODELS", "McGILL FUNCTIONS", "data-summary.R"))
      
      source(here("R MODELS", "McGILL FUNCTIONS", "fcts.R"))
      
      source(here("R MODELS", "McGILL FUNCTIONS", "model-Between-worker.R"))
      
      # Data interpretation
      
      source(here("RESULT INTERPRETATION", "BETWEEN WORKER ANALYSIS", "webexpo.between.interpretation.R"))
      
      source(here("RESULT INTERPRETATION", "BETWEEN WORKER ANALYSIS", "webexpo.between.summary.R"))
      
      source(here("RESULT INTERPRETATION", "SEG ANALYSIS", "webexpo.seg.interpretation.R"))
      
      source(here("RESULT INTERPRETATION", "SEG ANALYSIS", "webexpo.seg.summary.R"))


      
####### Illustration of low and high within worker correlation -  TABLE 8      
      
            
            ##### creation of random samples


            #Note : code disabled below. Find the actual samples in the associated EXCEL FILES or Report appendix
            
            ### lognormal sample low between worker variability : rho = 0.06 (25th percentile in Kromhout et al)
            
            #sample.1 <- webexpo.between.gener.LN(n.worker = 10,
            #                                                 n.days = rep(10,10), 
            #                                                 no.censoring = TRUE,
            #                                                 gm = 30,
            #                                                 gsd = 2.5,
            #                                                 rho = 0.06, 
            #                                                 error = "none")
            
            ### lognormal sample low between worker variability : rho = 0.66 (75th percentile in Kromhout et al)
            
            #sample.2 <- webexpo.between.gener.LN(n.worker = 10,
            #                                     n.days = rep(10,10), 
            #                                     no.censoring = TRUE,
            #                                     gm = 30,
            #                                     gsd = 2.5,
            #                                     rho = 0.66, 
            #                                     error = "none")    
            
  
          # Bayesian calculations      
      
          #####default webexpo analysis - JAGS
          
          mcmc.1 <- Webexpo.between.globalbayesian.jags(  data.sample = sample.1,
                                                         is.lognormal = TRUE,
                                                         error.type = "none", 
                                                         me.range = c(0.3,0.3), 
                                                         oel = 150,
                                                         prior.model = "informedvar",
                                                         mu.overall.lower = -20,
                                                         mu.overall.upper = 20,
                                                         log.sigma.between.mu=-0.8786,
                                                         log.sigma.between.prec=1.634,
                                                         log.sigma.within.mu=-0.4106,
                                                         log.sigma.within.prec=1.9002, 
                                                         n.iter = 50000, 
                                                         n.burnin = 5000, 
                                                         init.mu.overall = log(0.3),
                                                         init.sigma.within = log(2.5),
                                                         init.sigma.between = log(2.5),
                                                         init.log.sigma.within = log(log(2.5)),
                                                         init.log.sigma.between = log(log(2.5)))
          
          mcmc.2 <- Webexpo.between.globalbayesian.jags(  data.sample = sample.2,
                                                         is.lognormal = TRUE,
                                                         error.type = "none", 
                                                         me.range = c(0.3,0.3), 
                                                         oel = 150,
                                                         prior.model = "informedvar",
                                                         mu.overall.lower = -20,
                                                         mu.overall.upper = 20,
                                                         log.sigma.between.mu=-0.8786,
                                                         log.sigma.between.prec=1.634,
                                                         log.sigma.within.mu=-0.4106,
                                                         log.sigma.within.prec=1.9002, 
                                                         n.iter = 50000, 
                                                         n.burnin = 5000, 
                                                         init.mu.overall = log(0.3),
                                                         init.sigma.within = log(2.5),
                                                         init.sigma.between = log(2.5),
                                                         init.log.sigma.within = log(log(2.5)),
                                                         init.log.sigma.between = log(log(2.5)) )
          
          
          ###### numerical interpretation
          
          
          num.1 <- Webexpo.between.interpretation(  is.lognormal = TRUE , 
                                                    mu.overall = mcmc.1$mu.overall.chain,
                                                    sigma.between.chain = mcmc.1$sigma.between.chain ,
                                                    sigma.within.chain = mcmc.1$sigma.within.chain ,
                                                    probacred = 90, 
                                                    oel = 150,
                                                    target_perc = 95,
                                                    rappap_cover = 80, 
                                                    wwct = 0.2 , 
                                                    prob.ind.overexpo.threshold = 20)
          
          num.2 <- Webexpo.between.interpretation(  is.lognormal = TRUE , 
                                                    mu.overall = mcmc.2$mu.overall.chain,
                                                    sigma.between.chain = mcmc.2$sigma.between.chain ,
                                                    sigma.within.chain = mcmc.2$sigma.within.chain ,
                                                    probacred = 90, 
                                                    oel = 150,
                                                    target_perc = 95,
                                                    rappap_cover = 80, 
                                                    wwct = 0.2 , 
                                                    prob.ind.overexpo.threshold = 20)
          # Summaries
          
          sum.seg.T8 <- webexpo.between.summary( labels = c("Low rho" , "High rho") ,
                                                result.list = list(num.1,num.2) ,
                                                is.lognormal = TRUE )
          #Saving results
          
          write.csv2(sum.seg.T8,here("EXAMPLES", "output", "table8.csv"))
          
          
######## exposure metrics for the least and most exposed workers in sample 1 and sample 2 : TABLE9
      
      #### least and most exposed in sample 1 : #5 and #8   
      #### least and most exposed in sample 2 : #7 and #3
      
      
      # Numerical interpretation (worker specific)
              
      least.S1 <-Webexpo.between.interpretation.worker(worker.id  = "worker-5", 
                                                       is.lognormal = TRUE , 
                                                       mu.workers = mcmc.1$mu.workers,
                                                       worker.index = mcmc.1$worker.index, 
                                                       sigma.within.chain = mcmc.1$sigma.within.chain ,
                                                       probacred = 90, 
                                                       oel = 150,
                                                       frac_threshold = 5,
                                                       target_perc = 95 )
      most.S1 <-Webexpo.between.interpretation.worker(worker.id  = "worker-8", 
                                                      is.lognormal = TRUE , 
                                                      mu.workers = mcmc.1$mu.workers,
                                                      worker.index = mcmc.1$worker.index, 
                                                      sigma.within.chain = mcmc.1$sigma.within.chain ,
                                                      probacred = 90, 
                                                      oel = 150,
                                                      frac_threshold = 5,
                                                      target_perc = 95 ) 
      
      least.S2 <-Webexpo.between.interpretation.worker(worker.id  = "worker-7", 
                                                       is.lognormal = TRUE , 
                                                       mu.workers = mcmc.2$mu.workers,
                                                       worker.index = mcmc.2$worker.index, 
                                                       sigma.within.chain = mcmc.2$sigma.within.chain ,
                                                       probacred = 90, 
                                                       oel = 150,
                                                       frac_threshold = 5,
                                                       target_perc = 95 )
      
      most.S2 <-Webexpo.between.interpretation.worker(worker.id  = "worker-3", 
                                                      is.lognormal = TRUE , 
                                                      mu.workers = mcmc.2$mu.workers,
                                                      worker.index = mcmc.2$worker.index, 
                                                      sigma.within.chain = mcmc.2$sigma.within.chain ,
                                                      probacred = 90, 
                                                      oel = 150,
                                                      frac_threshold = 5,
                                                      target_perc = 95 )
      
      # Summaries
      
      sum.seg.T9 <- webexpo.seg.summary( labels = c("least.s1","most.s1","least.s2","most.s2") ,  
                                        result.list = list(least.S1,most.S1,least.S2,most.S2)
                                        , is.lognormal = TRUE )
      
      
      #Saving results
      
      write.csv2(sum.seg.T9,here("EXAMPLES", "output", "table9.csv")     )
      
      
      
### lognormal sample average between worker variability : rho = 0.22 (50th percentile in Kromhout et al), typical sample size TABLE 10
      
      #Note : code disabled below. Find the actual samples in the associated EXCEL FILES or Report appendix
      
      #sample.3 <- webexpo.between.gener.LN(n.worker = 3,
      #                                     n.days = rep(4,3), 
      #                                     no.censoring = TRUE,
      #                                     gm = 30,
      #                                     gsd = 2.5,
      #                                     rho = 0.22, 
      #                                     error = "none")      
      #
      
      
      #Bayesian models
      
      mcmc.3 <- Webexpo.between.globalbayesian.jags( data.sample = sample.3,
                                                     is.lognormal = TRUE,
                                                     error.type = "none", 
                                                     me.range = c(0.3,0.3), 
                                                     oel = 150,
                                                     prior.model = "informedvar",
                                                     mu.overall.lower = -20,
                                                     mu.overall.upper = 20,
                                                     log.sigma.between.mu=-0.8786,
                                                     log.sigma.between.prec=1.634,
                                                     log.sigma.within.mu=-0.4106,
                                                     log.sigma.within.prec=1.9002, 
                                                     n.iter = 50000, 
                                                     n.burnin = 5000, 
                                                     init.mu.overall = log(0.3),
                                                     init.sigma.within = log(2.5),
                                                     init.sigma.between = log(2.5),
                                                     init.log.sigma.within = log(log(2.5)),
                                                     init.log.sigma.between = log(log(2.5)))
      
      # Numerical interpretation
      
      num.3 <- Webexpo.between.interpretation(  is.lognormal = TRUE , 
                                                mu.overall = mcmc.3$mu.overall.chain,
                                                sigma.between.chain = mcmc.3$sigma.between.chain ,
                                                sigma.within.chain = mcmc.3$sigma.within.chain ,
                                                probacred = 90, 
                                                oel = 150,
                                                target_perc = 95,
                                                rappap_cover = 80, 
                                                wwct = 0.2 , 
                                                prob.ind.overexpo.threshold = 20)
      # Summaries
      
      sum.seg.T10 <- webexpo.between.summary( labels = c("realistic sample size") ,
                                             result.list = list(num.3) ,
                                             is.lognormal = TRUE )
      
      #Saving results
      
      write.csv(sum.seg.T10,here("EXAMPLES", "output", "table10.csv"))

         
      
############## the normal case - TABLE 11
      
      #Note : code disabled below. Find the actual samples in the associated EXCEL FILES or Report appendix
      
      #sample.4 <- webexpo.between.gener.N(n.worker = 10,
      #                                    n.days = rep(10,10), 
      #                                    no.censoring = TRUE,
      #                                    mu = 80,
      #                                    sd = 5,
      #                                    rho = 0.22,   
      #                                    error = "none")
      
      # Bayesian models
      
      mcmc.4 <- Webexpo.between.globalbayesian.jags( data.sample = sample.4 ,
                                                     is.lognormal = FALSE , 
                                                     error.type = "none" ,
                                                     oel = 85,
                                                     mu.overall.lower = 40,
                                                     mu.overall.upper = 125,
                                                     log.sigma.between.mu=1.099,
                                                     log.sigma.between.prec=1.191,
                                                     log.sigma.within.mu=1.099,
                                                     log.sigma.within.prec=1.191, 
                                                     init.mu.overall = 85,
                                                     init.sigma.within = 3,
                                                     init.sigma.between = 3,
                                                     init.log.sigma.within = 1.099,
                                                     init.log.sigma.between = 1.099)
      
      # Numerical interpretation
      
      num.4 <- Webexpo.between.interpretation(  is.lognormal = FALSE , 
                                                mu.overall = mcmc.4$mu.overall.chain,
                                                sigma.between.chain = mcmc.4$sigma.between.chain ,
                                                sigma.within.chain = mcmc.4$sigma.within.chain ,
                                                probacred = 90, 
                                                oel = 85,
                                                target_perc = 95,
                                                rappap_cover = 80, 
                                                wwct = 0.2 , 
                                                prob.ind.overexpo.threshold = 20)
      

      
      # Summary
      
      sum.seg.T11 <- webexpo.between.summary( labels = c("normal case") ,
                                             result.list = list(num.4) ,
                                             is.lognormal = FALSE )     
      
      #Saving results
      
      write.csv(sum.seg.T11,here("EXAMPLES", "output", "table11.csv"))
         

