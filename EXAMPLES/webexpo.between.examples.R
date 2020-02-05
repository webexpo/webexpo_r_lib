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


      
####### Illustration of low and high within worker correlation -  TABLE 6      
      
            
            ##### creation of random samples


            #Note : code disabled below. Find the actual samples in the associated EXCEL FILES or Report appendix
            
            ### lognormal sample low between worker variability : rho = 0.06 (25th percentile in Kromhout et al)
            
            #sample.3 <- webexpo.between.gener.LN(n.worker = 10,
            #                                                 n.days = rep(10,10), 
            #                                                 no.censoring = TRUE,
            #                                                 gm = 30,
            #                                                 gsd = 2.5,
            #                                                 rho = 0.06, 
            #                                                 error = "none")
            
            ### lognormal sample low between worker variability : rho = 0.66 (75th percentile in Kromhout et al)
            
            #sample.4 <- webexpo.between.gener.LN(n.worker = 10,
            #                                     n.days = rep(10,10), 
            #                                     no.censoring = TRUE,
            #                                     gm = 30,
            #                                     gsd = 2.5,
            #                                     rho = 0.66, 
            #                                     error = "none")    
            
  
          sample.3 <- data.frame(rbind(c("worker-01", "185"), c("worker-01", "34.8"), c("worker-01", "16.7"), c("worker-01", "12.4"), c("worker-01", "18.6"), c("worker-01", "47.4"), c("worker-01", "52.6"), c("worker-01", "15.3"), c("worker-01", "27.6"), c("worker-01", "26.3"), c("worker-02", "4.79"), c("worker-02", "23"), c("worker-02", "7.54"), c("worker-02", "62.3"), c("worker-02", "8.55"), c("worker-02", "9.28"), c("worker-02", "43.6"), c("worker-02", "94.2"), c("worker-02", "44.6"), c("worker-02", "66.6"), c("worker-03", "8.85"), c("worker-03", "31.7"), c("worker-03", "15.8"), c("worker-03", "89.6"), c("worker-03", "164"), c("worker-03", "40.5"), c("worker-03", "47.6"), c("worker-03", "75.5"), c("worker-03", "10.7"), c("worker-03", "62.3"), c("worker-04", "16.4"), c("worker-04", "6.91"), c("worker-04", "87.4"), c("worker-04", "20"), c("worker-04", "16.8"), c("worker-04", "7.12"), c("worker-04", "6.99"), c("worker-04", "16.4"), c("worker-04", "12.6"), c("worker-04", "63.9"), c("worker-05", "14.7"), c("worker-05", "59.6"), c("worker-05", "15"), c("worker-05", "21.8"), c("worker-05", "20.6"), c("worker-05", "96.1"), c("worker-05", "16.8"), c("worker-05", "15.8"), c("worker-05", "8.02"), c("worker-05", "26.7"), c("worker-06", "37.9"), c("worker-06", "96.9"), c("worker-06", "40.8"), c("worker-06", "106"), c("worker-06", "21.7"), c("worker-06", "25.8"), c("worker-06", "51.3"), c("worker-06", "23"), c("worker-06", "18.9"), c("worker-06", "20.2"), c("worker-07", "22"), c("worker-07", "44.8"), c("worker-07", "37.5"), c("worker-07", "16.6"), c("worker-07", "30.7"), c("worker-07", "7.07"), c("worker-07", "7.18"), c("worker-07", "80.9"), c("worker-07", "44.5"), c("worker-07", "135"), c("worker-08", "69.9"), c("worker-08", "30.5"), c("worker-08", "33.4"), c("worker-08", "53"), c("worker-08", "70.7"), c("worker-08", "78.3"), c("worker-08", "18"), c("worker-08", "45.2"), c("worker-08", "51.4"), c("worker-08", "33.7"), c("worker-09", "28.1"), c("worker-09", "7.49"), c("worker-09", "16"), c("worker-09", "23"), c("worker-09", "99.9"), c("worker-09", "12"), c("worker-09", "11.8"), c("worker-09", "57.4"), c("worker-09", "8.79"), c("worker-09", "24"), c("worker-10", "113"), c("worker-10", "7.68"), c("worker-10", "85.6"), c("worker-10", "196"), c("worker-10", "35"), c("worker-10", "17.6"), c("worker-10", "60.7"), c("worker-10", "15.5"), c("worker-10", "34.3"), c("worker-10", "12.1")), stringsAsFactors = F)
          sample.4 <- data.frame(rbind(c("worker-01", "66.8"), c("worker-01", "46"), c("worker-01", "61.1"), c("worker-01", "54.6"), c("worker-01", "31.7"), c("worker-01", "74.3"), c("worker-01", "60.9"), c("worker-01", "53.4"), c("worker-01", "38.9"), c("worker-01", "27.5"), c("worker-02", "14.2"), c("worker-02", "53.9"), c("worker-02", "21.8"), c("worker-02", "47.8"), c("worker-02", "48.8"), c("worker-02", "76.5"), c("worker-02", "41.3"), c("worker-02", "20.4"), c("worker-02", "31.9"), c("worker-02", "31.1"), c("worker-03", "186"), c("worker-03", "84.6"), c("worker-03", "94.4"), c("worker-03", "218"), c("worker-03", "189"), c("worker-03", "130"), c("worker-03", "107"), c("worker-03", "80.6"), c("worker-03", "288"), c("worker-03", "173"), c("worker-04", "23.5"), c("worker-04", "16.2"), c("worker-04", "40.2"), c("worker-04", "130"), c("worker-04", "42.2"), c("worker-04", "25.7"), c("worker-04", "35.4"), c("worker-04", "40.8"), c("worker-04", "109"), c("worker-04", "40.9"), c("worker-05", "43.8"), c("worker-05", "31.1"), c("worker-05", "13.1"), c("worker-05", "24.1"), c("worker-05", "27.7"), c("worker-05", "23.9"), c("worker-05", "40.2"), c("worker-05", "60.3"), c("worker-05", "29.8"), c("worker-05", "37.2"), c("worker-06", "41"), c("worker-06", "11.4"), c("worker-06", "4.44"), c("worker-06", "12.9"), c("worker-06", "22.7"), c("worker-06", "20.5"), c("worker-06", "12.6"), c("worker-06", "8.35"), c("worker-06", "13.6"), c("worker-06", "28.1"), c("worker-07", "6.56"), c("worker-07", "9.5"), c("worker-07", "6.97"), c("worker-07", "5.92"), c("worker-07", "2.42"), c("worker-07", "14"), c("worker-07", "12.3"), c("worker-07", "3.07"), c("worker-07", "7.01"), c("worker-07", "6.49"), c("worker-08", "9.21"), c("worker-08", "9.42"), c("worker-08", "28.7"), c("worker-08", "72.9"), c("worker-08", "35.6"), c("worker-08", "17.2"), c("worker-08", "20.2"), c("worker-08", "13.4"), c("worker-08", "10.5"), c("worker-08", "26.3"), c("worker-09", "19.6"), c("worker-09", "14.3"), c("worker-09", "22.8"), c("worker-09", "35.1"), c("worker-09", "28.9"), c("worker-09", "36.9"), c("worker-09", "13"), c("worker-09", "13.3"), c("worker-09", "13.6"), c("worker-09", "37"), c("worker-10", "78.7"), c("worker-10", "28.2"), c("worker-10", "41.3"), c("worker-10", "14.4"), c("worker-10", "72.9"), c("worker-10", "10.2"), c("worker-10", "16.2"), c("worker-10", "15.8"), c("worker-10", "42.2"), c("worker-10", "61")), stringsAsFactors = F)
          names(sample.3) <- c("worker", "x")
          names(sample.4) <- c("worker", "x")
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
          
          # Bayesian calculations      
      
          #####default webexpo analysis - JAGS
          
          mcmc.3 <- Webexpo.between.globalbayesian.jags(  data.sample = sample.3,
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
          
          mcmc.4 <- Webexpo.between.globalbayesian.jags(  data.sample = sample.4,
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
          
          num.4 <- Webexpo.between.interpretation(  is.lognormal = TRUE , 
                                                    mu.overall = mcmc.4$mu.overall.chain,
                                                    sigma.between.chain = mcmc.4$sigma.between.chain ,
                                                    sigma.within.chain = mcmc.4$sigma.within.chain ,
                                                    probacred = 90, 
                                                    oel = 150,
                                                    target_perc = 95,
                                                    rappap_cover = 80, 
                                                    wwct = 0.2 , 
                                                    prob.ind.overexpo.threshold = 20)
          # Summaries
          
          sum.seg.T6 <- webexpo.between.summary( labels = c("Low rho" , "High rho") ,
                                                result.list = list(num.3,num.4) ,
                                                is.lognormal = TRUE )
          #Saving results
          
          write.csv2(sum.seg.T6,here("EXAMPLES", "output", "table6.csv"))
          
          
######## exposure metrics for the least and most exposed workers in sample 1 and sample 2 : TABLE7
      
      #### least and most exposed in sample 1 : #5 and #8   
      #### least and most exposed in sample 2 : #7 and #3
      
      
      # Numerical interpretation (worker specific)
              
      least.S1 <-Webexpo.between.interpretation.worker(worker.id  = "worker-04", 
                                                       is.lognormal = TRUE , 
                                                       mu.workers = mcmc.3$mu.workers,
                                                       worker.index = mcmc.3$worker.index, 
                                                       sigma.within.chain = mcmc.3$sigma.within.chain ,
                                                       probacred = 90, 
                                                       oel = 150,
                                                       frac_threshold = 5,
                                                       target_perc = 95 )
      most.S1 <-Webexpo.between.interpretation.worker(worker.id  = "worker-08", 
                                                      is.lognormal = TRUE , 
                                                      mu.workers = mcmc.3$mu.workers,
                                                      worker.index = mcmc.3$worker.index, 
                                                      sigma.within.chain = mcmc.3$sigma.within.chain ,
                                                      probacred = 90, 
                                                      oel = 150,
                                                      frac_threshold = 5,
                                                      target_perc = 95 ) 
      
      least.S2 <-Webexpo.between.interpretation.worker(worker.id  = "worker-07", 
                                                       is.lognormal = TRUE , 
                                                       mu.workers = mcmc.4$mu.workers,
                                                       worker.index = mcmc.4$worker.index, 
                                                       sigma.within.chain = mcmc.4$sigma.within.chain ,
                                                       probacred = 90, 
                                                       oel = 150,
                                                       frac_threshold = 5,
                                                       target_perc = 95 )
      
      most.S2 <-Webexpo.between.interpretation.worker(worker.id  = "worker-03", 
                                                      is.lognormal = TRUE , 
                                                      mu.workers = mcmc.4$mu.workers,
                                                      worker.index = mcmc.4$worker.index, 
                                                      sigma.within.chain = mcmc.4$sigma.within.chain ,
                                                      probacred = 90, 
                                                      oel = 150,
                                                      frac_threshold = 5,
                                                      target_perc = 95 )
      
      # Summaries
      
      sum.seg.T7 <- webexpo.seg.summary( labels = c("least.s1","most.s1","least.s2","most.s2") ,  
                                        result.list = list(least.S1,most.S1,least.S2,most.S2)
                                        , is.lognormal = TRUE,
                                        show.bands = F)
      
      
      #Saving results
      
      write.csv2(sum.seg.T7,here("EXAMPLES", "output", "table7.csv"))
      
      
      
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
         

