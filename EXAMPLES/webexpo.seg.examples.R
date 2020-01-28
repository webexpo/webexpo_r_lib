##############################################################
#
#   WEBEXPO official R scripts
#
#   SEG ANALYSIS
#  
#   EXAMPLES FOUND IN THE FINAL REPORT
#
#   V1.0 Sept 2018  
#
#############################################################


###### library

library(rjags)

library(ggplot2)

require(gridExtra)

library(here)

#####  sourcing files

    ## data preparation / generation
    
    source(here("RANDOM SAMPLE GENERATION", "webexpo.seg.randomgeneration.R"))
    
    source(here("DATA PREPARATION", "SEG ANALYSIS", "webexpo.seg.dataprep.R"))
    
    # JAGS models
    
    source(here("JAGS MODELS", "SEG ANALYSIS", "webexpo.seg.mainbayesian.R"))
    
    source(here("JAGS MODELS", "SEG ANALYSIS", "webexpo.seg.informedvarbayesian.R"))
    
    source(here("JAGS MODELS", "SEG ANALYSIS", "webexpo.seg.informedvarbayesian.models.R"))
    
    source(here("JAGS MODELS", "SEG ANALYSIS", "webexpo.seg.uninformativebayesian.R"))
    
    source(here("JAGS MODELS", "SEG ANALYSIS", "webexpo.seg.uninformativebayesian.models.R"))
    
    source(here("JAGS MODELS", "SEG ANALYSIS", "webexpo.seg.riskbandbayesian.R"))
    
    source(here("JAGS MODELS", "SEG ANALYSIS", "webexpo.seg.riskbandbayesian.models_old.R"))
    
    # McGill models
    
    source(here("R MODELS", "WRAPPING FUNCTIONS", "webexpo.seg.mainbayesian.mcgill.R"))
    
    source(here("R MODELS", "McGILL FUNCTIONS", "data-summary.R"))
    
    source(here("R MODELS", "McGILL FUNCTIONS", "fcts.R"))
    
    source(here("R MODELS", "McGILL FUNCTIONS", "model-SEG-informedvar.R"))
    
    source(here("R MODELS", "McGILL FUNCTIONS", "model-SEG-uninformative.R"))
    
    source(here("R MODELS", "McGILL FUNCTIONS", "model-SEG-riskband.R"))
    
    # Data interpretation
    
    source(here("RESULT INTERPRETATION", "SEG ANALYSIS", "webexpo.seg.interpretation.R"))
    
    source(here("RESULT INTERPRETATION", "SEG ANALYSIS", "webexpo.seg.summary.R"))



##################  REPORT EXAMPLES ################

### SAMPLE 1 : lognormal sample, no measurement error, no censoring, section 4.2.1.6

          #Instructions disabled below were used to generate the sample 
          #sample.seg.1 <- webexpo.seg.gener.LN(n = 9,
          #                                 no.censoring = TRUE,
          #                                 gm = 30,
          #                                 gsd = 2,
          #                                 error = "none")
          
          sample.seg.1 <- c("24.7","64.1","13.8","43.7","19.9","133","32.1","15","53.7")
          
    
          #Bayesian calculation using JAGS models
    
          mcmc.seg.1 <- Webexpo.seg.globalbayesian.jags( data.sample = sample.seg.1 ,
                                                     is.lognormal = TRUE , 
                                                     error.type = "none" ,
                                                     oel = 100 )
          
          #Numerical interpretation of the MCMC chains
          
          num.seg.1 <- all.numeric(mu.chain = mcmc.seg.1$mu.chain,
                                sigma.chain = mcmc.seg.1$sigma.chain,
                                probacred = 90 ,
                                oel = 100,
                                frac_threshold =5 ,
                                target_perc = 95)
          
          #Summary of results into one table (TABLE 4)
          
          sum.seg.1 <- webexpo.seg.summary( labels=c("seg.1") ,
                                       result.list = list(num.seg.1) , 
                                       is.lognormal = TRUE) 
          
          
          # Graphical illustration of posteriors
          
          # mu and sigma
          
          data.to.plot <- data.frame(mu=as.numeric(mcmc.seg.1$mu.chain),sigma=as.numeric(mcmc.seg.1$sigma.chain))
          
          p1 <-ggplot(data=data.to.plot, aes(x=mu)) + geom_histogram(aes(y =..density.. )) 
          p1 <-p1 + labs(x="Posterior sample for mu")
          p1 <-p1 + geom_density(col=1)
          p1 <-p1 +theme(axis.title.x=element_text(vjust=-0.5,size=20))+
            theme(axis.title.y=element_text(size=20,angle=90))+
            theme(axis.text.x=element_text(size=20))+
            theme(axis.text.y=element_text(size=10,angle=90))
          p1 <- p1 + theme_bw()
          
          
          p2 <-ggplot(data=data.to.plot, aes(x=sigma)) + geom_histogram(aes(y =..density.. )) 
          p2 <-p2 + labs(x="Posterior sample for sigma")
          p2 <-p2 + geom_density(col=1)
          p2 <-p2 +theme(axis.title.x=element_text(vjust=-0.5,size=20))+
            theme(axis.title.y=element_text(size=20,angle=90))+
            theme(axis.text.x=element_text(size=20))+
            theme(axis.text.y=element_text(size=10,angle=90))
          p2 <- p2 + theme_bw()
          
          
          grid.arrange(p1, p2, ncol=2)
          
          
          #### posteriors for AM and the 95th percentile
          
          data.to.plot$P95 <- exp(data.to.plot$mu + qnorm(0.95)*data.to.plot$sigma)
          
          data.to.plot$am <- exp(data.to.plot$mu + 0.5*data.to.plot$sigma^2)
          
          p3 <-ggplot(data=data.to.plot, aes(x=P95)) + geom_histogram(aes(y =..density.. )) + xlim(c(0,500)) 
          p3 <-p3 + labs(x="Posterior sample for 95th percentile")
          p3 <-p3 + geom_density(col=1)
          p3 <-p3 +theme(axis.title.x=element_text(vjust=-0.5,size=20))+
            theme(axis.title.y=element_text(size=20,angle=90))+
            theme(axis.text.x=element_text(size=20))+
            theme(axis.text.y=element_text(size=10,angle=90))
          p3 <- p3 + theme_bw()
          p3 <- p3 + geom_segment(aes(x = 100 , y = 0, 
                                        xend = 100 , 
                                        yend = 0.012 ),size=1.3,color="black",linetype = "dashed")+
            annotate("text", x=140, y=0.012, label= "OEL")
          
          
          
          p4 <-ggplot(data=data.to.plot, aes(x=am)) + geom_histogram(aes(y =..density.. )) + xlim(c(0,150))
          p4 <-p4 + labs(x="Posterior sample for arithmetic mean")
          p4 <-p4 + geom_density(col=1)
          p4 <-p4 +theme(axis.title.x=element_text(vjust=-0.5,size=20))+
            theme(axis.title.y=element_text(size=20,angle=90))+
            theme(axis.text.x=element_text(size=20))+
            theme(axis.text.y=element_text(size=10,angle=90))
          p4 <- p4 + theme_bw()
          p4 <- p4 +   geom_segment(aes(x = 100 , y = 0, 
                                        xend = 100 , 
                                        yend = 0.033),size=1.3,color="black",linetype = "dashed")+
            annotate("text", x=110, y=0.033, label= "OEL")
          
          grid.arrange(p3, p4, ncol=2)
          

### TABLE 5 : Sensitivity to priors       

sample.seg.1 <- c("24.7","64.1","13.8","43.7","19.9","133","32.1","15","53.7")

        
        #informedvar JAGS          
        mcmc.seg.1 <- Webexpo.seg.globalbayesian.jags(  data.sample = sample.seg.1,
                                                        is.lognormal = TRUE , 
                                                        error.type = "none" ,
                                                        me.range = c(0.3,0.3) , 
                                                        oel = 100 ,
                                                        prior.model = "informedvar" ,
                                                        mu.lower = -20 ,
                                                        mu.upper = 20  ,
                                                        log.sigma.mu = -0.1744 ,
                                                        log.sigma.prec = 2.5523 , 
                                                        n.iter = 25000 ,
                                                        n.burnin = 5000 ,
                                                        init.mu = log(0.3) , 
                                                        init.sigma = log(2.5))
        
        
        #uninformative JAGS
        mcmc.seg.2 <- Webexpo.seg.globalbayesian.jags(data.sample = sample.seg.1 ,
                                                  is.lognormal = TRUE , 
                                                  error.type = "none" ,
                                                  me.range = c(0.3,0.3) ,
                                                  oel = 100 ,
                                                  prior.model = "uninformative" ,
                                                  mu.lower.uninf = -20 ,
                                                  mu.upper.uninf = 20 ,
                                                  sd.range = c(0, 2.3) ,
                                                  n.iter = 25000 ,
                                                  n.burnin = 5000 , 
                                                  init.mu = log(0.3) , 
                                                  init.sigma = log(2.5))
        
        
        #past data JAGS
        mcmc.seg.3 <- Webexpo.seg.globalbayesian.jags( data.sample = sample.seg.1,
                                                   is.lognormal = TRUE , 
                                                   error.type = "none" ,
                                                   me.range = c(0.3,0.3) , 
                                                   oel = 100 ,
                                                   prior.model = "informedvar" ,
                                                   mu.lower = -20 ,
                                                   mu.upper = 20  ,
                                                   log.sigma.mu = -0.1744 ,
                                                   log.sigma.prec = 2.5523 , 
                                                   n.iter = 25000 ,
                                                   n.burnin = 5000 ,
                                                   init.mu = log(0.3) , 
                                                   init.sigma = log(2.5),
                                                   past.data = list(mean = log(5), sd = log(2.4), n = 5 ))
        
        #riskband informed 
        
        # graphical assessment of the bands
        webexpo.riksbandcheck(is.lognormal = TRUE,
                              A = c(0.01,0.1,0.5,1) ,
                              gm.min = 0.001 ,
                              gm.max = 100 ,
                              gsd.min = 1.05 ,
                              gsd.max = 10 ,
                              target_perc = 95)
        #bayesian model
        mcmc.seg.4 <- Webexpo.seg.globalbayesian.jags( data.sample = sample.seg.1 ,
                                                         #distribution
                                                         is.lognormal = TRUE , 
                                                         error.type = "none" ,
                                                         me.range = c(0.3,0.3) , 
                                                         oel = 100 ,
                                                         prior.model = "riskband" ,
                                                         mu.lower.riskb = log(0.001) ,
                                                         mu.upper.riskb = log(100) ,
                                                         sigma.lower = log(1.05) ,
                                                         sigma.upper = log(10) ,
                                                         target_perc = 95 ,
                                                         A = c( 0.01 , 0.1 , 0.5 , 1 ) ,
                                                         region.prior.prob = c( 0.1 , 0.2 , 0.5 , 0.1 , 0.1 ) , 
                                                         n.iter = 25000 ,
                                                         n.burnin = 5000 )           
                  
                  
        # Numerical interpretation
        
        num.seg.1 <- all.numeric(mu.chain = mcmc.seg.1$mu.chain,
                              sigma.chain = mcmc.seg.1$sigma.chain,
                              probacred = 90 ,
                              oel = 100,
                              frac_threshold =5 ,
                              target_perc = 95)
        num.seg.2 <- all.numeric(mu.chain = mcmc.seg.2$mu.chain,
                              sigma.chain = mcmc.seg.2$sigma.chain,
                              probacred = 90 ,
                              oel = 100,
                              frac_threshold =5 ,
                              target_perc = 95)
        num.seg.3 <- all.numeric(mu.chain = mcmc.seg.3$mu.chain,
                              sigma.chain = mcmc.seg.3$sigma.chain,
                              probacred = 90 ,
                              oel = 100,
                              frac_threshold =5 ,
                              target_perc = 95)
        num.seg.4 <- all.numeric(mu.chain = mcmc.seg.4$mu.chain,
                              sigma.chain = mcmc.seg.4$sigma.chain,
                              probacred = 90 ,
                              oel = 100,
                              frac_threshold =5 ,
                              target_perc = 95)
        # summaries
        
        sum.seg.T5 <- webexpo.seg.summary( labels=c("Informedvar" , "Uninformative" , "PastData" , "Riskband") ,
                                     result.list = list(num.seg.1,num.seg.2,num.seg.3,num.seg.4) , 
                                     is.lognormal = TRUE)           
        # Saving the results
                  
        write.csv(sum.seg.T5,"C:/jerome/Dropbox/bureau/UdM/projets/sampling strats tool/rapport IRSST/data/table5.csv")          
                  
############# illustration of measurement error treatment - TABLE 6
        
        
        #code disabled below. Find the actual samples in the associated EXCEL FILES or Report appendix
        #sample.seg.2 <- webexpo.seg.gener.LN(n = 100,
        #                                 no.censoring = TRUE,
        #                                 gm = 60,
        #                                 gsd = 1.5,
        #                                error = "cv",
        #                                 me.cv = 0.30)
        
        
        # BayesiaN models
        mcmc.seg.1 <- Webexpo.seg.globalbayesian.jags(  data.sample = sample.seg.2,
                                                        is.lognormal = TRUE , 
                                                        error.type = "none" ,
                                                        me.range = c(0.3,0.3) , 
                                                        oel = 100 ,
                                                        prior.model = "informedvar" ,
                                                        mu.lower = -20 ,
                                                        mu.upper = 20  ,
                                                        log.sigma.mu = -0.1744 ,
                                                        log.sigma.prec = 2.5523 , 
                                                        n.iter = 25000 ,
                                                        n.burnin = 5000 ,
                                                        init.mu = log(0.3) , 
                                                        init.sigma = log(2.5) )
        
        
        mcmc.seg.2 <- Webexpo.seg.globalbayesian.jags(  data.sample = sample.seg.2,
                                                        is.lognormal = TRUE , 
                                                        error.type = "CV" ,
                                                        me.range = c(0.3,0.3) , 
                                                        oel = 100 ,
                                                        prior.model = "informedvar" ,
                                                        mu.lower = -20 ,
                                                        mu.upper = 20  ,
                                                        log.sigma.mu = -0.1744 ,
                                                        log.sigma.prec = 2.5523 , 
                                                        n.iter = 25000 ,
                                                        n.burnin = 5000 ,
                                                        init.mu = log(0.3) , 
                                                        init.sigma = log(2.5) )
        
        mcmc.seg.3 <- Webexpo.seg.globalbayesian.jags( data.sample = sample.seg.2,
                                                       is.lognormal = TRUE , 
                                                       error.type = "CV" ,
                                                       me.range = c(0.15,0.45) , 
                                                       oel = 100 ,
                                                       prior.model = "informedvar" ,
                                                       mu.lower = -20 ,
                                                       mu.upper = 20  ,
                                                       log.sigma.mu = -0.1744 ,
                                                       log.sigma.prec = 2.5523 , 
                                                       n.iter = 25000 ,
                                                       n.burnin = 5000 ,
                                                       init.mu = log(0.3) , 
                                                       init.sigma = log(2.5)  )
        # Numerical interpretation
        
        num.seg.1 <- all.numeric(mu.chain = mcmc.seg.1$mu.chain,
                              sigma.chain = mcmc.seg.1$sigma.chain,
                              probacred = 90 ,
                              oel = 100,
                              frac_threshold =5 ,
                              target_perc = 95)
        num.seg.2 <- all.numeric(mu.chain = mcmc.seg.2$mu.chain,
                              sigma.chain = mcmc.seg.2$sigma.chain,
                              probacred = 90 ,
                              oel = 100,
                              frac_threshold =5 ,
                              target_perc = 95)
        num.seg.3 <- all.numeric(mu.chain = mcmc.seg.3$mu.chain,
                              sigma.chain = mcmc.seg.3$sigma.chain,
                              probacred = 90 ,
                              oel = 100,
                              frac_threshold =5 ,
                              target_perc = 95)
        # Summaries
        
        sum.seg.T6 <- webexpo.seg.summary( labels=c("Naive","Known CV","Unknown CV") ,
                                     result.list = list( num.seg.1 , num.seg.2 , num.seg.3 ) , 
                                     is.lognormal = TRUE) 
        # Saving the results
        
        write.csv(sum.seg.T6,"C:/jerome/Dropbox/bureau/UdM/projets/sampling strats tool/rapport IRSST/data/table6.csv")              
 
##############  normal example -  TABLE 7
        
        
        # Sample
        sample.seg.3 <- c("79.6" , "78.6" , "78.7",  "83.6",  "74.6",  "79.1" , "76.9" , "75.8" , "80.8")
        
        # Bayesian models
        mcmc.seg.1 <- Webexpo.seg.globalbayesian.jags( data.sample = sample.seg.3 ,
                                                       is.lognormal = FALSE , 
                                                       error.type = "none" ,
                                                       me.range = c(0.3,0.3) ,
                                                       oel = 85 ,
                                                       prior.model = "uninformative" ,
                                                       mu.lower.uninf = 40 ,
                                                       mu.upper.uninf = 120 ,
                                                       sd.range = c(0, 100) ,
                                                       n.iter = 25000 ,
                                                       n.burnin = 5000 , 
                                                       init.mu = 70 , 
                                                       init.sigma = 10 )
        # Numerical interpretation
        num.seg.1 <- all.numeric.N(mu.chain = mcmc.seg.1$mu.chain,
                              sigma.chain = mcmc.seg.1$sigma.chain,
                              probacred = 90 ,
                              oel = 85,
                              frac_threshold = 5 ,
                              target_perc = 95 )
        
        # Summaries
        sum.seg.T7 <- webexpo.seg.summary( labels=c("A1") ,
                                     result.list = list(num.seg.1) , 
                                     is.lognormal = FALSE) 
        
        # Saving the results
        write.csv(sum.seg.T7,"C:/jerome/Dropbox/bureau/UdM/projets/sampling strats tool/rapport IRSST/data/table7.csv")
                 
          

           