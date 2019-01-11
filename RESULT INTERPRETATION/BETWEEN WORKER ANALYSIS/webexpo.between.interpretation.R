##############################################################
#
#   WEBEXPO official R scripts
#
#   BETWEEN WORKER ANALYSIS
#  
#   INTERPRETATION OF MCMC CHAINS 
#
#   V1.0 Sept 2018   
#
#############################################################


##
#
## INPUT FOR CREATING THE NUMERICAL OUTPUTS
#
# model : lognormal or normal
#
# MCMC chain for MU.OVERALL : mu.overall.chain
# MCMC chain for SIGMA.BETWEEN : sigma.between.chain
# MCMC chain for SIGMA.WITHIN : sigma.within.chain
# MCMC chains for each worker's mean

# OEL             :  occupational exposure limit      
# probacred       :  the level of confidence for credible intervals (%)
# frac_threshold  :  the acceptable exceedance value (%)
# target_perc     :  the percentile of interest (%)
# rappap_cover                        : the coverage (%) of the R ratio
# wwct                                : the threshold for the within worker correlation coefficient [0-1]
# prob.ind.overexpo.threshold         : the threshold for the probability of individual overexposure (%)
#
# worker.id                           : worker id for individual results
# worker.index                        : table which is an output from the bayesian analysis,
#                                       indicating the rank of each worker id in the results                                            


############# WRAPPING FUNCTION

#############  INDIVIDUAL WORKER RESULTS
            
 Webexpo.between.interpretation.worker <- function(worker.id, 
                                                        is.lognormal , 
                                                        mu.workers,
                                                        worker.index,
                                                        sigma.within.chain ,
                                                        probacred , 
                                                        oel,
                                                        frac_threshold,
                                                        target_perc
                                                       ) {
              

              mu.chain <- mu.workers[ worker.index$numberinoutput[worker.index$worker==worker.id] ,]
              
              sigma.chain <- sigma.within.chain
              
              if (is.lognormal) return( all.numeric( mu.chain = mu.chain,
                                                     sigma.chain = sigma.chain ,
                                                     probacred = probacred ,
                                                     oel = oel , 
                                                     frac_threshold = frac_threshold ,
                                                     target_perc = target_perc))
              
              
              if (!is.lognormal) return( all.numeric.N( mu.chain = mu.chain,
                                                     sigma.chain = sigma.chain ,
                                                     probacred = probacred ,
                                                     oel = oel , 
                                                     frac_threshold = frac_threshold ,
                                                     target_perc = target_perc))  
              
            }

            
            
#############  GROUP RESULTS
            
Webexpo.between.interpretation <- function(  is.lognormal , 
                                             mu.overall,
                                             sigma.between.chain ,
                                            sigma.within.chain ,
                                            probacred , 
                                            oel,
                                            target_perc,
                                            rappap_cover, 
                                            wwct , 
                                            prob.ind.overexpo.threshold  
            ) {
              

              
              if (is.lognormal) {return( all.numeric.between( mu.overall = mu.overall,
                                                             sigma.between = sigma.between.chain ,
                                                             sigma.within = sigma.within.chain ,
                                                     probacred = probacred ,
                                                     oel = oel , 
                                                     target_perc = target_perc,
                                                     rappap_cover = rappap_cover,
                                                     wwct = wwct,
                                                     prob.ind.overexpo.threshold = prob.ind.overexpo.threshold))}
              
              
              if (!is.lognormal) {return( all.numeric.between.N( mu.overall = mu.overall,
                                                                sigma.between = sigma.between.chain ,
                                                                sigma.within = sigma.within.chain ,
                                                                probacred = probacred ,
                                                                oel = oel , 
                                                                target_perc = target_perc,
                                                                rappap_cover = rappap_cover,
                                                                wwct = wwct,
                                                                prob.ind.overexpo.threshold = prob.ind.overexpo.threshold)) } 
              
            }
            
            
            
            
            
###############################    GROUP FUNCTIONS

############## A. LOGNORMAL MODEL

      #within worker correlation
      rho <- function( sigma.within , sigma.between , probacred) {  
        
        chain <- sigma.between^2 /( sigma.between^2 + sigma.within^2 )
        
        est <- median( chain )
        
        lcl <- quantile( chain , ( 100 - probacred ) / 200 )
        
        ucl <- quantile( chain , 1-( 100 - probacred ) / 200)
        
        return( list( est = est , lcl = lcl , ucl = ucl ) )
        
      }

      #proba that correlation >threshold
      prob.rho.overX <- function( sigma.within , sigma.between , wwct ) {  ## arithmetic mean
        
        chain <- sigma.between^2 /( sigma.between^2 + sigma.within^2 )
        
        risk <- 100 * length( chain[ chain > wwct ] ) / length( chain )
        
        return( risk )
        
      }

      ##rappaport ratio
      R.ratio <- function( sigma.between , probacred , rappap_cover ) {  ## arithmetic mean
        
        chain <- exp( 2 * qnorm( 1-( 100 -  rappap_cover ) / 200 ) * sigma.between )
        
        est <- median( chain )
        
        lcl <- quantile( chain , ( 100 - probacred ) / 200 )
        
        ucl <- quantile( chain , 1-( 100 - probacred ) / 200)
        
        return( list( est = est , lcl = lcl , ucl = ucl ) )
      }


      ##rappaport ratio >threshold
      R.ratio.greaterthanX <- function( sigma.between , rappap_cover , hetero.criterion ) {  ## arithmetic mean
        
        chain <- exp( 2 * qnorm( 1-( 100 -  rappap_cover ) / 200  ) * sigma.between )
        
        risk <- 100 * length( chain[ chain > hetero.criterion ] ) / length( chain )
        
        return(risk)
        
      }


      #probability of individual overexposure based on the critical percentile 
      prob.ind.overexpo.perc <- function( mu.overall, sigma.within , sigma.between , target_perc , oel , probacred ) {  
        
        qn <- qnorm( target_perc / 100 )
        
        chain <- 100*( 1- pnorm( ( log( ( oel ) ) - mu.overall - qn * sigma.within ) / sigma.between ) )
        
        est <- median( chain )
        
        lcl <- quantile( chain , ( 100 - probacred ) / 200 )
        
        ucl <- quantile( chain , 1-( 100 - probacred ) / 200)
        
        return( list( est = est , lcl = lcl , ucl = ucl ) )
        
      }

        #chances that probability of individual overexposure based on the critical percentile is overX
        prob.ind.overexpo.perc.proboverX <- function( mu.overall, sigma.within , sigma.between ,  target_perc , oel , X ) {  
          
          qn <- qnorm( target_perc / 100 )
          
          chain <- 100*( 1- pnorm( ( log( ( oel ) ) - mu.overall - qn * sigma.within ) / sigma.between ) )
          
          risk <- 100 * length( chain[ chain > X ] ) / length( chain )
          
          return( risk )
          
        }


        #probability of individual overexposure based on the arithmetic mean
        prob.ind.overexpo.am <- function( mu.overall, sigma.within , sigma.between , oel , probacred) {  
          
          chain <- 100 * ( 1 - pnorm( ( log( ( oel ) ) - mu.overall - 0.5 * sigma.within^2 ) / sigma.between ) )
          
          est <- median( chain )
          
          lcl <- quantile( chain , ( 100 - probacred ) / 200 )
          
          ucl <- quantile( chain , 1-( 100 - probacred ) / 200)
          
          return( list( est = est , lcl = lcl , ucl = ucl ) )
          
        }

        #chances that probability of individual overexposure based on the arithmetic mean is overX
        prob.ind.overexpo.am.proboverX <- function( mu.overall, sigma.within , sigma.between , oel , X) {  
          
          chain <- 100 * ( 1 - pnorm( ( log( ( oel ) ) - mu.overall - 0.5 * sigma.within^2 ) / sigma.between ) )
          
          risk <- 100 * length( chain[ chain > X ] ) / length( chain )
          
          return( risk )
          
        }

          ####### creating a list containing everything 
          
          all.numeric.between <-function(mu.overall, sigma.within , sigma.between,
                                         oel ,
                                         probacred ,
                                         target_perc ,
                                         rappap_cover , 
                                         wwct , 
                                         prob.ind.overexpo.threshold  ){
            
            return( list(
              
              group.gm = gm( mu.chain = mu.overall, probacred = probacred ),
              
              gsdb = gsd( sigma.chain = sigma.between , 
                          probacred = probacred ),
              
              gsdw = gsd( sigma.chain = sigma.within , 
                          probacred = probacred ),
              
              rho = rho( sigma.within = sigma.within , 
                         sigma.between = sigma.between ,
                         probacred = probacred ),
              
              prob.rho.overX = prob.rho.overX( sigma.within = sigma.within ,
                                               sigma.between = sigma.between , 
                                               wwct = wwct ),
              
              R.ratio = R.ratio( sigma.between = sigma.between  , 
                                 probacred = probacred , 
                                 rappap_cover = rappap_cover ),
              
              prob.R.over2 = R.ratio.greaterthanX( sigma.between = sigma.between  ,
                                                   rappap_cover =rappap_cover ,
                                                   hetero.criterion = 2 ),
              
              prob.R.over10 = R.ratio.greaterthanX( sigma.between = sigma.between  , 
                                                    rappap_cover = rappap_cover ,
                                                    hetero.criterion = 10),
              
              
              prob.ind.overexpo.perc = prob.ind.overexpo.perc( mu.overall = mu.overall, 
                                                               sigma.within = sigma.within , 
                                                               sigma.between = sigma.between , 
                                                               target_perc = target_perc , 
                                                               oel = oel ,
                                                               probacred = probacred),
              
              prob.ind.overexpo.perc.proboverX = prob.ind.overexpo.perc.proboverX( mu.overall = mu.overall,
                                                                                   sigma.within = sigma.within ,
                                                                                   sigma.between = sigma.between , 
                                                                                   target_perc = target_perc  ,
                                                                                   oel = oel ,
                                                                                   X = prob.ind.overexpo.threshold),
              
              prob.ind.overexpo.am = prob.ind.overexpo.am( mu.overall = mu.overall, 
                                                           sigma.within = sigma.within , 
                                                           sigma.between = sigma.between , 
                                                           oel=oel , 
                                                           probacred = probacred),
              
              prob.ind.overexpo.am.proboverX = prob.ind.overexpo.am.proboverX(mu.overall = mu.overall,
                                                                              sigma.within = sigma.within , 
                                                                              sigma.between = sigma.between , 
                                                                              oel = oel ,
                                                                              X = prob.ind.overexpo.threshold)
            ))
            
          }
            
            
            
###################### B. NORMAL MODEL
            
          ##R difference : difference between the "worst" and "best" worker relative to group mean
          R.diff <- function( mu.overall , sigma.between , probacred , rappap_cover ) {  ## arithmetic mean
            
            group.mean <- median( mu.overall )
            
            chain <- 100 * 2 * qnorm( 1-( 100 -  rappap_cover ) / 200  ) * sigma.between / group.mean
            
            est <- median( chain )
            
            lcl <- quantile( chain , ( 100 - probacred ) / 200 )
            
            ucl <- quantile( chain , 1-( 100 - probacred ) / 200)
            
            return( list( est = est , lcl = lcl , ucl = ucl ) )
          }
          
          
          #probability of individual overexposure based on the critical percentile 
          prob.ind.overexpo.perc.N <- function( mu.overall, sigma.within , sigma.between , target_perc , oel , probacred ) {  
            
            qn <- qnorm( target_perc / 100 )
            
            chain <- 100*( 1- pnorm( ( oel  - mu.overall - qn * sigma.within ) / sigma.between ) )
            
            est <- median( chain )
            
            lcl <- quantile( chain , ( 100 - probacred ) / 200 )
            
            ucl <- quantile( chain , 1-( 100 - probacred ) / 200)
            
            return( list( est = est , lcl = lcl , ucl = ucl ) )
            
          }
          
          #chances that probability of individual overexposure based on the critical percentile is overX
          prob.ind.overexpo.perc.proboverX.N <- function( mu.overall, sigma.within , sigma.between ,  target_perc , oel , X ) {  
            
            qn <- qnorm( target_perc / 100 )
            
            chain <- 100*( 1- pnorm( ( oel  - mu.overall - qn * sigma.within ) / sigma.between ) )
            
            risk <- 100 * length( chain[ chain > X ] ) / length( chain )
            
            return( risk )
            
          }
          
          
          #probability of individual overexposure based on the arithmetic mean
          prob.ind.overexpo.am.N <- function( mu.overall, sigma.between , oel , probacred) {  
            
            chain <- 100 * ( 1 - pnorm( ( oel - mu.overall ) / sigma.between ) )
            
            est <- median( chain )
            
            lcl <- quantile( chain , ( 100 - probacred ) / 200 )
            
            ucl <- quantile( chain , 1-( 100 - probacred ) / 200)
            
            return( list( est = est , lcl = lcl , ucl = ucl ) )
            
          }
          
          #chances that probability of individual overexposure based on the arithmetic mean is overX
          prob.ind.overexpo.am.proboverX.N <- function( mu.overall, sigma.between , oel , X) {  
            
            chain <- 100 * ( 1 - pnorm( ( oel - mu.overall ) / sigma.between ) )
            
            risk <- 100 * length( chain[ chain > X ] ) / length( chain )
            
            return( risk )
            
          }
          
          ####### creating a list containing everything 
          all.numeric.between.N <-function(mu.overall, sigma.within , sigma.between,
                                           oel ,
                                           probacred ,
                                           target_perc ,
                                           rappap_cover , 
                                           wwct , 
                                           prob.ind.overexpo.threshold  ) {
            
            return( list(
              
              group.mean = am.N( mu.chain = mu.overall, probacred = probacred ),
              
              asdb = asd.N( sigma.chain = sigma.between , 
                            probacred = probacred ),
              
              asdw = asd.N( sigma.chain = sigma.within , 
                            probacred = probacred ),
              
              rho = rho( sigma.within = sigma.within , 
                         sigma.between = sigma.between ,
                         probacred = probacred ),
              
              prob.rho.overX = prob.rho.overX( sigma.within = sigma.within ,
                                               sigma.between = sigma.between ,
                                               wwct = wwct ),
              
              R.diff = R.diff(  mu.overall = mu.overall , 
                                sigma.between = sigma.between  , 
                                probacred = probacred , 
                                rappap_cover = rappap_cover ),
              
              
              prob.ind.overexpo.perc = prob.ind.overexpo.perc.N( mu.overall = mu.overall, 
                                                               sigma.within = sigma.within , 
                                                               sigma.between = sigma.between , 
                                                               target_perc = target_perc , 
                                                               oel = oel ,
                                                               probacred = probacred),
              
              prob.ind.overexpo.perc.proboverX = prob.ind.overexpo.perc.proboverX.N( mu.overall = mu.overall,
                                                                                   sigma.within = sigma.within ,
                                                                                   sigma.between = sigma.between , 
                                                                                   target_perc = target_perc  ,
                                                                                   oel = oel ,
                                                                                   X = prob.ind.overexpo.threshold),
              
              prob.ind.overexpo.am = prob.ind.overexpo.am.N( mu.overall = mu.overall, 
                                                           sigma.between = sigma.between , 
                                                           oel=oel , 
                                                           probacred = probacred),
              
              prob.ind.overexpo.am.proboverX = prob.ind.overexpo.am.proboverX.N(mu.overall = mu.overall,
                                                                              
                                                                              sigma.between = sigma.between , 
                                                                              oel = oel ,
                                                                              X = prob.ind.overexpo.threshold)
            ))
            
            
          }
          
          