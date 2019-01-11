##############################################################
#
#   WEBEXPO official R scripts
#
#   SEG ANALYSIS
#  
#   INTERPRETATION OF MCMC CHAINS 
#
#   V1.0 Sept 2018  
#   
#
#############################################################



##
#
## INPUT FOR CREATING THE NUMERICAL OUTPUTS
#
# model : lognormal or normal
#
# MCMC chain for MU : mu.chain
# MCMC chain for SIGMA : sigma.chain

# OEL             :  occupational exposure limit      
# probacred       :  the level of confidence for credible intervals (%)
# frac_threshold  :  the acceptable exceedance value (%)
# target_perc     :  the percentile of interest (%)



############# WRAPPING FUNCTION

Webexpo.seg.interpretation <- function( is.lognormal , 
                                        mu.chain ,
                                        sigma.chain ,
                                        probacred , 
                                        oel,
                                        frac_threshold,
                                        target_perc) {
  
  
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


############## A. LOGNORMAL MODEL

#### parameter calculations - functions
            
            #geometric mean
            gm <- function( mu.chain , probacred ) {
              
              chain <- exp( mu.chain )
              
              est <- median( chain )
              
              lcl <- quantile( chain , (100 - probacred ) / 200 )
              
              ucl <- quantile( chain , 1 - (100 - probacred )/ 200 )
              
              return( list( est = est , lcl = lcl , ucl = ucl ) )
              
            }  
            
            #geometric standard deviation
            gsd <- function( sigma.chain , probacred ) {
              
              chain <- exp( sigma.chain )
              
              est <- median( chain )
              
              lcl <- quantile( chain , ( 100 - probacred ) / 200 )
              
              ucl <- quantile( chain , 1 - ( 100 -probacred ) / 200 ) 
              
              return( list( est = est , lcl = lcl , ucl = ucl ) )
              
            } 
            
            
            #exceedance fraction
            frac <- function( mu.chain , sigma.chain , probacred , oel ) {  
              
              chain <- 100 * ( 1 - pnorm( ( log( ( oel ) ) - mu.chain ) / sigma.chain ) )
              
              est <- median( chain )
              
              lcl <- quantile( chain , ( 100 - probacred ) / 200 )
              
              ucl <- quantile( chain , 1 - ( 100 - probacred ) / 200 )
              
              return( list( est = est , lcl = lcl , ucl = ucl ) )
              
            }
            
            #percentile of interest
            perc <- function( mu.chain , sigma.chain , target_perc , probacred) {
              
              chain <- exp( mu.chain + qnorm( target_perc / 100 ) * sigma.chain )
              
              est <- median( chain )
              
              lcl <- quantile( chain , ( 100 - probacred ) / 200 )
              
              ucl <- quantile( chain , 1 - ( 100 - probacred ) / 200 )
              
              return( list( est = est , lcl = lcl , ucl = ucl ) )
              
            }
            
            ## arithmetic mean
            am <- function( mu.chain , sigma.chain , probacred ) {  
              
              chain <- exp( mu.chain + 0.5 * sigma.chain^2 )
              
              est <- median( chain )
              
              lcl <- quantile( chain , ( 100 - probacred ) / 200 )
              
              ucl <- quantile( chain , 1 - ( 100 - probacred ) / 200 )
              
              return( list( est = est , lcl = lcl , ucl = ucl ) )
              
            }
            
            
            #### Overexposure risk based on exceedance fraction
            frac.risk <- function( mu.chain , sigma.chain , frac_threshold , oel) {  
              
              chain <- 100 * ( 1 - pnorm( ( log( ( oel ) ) - mu.chain ) / sigma.chain ) )
              
              risk <- 100 * length( chain[ chain > frac_threshold ] ) / length( chain )
              
              return( risk )
              
            }
            
            
            #### Overexposure risk based on critical percentile
            perc.risk <- function( mu.chain , sigma.chain , target_perc , oel ) { 
              
              chain <- exp( mu.chain + qnorm( target_perc / 100 ) * sigma.chain )
              
              risk <- 100 * length( chain[ chain > oel ] ) / length( chain )
              
              return( risk )
              
            }
            
            #### Overexposure risk based on arithmetic mean
            am.risk <- function( mu.chain, sigma.chain , oel ) {  ## arithmetic mean
              
              chain <- exp( mu.chain + 0.5 * sigma.chain^2 )
              
              risk <- 100 * length( chain[ chain > oel ] ) / length( chain )
              
              return( risk )
              
            }
            
            
            #### probability of the critical percentile in each AIHA riskband 
            perc.riskbands <- function( mu.chain , sigma.chain , target_perc , oel ) { 
              
              chain <- exp( mu.chain + qnorm( target_perc / 100 ) * sigma.chain )
              
              riskbands <- 100*table( cut( chain , c(0.000001 , 0.01 , 0.1 , 0.5 , 1, 1000000)*oel)) / length(chain)
              
              names(riskbands) <- c("<0.01*OEL" , "[0.01-0.1]*OEL" , "[0.1-0.5]*OEL" , "[0.5-1]*OEL" , ">OEL")
              
              return( riskbands )
              
            }
            
            
            
            #### probability of the critical percentile in each AIHA riskband 
            am.riskbands <- function( mu.chain , sigma.chain , oel ) { 
              
              chain <- exp( mu.chain + 0.5 * sigma.chain^2 )
              
              riskbands <- 100*table( cut( chain , c(0.000001 , 0.01 , 0.1 , 0.5 , 1, 1000000)*oel)) / length(chain)
              
              names(riskbands) <- c("<0.01*OEL" , "[0.01-0.1]*OEL" , "[0.1-0.5]*OEL" , "[0.5-1]*OEL" , ">OEL")
              
              return( riskbands )
              
            }
            
            
            ####### creating a global list containing all results
            
            all.numeric <-function( mu.chain , sigma.chain , probacred ,
                                   oel , frac_threshold , target_perc) {
              
              return(list(
                
                gm = gm( mu.chain = mu.chain, probacred = probacred),
                
                gsd = gsd( sigma.chain = sigma.chain , probacred = probacred ),
                
                frac = frac( mu.chain = mu.chain , sigma.chain = sigma.chain , probacred = probacred, oel = oel ),
                
                perc = perc( mu.chain = mu.chain , sigma.chain = sigma.chain, target_perc = target_perc , probacred = probacred ),
                
                am = am( mu.chain = mu.chain, sigma.chain =sigma.chain , probacred = probacred),
                
                frac.risk = frac.risk( mu.chain =mu.chain , sigma.chain = sigma.chain, frac_threshold = frac_threshold , oel = oel),
                
                perc.risk = perc.risk( mu.chain = mu.chain , sigma.chain = sigma.chain , target_perc = target_perc , oel = oel),
                
                am.risk = am.risk( mu.chain = mu.chain , sigma.chain = sigma.chain , oel = oel),
                
                perc.riskbands = perc.riskbands( mu.chain = mu.chain , sigma.chain = sigma.chain , target_perc = target_perc , oel = oel),
                
                am.riskbands = am.riskbands( mu.chain = mu.chain , sigma.chain = sigma.chain , oel = oel)
                ))
              
              
            }
            
###################### B. NORMAL MODEL
            

            #arithmetic mean
            am.N <- function( mu.chain , probacred ) {
              
              chain <- mu.chain
              
              est <- median( chain )
              
              lcl <- quantile( chain , (100 - probacred ) / 200 )
              
              ucl <- quantile( chain , 1 - (100 - probacred )/ 200 )
              
              return( list( est = est , lcl = lcl , ucl = ucl ) )
              
            }  
            
            #arithmetic standard deviation
            asd.N <- function( sigma.chain , probacred ) {
              
              chain <- sigma.chain
              
              est <- median( chain )
              
              lcl <- quantile( chain , (100 - probacred ) / 200 )
              
              ucl <- quantile( chain , 1 - (100 - probacred )/ 200 )
              
              return( list( est = est , lcl = lcl , ucl = ucl ) )
              
            } 
            
            #exceedance fraction
            frac.N <- function( mu.chain , sigma.chain , probacred , oel) {  
              
              chain <- 100 * ( 1- pnorm( ( oel - mu.chain ) / sigma.chain ) )
              
              est <- median( chain )
              
              lcl <- quantile( chain , (100 - probacred ) / 200 )
              
              ucl <- quantile( chain , 1 - (100 - probacred )/ 200 )
              
              return( list( est = est , lcl = lcl , ucl = ucl ) )
              
            }
            
            #percentile of interest
            perc.N <- function( mu.chain , sigma.chain, target_perc , probacred) { 
              
              chain <- mu.chain + qnorm( target_perc / 100 ) * sigma.chain
              
              est <- median( chain )
              
              lcl <- quantile( chain , (100 - probacred ) / 200 )
              
              ucl <- quantile( chain , 1 - (100 - probacred )/ 200 )
              
              return( list( est = est , lcl = lcl , ucl = ucl ) )
              
            }
            
            
            #### Overexposure risk based on exceedance fraction
            frac.risk.N <- function( mu.chain , sigma.chain , frac_threshold , oel ) { 
              
              chain <-100 * ( 1 - pnorm( ( oel - mu.chain ) / sigma.chain ) )
              
              risk <- 100 * length( chain[ chain > frac_threshold ] ) / length( chain )
              
              return( risk )
              
            }
            
            #### Overexposure risk based on critical percentile
            perc.risk.N <- function( mu.chain , sigma.chain , target_perc , oel ) { 
              
              chain <-mu.chain + qnorm( target_perc / 100 ) * sigma.chain
              
              risk <- 100 * length( chain[ chain > oel ] ) / length( chain )
              
              return( risk )
              
            }
            
            #### Overexposure risk based on arithmetic mean
            am.risk.N <- function( mu.chain , oel) {  
              
              chain <- mu.chain
              
              risk <- 100 * length( chain[ chain > oel ] ) / length( chain )
              
              return( risk )
              
            }
            
            
            
            ####### creating a list containing everything
            
            all.numeric.N <-function(mu.chain , sigma.chain , probacred ,
                                     oel , frac_threshold , target_perc) {
              
              return(list(
                
                am = am.N( mu.chain = mu.chain , probacred = probacred ),
                
                asd = asd.N( sigma.chain = sigma.chain , probacred = probacred ),
                
                am.risk = am.risk.N( mu.chain = mu.chain , oel = oel ),
                
                frac = frac.N( mu.chain = mu.chain , sigma.chain = sigma.chain , probacred = probacred , oel = oel  ),
                
                frac.risk = frac.risk.N( mu.chain = mu.chain , sigma.chain = sigma.chain , frac_threshold = frac_threshold , oel = oel  ),
                
                perc = perc.N( mu.chain = mu.chain , sigma.chain = sigma.chain , target_perc = target_perc , probacred = probacred),
                
                perc.risk = perc.risk.N( mu.chain = mu.chain , sigma.chain = sigma.chain , target_perc , oel = oel )
              ))
              
            }
            
            
