##############################################################
#
#   WEBEXPO official R scripts
#
#   BETWEEN WORKER ANALYSIS
#  
#   Creation of summary result tables 
#
#   V1.0 Sept 2018  
#
#############################################################


##
#
## INPUT 
#
# distributional option (lognormal , normal)
#
# list of output of Webexpo.between.interpretation
#

### WRPPING FUNCTION

webexpo.between.summary <- function( labels ,  result.list , is.lognormal ) {
  
  # distribution specific table of results
  
  if (is.lognormal) results <- data.frame(names=c("group.gm" , "gsd.between" , "gsd.within" ,  "rho" , "R.ratio", " prob.ind.perc", "prob.ind.am" ), stringsAsFactors = FALSE)
 
  if (!is.lognormal) results <- data.frame(names=c("group.am" , "sd.between" , "sd.within" ,  "rho" , "R.diff", " prob.ind.perc", "prob.ind.am" ), stringsAsFactors = FALSE)
 
  n <-length(result.list)
  
  for (i in 1:n) {
    
    if (is.lognormal) results <- cbind( results , fun.summary.between.LN( result.list[[i] ]))
    
    if (!is.lognormal) results <- cbind( results , fun.summary.between.N( result.list[[i] ]))
    
    
  }
  
  names(results) <- c("names" , labels) 
  
  return(results)
  
}





###################  SUPPORT FUNCTIONS

## LOGNORMAL

#input is one object from Webexpo.between.interpretation

fun.summary.between.LN <- function( num.res ) {
  
  result <-character(7)
  
  result[1] <- paste(signif(num.res$group.gm$est,3),
                     " [ ",
                     signif(num.res$group.gm$lcl,3),
                     " - ",
                     signif(num.res$group.gm$ucl,3),
                     " ]", sep="")
  
  result[2] <- paste(signif(num.res$gsdb$est,3),
                     " [ ",
                     signif(num.res$gsdb$lcl,3),
                     " - ",
                     signif(num.res$gsdb$ucl,3),
                     " ]", sep="")
  
  result[3] <- paste(signif(num.res$gsdw$est,3),
                     " [ ",
                     signif(num.res$gsdw$lcl,3),
                     " - ",
                     signif(num.res$gsdw$ucl,3),
                     " ]", sep="")
  
  result[4] <- paste(signif(num.res$rho$est,3),
                     " [ ",
                     signif(num.res$rho$lcl,3),
                     " - ",
                     signif(num.res$rho$ucl,3),
                     " ] Risk: ",
                     signif(num.res$prob.rho.overX,2), 
                      "%", sep="")
  
  result[5] <- paste(signif(num.res$R.ratio$est,3),
                     " [ ",
                     signif(num.res$R.ratio$lcl,3),
                     " - ",
                     signif(num.res$R.ratio$ucl,3),
                     " ] Risk: ",
                     signif(num.res$prob.R.over2,2), 
                     "%",
                     signif(num.res$prob.R.over10,2), 
                     "%",sep="") 
  
  result[6] <- paste(signif(num.res$prob.ind.overexpo.perc$est,3),
                     " [ ",
                     signif(num.res$prob.ind.overexpo.perc$lcl,3),
                     " - ",
                     signif(num.res$prob.ind.overexpo.perc$ucl,3),
                     " ] Risk: ",
                     signif(num.res$prob.ind.overexpo.perc.proboverX,2), 
                     "%",sep="") 
  
  
  result[7] <- paste(signif(num.res$prob.ind.overexpo.am$est,3),
                     " [ ",
                     signif(num.res$prob.ind.overexpo.am$lcl,3),
                     " - ",
                     signif(num.res$prob.ind.overexpo.am$ucl,3),
                     " ] Risk: ",
                     signif(num.res$prob.ind.overexpo.am.proboverX,2), 
                     "%",sep="") 
  

  return(result)
}


## NORMAL

#input is one object from Webexpo.between.interpretation

fun.summary.between.N <- function( num.res ) {
  
  result <-character(7)
  
  result[1] <- paste(signif(num.res$group.mean$est,3),
                     " [ ",
                     signif(num.res$group.mean$lcl,3),
                     " - ",
                     signif(num.res$group.mean$ucl,3),
                     " ]", sep="")
  
  result[2] <- paste(signif(num.res$asdb$est,3),
                     " [ ",
                     signif(num.res$asdb$lcl,3),
                     " - ",
                     signif(num.res$asdb$ucl,3),
                     " ]", sep="")
  
  result[3] <- paste(signif(num.res$asdw$est,3),
                     " [ ",
                     signif(num.res$asdw$lcl,3),
                     " - ",
                     signif(num.res$asdw$ucl,3),
                     " ]", sep="")
  
  result[4] <- paste(signif(num.res$rho$est,3),
                     " [ ",
                     signif(num.res$rho$lcl,3),
                     " - ",
                     signif(num.res$rho$ucl,3),
                     " ] Risk: ",
                     signif(num.res$prob.rho.overX,2), 
                     "%", sep="")
  
  result[5] <- paste(signif(num.res$R.diff$est,3),
                     " [ ",
                     signif(num.res$R.diff$lcl,3),
                     " - ",
                     signif(num.res$R.diff$ucl,3),
                     " ]" ,sep="") 
  
  result[6] <- paste(signif(num.res$prob.ind.overexpo.perc$est,3),
                     " [ ",
                     signif(num.res$prob.ind.overexpo.perc$lcl,3),
                     " - ",
                     signif(num.res$prob.ind.overexpo.perc$ucl,3),
                     " ] Risk: ",
                     signif(num.res$prob.ind.overexpo.perc.proboverX,2), 
                     "%",sep="") 
  
  
  result[7] <- paste(signif(num.res$prob.ind.overexpo.am$est,3),
                     " [ ",
                     signif(num.res$prob.ind.overexpo.am$lcl,3),
                     " - ",
                     signif(num.res$prob.ind.overexpo.am$ucl,3),
                     " ] Risk: ",
                     signif(num.res$prob.ind.overexpo.am.proboverX,2), 
                     "%",sep="") 
  
  
  return(result)
}