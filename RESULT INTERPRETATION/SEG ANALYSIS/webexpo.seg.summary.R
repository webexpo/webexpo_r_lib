##############################################################
#
#   WEBEXPO official R scripts
#
#   SEG ANALYSIS
#  
#   Creation of summary result tables to compare several analyses 
#
#   V1.0 Sept 2018 
#   
#
#############################################################


##
#
## INPUT 
#
# distributional option (lognormal , normal)
#
# list of outputs from the Webexpo.seg.interpretation function
#



#### WRAPPING FUNCTION

webexpo.seg.summary <- function( labels ,  result.list , is.lognormal ) {
  
  
# table of results specific to the distribution  
  
 if (is.lognormal) results <- data.frame(names=c("gm" , "gsd" , "frac" ,  "perc" , "perc.band", "am" , "am.band" ), stringsAsFactors = FALSE)
 
 if (!is.lognormal) results <- data.frame(names=c("am" , "asd" , "frac" ,  "perc" , "am" ), stringsAsFactors = FALSE)
 
 
  n <-length(result.list)
  
  for (i in 1:n) {
    
    if (is.lognormal) results <- cbind( results , fun.summary.seg.LN( result.list[[i] ]))
    
    if (!is.lognormal) results <- cbind( results , fun.summary.seg.N( result.list[[i] ]))
    
    
  }
  
  names(results) <- c("parameter", labels)
  
  return(results)
  
}





###################  SUPPORT FUNCTIONS

## LOGNORMAL DISTRIBUTION

#input is one object from Webexpo.seg.interpretation

fun.summary.seg.LN <- function( num.res ) {
  
  result <-character(7)
  
  result[1] <- paste(signif(num.res$gm$est,3),
                     " [ ",
                     signif(num.res$gm$lcl,3),
                     " - ",
                     signif(num.res$gm$ucl,3),
                     " ]", sep="")
  
  result[2] <- paste(signif(num.res$gsd$est,3),
                     " [ ",
                     signif(num.res$gsd$lcl,3),
                     " - ",
                     signif(num.res$gsd$ucl,3),
                     " ]", sep="")
  
  result[3] <- paste(signif(num.res$frac$est,3),
                     " [ ",
                     signif(num.res$frac$lcl,3),
                     " - ",
                     signif(num.res$frac$ucl,3),
                     " ] Risk: ",
                     signif(num.res$frac.risk,2), 
                      "%", sep="")
  
  result[4] <- paste(signif(num.res$perc$est,3),
                     " [ ",
                     signif(num.res$perc$lcl,3),
                     " - ",
                     signif(num.res$perc$ucl,3),
                     " ] Risk: ",
                     signif(num.res$perc.risk,2), 
                     "%", sep="") 
  
  result[5] <- paste(signif(num.res$perc.riskbands,2), 
                     collapse = " / ") 
  
  
  result[6] <- paste(signif(num.res$am$est,3),
                     " [ ",
                     signif(num.res$am$lcl,3),
                     " - ",
                     signif(num.res$am$ucl,3),
                     " ] Risk: ",
                     signif(num.res$am.risk,2), 
                     "%", sep="")
  
  result[7] <- paste(signif(num.res$am.riskbands,2), 
                     collapse = " / ") 
  
  return(result)
}


## NORMAL DISTRIBUTION

#input is one object from Webexpo.seg.interpretation


fun.summary.seg.N <- function( num.res ) {
  
  result <-character(5)
  
  result[1] <- paste(signif(num.res$am$est,3),
                     " [ ",
                     signif(num.res$am$lcl,3),
                     " - ",
                     signif(num.res$am$ucl,3),
                     " ]", sep="")
  
  result[2] <- paste(signif(num.res$asd$est,3),
                     " [ ",
                     signif(num.res$asd$lcl,3),
                     " - ",
                     signif(num.res$asd$ucl,3),
                     " ]", sep="")
  

  result[3] <- paste(signif(num.res$frac$est,3),
                     " [ ",
                     signif(num.res$frac$lcl,3),
                     " - ",
                     signif(num.res$frac$ucl,3),
                     " ] Risk: ",
                     signif(num.res$frac.risk,2), 
                     "%", sep="")
  
  result[4] <- paste(signif(num.res$perc$est,3),
                     " [ ",
                     signif(num.res$perc$lcl,3),
                     " - ",
                     signif(num.res$perc$ucl,3),
                     " ] Risk: ",
                     signif(num.res$perc.risk,2), 
                     "%", sep="") 
  
  result[5] <- paste(signif(num.res$am$est,3),
                     " [ ",
                     signif(num.res$am$lcl,3),
                     " - ",
                     signif(num.res$am$ucl,3),
                     " ] Risk: ",
                     signif(num.res$am.risk,2), 
                     "%", sep="")

  
  
  return(result)
  
}