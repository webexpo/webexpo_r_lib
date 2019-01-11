##############################################################
#
#   WEBEXPO official R scripts
#
#   BETWEEN WORKER ANALYSIS
#  
#   DATA PREPARATION BEFORE BAYESIAN ANALYSIS 
#
#   V1.0  Sept 2018
#
#############################################################

##
#
#  INPUT : 
#
# data.frame of  observations 
# should be data.frame with 2 variables:
# x : a vector of strings, with observations as e.g. "2.1" , "<5" , "[2.3-5]" , ">23.1"
# worker : a vector of worker identifiers
#
#  
#
##


webexpo.between.datapreparation <-function(data.in) {
  
  
  
  result <-list()
  
  result$data <-data.in

  result$leftcensored <-grepl('<' , data.in$x, fixed = TRUE) 
  result$rightcensored <-grepl('>' , data.in$x , fixed = TRUE)
  result$intcensored <-grepl('[' , data.in$x , fixed = TRUE)
  result$notcensored <-!(result$leftcensored | result$rightcensored | result$intcensored)
  
  return(result)
  
}
