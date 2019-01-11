##############################################################
#
#   WEBEXPO official R scripts
#
#   SEG ANALYSIS
#  
#   DATA PREPARATION BEFORE BAYESIAN ANALYSIS 
#
#   V1.0    Sept. 2018
#
#############################################################

##
#
#  INPUT : 
#
# vector of  observations 
# should be a vector of strings, with observations as e.g. "2.1" , "<5" , "[2.3-5]" , ">23.1"
#
#
#
##


webexpo.seg.datapreparation <-function(data.in) {
  
  
  result <-list()
  
  result$data <-data.in
  
  result$leftcensored <-grepl('<' , data.in, fixed = TRUE) 
  result$rightcensored <-grepl('>' , data.in , fixed = TRUE)
  result$intcensored <-grepl('[' , data.in , fixed = TRUE)
  result$notcensored <-!(result$leftcensored | result$rightcensored | result$intcensored)
  
  return(result)
  
}
