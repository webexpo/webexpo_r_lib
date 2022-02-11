################################################################################################################################
#
#
#        NDEXPO R FUNCTION
#
#
#  The following function allows to process left-censored lognormal data in the same manner as NDexpo
#
#  Visit www.expostats.ca for complete documentation on NDEXPO and the Regression on order statistics (ROS) approach to 
#  censored data treatment
#
#  Also consider using the Expostats / Webexpo instread of NDexpo : Bayesian statistics allow the optimal treatement of NDs along
#  with other advantages. Expostats allow direct online calculations, while Webexpo includes a library of R code.
#
#  In the minds of its author, NDexpo is "deprecated"
#
#
#  See : https://github.com/webexpo/webexpo_r_lib
#
#  
#
#  FUNCTION PARAMETERS
#
#
##input : vector of measurement values formatted as in NDexpo : character string including numerical values and <x values (see example)

##output :
#
# A list with 3 elements : 
#
#$data is a data.frame that contains (in the same order as initial data vector)
	    #$xfin : the final predicted values
	    #$yfin : the final predicted log-transformed values
	    #$is.censored : index of which points are censored
	    #$pp : plotting positions
	    #$order.index : index for ordering the value in increasing order
	    #$x : initial data (censored value are identified using $is.censored)
#$alpha : intercept for the QQ regression line
#$beta : slope for the QQ regression line



##### example (The default example in Expostats TOOL 1)

mydata <-c("28.9","19.4","<5.5","149.9","26.42","56.1")


fun.NdExpo.lognorm( mydata )$data$xfin


#START DEFINING FUNCTION

fun.NdExpo.lognorm <-function(x){

        
        #is.censored is 0 for censored, 1 for observed
        
        is.censored <- as.integer(!grepl( '<',x))     #0 means censored
        
        #eliminating the  '<' sign
        
        x[regexpr('<',x)==1] <-as.numeric(substring(x[regexpr('<',x)==1],2))
        
        x <-as.numeric(x)
        
        ###ordering X and keeping the order
        
        x.order <-rank(x, ties.method='first')
        
        x.ordered <-numeric(length(x))
        x.ordered[x.order] <-x
        
        is.censored[x.order] <- is.censored
        
######################IF NO NDs #################################

        if (length(x[is.censored==0])==0) {
        
        y <- log(as.numeric(x.ordered))
        
        pp <-((1:length(x))-0.375)/(length(x)+0.25)
        
        
        
        res <-list(data=data.frame(xfin=x,yfin=y[x.order],is.censored=is.censored[x.order],pp=pp[x.order],order.index=x.order,x=x.ordered[x.order]),
                      alpha=0,beta=0)
        
        
        
        }

#########when there are NDs

      else {
      
      
      ####correcting adjacency : The NDExpo approach (see documentation on the NDexpo website)
      
      ##first ordering so that no case of 5 <5 : to become <5 5
      
      for (i in (length(x)):2) {
      
              if (is.censored[i]==0 & is.censored[i-1]==1 & x.ordered[i]==x.ordered[i-1]) {
      
                                      is.censored[i] <-1
                                      is.censored[i-1] <-0
      
                                      temp <-x.order[i-1]
                                      x.order[i-1] <- x.order[i]
                                      x.order[i] <- temp
      
                                  }
      
      }


      #####correcting adjacency (see documentation on the NDexpo website)

      for (i in 1:(length(x)-1) ) {
      
                          if (is.censored[i]==0 & is.censored[i+1]==0 & x.ordered[i]!=x.ordered[i+1]) {
      
                          x.ordered[i+1] <-  x.ordered[i]
      
                          }
                          }
      

      #log transformation
      
      y <- log(as.numeric(x.ordered))
      
      
      #nb detection limits
      
      NL <-length(unique(y[is.censored==0]))
      
      
      
      #######the case where NL is 1
      
      if (NL==1) {
      
      NbA <- length(y[is.censored==1 & y>=min(y[is.censored==0])])
      NbB <- length(y[is.censored==1 & y<min(y[is.censored==0])])
      NbC <-length(y[is.censored==0])
      
      index.A <-1:NbA
      index.B <-1:NbB
      index.C <-1:NbC
      
      Pdep <-NbA/(NbB+NbC+NbA)
      
      if (NbA!=0) PdA <-(1-Pdep)+Pdep*index.A/(NbA+1) else PdA <-numeric(0)
      
      if (NbB!=0)PdB <-(1-Pdep)*index.B/(NbB+1) else PdB <-numeric(0)
      
      if (NbC!=0)PdC <-(1-Pdep)*index.C/(NbC+1) else PdC <-numeric(0)
      
      pp <-c(PdB,PdC,PdA)
      
                  }
      
      
      else {
      
      #vector of limits
      
      L <-unique(y[is.censored==0][order(y[is.censored==0])])
      
      L <-c(-100000,L,10000)
      
      #vector of detected values between each LD
      
      NbA <-numeric(NL+1)
      
      for (i in 1:(NL+1)) NbA[i] <-length(y[is.censored==1 & y>=L[i] & y<L[i+1]])
      
      #vector of all values smaller than Li
      
      NbB.ND <-numeric(NL+1)    #inequality  is lenient
      
      NbB.D  <-numeric(NL+1)    #inequality  is strict
      
      
      for (i in 1:(NL+1)) NbB.ND[i] <-length(y[is.censored==0 & y<=L[i]])
      for (i in 1:(NL+1)) NbB.D[i] <-length(y[is.censored==1 & y<L[i]])
      
      NbB <-NbB.ND+NbB.D
      
      #vector of ND values at Lj
      
      NbC <-numeric(NL)
      
      for (i in 1:(NL)) NbC[i] <-NbB[i+1]-NbB[i]-NbA[i]
      
      NbC <-c(NbC,0)
      
      ##probability of exceedance for limit i
      
      Pdep <-  numeric(NL+2)
      Pdep[1] <-1
      Pdep[NL+2] <-0
      
      for (i in (NL+1):2) Pdep[i] <-Pdep[i+1] + (1-Pdep[i+1])*NbA[i]/(NbA[i]+NbB[i])
      
      
      ######vector of plotting positions
      
      pp <-numeric(length(y))
      
      
      ###forming the detected groups
      
      start <-1
      
      for (j in 1:(NL+1)) {
      
      
      
          if (NbA[j]!=0) {
      
          Pdet <-numeric(NbA[j])
      
          for (i in 1:NbA[j]) Pdet[i] <-(1-Pdep[j])+i*(Pdep[j]-Pdep[j+1])/(NbA[j]+1)
      
          index <-start:(start+NbA[j]-1)
      
          pp[index] <-Pdet
      
                          }
      
          start <-start+NbA[j]+NbC[j]
      
      
                          }
      
      
       ### forming the not detected groups
      
       start <-1+NbA[1]
      
      for (j in 1:(NL)) {
      
      
          if (NbC[j]!=0) {
      
          Pdet <-numeric(NbC[j])
      
          for (i in 1:NbC[j]) Pdet[i] <-i*(1-Pdep[j+1])/(NbC[j]+1)
      
          index <-start:(start+NbC[j]-1)
      
          pp[index] <-Pdet
      
                          }
          start <-start+NbA[j+1]+NbC[j]
      
      
                          }
      
      }
      #########predicting the NDs
      
      normal.scores <-qnorm(pp)
      
      regression <-lm(y[is.censored==1]~normal.scores[is.censored==1])
      
      y[is.censored==0] <- regression$coefficients[1]+regression$coefficients[2]*normal.scores[is.censored==0]

#########final result preparation
      
      #####
      
      #$data is a data.frame that contains (in the same order as initial)
      	#$xfin : the final predicted values
      	#$yfin : the final predicted log-transformed values
      	#$is.censored : index of which points are censored
      	#$pp : plotting positions
      	#$order.index : index for ordering the value in increasing order
      	#$x : initial data (censored value are identified using $is.censored)
      #$alpha : intercept for the QQ regression line
      #$beta : slope for the QQ regression line
      
      res <-list(data=data.frame(xfin=exp(y[x.order]),yfin=y[x.order],is.censored=is.censored[x.order],pp=pp[x.order],order.index=x.order,x=x.ordered[x.order]),
                    alpha=regression$coefficients[1],beta=regression$coefficients[2])
      
                    }

return(res)
}

#END DEFINING FUNCTION



#########################################################################################################################

