#  Version 0.1 (Fev 2022)
                                                                                 

# ==============================================================================
#  This code was first written for the Webexpo project, more precisely
#  for the elication of a continuous-piecewise linear [CPL] prior distrn.
#
#  However, it is generic and could be used in any bisectional search
#  for the MINIMUM value (or its max value) on a specified domain [range] 
#  --- which could be infinite --- if the function is monotonic on both sides 
#  of the solution (its min/max location);
#
#    f(x1) < f(x2) is assessed through the fct better.result(f1, f2) [logical]
#
#    where the function better.result is defined beforehand.
#
#    Hence the function bs.update() could be used to find the maximum value of
#    a function f if better.result(f1,f2) = T  iff  f(x1) >= f(x2)
#  
#    It could also be used to find a root instead, 
#    e.g. find x such that f(x) = f0 
#      if better.result(f1, f2) = T  iff |f(x1)-f0| <= |f(x2)-f0|
#      *if* the function |f(x)-f0| is monotonic on both sides of the solution.  


# Change Log *******************************************************************



bs.init <- function(range, epsilon, 
                    init=start$init, margin=start$margin,
                    start.lo=F, start.hi=F,
                    start=list(init=Inf, margin=0),
                    untouchable.bounds=F, force.1pass=F)
{        
  # Call with force.1pass=T if you want to force continue=T, that is,
  # you want to go through the while-loop (that will inevitable follow a call to bs.init())
  # at least once.
  
  # In the absence of a finite initial value [init],
  # use the - upper value in range if start.hi=T 
  #         - lower value in range if start.lo=T
  #         - middle-range value if both start.lo & start.hi are F [the default]
  #
  # margin: If a null value (0) is specified for margin, then the bisectional search will aim for
  # the closer boundary to the initial value as a second step.   
  
  range <- sort(range)

  if (init < range[1] | init > range[2])  init <- Inf
  
  bounded <- is.finite(range)
  
  if (any(!bounded) & margin == 0 & (is.infinite(init) | any(init == range)))
           stop('margin must be specified in presence of an infinite range.')
  
  if (all(!bounded) & is.infinite(init)) stop('init must be specified in presence of an infinite (both sides) range.')

  started.on.bound <- 0
  

  if (is.infinite(init)) 
  {
    if (any(!bounded))
    {  
      init <- range[bounded]
      started.on.bound <- which(bounded)
    }
    else if (start.lo)
    {  
      init <- range[1]
      started.on.bound <- 1
    }
    else if (start.hi)
    {  
      init <- range[2]
      started.on.bound <- 2
    }
    else init <- mean(range)
  }
  else
  {
    starting.on.bound <- init %in% range
    if (starting.on.bound) started.on.bound <- match(init, range)
  }



  continue <- force.1pass | diff(range) >= epsilon
    # due to numerical imprecision, diff() happens to be negative even when values in range are sorted
    # [and are extremely close]


  list(continue=continue, range=range, value=init, margin=margin, epsilon=epsilon,
       started.on.bound=started.on.bound, untouchable.bounds=untouchable.bounds)
} # end of bs.init


bs.update <- function(bs, score, monitor)
{ 
  next.value.within.x <- function(bs.x, new.score.was.the.best, explored.right.side=T)
  {
    # Used when bs.x = bs$x is of length 3
    
    if (new.score.was.the.best)
    {
      midpoint <- (bs.x[1] + bs.x[3]) / 2
      side2explore <- ifelse(bs.x[2] < midpoint, 1, 3) # shorter side
    }
    else
    {
      side2explore <- ifelse(explored.right.side, 1, 3) # side opposite to last explored
    }
    
    (bs.x[2] + bs.x[side2explore]) / 2
  } # end of next.value.within.x
  
  
  stay.within.range <- function(bs, next.value.proposal)
  { 
    # Used on first entry only (bs$x is of length 1)
    
    w <- 0
    
    if      (next.value.proposal <= bs$range[1])  w <- 1
    else if (next.value.proposal >= bs$range[2])  w <- 2 

    if (w > 0)
    { 
      next.value.proposal <- ifelse(bs$untouchable.bounds, (bs$value + bs$range[w]) / 2, bs$range[w])
      bs$touched.bound.on.2nd.step <- !bs$untouchable.bounds
    }
    
    bs$value <- next.value.proposal    
    
    bs
  } # end of stay.within.range
  
  
  update.info.going.outward <- function(bs, going.east, declare.found.bounds.when.touching.bound=T)
  {         
    # To be used when bs$x is of length 2    
    
    exploring.side <- ifelse(going.east, 2,  1)
    lim <- bs$range[exploring.side]

    if (bs$margin == 0)
    {
      going4bounds <- T 
    }
    else
    {
      bs$margin <- 2 * bs$margin
      direction <- ifelse(going.east, 1, -1)
      next.value <- bs$value + direction * bs$margin
      s <- sign(next.value - lim)
      going4bounds <- s == direction | s == 0
    }
    
    
    if (going4bounds)
    {  
      if (declare.found.bounds.when.touching.bound) bs$found.bounds <- !bs$untouchable.bounds
      bs$value <- ifelse(!bs$untouchable.bounds, lim, (bs$value + lim) / 2)
    }
    else
    {
      bs$value <- next.value
    }
    
       
    bs     
  } # end of update.info.going.outward
 
 
  # ----------------------------------------------------------------------------
        
  if (length(bs$x) == 0)
  {          
    new.score.is.the.best <- T
    
    bs$monitor <- list() # declare it now to have it appearing on top of list
    
    bs$touched.bound.on.2nd.step <- F
    bs$found.bounds <- F

    bs$x <- bs$value
    bs$scores <- matrix(score, ncol=1)
    bs$which.best <- 1 # declare it now to have it appearing just after 'scores'
    

    if (bs$margin == 0)
    {
      if (bs$started.on.bound > 0)  
      {
        w <- 3 - bs$started.on.bound # other bound
        
        if (!bs$untouchable.bounds)
        {
          bs$value <- bs$range[w]
          bs$touched.bound.on.2nd.step <- T
          bs$found.bounds <- T
        }
        else
        {
          bs$value <- mean(bs$range)
        }
      }
      else
      {
        closer.bound <- which.min(abs(bs$value-bs$range))
        lim <- bs$range[closer.bound]
        bs$margin <- abs(lim - bs$value)
        
        bs$touched.bound.on.2nd.step <- !bs$untouchable.bounds
        bs$value <- ifelse(!bs$untouchable.bounds, lim, (bs$value + lim) / 2)
      }
    }
    else
    {
      if (bs$started.on.bound > 0)
      {
        direction <- ifelse(bs$started.on.bound == 1, 1, -1)
        next.value <- bs$value + direction * bs$margin
        bs <- stay.within.range(bs, next.value)
      }
      else
      {
        w <- which.max(abs(bs$range - bs$value))
        direction <- ifelse(w == 2, 1, -1)
        next.value <- bs$value + direction * bs$margin
        bs <- stay.within.range(bs, next.value)
      }
    }
    
    
    bs$continue <- NULL # to remove it from its initial location in the list
    bs$continue <- (bs$range[2] - bs$range[1]) >= bs$epsilon
    
    # Change also the position of 'value' in the list
    bs.value <- bs$value
    bs$value <- NULL
    bs$value <- bs.value
  }
  else
  {
    # length(bs$x) > 0
    new.score.is.the.best <- better.result(score, bs$scores[,bs$which.best])


    if (length(bs$x) == 3)
    {
      explored.right.side <- bs$value > bs$x[2]

      if (new.score.is.the.best)
      {
        j <- ifelse(explored.right.side, 2, 1)

        bs$x <- c(bs$x[j], bs$value, bs$x[j+1])
        bs$scores <- matrix(c(bs$scores[,j], score, bs$scores[,j+1]), ncol=3)
      }
      else
      {
        explored.side <- ifelse(explored.right.side, 3, 1)
        bs$x[explored.side] <- bs$value
        bs$scores[,explored.side] <- score
      }      
      
      bs$value <- next.value.within.x(bs$x, new.score.is.the.best, explored.right.side)
      bs$continue <- (bs$x[3] - bs$x[1]) >= bs$epsilon
    }
    else if (length(bs$x) == 2)
    {        
      within <- bs$value > bs$x[1] & bs$value < bs$x[2]
         
      if (new.score.is.the.best)
      {
        if (within)
        {
          bs$x <- c(bs$x[1], bs$value, bs$x[2])
          bs$scores <- matrix(c(bs$scores[,1], score, bs$scores[,2]), ncol=3)
          bs$which.best <- 2 
          bs$value <- next.value.within.x(bs$x, T)
        }
        else
        {
          going.east <- bs$value > bs$x[2]
          
          if (going.east)
          {
            bs$x <- c(bs$x[2], bs$value)
            bs$scores <- matrix(c(bs$scores[,2], score), ncol=2)
          }
          else
          {
            bs$x <- c(bs$value, bs$x[1])
            bs$scores <- matrix(c(score, bs$scores[,1]), ncol=2)
          }
          
          
          if (bs$found.bounds)  bs$value <- mean(bs$x)          
          else                  bs       <- update.info.going.outward(bs, going.east)
        }
      }
      else
      {
        # Not a new best score
        
        if (within)
        {
          # The new value & score replace those which were previously in 2nd place
          j <- 3 - bs$which.best  # j = 2nd place index
          bs$x[j] <- bs$value
          bs$scores[,j] <- score
          bs$value <- mean(bs$x)
        }
        else
        {
          going.east <- bs$value > bs$x[2]
          
          if (going.east)
          {
            bs$x <- c(bs$x, bs$value)
            bs$scores <- cbind(bs$scores, score)
          }
          else
          {
            bs$x <- c(bs$value, bs$x)
            bs$scores <- cbind(score, bs$scores)
          }

          
          # bs$x is now of length 3: bs$which.best & bs$found.bounds are necessary anymore (they are known --- 2 & T respectively)
          bs$which.best <- 2 # even if it is known, I set it to its (known) value here, in case fcts using the output of bs.update() are based on it
          
          bs$value <- next.value.within.x(bs$x, T) # do as if the last calculation was for the point in the middle of bs$x [hence with the best score so far]
        }  
      }
      
      
      if (length(bs$x) == 3)  bs$continue <- (bs$x[3] - bs$x[1]) >= bs$epsilon
      else                    bs$continue <- diff(bs$x)          >= bs$epsilon
    }
    else
    {
      # length(bs$x) == 1

      going.east <- bs$value > bs$x
      
      
      if (going.east)
      {
        bs$x <- c(bs$x, bs$value)
        bs$scores <- cbind(bs$scores, score)
      }
      else
      {
        bs$x <- c(bs$value, bs$x)
        bs$scores <- cbind(score, bs$scores)
      }


      if (new.score.is.the.best)
      {
        if (bs$touched.bound.on.2nd.step)
        {
          bs$found.bounds <- T
          bs$value <- mean(bs$x)
          bs$continue <- diff(bs$x) >= bs$epsilon
        }
        else
        {
          bs <- update.info.going.outward(bs, going.east, F)
        }
        
        bs$which.best <- ifelse(going.east, 2, 1)
      }
      else if (bs$started.on.bound > 0)
      {
        bs$value <- mean(bs$x)
        bs$which.best <- ifelse(going.east, 1, 2) 
      }
      else
      {
        # we change direction
        bs$value <- ifelse(going.east, bs$x[1], bs$x[2])    # do as if this was the last step performed
        bs <- update.info.going.outward(bs, !going.east)
        
        bs$which.best <- ifelse(going.east, 1, 2)
      }
    }
  }


  if (new.score.is.the.best)  bs$monitor <- monitor

  bs
} # end of bs.update
