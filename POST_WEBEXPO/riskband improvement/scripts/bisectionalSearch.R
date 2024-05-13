
# Version 1.1 (Jun 2023)
#		Shared/distributed: yes


# -----------------------------------
                                                                                 
                                                                                 
# ==============================================================================
#  This code was first written for the Webexpo project, more precisely
#  for the elication of a continuous-piecewise linear [CPL] prior distrn.
#
#  However, it is generic and could be used in any bisectional search
#  for the MINIMUM value (or its max value) on a specified domain [range] 
#  --- which could be infinite --- if the function is monotonic on both sides 
#  of the solution (its min/max location);
#
#    f(x1) < f(x2) is assessed through the fct better.score(f1, f2) [logical]
#
#    where the function better.score is defined through an argument of
#    the function bs.init
#
#    Hence the function bs.update() could be used to find the maximum value of
#    a function f if bs$better.score(f1,f2) = TRUE  iff  f(x1) >= f(x2)
#  
#    It could also be used to find a root instead, 
#    e.g. find x such that f(x) = f0 
#      if bs$better.score(f1, f2) = TRUE  iff |f(x1)-f0| <= |f(x2)-f0|
#      *if* the function |f(x)-f0| is monotonic on both sides of the solution.
#
#
#     Since version 0.8,
#     we allow the user to indicate in what direction (starting from current value $x)
#     the optimal solution is known to be found.
#
#
#     dir: the direction in which the optimal solution is KNOWN to sit, relatively to the current visited point ($x).
#          In other words, dir = the opposite sign of the gradiant [of the scoring fct to minimize] at the visited point
#          Leave dir to its default value [0] when it is not known.
#
#          dir =  1 -> look for a value > bs$x [current value]
#              = -1 -> look for a value < bs$x
#              =  0 -> solution is on either side of bs$x [leave the algorithm search on both sides of it]


# Change Log *******************************************************************
#
#
# Version 1.1 (Jun 2023)
# ----------------------
#
#
# Version 1.0 (May 2023)
# -----------------------
#  First official version.
#
#  We have changed the parameter 'margin' for 'change'
#  and changed its default value (now Inf).
#  We have also modified the first steps to avoid falling in an infinite loop
#  when range is infinite and 'change' left to its default infinite value.
#  We also prevent starting on a range limit when either of start.hi or start.lo is TRUE
#  but the corresponding boundary is 'untouchable' (as indicated by untouchable.bounds).
#
#
# Version 0.9 (Apr 2023)
# ----------------------
#   We have removed the limitation of the previous version regarding the (before impossible)
#   alternation between known & unknown direction (towards best solution) --- it finally proves useful
#   in the context of CPL.prior.
#   That (almost!) involved a rewriting of bs.update from scratch.
#
#
# Version 0.8 (Mar 2023)
# ----------------------
#   Added code (in fct bs.update) allowing the user to indicate in which direction
#   the optimal solution is known to be sitting.
#     (IMPORTANT: Note however that this algorithm must ALWAYS be called with either dir = 0 or dir != 0
#                [that is, you cannot alternate between known & unknown directions --- if you need to do so, you will
#                need to modify bs.update code])
#   This version was known to work perfectly (albeit the limitation noted above).
#
#
# Version 0.7 (Mar 2023)
# ----------------------
#   Removed touched.bound.on.2nd.step from the code
#   Removed the argument declare.found.bounds.when.touching.bound in function update.info.going.outward
#     as it made the output odd-looking, indicating found.bounds=TRUE when it is not true *yet*,
#     but will only be true at next step.
#   This version was known to work perfectly.
#
#
# Version 0.6 (Mar 2023)
# ----------------------
#   Added the possibility to specify the second step visited 'x'
#     [through the argument suggested.2nd.step]
#
#
# Version 0.5 (Mar 2023)
# ----------------------
#   Added a protection against an initial margin < epsilon
#   Added an automatic incrementation of margin when #steps grows and !found.bounds
#   Added a epsilon.safety.mult parameter
#   Added last.score.was.the.best to the output of bs.update()
#   Added the fct bs.best_score()
#
#
# Version 0.4 (Sep 2022)
# ----------------------
#   Added monitoring & use of x.on.bound & newx.on.bound
#
#
# Version 0.3 (Jun 2022)
# ----------------------
#   'monitor' was renamed this.parms & best.parms
#   Added first.step.was.going.east
#   Modified next side to explore in fct next.x.within.visited
#
#
# Version 0.2 (Feb 2022)
# ----------------------
#  The fct used to compare two outcomes must now be passed as an argument
#  to the function bs.init()
#
#                                                            (End of Change Log)


bs.best_score <- function(bs)
{
  return(bs$scores[,bs$which.best])
} # end of bs.best_score



bs.init <- function(range, epsilon,
                    init=start$init, change=start$change,
                    start.lo=FALSE, start.hi=FALSE,
                    start=list(init=Inf, change=Inf),
                    untouchable.bounds=rep(FALSE, 2), force.1pass=FALSE,
                    better.score=function(out1, out2){out1<out2},
                    disallowed.proximity=5*epsilon,
                    suggested.2nd.step=NULL,
                    epsilon.safety.mult=100)
{
  # Call with force.1pass=TRUE if you want to force continue=TRUE, that is,
  # you want to go through the while-loop (that will inevitably follow a call to bs.init())
  # at least once.

  # In the absence of a finite initial value [init],
  # use the - upper value in range if start.hi=TRUE
  #         - lower value in range if start.lo=TRUE
  #         - middle-range value if both start.lo & start.hi are FALSE [the default]
  #
  # change: If infinity (Inf) is specified for change, then the bisectional search will aim for
  # the closer boundary to the initial value as a second step.
  #
  # better.score: a boolean function comparing two outcomes, out1 & out2, returning TRUE if out1 is better than out2
  #
  # disallowed.proximity: if a starting point not falling on a BORDER of the variable range
  #                       is within that distance from a border, it will not be allowed
  #                       (and the algorithm will set the initial value to be equal to the default starting point)
  #
  # suggested.2nd.step:   The value of 'x' on the second step of the bisectional search.
  #                       For example, the initial value [init] may be calculated from a regression model as a best-guess
  #                       to the optimal solution; however, sometimes we may not be sure of the appropriateness of our
  #                       initial value and want to suggest a value for 'x' on the 2nd step which is mathematically sensible
  #                       hence making sure to return a mathematically sensible optimal value.


  range <- sort(range)


  if (init < range[1] | init > range[2])
  {
    init <- Inf
  }
  else
  {
    side <- 0

    if      (abs(init - range[1]) < disallowed.proximity)   side <- 1
    else if (abs(init - range[2]) < disallowed.proximity)   side <- 2

    if (side > 0)
    {
      if (untouchable.bounds[side])
      {
        init <- Inf
      }
      else
      {
        init <- range[side]
        started.on.boundj <- side
      }
    }
  }


  bounded <- is.finite(range)


  if (all(!bounded) & is.infinite(init)) stop('init must be specified in presence of an infinite (both sides) range.')

  started.on.boundj <- 0


  if (is.infinite(init))
  {
    if (any(!bounded))
    {
      init <- range[bounded]
      started.on.boundj <- which(bounded)
    }
    else if (start.lo)
    {
      if (!untouchable.bounds[1])
      {
        init <- range[1]
        started.on.boundj <- 1
      }
      else if (!bounded[2])
      {
        init <- range[1] + epsilon * epsilon.safety.mult
      }
      else
      {
        init <- range[1] + diff(range)/1000
      }
    }
    else if (start.hi)
    {
      if (!untouchable.bounds[2])
      {
        init <- range[2]
        started.on.boundj <- 2
      }
      else if (!bounded[1])
      {
        init <- range[2] - epsilon * epsilon.safety.mult
      }
      else
      {
        init <- range[2] - diff(range)/1000
      }
    }
    else init <- mean(range)
  }
  else
  {
    starting.on.bound <- init %in% range
    if (starting.on.bound) started.on.boundj <- match(init, range)
  }


  continue <- force.1pass | diff(range) >= epsilon
    # due to numerical imprecision, diff() happens to be negative even when values in range are sorted
    # [and are extremely close]

  found.bounds <- rep(FALSE, 2)
  if (started.on.boundj > 0)  found.bounds[started.on.boundj] <- TRUE


  list(continue=continue, range=range, x=init, change=change, epsilon=epsilon,
       untouchable.bounds=untouchable.bounds,
       better.score=better.score, suggested.2nd.step=suggested.2nd.step,
       found.Bounds=FALSE, found.bounds=found.bounds, x.on.bound=any(found.bounds),
       epsilon.safety.mult=epsilon.safety.mult, n.steps=0)
} # end of bs.init


bs.update <- function(bs, score, this.parms, dir=0)
{
  # bs: a bs.init()- or bs.update()-returned object
  # score: can be an array or a scalar -- normally the result of a scoring function evaluated at bs$x
  # this.parms: can be just about anything -- bs$best.parms will be = to this.parms in the returned object iff 'score' was the best score so far
  # epsilon.safety.mult: a multiplication factor for epsilon as a safe minimal initial search change (around the target parm)
  #
  #   dir: the direction in which the optimal solution is KNOWN to sit, relatively to the current visited point ($x)$
  #        In other words, dir = the opposite sign of the gradiant [of the scoring fct to minimize] at the visited point
  #        Leave dir to its default value [0] when it is not known
  #
  #   dir =  1 -> look for a value > bs$x [current value]
  #       = -1 -> look for a value < bs$x
  #       =  0 -> solution is on either side of bs$x [leave the algorithm search on both sides of it]


  Clean.visits <- function(bs, clean.right.side, n=1)
  {
    if (clean.right.side)
    {
      bs$visited <- pop(bs$visited, n)
      bs$scores  <- pop.matrix(bs$scores, n)
    }
    else
    {
      bs$visited <- bs$visited[-seq(n)]
      bs$scores  <- bs$scores[, -seq(n), drop=FALSE]
    }

    return(bs)
  } # end of Clean.visits


  Insert.new.best <- function(bs, score)
  {
    bs$visited <- c(bs$visited[1], bs$x, bs$visited[2])
    bs$scores  <- matrix(c(bs$scores[,1], score, bs$scores[,2]), ncol=3)
    bs$which.best <- 2

    return(bs)
  } # end of Insert.new.best


  playSafe.change <- function(bs)
  {
    min.change <- bs$epsilon.safety.mult * bs$epsilon
    if (bs$change < min.change)  bs$change <- min.change

    return(bs)
  } # end of playSafe.change


  pop <- function(x, n=1)
  {
    j <- seq(to=length(x), length=n)
    return(x[-j])
  } # end of pop


  pop.matrix <- function(m, n=1)
  {
    j <- seq(to=ncol(m), length=n)
    return(m[,-j,drop=FALSE])
  } # end of pop.matrix


  Register.initial.visit <- function(bs, score)
  {
    bs$visited <- bs$x
    bs$scores <- matrix(score, ncol=1)
    bs$which.best <- 1

    if (bs$x.on.bound)  bs <- Updated.bound.info(bs, bs$x == bs$range[2])
    return(bs)
  } # end of Register.initial.visit


  Register.new.visit <- function(bs, score, right.side, new.score.is.the.best=TRUE)
  {
    # To register a new visit on either side of the domain already visited

    if (right.side)
    {
      bs$scores <- cbind(bs$scores, score)
      bs$visited <- c(bs$visited, bs$x)
    }
    else
    {
      bs$scores <- cbind(score, bs$scores)
      bs$visited <- c(bs$x, bs$visited)
    }


    if (new.score.is.the.best)
    {
      bs$which.best <- ifelse(right.side, length(bs$visited), 1)
    }
    else if (!right.side)  bs$which.best <- bs$which.best + 1


    if (bs$x.on.bound)  bs <- Updated.bound.info(bs, right.side)
    return(bs)
  } # end of Register.new.visit


  Replace.visit <- function(bs, score, right.side=TRUE)
  {
    j <- ifelse(right.side, length(bs$visited), 1)

    bs$visited[j] <- bs$x
    bs$scores[,j] <- score

    if (bs$x.on.bound)  bs <- Updated.bound.info(bs, right.side)

    return(bs)
  } # end of Replace.visit


  Stay.within2visits <- function(bs)
  {
    bs$continue <- (bs$visited[2] - bs$visited[1]) >= bs$epsilon

    if (bs$continue)
    {
      bs$x <- mean(bs$visited)
      bs$x.on.bound <- FALSE
    }

    return(bs)
  } # end of Stay.within2visits


  Stay.within.range <- function(bs, go.east)
  {
    if      ( go.east && !bs$found.bounds[2] && bs$untouchable.bounds[2] && (bs$range[2]-bs$x) < bs$epsilon)
    {
      bs$continue <- FALSE
      return(bs)
    }
    else if (!go.east && !bs$found.bounds[1] && bs$untouchable.bounds[1] && (bs$x-bs$range[1]) < bs$epsilon)
    {
      bs$continue <- FALSE
      return(bs)
    }

    if (!bs$found.Bounds && bs$n.steps > 5 && bs$change < diff(bs$range)/2)  bs$change <- 2 * bs$change

    if      ( go.east && bs$found.bounds[2])  go.east <- FALSE
    else if (!go.east && bs$found.bounds[1])  go.east <- TRUE


    direction <- ifelse(go.east, 1, -1)
    next.x.proposal <- bs$x + direction * bs$change
    w <- 0


    if      (next.x.proposal <= bs$range[1])  w <- 1
    else if (next.x.proposal >= bs$range[2])  w <- 2

    if (w > 0)
    {
      next.x.proposal <- ifelse(bs$untouchable.bounds[w],  (bs$x + bs$range[w]) / 2, bs$range[w])
      bs$x.on.bound   <- !bs$untouchable.bounds[w]
    }
    else bs$x.on.bound <- FALSE


    # Prevent going to Infinity (when change=Inf) when Infinity is not 'touchable'

    if (is.infinite(next.x.proposal))
    {
      j <- ifelse(go.east, 2, 1)

      if (bs$untouchable.bounds[j])
      {
        bs$change <- ifelse(length(bs$visited) > 1, bs$visited[2] - bs$visited[1], 0)
        bs <- playSafe.change(bs)

        next.x.proposal <- bs$x + direction * bs$change
      }
    }



    bs$x <- next.x.proposal

    return(bs)
  } # end of Stay.within.range


  Step2Longer.side <- function(bs)
  {
    # Used when bs$visited is of length 3

    midpoint     <- (bs$visited[1] + bs$visited[3]) / 2
    side2explore <- ifelse(bs$visited[2] < midpoint, 3, 1) # larger side

    bs$x <- (bs$visited[2] + bs$visited[side2explore]) / 2
    bs$x.on.bound <- FALSE

    return(bs)
  } # end of Step2Longer.side


  trim <- function(x)
  {
    return(x[-c(1, length(x))])
  } # end of trim


  trim.matrix <- function(m)
  {
    return(m[,-c(1, ncol(m)), drop=FALSE])
  } # end of trim.matrix


  Updated.bound.info <- function(bs, rightside.bound)
  {
    # rightside.bound: a boolean value

    side <- ifelse(rightside.bound, 2, 1)
    bs$found.bounds[side] <- TRUE

    opp.side <- 3 - side
    if (bs$found.bounds[opp.side])  bs$found.Bounds <- TRUE

    return(bs)
  } # end of Updated.bound.info


  # ----------------------------------------------------------------------------------------------------
  bs$n.steps <- bs$n.steps + 1


  if (length(bs$visited) == 0)
  {
    bs$last.score.was.the.best <- TRUE

    bs$best.parms <- this.parms
    bs <- Register.initial.visit(bs, score)


    if ((dir > 0 && bs$found.bounds[2]) || (dir < 0 && bs$found.bounds[1]))
    {
      bs$continue <- FALSE
      return(bs)
    }


    if (dir != 0)  bs <- Updated.bound.info(bs, dir < 0)


    if (!is.null(bs$suggested.2nd.step) && bs$suggested.2nd.step != bs$x)
    {
      # We consider the suggested.2nd.step value
      # (and validate it to make sure it is within feasible & potentially optimal range)

      next.x.defined <- FALSE
      current.x <- bs$x
      bs$x.on.bound <- FALSE


      if (dir >= 0 && bs$suggested.2nd.step > bs$x)
      {
        next.x.defined <- TRUE

        if (bs$suggested.2nd.step >= bs$range[2])
        {
          if (bs$untouchable.bounds[2])  bs$x <- (bs$x + bs$range[2]) / 2
          else
          {
            bs$x <- bs$range[2]
            bs$x.on.bound <- TRUE
          }
        }
        else bs$x <- bs$suggested.2nd.step
      }
      else if (dir <= 0 && bs$suggested.2nd.step < bs$x)
      {
        # dir < 0 & suggested.2nd.step goes in the expected direction

        next.x.defined <- TRUE

        if (bs$suggested.2nd.step <= bs$range[1])
        {
          if (bs$untouchable.bounds[1])  bs$x <- (bs$x + bs$range[1]) / 2
          else
          {
            bs$x <- bs$range[1]
            bs$x.on.bound <- TRUE
          }
        }
        else bs$x <- bs$suggested.2nd.step
      }


      if (next.x.defined)
      {
        if (bs$change == 0)  bs$change <- abs(bs$x - current.x)
        return(bs)
      }
    }


    bs <- playSafe.change(bs)
    if (dir == 0)  dir <- ifelse((bs$range[2] - bs$x) <= (bs$x - bs$range[1]) & bs$x != bs$range[2], 1, -1)

    bs <- Stay.within.range(bs, dir > 0)
    return(bs)
  }


  # length(bs$visited) > 0


  new.score.is.the.best <- bs$better.score(score, bs$scores[,bs$which.best])
  if (new.score.is.the.best) bs$best.parms <- this.parms
    bs$last.score.was.the.best <- new.score.is.the.best # useful when debuggin (your own code calling this fct)



  if (bs$n.steps == 2)
  {
    going.east <- bs$x > bs$visited

    if (new.score.is.the.best)  bs <- Updated.bound.info(bs, !going.east)
    else                        bs <- Updated.bound.info(bs,  going.east)
  }



  if (new.score.is.the.best && bs$x.on.bound && dir != 0)
  {
    right.side <- bs$x == bs$range[2]

    if (length(bs$visited) > 1)  bs <- Clean.visits(bs, !right.side)
    bs <- Register.new.visit(bs, score, right.side)

    is.optimal  <- (dir < 0 && !right.side) || (dir > 0 && right.side)

    bs$continue <- !is.optimal
    if (!bs$continue)  return(bs)


    # bs$continue = TRUE (for the moment)

    bs$continue <- diff(bs$visited) >= bs$epsilon

    if (bs$continue)  bs <- Stay.within2visits(bs)
    return(bs)
  }



  if (length(bs$visited) == 1)
  {
    right.side <- bs$x > bs$visited
    last.dir   <- ifelse(right.side, 1, -1)
    side       <- ifelse(right.side, 2, 1)
    opp.side   <- 3 - side


    if (bs$x.on.bound)
    {
      bs <- Register.new.visit(bs, score, right.side, new.score.is.the.best)

      if (bs$found.Bounds)
      {
        bs <- Stay.within2visits(bs)
      }
      else
      {
        bs$x <- bs$visited[opp.side]
        bs <- Stay.within.range(bs, !right.side)
      }

      return(bs)
    }


    # !bs$x.on.bound

    if (new.score.is.the.best)
    {
      if (dir == last.dir)
      {
        bs <- Replace.visit(bs, score, right.side)
        bs <- Stay.within.range(bs, right.side)

        if (bs$untouchable.bounds[side])  bs$continue <- abs(bs$x - bs$range[side]) >= bs$epsilon
        return(bs)
      }


      bs <- Register.new.visit(bs, score, right.side)

      if (dir != 0)
      {
        bs <- Stay.within2visits(bs)
        bs$found.bounds <- rep(TRUE, 2)
        bs$found.Bounds <-     TRUE
      }
      else  bs <- Stay.within.range(bs, right.side)
    }
    else
    {
      # !new.score.is.the.best

      bs <- Register.new.visit(bs, score, right.side, FALSE)

      if (!bs$found.Bounds)  bs <- Updated.bound.info(bs, right.side)


      if (bs$found.Bounds)
      {
        bs <- Stay.within2visits(bs)
      }
      else
      {
        bs$x <- bs$visited[opp.side]
        bs <- Stay.within.range(bs, !right.side)
      }
    }

    return(bs)
  }


  if (length(bs$visited) == 3)
  {
    right.side <- bs$x > bs$visited[2]


    if (new.score.is.the.best)
    {
      dir2center <- ifelse(right.side, -1, 1)

      if (dir == dir2center)
      {
        bs$visited <- trim(bs$visited)
        bs$scores  <- trim.matrix(bs$scores)

        bs <- Register.new.visit(bs, score, right.side)
        bs <- Stay.within2visits(bs)
      }
      else if (dir != 0)
      {
        bs <- Clean.visits(bs, !right.side, 2)
        bs <- Register.new.visit(bs, score, !right.side)
        bs <- Stay.within2visits(bs)
      }
      else
      {
        bs <- Clean.visits(bs, !right.side)
        bs <- Insert.new.best(bs, score)
        bs <- Step2Longer.side(bs)
      }
    }
    else
    {
      bs <- Replace.visit(bs, score, right.side)
      bs <- Step2Longer.side(bs)
    }


    bs$continue <- (bs$visited[length(bs$visited)] - bs$visited[1]) >= bs$epsilon
    return(bs)
  }


  # length(bs$visited) == 2

  within <- bs$x > bs$visited[1] && bs$x < bs$visited[2]


  if (within)
  {
    if (new.score.is.the.best)
    {
      if (dir == 0)
      {
        bs <- Insert.new.best(bs, score)
        bs <- Step2Longer.side(bs)
      }
      else
      {
        bs <- Replace.visit(bs, score, dir < 0)
        bs <- Stay.within2visits(bs)
      }
    }
    else
    {
      bs <- Replace.visit(bs, score, bs$which.best == 1)
      bs <- Stay.within2visits(bs)
    }
  }
  else
  {
    # !within

    right.side <- bs$x > bs$visited[1]
    last.dir <- ifelse(right.side, 1, -1)

    if (new.score.is.the.best)
    {
      if (dir == last.dir)
      {
        bs <- Register.initial.visit(bs, score)
        bs <- Stay.within.range(bs, dir > 0)
      }
      else
      {
        bs <- Clean.visits(bs, !right.side)

        bs <- Register.new.visit(bs, score, right.side)


        if (dir != 0)
        {
          bs <- Updated.bound.info(bs, right.side)
          bs <- Stay.within2visits(bs)
        }
        else
        {
          if (bs$found.Bounds)  bs <- Stay.within2visits(bs)
          else                  bs <- Stay.within.range(bs, last.dir > 0)
        }
      }
    }
    else
    {
      # !new.score.is.the.best

      bs <- Updated.bound.info(bs, right.side)
      bs <- Register.new.visit(bs, score, right.side, FALSE)
      bs <- Step2Longer.side(bs)
    }
  }

  return(bs)
} # end of bs.update
