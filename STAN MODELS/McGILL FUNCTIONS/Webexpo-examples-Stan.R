# Exemples d'appels de la version Stan de quelques algorithmes de WebExpo

# steps to running an example 


# 1) source the stan functions from github
# 2) compile the selected models and create the model list as describe in compile-stan-models.R
# 4) run the stan function with the data



# --- sourcing the stan functions -------------------------------------------------------------------------

library(rstan)

source('https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/McGILL%20FUNCTIONS/stan-fcts.R')
source('https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/McGILL%20FUNCTIONS/SEG-informedVar-stan.R')
source('https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/McGILL%20FUNCTIONS/SEG-uninformative-stan.R')
source('https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/McGILL%20FUNCTIONS/SEG-informedMean-stan.R')

# --- compiling the models -------------------------------------------------------------------------

stan.models.list <- list()

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_uninformative.stan'

stan.models.list <- augment.stan.models.list(stan.models.list, stan.file = f)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedvar.stan'

stan.models.list <- augment.stan.models.list(stan.models.list, f)

compiled.models.list(stan.models.list)


# --- Uninformative / normal / no ME -----------------------------------------------------------------------

# compiling and saving the model


y <- c(6, 7.2, 5.4, 6.08)
o <- SEG.uninformative.stan(y, outcome.is.logNormally.distributed=FALSE,
                            models.list=stan.models.list)

plot(o$mu, o$sigma, type='p', pch='.')


# --- InformedVar -------------------------------------------------------------------------


y <- c(6, 7.2, 5.4, 6.08)
o <- SEG.informedvar.stan(y, outcome.is.logNormally.distributed=FALSE,
                          models.list=stan.models.list)

plot(o$mu, o$sigma, type='p', pch='.')


