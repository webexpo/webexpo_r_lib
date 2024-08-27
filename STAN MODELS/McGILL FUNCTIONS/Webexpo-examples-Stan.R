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


# --- Uninformative / normal / no ME -----------------------------------------------------------------------

# compiling and saving the model

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_uninformative.stan'
code <- readLines(f)
stan.model.uninformative <- stan_model(model_code=code)
saveRDS(stan.model.uninformative,"SEG_stan_model_uninformative.RDS")

y <- c(6, 7.2, 5.4, 6.08)
o <- SEG.uninformative.stan(y, outcome.is.logNormally.distributed=FALSE,
                            models.folder=stan.folder)

plot(o$mu, o$sigma, type='p', pch='.')


# --- InformedVar -------------------------------------------------------------------------

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedVar.stan'
code <- readLines(f)
stan.model.informedVar <- stan_model(model_code=code)
saveRDS(stan.model.informedVar,"SEG_stan_model_informedVar.RDS")



y <- c(6, 7.2, 5.4, 6.08)
o <- SEG.informedvar.stan(y, outcome.is.logNormally.distributed=FALSE,
                          models.folder=stan.folder)

plot(o$mu, o$sigma, type='p', pch='.')


