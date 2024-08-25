# Exemples d'appels de la version Stan de quelques algorithmes de WebExpo

# Attention: il y a des paths a modifier avant la premiere utilisation des fichiers suivants:
#   (cf. les lignes où on trouve un commentaire 'ICI' dans ces fichiers)
#            1) SEG*-stan.R
#            2) models/compile*.R
#            3) stan-fcts.R
#            4) dans le fichier que vous etes en train de lire, ci-dessous
#


# steps to running an example 

# 1) set the working directory to the folder where the local stan model files will be located
# 2) source the stan functions from github
# 3) compile de selected model and store it in tghe above folder
# 4) run the stan function with the data

stan.folder <- 'F:\Dropbox\temp\stanmodels'

setwd(stan.folder)

# --- preparing the compiled models -------------------------------------------------------------------------

f <- 'SEG-uniformative.stan'
code <- readLines(f)
stan.model.Uninformative <- stan_model(model_code=code)
saveRDS(stan.model.Uninformative, file='stan_model_Uninformative.RDS')


# --- Uninformative / normal / no ME -----------------------------------------------------------------------

source('SEG-uninformative-stan.R') # lit fct et pre-compile le modele

y <- c(6, 7.2, 5.4, 6.08)
o <- SEG.uninformative.stan(y, outcome.is.logNormally.distributed=FALSE)

plot(o$mu, o$sigma, type='p', pch='.')


# --- InformedVar -------------------------------------------------------------------------

source('SEG-informedVar-stan.R') # lit fct et pre-compile le modele

y <- c(6, 7.2, 5.4, 6.08)
o <- SEG.informedvar.stan(y, outcome.is.logNormally.distributed=FALSE)

plot(o$mu, o$sigma, type='p', pch='.')


obs <- c(0.123755622, 0.003125879, 0.095770970, 0.047123777, 0.040842749, 0.001553962, 0.002537681, 0.003694534)
lt <- 0.0008851732

res <- SEG.informedvar.stan(y=obs, lt=lt, mu.lower = -20, mu.upper=20)




# --- InformedVar & InformedMean ----------------------------------------------------------

source('SEG-informedMean-stan.R') # lit fct et pre-compile le modele

y <- c(6, 7.2, 5.4, 6.08)
o <- SEG.informedmean.stan(y, mu.mean=3, mu.sd=1.2, outcome.is.logNormally.distributed=FALSE)

plot(o$mu, o$sigma, type='p', pch='.')
