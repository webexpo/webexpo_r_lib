############# WEBEXPO STAN LIBRARY #######################

##### SCRIPT FOR CREATING THE STAN MODEL OBJECTS 

# the compiled objects need to be saved in order as RDS objects to be used in the analysis

# Requires the RSTAN library to be active

#### SEG ANALYSIS #####



##### Uninformative models ---------------------------


f <- 'SEG-uniformative.stan'
code <- readLines(f)
stan.model.Uninformative <- stan_model(model_code=code)
#saveRDS(stan.model.Uninformative, file='stan_model_Uninformative.RDS')


f <- 'SEG-uniformative-logNormal+me=classic.stan'
code <- readLines(f)
stan.model.Uninformative.logNormal.mesd <- stan_model(model_code=code)
#saveRDS(stan.model.Uninformative.logNormal.mesd, file='stan_model_Uninformative_logNormal_mesd.RDS')

f <- 'SEG-uniformative-Normal+me=classic.stan'
code <- readLines(f)
stan.model.Uninformative.Normal.mesd <- stan_model(model_code=code)
#saveRDS(stan.model.Uninformative.Normal.mesd, file='stan_model_Uninformative_Normal_mesd.RDS')

f <- 'SEG-uniformative-logNormal+me=CV.stan'
code <- readLines(f)
stan.model.Uninformative.logNormal.mecv <- stan_model(model_code=code)
#saveRDS(stan.model.Uninformative.logNormal.mecv, file='stan_model_Uninformative_logNormal_mecv.RDS')

f <- 'SEG-uniformative-Normal+me=CV.stan'
code <- readLines(f)
stan.model.Uninformative.Normal.mecv <- stan_model(model_code=code)
#saveRDS(stan.model.Uninformative.Normal.mecv, file='stan_model_Uninformative_Normal_mecv.RDS')


##### InformedVar models ---------------------------


f <- 'SEG-informedVar.stan'
code <- readLines(f)
stan.model.InformedVar <- stan_model(model_code=code)
#saveRDS(stan.model.InformedVar, file='stan_model_InformedVar.RDS')


f <- 'SEG-informedVar-logNormal+me=classic.stan'
code <- readLines(f)
stan.model.InformedVar.logNormal.mesd <- stan_model(model_code=code)
#saveRDS(stan.model.InformedVar.logNormal.mesd, file='stan_model_InformedVar_logNormal_mesd.RDS')

f <- 'SEG-informedVar-Normal+me=classic.stan'
code <- readLines(f)
stan.model.InformedVar.Normal.mesd <- stan_model(model_code=code)
#saveRDS(stan.model.InformedVar.Normal.mesd, file='stan_model_InformedVar_Normal_mesd.RDS')

f <- 'SEG-informedVar-logNormal+me=CV.stan'
code <- readLines(f)
stan.model.InformedVar.logNormal.mecv <- stan_model(model_code=code)
#saveRDS(stan.model.InformedVar.logNormal.mecv, file='stan_model_InformedVar_logNormal_mecv.RDS')

f <- 'SEG-informedVar-Normal+me=CV.stan'
code <- readLines(f)
stan.model.InformedVar.Normal.mecv <- stan_model(model_code=code)
#saveRDS(stan.model.InformedVar.Normal.mecv, file='stan_model_InformedVar_Normal_mecv.RDS')


##### InformedVar+Mean models ---------------------------


f <- 'SEG-informedMean.stan'
code <- readLines(f)
stan.model.InformedMean <- stan_model(model_code=code)
#saveRDS(stan.model.InformedMean, file='stan_model_InformedMean.RDS')


f <- 'SEG-informedMean-logNormal+me=classic.stan'
code <- readLines(f)
stan.model.InformedMean.logNormal.mesd <- stan_model(model_code=code)
#saveRDS(stan.model.InformedMean.logNormal.mesd, file='stan_model_InformedMean_logNormal_mesd.RDS')

f <- 'SEG-informedMean-Normal+me=classic.stan'
code <- readLines(f)
stan.model.InformedMean.Normal.mesd <- stan_model(model_code=code)
#saveRDS(stan.model.InformedMean.Normal.mesd, file='stan_model_InformedMean_Normal_mesd.RDS')

f <- 'SEG-informedMean-logNormal+me=CV.stan'
code <- readLines(f)
stan.model.InformedMean.logNormal.mecv <- stan_model(model_code=code)
#saveRDS(stan.model.InformedMean.logNormal.mecv, file='stan_model_InformedMean_logNormal_mecv.RDS')

f <- 'SEG-informedMean-Normal+me=CV.stan'
code <- readLines(f)
stan.model.InformedMean.Normal.mecv <- stan_model(model_code=code)
#saveRDS(stan.model.InformedMean.Normal.mecv, file='stan_model_InformedMean_Normal_mecv.RDS')
