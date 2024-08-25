############# WEBEXPO STAN LIBRARY #######################

##### SCRIPT FOR CREATING THE STAN MODEL OBJECTS 

# the compiled objects need to be saved in order as RDS objects to be used in the analysis

# Requires the RSTAN library to be active

#### SEG ANALYSIS #####



##### Uninformative models ---------------------------


f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_uninformative.stan'
code <- readLines(f)
stan.model.uninformative <- stan_model(model_code=code)


f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_uninformative_logNormal_mesd.stan'
code <- readLines(f)
stan.model.uninformative.logNormal.mesd <- stan_model(model_code=code)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_uninformative_Normal_mesd.stan'
code <- readLines(f)
stan.model.uninformative.Normal.mesd <- stan_model(model_code=code)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_uninformative_logNormal_mecv.stan'
code <- readLines(f)
stan.model.uninformative.logNormal.mecv <- stan_model(model_code=code)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_uninformative_Normal_mecv.stan'
code <- readLines(f)
stan.model.uninformative.Normal.mecv <- stan_model(model_code=code)




##### InformedVar models ---------------------------




f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedVar.stan'
code <- readLines(f)
stan.model.informedVar <- stan_model(model_code=code)


f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedVar_logNormal_mesd.stan'
code <- readLines(f)
stan.model.informedVar.logNormal.mesd <- stan_model(model_code=code)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedVar_Normal_mesd.stan'
code <- readLines(f)
stan.model.informedVar.Normal.mesd <- stan_model(model_code=code)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedVar_logNormal_mecv.stan'
code <- readLines(f)
stan.model.informedVar.logNormal.mecv <- stan_model(model_code=code)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedVar_Normal_mecv.stan'
code <- readLines(f)
stan.model.informedVar.Normal.mecv <- stan_model(model_code=code)


# Versions with use of past data

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedVar_pastData.stan'
code <- readLines(f)
stan.model.informedVar.pastData <- stan_model(model_code=code)


f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedVar_logNormal_mesd_pastData.stan'
code <- readLines(f)
stan.model.informedVar.logNormal.mesd.pastData <- stan_model(model_code=code)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedVar_Normal_mesd_pastData.stan'
code <- readLines(f)
stan.model.informedVar.Normal.mesd.pastData <- stan_model(model_code=code)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedVar_logNormal_mecv_pastData.stan'
code <- readLines(f)
stan.model.informedVar.logNormal.mecv.pastData <- stan_model(model_code=code)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedVar_Normal_mecv_pastData.stan'
code <- readLines(f)
stan.model.informedVar.Normal.mecv.pastData <- stan_model(model_code=code)

##### InformedVar+Mean models ---------------------------


f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedMean.stan'
code <- readLines(f)
stan.model.InformedMean <- stan_model(model_code=code)


f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedMean_logNormal_mesd.stan'
code <- readLines(f)
stan.model.InformedMean.logNormal.mesd <- stan_model(model_code=code)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedMean_Normal_mesd.stan'
code <- readLines(f)
stan.model.InformedMean.Normal.mesd <- stan_model(model_code=code)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedMean_logNormal_mecv.stan'
code <- readLines(f)
stan.model.InformedMean.logNormal.mecv <- stan_model(model_code=code)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedMean_Normal_mecv.stan'
code <- readLines(f)
stan.model.InformedMean.Normal.mecv <- stan_model(model_code=code)


# Versions with use of past data

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedMean_pastData.stan'
code <- readLines(f)
stan.model.InformedMean.pastData <- stan_model(model_code=code)


f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedMean_logNormal_mesd_pastData.stan'
code <- readLines(f)
stan.model.InformedMean.logNormal.mesd.pastData <- stan_model(model_code=code)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedMean_Normal_mesd_pastData.stan'
code <- readLines(f)
stan.model.InformedMean.Normal.mesd.pastData <- stan_model(model_code=code)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedMean_logNormal_mecv_pastData.stan'
code <- readLines(f)
stan.model.InformedMean.logNormal.mecv.pastData <- stan_model(model_code=code)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/SEG_informedMean_Normal_mecv_pastData.stan'
code <- readLines(f)
stan.model.InformedMean.Normal.mecv.pastData <- stan_model(model_code=code)

