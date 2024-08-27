############# WEBEXPO STAN LIBRARY #######################

##### SCRIPT FOR CREATING THE STAN MODEL OBJECTS 

# The scripts below permit to read the model code from the webexpo github repository and create stan models object

# Once the desired models objects are created, they should be assigned to a list using the augment.stan.models.list function, as shown below

# this list is an input for the calculation functions

# Requires the RSTAN library to be active

#### seg ANALYSIS #####

stan.models.list <- list()


#### names the models in the created list : 

compiled.models.list(stan.models.list)

##### Uninformative models ---------------------------


f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_uninformative.stan'

stan.models.list <- augment.stan.models.list(stan.models.list, stan.file = f)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_uninformative_lognormal_mesd.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_uninformative_normal_mesd.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_uninformative_lognormal_mecv.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_uninformative_normal_mecv.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)




##### Informedvar models ---------------------------




f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedvar.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)


f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedvar_lognormal_mesd.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedvar_normal_mesd.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedvar_lognormal_mecv.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedvar_normal_mecv.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)


# Versions with use of past data

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedvar_pastdata.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)


f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedvar_lognormal_mesd_pastdata.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedvar_normal_mesd_pastdata.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedvar_lognormal_mecv_pastdata.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedvar_normal_mecv_pastdata.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)

##### Informedvar+mean models ---------------------------


f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedmean.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)


f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedmean_lognormal_mesd.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedmean_normal_mesd.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedmean_lognormal_mecv.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedmean_normal_mecv.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)


# Versions with use of past data

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedmean_pastdata.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)


f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedmean_lognormal_mesd_pastdata.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedmean_normal_mesd_pastdata.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedmean_lognormal_mecv_pastdata.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)

f <- 'https://raw.githubusercontent.com/webexpo/webexpo_r_lib/master/STAN%20MODELS/MODELS/seg_informedmean_normal_mecv_pastdata.stan'
stan.models.list <- augment.stan.models.list(stan.models.list, f)

