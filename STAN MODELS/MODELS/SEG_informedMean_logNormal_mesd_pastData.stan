// Version 0.1 (Aug 2024)
//         Shared/distributed: yes
// -------------------------------


data
{
  int N;
  vector[N] y;    // Data points known exactly

  int L;
  vector[L] lt;   // Right-censored Data points  (x < some value)

  int G;
  vector[G] gt;   // Left-censored Data points  (x > some value)

  int I;
  matrix[I, 2] interval;  // Interval-censored Data points  (some value < x < some value)


  int pastData_n;      // Past data summary stats
  real pastData_mean;
  real pastData_sd;


  // Hyperparameters

  real MU_MEAN;
  real MU_SD;

  real LOGSIGMA_MEAN;
  real LOGSIGMA_SD;

  real ME_SD_LO;
  real ME_SD_HI;
}

transformed data
{
  real pastData_a;
  real pastData_sqrtN;
  real pastData_s2;

  pastData_a = (pastData_n - 1) / 2;
  pastData_sqrtN = sqrt(pastData_n);
  pastData_s2 = pow(pastData_sd, 2);
}

parameters
{
  real mu;              // Global mean
  real<lower=0> sigma;  // Std Deviation
  real me_sd;


  // (latent) True values

  vector[N] true_y;
  vector[L] true_lt;
  vector[G] true_gt;
  vector[I] true_int;
}

transformed parameters
{
  real pastData_b;
  real pastData_sigma = sigma / pastData_sqrtN;

  pastData_b = pastData_a / pow(sigma, 2);
}

model
{
  vector[G] loglik_gt;
  vector[L] loglik_lt;
  vector[I] loglik_interval;


  // Likelihood

  for (i in 1:N)
  {
    true_y[i] ~ lognormal(mu, sigma);

    target += normal_lpdf(y[i] | true_y[i], me_sd);
  }



  // Right-censored data values

  for (i in 1:L)
  {
    true_lt[i] ~ lognormal(mu, sigma);
    loglik_lt[i] = normal_lcdf(lt[i] | true_lt[i], me_sd);
  }


  // Left-censored data values

  for (i in 1:G)
  {
    true_gt[i] ~ lognormal(mu, sigma);
    loglik_gt[i] = normal_lccdf(gt[i] | true_gt[i], me_sd);
  }


  // Interval-censored data values

  for (i in 1:I)
  {
    true_int[i] ~ lognormal(mu, sigma);
    loglik_interval[i] = log(exp(normal_lcdf(interval[i,2] | true_int[i], me_sd)) - exp(normal_lcdf(interval[i,1] | true_int[i], me_sd)));
  }


  target += sum(loglik_lt) + sum(loglik_gt) + sum(loglik_interval);


  // Past data likelihood

  pastData_mean ~ normal(mu, pastData_sigma);
  pastData_s2   ~ gamma(pastData_a, pastData_b);


  // Prior distrns

  mu    ~ normal(MU_MEAN, MU_SD);
  sigma ~ lognormal(LOGSIGMA_MEAN, LOGSIGMA_SD);

  me_sd ~ uniform(ME_SD_LO, ME_SD_HI);
}
// label: seg_informedmean_lognormal_mesd_pastdata
