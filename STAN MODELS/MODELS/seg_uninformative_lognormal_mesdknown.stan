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


  // Hyperparameters

  real MU_LO;
  real MU_HI;

  real SIGMA_LO;
  real SIGMA_HI;

  real me_sd;
}

parameters
{
  real mu;              // Global mean
  real<lower=0> sigma;  // Std Deviation


  // (latent) True values

  vector[N] true_y;
  vector[L] true_lt;
  vector[G] true_gt;
  vector[I] true_int;
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


  // Prior distrns

  mu    ~ uniform(MU_LO, MU_HI);
  sigma ~ uniform(SIGMA_LO, SIGMA_HI);
}
// label: seg_uninformative_lognormal_mesdknown
