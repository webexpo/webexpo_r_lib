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

  real MU_MEAN;
  real MU_SD;

  real LOGSIGMA_MEAN;
  real LOGSIGMA_SD;
}

parameters
{
  real mu;              // Global mean
  real<lower=0> sigma;  // Std Deviation
}

model
{
  vector[G] loglik_gt;
  vector[L] loglik_lt;
  vector[I] loglik_interval;


  // Likelihood

  y ~ normal(mu, sigma);

    // Right-censored data values

    for (i in 1:L)
    {
      loglik_lt[i] = normal_lcdf(lt[i] | mu, sigma);
    }

    // Left-censored data values

    for (i in 1:G)
    {
      loglik_gt[i] = normal_lccdf(gt[i] | mu, sigma);
    }

    // Interval-censored data values

    for (i in 1:I)
    {
      loglik_interval[i] = log(exp(normal_lcdf(interval[i,2] | mu, sigma)) - exp(normal_lcdf(interval[i,1] | mu, sigma)));
    }


  target += sum(loglik_lt) + sum(loglik_gt) + sum(loglik_interval);


  // Prior distrns

  mu    ~ normal(MU_MEAN, MU_SD);
  sigma ~ lognormal(LOGSIGMA_MEAN, LOGSIGMA_SD);
}
