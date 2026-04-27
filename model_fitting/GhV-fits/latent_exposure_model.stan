data {
  int<lower=1> N;                            // total number of observations
  int<lower=1> N_bats;                       // number of bats
  int<lower=1,upper=N_bats> bat_id[N];       // mapped barcode_id
  vector[N] t_obs;                            // observation times (days)
  vector[N] MFI;                              // observed antibody titers
  vector[N] t_exp_min;                        // lower bound for latent exposure
  vector[N] t_exp_max;                        // upper bound for latent exposure
}

parameters {
  // Population-level (hyper)parameters
  real<lower=0> mu_A;
  real<lower=0> sigma_A;
  real<lower=0> mu_B;
  real<lower=0> sigma_B;

  // Bat-level random effects (standardized)
  vector[N_bats] A_raw;
  vector[N_bats] B_raw;

  // Latent exposure timing per bat
  vector<lower=0,upper=1>[N_bats] t_exp_raw;

  // Observation noise
  real<lower=0> sigma;
}

transformed parameters {
  // De-standardize random effects to get bat-level parameters
  vector<lower=0>[N_bats] A;
  vector<lower=0>[N_bats] B;
  vector[N_bats] t_exp;

  for (i in 1:N_bats) {
    A[i] = fmax(1e-3, mu_A + sigma_A * A_raw[i]);  // enforce positivity
    B[i] = fmax(1e-6, mu_B + sigma_B * B_raw[i]);
    t_exp[i] = t_exp_min[i] + t_exp_raw[i] * (t_exp_max[i] - t_exp_min[i]);
  }
}

model {
  // Hyperpriors
  mu_A ~ normal(400, 200);           // broad prior, center near your empirical scale
  sigma_A ~ normal(0, 200);          // variability across bats
  mu_B ~ normal(1/42, 0.02);         // prior for decay rate
  sigma_B ~ normal(0, 0.02);
  sigma ~ normal(0, 50);

  // Standardized bat-level effects
  A_raw ~ normal(0, 1);
  B_raw ~ normal(0, 1);

  // Likelihood
  for (n in 1:N) {
    real t_since = t_obs[n] - t_exp[bat_id[n]];
    real mu = (t_since > 0)
      ? A[bat_id[n]] * t_since * exp(-B[bat_id[n]] * t_since)
      : 0;
    MFI[n] ~ normal(mu, sigma);
  }
}

generated quantities {
  // Posterior predictive draws for validation or plotting
  vector[N] MFI_pred;
  for (n in 1:N) {
    real t_since = t_obs[n] - t_exp[bat_id[n]];
    real mu = (t_since > 0)
      ? A[bat_id[n]] * t_since * exp(-B[bat_id[n]] * t_since)
      : 0;
    MFI_pred[n] = normal_rng(mu, sigma);
  }
}
