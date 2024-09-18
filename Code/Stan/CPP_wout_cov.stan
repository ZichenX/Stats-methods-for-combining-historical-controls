data {
  int<lower=0> N0;
  int<lower=0> P;
  matrix[N0, P] X0;
  int<lower=0, upper=1> y0[N0];
  int<lower=0> N1;
  matrix[N1, P] X1;
  int<lower=0, upper=1> y1[N1];
}
parameters {
  real alpha0;
  real alpha1;
  vector[P] beta0;
  vector[P] beta1;
  real<lower=0, upper=1> a_0;
  real<lower=0, upper=1> tau;
}
model {
  /* Power prior */
  target += normal_lpdf(alpha0 | 0, 100);
  target += normal_lpdf(beta0 | 0, 100);
  target += normal_lpdf(alpha1 | alpha0, tau);
  target += normal_lpdf(beta1 | beta0, tau); 
  // here we assume that the heterogeneity par. is the same.
  target += a_0 * bernoulli_logit_lpmf(y0 | alpha0 + X0*beta0);
  target += beta_lpdf(a_0 | 1/square(tau), 1); 
  target += beta_lpdf(tau | 1, 1);
  /* Likelihood */
  target += bernoulli_logit_lpmf(y1 | alpha1 + X1*beta1);
}
