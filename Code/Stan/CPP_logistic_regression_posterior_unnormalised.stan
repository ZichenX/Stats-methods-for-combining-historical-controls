data {
  int<lower=0> N0;
  int<lower=0> P;
  matrix[N0, P] X0;
  int<lower=0, upper=1> y0[N0];
  int<lower=0> N;
  matrix[N, P] X;
  int<lower=0, upper=1> y[N];
  real<lower=0> eta;
  real<lower=0> nu;
}
parameters {
  real alpha0;
  real alpha1;
  vector[P] beta0;
  vector[P] beta1;
  real<lower=0,upper=1> a_0;
}
model {
  /* Commensurate Power prior */
  target += normal_lpdf(alpha0 | 0, 1);
  target += normal_lpdf(beta0 | 0, 1);
  target += normal_lpdf(alpha1 | 0, 1);
  target += normal_lpdf(beta1 | 0, 1);
  target += a_0 * bernoulli_logit_lpmf(y0 | alpha0 + X0*beta0);
  target += beta_lpdf(a_0 | eta, nu);
  /* Likelihood */
  target += bernoulli_logit_lpmf(y | alpha1 + X*beta1);
}
