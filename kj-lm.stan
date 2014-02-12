// Linear regression in Stan
data {
  int<lower=1> N; real y[N]; }
parameters {
  real baseline; real trend; real<lower=0> sigma; }
model {
  // Priors should be added here, for sigma, baseline and trend.
  for (t in 1:N) y[t] ~ normal(baseline + trend*t, sigma); }
