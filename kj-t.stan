// Linear regression with t residuals.
data {
  int<lower=1> N; real y[N]; }
parameters {
  real baseline; real trend; real<lower=0> sigma; real<lower=0> df; }
model {
  // Add priors for baseline, trend, and sigma here.
  df ~ normal(0, 20);
  for (t in 1:N) y[t] ~ student_t(df+1, baseline + trend*t, sigma); }
