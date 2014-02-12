// Time series model with pooled MA coefficients and 
// linear baseline.

data {
  int<lower=1> N; real y[N]; 
  int<lower=1> D; // Maximum lag. 
}
parameters {
  real baseline; real trend; 
  real<lower=0> sigma; real<lower=0> theta_sd;
  vector[D] theta; 
}
transformed parameters {
  vector[N] eta; 
  { vector[N] err;
    for (t in 1:N) { 
    eta[t] <- trend*t + baseline; 
    for (lag in 1:D) { if (t>lag)  eta[t] <- eta[t] + tanh(theta[lag]) * err[t-lag]; }
    err[t] <- y[t] - eta[t]; }}
}
model { 
  // Priors should be added here for sigma, baseline and trend.
  theta_sd ~ normal(0, .5); theta ~ normal(0, theta_sd);
  y ~ normal(eta, sigma); 
}
