// MA time series model with linear baseline.
// May not work well for larger D because the trend 
// may not be identifiable.
data {
  int<lower=1> N; real y[N]; 
  int<lower=1> D; // Maximum lag. 
}
parameters {
  real baseline; real trend; 
  real<lower=0> sigma;
  vector<lower=-1, upper=1>[D] theta; 
}
transformed parameters {
  vector[N] eta; 
  { vector[N] err;
    for (t in 1:N) { 
    eta[t] <- trend*t + baseline; 
    for (lag in 1:D) { if (t>lag)  eta[t] <- eta[t] + theta[lag] * err[t-lag]; }
    err[t] <- y[t] - eta[t]; }}
}
// The model is missing priors for baseline, trend, and sigma.
model { y ~ normal(eta, sigma); }
