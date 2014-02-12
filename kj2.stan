// Time series model with pooled, possiby decreasing MA coefficients, 
// t residuals and quadratic baseline. 

data {
  int<lower=1> N; int<lower=1> D; real y[N]; 
}
parameters {
  real baseline; real trend; real k2;
  real<lower=0> sigma; 
  real theta_lsd1; real theta_lsd2;
  real<lower=0> df;
  vector[D] theta_nrm; 
}
transformed parameters {
  vector[N] eta; 
  vector[D] theta;
  { vector[N] err;
    for (i in 1:D) theta[i] <- exp((theta_lsd1*(D-i) + theta_lsd2*i)/D)*theta_nrm[i];
    for (t in 1:N) { 
    eta[t] <- k2*(t-N/2)*(t-N/2) + trend*t + baseline; 
    for (lag in 1:D) { if (t>lag)  eta[t] <- eta[t] + tanh(theta[lag]) * err[t-lag]; }
    err[t] <- y[t] - eta[t]; }}
}
model { 
  // Add priors for k2, baseline, trend, and sigma.
  df ~ normal(0, 20);
  theta_lsd1 ~ normal(-2, 2); theta_lsd2 ~ normal(-2, 2);
  theta_nrm ~ normal(0, 1);
  y ~ student_t(df+3, eta, sigma); 
}

