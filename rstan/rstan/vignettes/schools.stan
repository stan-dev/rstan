data {
  int<lower=0> J; // number of schools 
  real y[J];      // estimated treatment effects
  real<lower=0> sigma[J]; // s.e. of effect estimates 
}
parameters {
  real mu; 
  real<lower=0> tau;
  vector[J] eta;
}
transformed parameters {
  vector[J] theta;
  theta <- mu + tau * eta;
}
model {
  eta ~ normal(0, 1);
  y ~ normal(theta, sigma);
}
