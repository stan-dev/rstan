data {
  int<lower=0> J;          // number of schools
  array[N] real y;               // estimated treatment effects
  array[J] real<lower=0> sigma;  // s.e. of effect estimates
}
parameters {
  real mu;
  real<lower=0> tau;
  vector[J] eta;
}
transformed parameters {
  vector[J] theta;
  theta = mu + tau * eta;
}
model {
  target += normal_lpdf(eta | 0, 1);
  target += normal_lpdf(y | theta, sigma);
}
