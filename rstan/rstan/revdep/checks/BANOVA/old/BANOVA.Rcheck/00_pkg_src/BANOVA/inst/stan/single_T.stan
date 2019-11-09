data {
  int<lower=0> N;
  int<lower=0> J;
  matrix[N, J] X;
  vector[N] y;
}

parameters {
  vector[J] beta1;
  real<lower=0> tau_yS;
  real<lower=1, upper=10> df;
} 

transformed parameters {
  vector[N] y_hat;
  y_hat = X*beta1;
}

model {
  real tau_y;
  tau_y = sqrt(tau_yS);
  y ~ student_t(df, y_hat, tau_y);
  df ~ normal(5,5);
  tau_yS ~ inv_gamma(1, 1);
  for (i in 1:J){
    beta1[i] ~ normal(0, 1);
  }
}

generated quantities {
  real var_f;
  real r_2;
  real tau_ySq;
  var_f = variance(y_hat);
  r_2 = var_f/variance(y);
  tau_ySq = variance(y - y_hat);
}

