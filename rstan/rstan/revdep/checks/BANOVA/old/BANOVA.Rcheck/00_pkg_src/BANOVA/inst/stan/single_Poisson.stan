data {
  int<lower=0> N;
  int<lower=0> J;
  matrix[N, J] X;
  int y[N];
}

parameters {
  vector[J] beta1;
} 

transformed parameters {
  vector[N] y_hat;
  y_hat = X*beta1;
}

model {
  y ~ poisson_log(y_hat);
  for (i in 1:J){
    beta1[i] ~ normal(0, 1);
  }
}

generated quantities {
  real var_f;
  real r_2;
  var_f = variance(y_hat);
  r_2 = var_f/(var_f + log(1/exp(beta1[1]) + 1));
}