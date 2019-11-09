data {
  int<lower=0> N;
  int<lower=0> J;
  int<lower=1> trials[N];
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
  y ~ binomial_logit(trials, y_hat);
  for (i in 1:J){
    beta1[i] ~ normal(0, 1);
  }
}

generated quantities {
  real var_f;
  real r_2;
  var_f = variance(y_hat);
  r_2 = var_f/(var_f + pi()^2/3);
}
