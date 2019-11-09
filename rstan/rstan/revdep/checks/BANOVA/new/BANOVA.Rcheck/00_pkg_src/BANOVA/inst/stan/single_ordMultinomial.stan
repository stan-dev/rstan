data {
  int<lower=3> cat;
  int<lower=0> N;
  int<lower=0> J;
  matrix[N, J] X;
  int<lower=1, upper=cat> y[N];
}

parameters {
  vector[J] beta1;
  ordered[cat-1] c;  
}
transformed parameters {
  ordered[cat-1] c_trans;
  vector[N] y_hat;
  y_hat = X*beta1;
  for (i in 1:(cat-1))
    c_trans[i] = c[i] - c[1];
}

model {
  for (i in 1:N){
      y[i] ~ ordered_logistic(y_hat[i], c_trans);
  }
  c ~ normal(0, 10);
  for (i in 1:J){
    beta1[i] ~ normal(0, 10);
  }
}

generated quantities {
  real var_f;
  real r_2;
  var_f = variance(y_hat);
  r_2 = var_f/(var_f + pi()^2/6);
}