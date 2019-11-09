data {
  int<lower=0> N;
  int<lower=0> J;
  int<lower=3> n_cat;
  vector[J] X[N, n_cat];
  int<lower=1, upper=n_cat> y[N];
}

parameters {
  vector[J] beta1;
}
transformed parameters {
  matrix[N, n_cat] mu;
  for (i in 1:N){
    for (j in 1:n_cat){
      mu[i, j] = X[i,j]'*(beta1);
    }
  }
}
model {
  for (i in 1:N){
      y[i] ~ categorical_logit(mu[i, ]');
  }
  for (i in 1:J){
    beta1[i] ~ normal(0, 100);
  }
}

