data {
  int<lower=0> N;
  int<lower=0> J;
  int<lower=3> n_cat;
  int<lower=0> M;
  int<lower=0> K;
  vector[J] X[N, n_cat];
  matrix[M, K] Z;
  int<lower=0> id[N];
  int<lower=1, upper=n_cat> y[N];
}

parameters {
  matrix[J, M] beta1;
  matrix[K, J] beta2; 
  vector<lower=0>[J] tau_beta1Sq;
}
transformed parameters {
  matrix[N, n_cat] mu;
  for (i in 1:N){
    for (j in 1:n_cat){
      mu[i, j] = X[i,j]'*(beta1[,id[i]]);
    }
  }
}
model {
  matrix[M, J] mu_beta1;
  vector[J] tau_beta1;
  tau_beta1 = sqrt(tau_beta1Sq);
  for (i in 1:N){
      y[i] ~ categorical_logit(mu[i, ]');
  }
  mu_beta1 = Z*beta2;
  for (i in 1:J){
    beta1[i,] ~ normal(mu_beta1[,i], tau_beta1[i]);
  }
  tau_beta1Sq ~ inv_gamma(1, 1);
  for (i in 1:J){
    beta2[,i] ~ normal(0, 100);
  }
}

