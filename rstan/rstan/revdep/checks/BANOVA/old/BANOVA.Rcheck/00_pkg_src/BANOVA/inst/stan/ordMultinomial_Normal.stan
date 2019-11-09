data {
  int<lower=3> cat;
  int<lower=0> N;
  int<lower=0> J;
  int<lower=0> M;
  int<lower=0> K;
  matrix[N, J] X;
  matrix[M, K] Z;
  int<lower=0> id[N];
  int<lower=1, upper=cat> y[N];
}

parameters {
  ordered[cat-1] c;
  matrix[J, M] beta1;
  matrix[K, J] beta2; 
  vector<lower=0>[J] tau_beta1Sq;
}
transformed parameters {
  ordered[cat-1] c_trans;
  vector[N] y_pred;
  vector[N] y_hat;
  matrix[M, J] mu_beta1;
  mu_beta1 = Z*beta2;
  for (i in 1:(cat-1))
    c_trans[i] = c[i] - c[1]; 
  for (i in 1:N){
    y_hat[i] = X[i,]*beta1[,id[i]];
  }
  for (i in 1:N){
    y_pred[i] = X[i,]*mu_beta1[id[i],]';
  }
}

model {
  vector[J] tau_beta1;
  tau_beta1 = sqrt(tau_beta1Sq);
  for (i in 1:N){
      y[i] ~ ordered_logistic(y_hat[i], c_trans);
  }
  c ~ normal(0, 10);
  for (i in 1:J){
    beta1[i,] ~ normal(mu_beta1[,i], tau_beta1[i]);
  }
  tau_beta1Sq ~ inv_gamma(1, 1);
  for (i in 1:J){
    beta2[,i] ~ normal(0, 100);
  }
}

generated quantities {
  real r_2;
  r_2 = variance(y_pred)/(variance(y_hat) + pi()^2/6);
}