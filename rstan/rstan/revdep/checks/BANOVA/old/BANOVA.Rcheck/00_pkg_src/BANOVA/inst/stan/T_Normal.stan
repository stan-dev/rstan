data {
  int<lower=0> N;
  int<lower=0> J;
  int<lower=0> M;
  int<lower=0> K;
  matrix[N, J] X;
  matrix[M, K] Z;
  int<lower=0> id[N];
  vector[N] y;
}

parameters {
  matrix[J, M] beta1;
  matrix[K, J] beta2; 
  real<lower=0> tau_yS;
  vector<lower=0>[J] tau_beta1Sq;
  real<lower=1, upper=10> df;
} 

transformed parameters {
  vector[N] y_hat;
  vector[N] y_pred;
  matrix[M, J] mu_beta1;
  mu_beta1 = Z*beta2;
  for (i in 1:N){
    y_hat[i] = X[i,]*beta1[,id[i]];
  }
  for (i in 1:N){
    y_pred[i] = X[i,]*mu_beta1[id[i],]';
  }
}

model {
  real tau_y;
  vector[J] tau_beta1;
  tau_y = sqrt(tau_yS);
  tau_beta1 = sqrt(tau_beta1Sq);
  y ~ student_t(df, y_hat, tau_y);
  df ~ normal(5,5);
  tau_yS ~ inv_gamma(1, 1);
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
  real tau_ySq;
  r_2 = variance(y_pred)/variance(y);
  tau_ySq = variance(y - y_hat);
}


