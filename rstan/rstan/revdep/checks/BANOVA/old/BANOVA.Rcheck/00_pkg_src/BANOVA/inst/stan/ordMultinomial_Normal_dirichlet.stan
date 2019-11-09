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
  simplex[cat-2] theta; 
  real<lower=0> scale;
  vector<lower=0>[cat-2] alpha;
  matrix[J, M] beta1;
  matrix[K, J] beta2; 
  vector<lower=0>[J] tau_beta1Sq;
}
transformed parameters {
  ordered[cat-1] c_trans;
  vector[cat-2] theta_sum;
  c_trans[1] = 0;
  theta_sum = cumulative_sum(theta);
  for (i in 2:(cat-1))
    c_trans[i] = theta_sum[i-1] * scale; 
}

model {
  vector[N] y_hat;
  matrix[M, J] mu_beta1;
  vector[J] tau_beta1;
  tau_beta1 = sqrt(tau_beta1Sq);
  for (i in 1:N){
    y_hat[i] = X[i,]*beta1[,id[i]];
  }
  for (i in 1:N){
      y[i] ~ ordered_logistic(y_hat[i], c_trans);
  }
  theta ~ dirichlet(alpha);
  scale ~ uniform(0, 100);
  alpha ~ normal(0, 10);
  mu_beta1 = Z*beta2;
  for (i in 1:J){
    beta1[i,] ~ normal(mu_beta1[,i], tau_beta1[i]);
  }
  tau_beta1Sq ~ inv_gamma(1, 1);
  for (i in 1:J){
    beta2[,i] ~ normal(0, 100);
  }
}