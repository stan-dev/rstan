\name{BANOVA.build}
\alias{BANOVA.build}
\title{Build BANOVA models}
\description{
\code{BANOVA.build} builds(compiles) BANOVA models.
}
\usage{
BANOVA.build(BANOVA_model)
}
\arguments{
  \item{BANOVA_model}{an object of class \code{"BANOVA.model"}}
}

\value{
  \code{BANOVA.build} returns an object of class \code{"BANOVA.build"}. The returned object is a list containing:
  
  \item{stanmodel}{the compiled 'Stan' model}
  \item{model_name}{the model name}
  \item{single_level}{if the model is a single level model}

}

\examples{
\donttest{
model <- BANOVA.model('Poisson', single_level = FALSE) 
Poisson_model <- BANOVA.build(model)
# visualize the model
cat(model$model_code)
# modify the model code and rebuild
# be careful to change any parameters
model$model_code <-"
data {
  int<lower=0> N;
  int<lower=0> J;
  int<lower=0> M;
  int<lower=0> K;
  matrix[N, J] X;
  matrix[M, K] Z;
  int<lower=0> id[N];
  int y[N];
}

parameters {
  matrix[J, M] beta1;
  matrix[K, J] beta2; 
  vector<lower=0>[J] tau_beta1Sq;
} 

model {
  vector[N] y_hat;
  matrix[M, J] mu_beta1;
  vector[J] tau_beta1;
  tau_beta1 = sqrt(tau_beta1Sq);
  for (i in 1:N){
    y_hat[i] = X[i,]*beta1[,id[i]];
  }
  y ~ poisson_log(y_hat);
  mu_beta1 = Z*beta2;
  for (i in 1:J){
    beta1[i,] ~ normal(mu_beta1[,i], tau_beta1[i]);
  }
  tau_beta1Sq ~ inv_gamma(1, 1);
  for (i in 1:J){
    beta2[,i] ~ normal(0, 10);
  }
}
"
Poisson_model_new <- BANOVA.build(model)

}
}



