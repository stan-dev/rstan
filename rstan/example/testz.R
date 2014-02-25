library(rstan)

# an example for test issue 51
stanmodelcode <- '
data {
  int<lower=0> n;
}
parameters {
  real a0;
  real a[n];
}
model {
   a0 ~ normal(0,1);
   a ~ normal(a0,1);
}
'

dat <- list(n = as.integer(0))
model_name <- "normal1"; 


stan(model_code = stanmodelcode, model_name = model_name, 
     data = dat, sample_file = 'tz.csv')

