
library(rstan)
mc <- '
 parameters {
   simplex[2] a;
   real y;
 }
 model {
   y ~ normal(0,1);
   increment_log_prob(-square(a[1])-square(a[2]));
 } 
'

sm <- stan_model(model_code = mc)
op3 <- optimizing(sm, hessian = !TRUE)
op3 
op4 <- optimizing(sm, hessian = TRUE)
op4



 
