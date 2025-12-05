# an example of using Fixed_param algorithm 
code <- '
model {
} 

generated quantities {
  real y;
  y <- normal_rng(0, 1);
}
'

library(rstan)

fit <- stan(model_code = code)

fit2 <- stan(fit = fit, algorithm = 'Fixed_param')
print(fit2)
