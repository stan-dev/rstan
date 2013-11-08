library(rstan)

y <- c(-0.01086, -0.08631,  1.50312,  2.06884,  0.12150, -0.41584,  1.23002,  0.06848,
       1.05639, -0.11791,  0.17318,  1.62376,  0.03792, -0.58818,  0.45146, -2.03457,
       -0.40714,  1.25568,  1.49409,  0.23507)

sd(y)
sigma <- sd(y)

f1 <- function(sigma) {
  sum(dnorm(y, sd = sigma, log = TRUE))
} 


# op <- optim(1, f1, hessian = TRUE, lower = 0, method = 'L-BFGS-B', control = list(fnscale = -1))
op <- optim(1, f1, hessian = TRUE, method = "BFGS", control = list(fnscale = -1))
op

mc <- '
  data { real y[20]; }
  parameters { real<lower=0> sigma; }
  model { y ~ normal(0, sigma); }
'

op2 <- optimizing(stan_model(model_code = mc), data = list(y = y), hessian = TRUE)
op2 

