library(rstan)
stdnorm <- '
data {
  int N;
  real y[N];
}

parameters {
  real mu;
  real<lower=0> sigma;
}

model {
  y ~ normal(mu, sigma);
}
'

N <- 30
y <- rnorm(N)
cat("mean(y)=", mean(y), ", and sd(y)=", sd(y), "\n", sep = '')
dat <- list(N = N, y = y)
dump(c("N", "y"), file = 'optim.data.R')
sm <- stan_model(model_code = stdnorm)

optim1 <- rstan:::optimizing(sm, data = dat, algorithm = "Newton")
print(optim1)

optim2 <- rstan:::optimizing(sm, data = dat, algorithm = 'BFGS', 
                             sample_file = 'opt.csv', init_alpha = 0.02,
                             tol_obj = 1e-7, tol_grad=1e-9, tol_param=1e-7)
print(optim2)

optim3 <- rstan:::optimizing(sm, data = dat, algorithm = 'Nesterov', stepsize = 0.5)
print(optim3)

optim4 <- rstan:::optimizing(sm, data = dat, algorithm = 'Nesterov', stepsize = 0.5, as_vector = FALSE)
print(optim4)





