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

optim2 <- rstan:::optimizing(sm, data = dat, algorithm = 'BFGS')
print(optim2)

optim3 <- rstan:::optimizing(sm, data = dat, algorithm = 'Nesterov')
print(optim3)


