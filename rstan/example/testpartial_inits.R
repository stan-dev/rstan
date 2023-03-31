
library(rstan)
# example(stanc)


stanmodelcode <- '
data {
  int<lower=0> N;
  array[N] real;
}

parameters {
  real mu;
  real alpha;
  simplex[3] delta;
  real eta;
  real<lower=0> sigma;
}

transformed parameters {
  real beta;
  beta <- -alpha;
}

model {
  mu ~ normal(0, 10);
  y ~ normal(mu, 1);
  alpha ~ normal(0, 1);
  eta ~ normal(0, 1);
  sigma ~ exponential(2);
}

generated quantities {
  real gamma;
  gamma <- mu + alpha;
}
'
model_name <- "normal1";

rr <- stan_model(model_code = stanmodelcode, model_name = model_name,
                 verbose = TRUE)

y <- rnorm(20)
mean(y)
sd(y)
dat <- list(N = 20, y = y)
f <- sampling(rr, data = dat, init = 0, iter = 10, sample_file = 'norm1.csv')

fit <- sampling(rr, data = dat, iter = 10, chains = 1,
                init = list(list(mu2 = 2)), seed = 3, thin = 1,
                sample_file = 'norm1.csv')

print(get_inits(fit))

fit <- sampling(rr, data = dat, iter = 10, chains = 1,
                init = list(list(mu2 = 2)), seed = 3, thin = 1,
                sample_file = 'norm1.csv', enable_random_init = FALSE)

print(get_inits(fit))

## initialization with dimension not matching
fit <- sampling(rr, data = dat, iter = 10, chains = 1,
                init = list(list(delta = 2)), seed = 3, thin = 1,
                sample_file = 'norm1.csv', enable_random_init = !FALSE)

print(get_inits(fit))

