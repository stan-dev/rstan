
library(rstan)
# example(stanc)


stanmodelcode <- '
data {
  int<lower=0> N;
  real y[N];
}

parameters {
  real mu;
  real alpha;
  real eta;
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
                init = list(list(mu = 2)), seed = 3, thin = 1,
                sample_file = 'norm1.csv')

print(get_inits(fit))

fit <- sampling(rr, data = dat, iter = 10, chains = 1,
                init = "random", init_r = .5, seed = 3, thin = 1,
                sample_file = 'norm1.csv')

print(get_inits(fit))

fit <- sampling(rr, data = dat, iter = 10, chains = 1,
                init = "0", seed = 3, thin = 1,
                sample_file = 'norm1.csv')

print(get_inits(fit))
