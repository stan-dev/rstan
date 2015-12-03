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
sm <- stan_model(model_code = stdnorm, verbose = TRUE)


# default algorithm 
vbf <- rstan:::vb(sm, data = dat, sample_file = 'vb.csv')
print(vbf)

vbf2 <- rstan:::vb(sm, data = dat, sample_file = 'vb2.csv', 
                   algorithm = "fullrank")
vbf2

vbf3 <- rstan:::vb(sm, data = dat, algorithm = "fullrank")
vbf3

vbf4 <- rstan:::vb(sm, data = dat, iter = 10001, seed = 12354, 
                   algorithm = 'fullrank', grad_samples = 2,
                   elbo_samples = 50, eval_elbo = 48, 
                   output_samples = 500, iter = 50,
                   eta = 1.0, tol_rel_obj = 0.001)

vbf4
