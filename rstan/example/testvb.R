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
vbf <- rstan:::variational(sm, data = dat, sample_file = 'vb.csv')
print(vbf)

vbf2 <- rstan:::variational(sm, data = dat, sample_file = 'vb2.csv', algorithm = "fullrank")
vbf2

