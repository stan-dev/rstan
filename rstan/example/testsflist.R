library(rstan)
scode <- "
parameters {
  real y[2];
}
model {
  y[1] ~ normal(0, 1);
  y[2] ~ double_exponential(0, 2);
}
"
seed <- 123 # or any other integer
f1 <- stan(model_code = scode, chains = 1, seed = seed, chain_id = 1)
f2 <- stan(fit = f1, chains = 2, seed = seed, chain_id = 2:3)
f12 <- sflist2stanfit(list(f1, f2))
f3 <- stan(fit = f1, test_grad = TRUE)
f123 <- sflist2stanfit(list(f1, f2, f3))

## parallel stan call
library(parallel)

## this may not work on Windows
sflist1 <-
  mclapply(1:4, mc.cores = 4,
           function(i) stan(fit = f1, seed = seed, chains = 1, chain_id = i, refresh = -1))

sf2 <- sflist2stanfit(sflist1)
