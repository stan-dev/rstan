test_loo_moment_matching <- function() {



  stanmodelcode <- '
  data {
  int<lower=0> N;
  vector[N] x;
  }
  parameters {
  real mu;
  real<lower=0> sigma;
  }
  model {
  x ~ normal(mu,sigma);
  }
  generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
  log_lik[n] = normal_lpdf(x[n] | mu, sigma);
  }
  '

  rr <- stan_model(model_code = stanmodelcode, model_name = "m1",
                   verbose = !TRUE)

  CHAINS <- 4
  ITER <- 2000
  SEED <- 12
  set.seed(SEED)


  # generate toy data
  data_sd <- 1
  n <- as.integer(30)
  x <- rnorm(n = n, mean = 0, sd = data_sd)
  x[1] <- 11

  standata <- list(N = n, x = x, loo = 0)

  fit_rstan <- sampling(rr, data = standata,
                        chains = CHAINS, seed = SEED, iter = ITER)

  checkException(loo(fit_rstan,moment_match = TRUE, k_threshold = "error"))
  checkException(loo(fit_rstan,moment_match = "error"))

  loo3 <- loo(fit_rstan,moment_match = TRUE, save_psis = TRUE)
  checkEquals(loo3$diagnostics,loo3$psis_object$diagnostics)

}


