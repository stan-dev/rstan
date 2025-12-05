test_that("loo moment matching works", {
  skip("Backwards compatibility")

  stan_model_code <- "
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
      for (n in 1:N) {
        log_lik[n] = normal_lpdf(x[n] | mu, sigma);
      }
   }
  "

  rr <- stan_model(
    model_code = stan_model_code,
    model_name = "m1",
    verbose = FALSE
  )

  CHAINS <- 4
  ITER <- 2000
  SEED <- 12
  set.seed(SEED)

  # Generate toy data.
  data_sd <- 1
  n <- as.integer(30)
  x <- rnorm(n = n, mean = 0, sd = data_sd)
  x[[1]] <- 11

  stan_data <- list(N = n, x = x, loo = 0)

  fit_rstan <- sampling(
    rr,
    data = stan_data,
    chains = CHAINS,
    seed = SEED,
    iter = ITER
  )

  expect_error(loo(fit_rstan, moment_match = TRUE, k_threshold = "error"))
  expect_error(loo(fit_rstan, moment_match = "error"))

  loo_3 <- loo(fit_rstan, moment_match = TRUE, save_psis = TRUE)
  expect_equal(loo_3$diagnostics, loo_3$psis_object$diagnostics)
})
