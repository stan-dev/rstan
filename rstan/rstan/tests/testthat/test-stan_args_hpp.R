# FIXME: This should be uncommented (and possibly refactored) when the tests in
# this file are unskipped.

# Needed to test part of stan_args.hpp.
# src <- "
#   BEGIN_RCPP
#   Rcpp::List lst(x);
#   rstan::stan_args args(lst);
#   return args.stan_args_to_rlist();
#   END_RCPP
# "
#
# fx <- inline::cxxfunction(signature(x = "list"),
#   body = src,
#   includes = "#include <rstan/stan_args.hpp>",
#   plugin = "rstan", verbose = TRUE
# )

test_that("stan_args_hppb works", {
  skip("Backwards compatibility")

  b1 <- fx(list(iter = 100, seed = 12354, method = "optim"))
  expect_equal(b1$random_seed, "12354")
  expect_equal(b1$algorithm, "LBFGS")
  b2 <- fx(list(iter = 100, seed = 12354, method = "optim", algorithm = "BFGS"))
  b3 <- fx(list(iter = 100, seed = 12354, method = "optim", algorithm = "LBFGS"))
  b4 <- fx(list(iter = 100, seed = 12354, method = "optim", algorithm = "LBFGS", history_size = 6))
  expect_equal(b1$random_seed, "12354")
  expect_true(is.null(b2$history_size))
  expect_equal(b3$history_size, 5)
  expect_equal(b4$history_size, 6)
})

test_that("stan_args_hppvb works", {
  skip("Backwards compatibility")

  b1 <- fx(list(iter = 100, seed = 12354, method = "variational"))
  expect_equal(b1$random_seed, "12354")
  expect_equal(b1$algorithm, "meanfield")
  expect_true(b1$adapt_engaged)
  b2 <- fx(list(iter = 100, seed = 12354, method = "variational", algorithm = "fullrank"))
  expect_equal(b2$algorithm, "fullrank")
  expect_equal(b2$iter, 100)
  expect_equal(b2$grad_samples, 1)
  expect_equal(b2$eta, 1)
  b3 <- fx(list(
    iter = 101, seed = 12354, method = "variational",
    algorithm = "fullrank", grad_samples = 2,
    adapt_iter = 102,
    elbo_samples = 50, eval_elbo = 48, output_samples = 500,
    eta = .5, tol_rel_obj = 0.001
  ))
  expect_equal(b3$iter, 101)
  expect_equal(b3$adapt_iter, 102)
  expect_equal(b3$random_seed, "12354")
  expect_equal(b3$grad_samples, 2)
  expect_equal(b3$eval_elbo, 48)
  expect_equal(b3$tol_rel_obj, 0.001)
  expect_equal(b3$eta, 0.5)
  expect_equal(b3$elbo_samples, 50)
  expect_equal(b3$output_samples, 500)
})

test_that("stan_args_hpp works", {
  skip("Backwards compatibility")

  a1 <- fx(list(iter = 100, thin = 100))
  a2 <- fx(list(iter = 100, thin = 3))
  a3 <- fx(list(iter = 5, thin = 3, refresh = -1, seed = "12345"))
  a4 <- fx(list(iter = 5, thin = 3, refresh = -1, method = "test_grad"))
  a4b <- fx(list(
    iter = 5, thin = 3, refresh = -1, method = "test_grad",
    control = list(epsilon = 1.3, error = 0.1)
  ))
  a4c <- fx(list(
    iter = 5, thin = 3, refresh = -1, method = "test_grad",
    control = list(epsilon = 1.3)
  ))
  a4d <- fx(list(
    iter = 5, thin = 3, refresh = -1, method = "test_grad",
    control = list(error = 1.3)
  ))
  a5 <- fx(list(
    iter = 5, thin = 3, algorithm = "HMC",
    control = list(stepsize = .1)
  ))
  a6 <- fx(list(
    iter = 5, thin = 3, algorithm = "NUTS",
    control = list(stepsize = .1, metric = "unit_e")
  ))
  a7 <- fx(list(
    iter = 5, thin = 3, algorithm = "NUTS",
    control = list(stepsize = .1, metric = "diag_e")
  ))
  a7b <- fx(list(
    iter = 5, thin = 3, algorithm = "NUTS",
    control = list(
      stepsize = .1, metric = "diag_e",
      adapt_term_buffer = 4,
      adapt_window = 30,
      adapt_init_buffer = 40
    )
  ))
  a8 <- fx(list(
    iter = 5, thin = 3, algorithm = "NUTS",
    control = list(stepsize = .1, metric = "dense_e")
  ))
  a9 <- fx(list(
    iter = 5, thin = 3, algorithm = "HMC",
    control = list(stepsize = .1, metric = "unit_e")
  ))
  a10 <- fx(list(
    iter = 5, thin = 3, algorithm = "HMC",
    control = list(stepsize = .1, metric = "diag_e")
  ))
  a11 <- fx(list(
    iter = 5, thin = 3, algorithm = "HMC",
    control = list(stepsize = .1, metric = "dense_e")
  ))

  a12 <- fx(list(
    iter = 5, thin = 3, algorithm = "Fixed_param",
    control = list(adapt_engaged = TRUE)
  ))

  expect_equal(a1$iter, 100)
  expect_equal(a1$thin, 100)
  expect_equal(a2$iter, 100)
  expect_equal(a2$thin, 3)
  expect_equal(a2$control$adapt_init_buffer, 75)
  expect_equal(a2$control$adapt_term_buffer, 50)
  expect_equal(a2$control$adapt_window, 25)
  expect_equal(a3$iter, 5)
  expect_equal(a3$random_seed, "12345")
  expect_equal(a3$test_grad, FALSE)
  expect_equal(a3$refresh, -1)
  expect_equal(a3$control$adapt_init_buffer, 75)
  expect_equal(a3$control$adapt_term_buffer, 50)
  expect_equal(a3$control$adapt_window, 25)
  expect_equal(a4$test_grad, TRUE)
  expect_equal(a4$method, "test_grad")
  expect_equal(a4$control$epsilon, 1e-6)
  expect_equal(a4$control$error, 1e-6)
  expect_equal(a4b$control$epsilon, 1.3)
  expect_equal(a4b$control$error, 1e-1)
  expect_equal(a4c$control$epsilon, 1.3)
  expect_equal(a4c$control$error, 1e-6)
  expect_equal(a4c$control$error, 1e-6)
  expect_equal(a4c$control$epsilon, 1.3)

  expect_error(fx(list(iter = 5, thin = 3, refresh = -1, test_grad = TRUE, seed = "111111111111111111111")))

  expect_equal(a6$control$metric, "unit_e")
  expect_equal(a7$control$metric, "diag_e")
  expect_equal(a7b$control$adapt_init_buffer, 40)
  expect_equal(a7b$control$adapt_term_buffer, 4)
  expect_equal(a7b$control$adapt_window, 30)
  expect_equal(a8$control$metric, "dense_e")
  expect_equal(a9$control$metric, "unit_e")
  expect_equal(a10$control$metric, "diag_e")
  expect_equal(a11$control$metric, "dense_e")
  expect_equal(a11$control$adapt_window, 25)
})
