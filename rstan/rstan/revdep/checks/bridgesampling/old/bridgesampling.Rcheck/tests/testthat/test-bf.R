
context('bridge sampling bf function')

test_that("bf various basic checks", {

  # library(bridgesampling)
  library(mvtnorm)

  x <- rmvnorm(1e4, mean = rep(0, 2), sigma = diag(2))
  colnames(x) <- c("x1", "x2")
  log_density <- function(s, data) {
    -.5*t(s)%*%s
  }

  lb <- rep(-Inf, 2)
  ub <- rep(Inf, 2)
  names(lb) <- names(ub) <- colnames(x)

  # repetitions = 1
  bridge_normal <- bridge_sampler(samples = x, log_posterior = log_density,
                                  data = NULL, lb = lb, ub = ub,
                                  method = "normal", silent = TRUE)
  bridge_warp3 <- bridge_sampler(samples = x, log_posterior = log_density,
                                  data = NULL, lb = lb, ub = ub,
                                 method = "warp3", silent = TRUE)
  expect_error(bf(bridge_normal, 4), "class 'bridge' or 'bridge_list'")

  BF <- bf(bridge_normal, bridge_warp3)
  log_BF <- bf(bridge_normal, bridge_warp3, log = TRUE)

  expect_output(print(BF), "Estimated Bayes factor")
  expect_output(print(log_BF), "Estimated log Bayes factor")

  BF2 <- bayes_factor(bridge_normal, bridge_warp3)
  log_BF2 <- bayes_factor(bridge_normal, bridge_warp3, log = TRUE)

  expect_output(print(BF2), "Estimated Bayes factor")
  expect_output(print(log_BF2), "Estimated log Bayes factor")


  # repetitions > 1
  bridge_normal_mult <- bridge_sampler(samples = x, log_posterior = log_density,
                                  data = NULL, lb = lb, ub = ub,
                                  method = "normal", silent = TRUE, repetitions = 2)
  bridge_warp3_mult <- bridge_sampler(samples = x, log_posterior = log_density,
                                 data = NULL, lb = lb, ub = ub,
                                 method = "warp3", silent = TRUE, repetitions = 2)

  BF_mult <- bf(bridge_normal_mult, bridge_warp3_mult)
  log_BF_mult <- bf(bridge_normal_mult, bridge_warp3_mult, log = TRUE)

  expect_output(print(BF_mult), "based on medians")
  expect_output(print(log_BF_mult), "based on medians")

  ## bf with multi and singular objects
  expect_is(suppressWarnings(bf(bridge_normal_mult, bridge_normal)), "bf_bridge_list")
  expect_is(bf(bridge_normal, bridge_normal_mult), "bf_bridge")
  expect_error(bf(bridge_normal_mult, 4), "class 'bridge' or 'bridge_list'")

  # default
  BF <- bf(1, 2)
  log_BF <- bf(1, 2, log = TRUE)

  expect_output(print(BF), "Bayes factor")
  expect_output(print(log_BF), "Log Bayes factor")

})
