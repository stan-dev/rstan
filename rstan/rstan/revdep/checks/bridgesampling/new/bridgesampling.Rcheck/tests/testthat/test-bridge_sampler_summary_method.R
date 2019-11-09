
context('bridge sampling summary method')

test_that("bridge sampler summary method correctly displayed", {

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
  s_normal <- summary(bridge_normal)
  s_warp3 <- summary(bridge_warp3)

  expect_equal(names(s_normal), c("Logml_Estimate",
                                  "Relative_Mean_Squared_Error",
                                  "Coefficient_of_Variation",
                                  "Percentage_Error",
                                  "Method", "Repetitions"))
  expect_equal(names(s_warp3), c("Logml_Estimate", "Method", "Repetitions"))

  expect_output(print(s_normal), 'All error measures are approximate.')
  expect_output(print(s_warp3), 'No error measures are available for method = "warp3"')


  # repetitions > 1
  bridge_normal_2 <- bridge_sampler(samples = x, log_posterior = log_density,
                                  data = NULL, lb = lb, ub = ub,
                                  method = "normal", silent = TRUE,
                                  repetitions = 2)
  bridge_warp3_2 <- bridge_sampler(samples = x, log_posterior = log_density,
                                 data = NULL, lb = lb, ub = ub,
                                 method = "warp3", silent = TRUE,
                                 repetitions = 2)
  s_normal_2 <- summary(bridge_normal_2)
  s_warp3_2 <- summary(bridge_warp3_2)

  expect_equal(names(s_normal_2), c("Logml_Estimate",
                                    "Min",
                                    "Max",
                                    "Interquartile_Range",
                                    "Method", "Repetitions"))
  expect_equal(names(s_warp3_2), c("Logml_Estimate",
                                   "Min",
                                   "Max",
                                   "Interquartile_Range",
                                   "Method", "Repetitions"))

  expect_output(print(s_normal_2), 'All error measures are based on 2 estimates.')
  expect_output(print(s_warp3_2), 'All error measures are based on 2 estimates.')

})
