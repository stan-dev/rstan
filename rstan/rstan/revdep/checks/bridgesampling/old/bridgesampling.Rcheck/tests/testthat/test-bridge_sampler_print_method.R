
context('bridge sampling print method')

test_that("bridge sampler print method correctly displayed", {

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

  expect_output(print(bridge_normal), "Bridge sampling estimate of the log marginal likelihood")
  expect_output(print(bridge_warp3), "Bridge sampling estimate of the log marginal likelihood")

  # repetitions > 1
  bridge_normal <- bridge_sampler(samples = x, log_posterior = log_density,
                                  data = NULL, lb = lb, ub = ub,
                                  method = "normal", silent = TRUE, repetitions = 2)
  bridge_warp3 <- bridge_sampler(samples = x, log_posterior = log_density,
                                 data = NULL, lb = lb, ub = ub,
                                 method = "warp3", silent = TRUE, repetitions = 2)

  expect_output(print(bridge_normal), "Median of")
  expect_output(print(bridge_warp3), "Median of")

})

test_that("prints with NAs with warning.", {
  bridge_o <- structure(list(logml = c(4291.14352476047, 4293.29076119542,
4291.96372581169, 4293.02187182362, NA, NA, 4290.9761730488,
4293.32075269401, 4293.5762219227, 4294.02761288449), niter = c(104,
16, 52, 8, 1000, 1000, 167, 16, 21, 44), method = "normal", repetitions = 10), .Names = c("logml",
"niter", "method", "repetitions"), class = "bridge_list")

  expect_warning(print(bridge_o), "NA")
})
