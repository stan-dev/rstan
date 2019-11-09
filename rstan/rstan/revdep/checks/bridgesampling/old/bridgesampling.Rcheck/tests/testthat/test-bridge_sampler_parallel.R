
context('basic bridge sampling behavior normal parallel')

test_that("bridge sampler matches anlytical value normal example", {

  testthat::skip_on_cran()
  testthat::skip_on_travis()

  # library(bridgesampling)
  library(mvtnorm)

  x <- rmvnorm(1e4, mean = rep(0, 2), sigma = diag(2))
  colnames(x) <- c("x1", "x2")
  log_density <- function(s, data) {
    -.5*t(s)%*%s
  }
  assign("log_density", log_density, envir = .GlobalEnv)
  lb <- rep(-Inf, 2)
  ub <- rep(Inf, 2)
  names(lb) <- names(ub) <- colnames(x)
  bridge_normal <- bridge_sampler(samples = x, log_posterior = log_density,
                                  data = NULL, lb = lb, ub = ub,
                                  method = "normal", cores = 2, silent = TRUE)
  bridge_warp3 <- bridge_sampler(samples = x, log_posterior = log_density,
                                 data = NULL, lb = lb, ub = ub,
                                 method = "warp3", cores = 2, silent = TRUE)
  bridge_normal_c <- bridge_sampler(samples = x, log_posterior = "log_density",
                                    data = NULL, lb = lb, ub = ub,
                                    method = "normal", cores = 2, silent = TRUE,
                                    envir = sys.frame(sys.nframe()))
  bridge_warp3_c <- bridge_sampler(samples = x, log_posterior = "log_density",
                                   data = NULL, lb = lb, ub = ub,
                                   method = "warp3", cores = 2, silent = TRUE,
                                   envir = sys.frame(sys.nframe()))

  expect_equal(bridge_normal$logml, expected = log(2*pi), tolerance = 0.01)
  expect_equal(bridge_warp3$logml, expected = log(2*pi), tolerance = 0.01)
  expect_equal(bridge_normal_c$logml, expected = log(2*pi), tolerance = 0.01)
  expect_equal(bridge_warp3_c$logml, expected = log(2*pi), tolerance = 0.01)

  # test dots argument
  mu <- c(1, 2)
  x <- rmvnorm(1e4, mean = mu, sigma = diag(2))
  colnames(x) <- c("x1", "x2")
  log_density <- function(s, data, ...) {
    -.5*t(s - mu) %*% (s - mu)
  }
  assign("log_density", log_density, envir = .GlobalEnv)

  lb <- rep(-Inf, 2)
  ub <- rep(Inf, 2)
  names(lb) <- names(ub) <- colnames(x)

  bridge_normal_dots <- bridge_sampler(samples = x, log_posterior = log_density,
                                       mu, data = NULL, lb = lb, ub = ub,
                                       method = "normal", cores = 2, silent = TRUE)
  bridge_warp3_dots <- bridge_sampler(samples = x, log_posterior = log_density,
                                      mu, data = NULL, lb = lb, ub = ub,
                                      method = "warp3", cores = 2, silent = TRUE)
  bridge_normal_c_dots <- bridge_sampler(samples = x, log_posterior = "log_density",
                                         mu, data = NULL, lb = lb, ub = ub,
                                         method = "normal", cores = 2, silent = TRUE,
                                         envir = sys.frame(sys.nframe()))
  # ls.str(envir = sys.frame(sys.nframe()))
  bridge_warp3_c_dots <- bridge_sampler(samples = x, log_posterior = "log_density",
                                        mu, data = NULL, lb = lb, ub = ub,
                                        method = "warp3", cores = 2, silent = TRUE,
                                        envir = sys.frame(sys.nframe()))

  expect_equal(bridge_normal_dots$logml, expected = log(2*pi), tolerance = 0.01)
  expect_equal(bridge_warp3_dots$logml, expected = log(2*pi), tolerance = 0.01)
  expect_equal(bridge_normal_c_dots$logml, expected = log(2*pi), tolerance = 0.01)
  expect_equal(bridge_warp3_c_dots$logml, expected = log(2*pi), tolerance = 0.01)

})
