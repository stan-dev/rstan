
context('basic bridge sampling behavior normal Rcpp')

test_that("bridge sampler matches anlytical value normal example", {

  testthat::skip_on_cran()
  testthat::skip_on_travis()

  # library(bridgesampling)
  library(mvtnorm)
  if(require(RcppEigen)) {

    x <- rmvnorm(1e4, mean = rep(0, 2), sigma = diag(2))
    colnames(x) <- c("x1", "x2")

    lb <- rep(-Inf, 2)
    ub <- rep(Inf, 2)
    names(lb) <- names(ub) <- colnames(x)

    Rcpp::sourceCpp("unnormalized_normal_density.cpp")

    bridge_normal <- bridge_sampler(samples = x, log_posterior = log_densityRcpp,
                                    data = NULL, lb = lb, ub = ub,
                                    method = "normal",
                                    silent = TRUE)
    bridge_warp3 <- bridge_sampler(samples = x, log_posterior = log_densityRcpp,
                                   data = NULL, lb = lb, ub = ub,
                                   method = "warp3",
                                   silent = TRUE)

    expect_equal(bridge_normal$logml, expected = log(2*pi), tolerance = 0.01)
    expect_equal(bridge_warp3$logml, expected = log(2*pi), tolerance = 0.01)

    # test dots argument
    mu <- c(1, 2)
    x <- rmvnorm(1e4, mean = mu, sigma = diag(2))
    colnames(x) <- c("x1", "x2")

    lb <- rep(-Inf, 2)
    ub <- rep(Inf, 2)
    names(lb) <- names(ub) <- colnames(x)

    Rcpp::sourceCpp("unnormalized_normal_density_mu.cpp")

    bridge_normal_dots <- bridge_sampler(samples = x, log_posterior = log_densityRcpp_mu,
                                         mu, data = NULL, lb = lb, ub = ub,
                                         method = "normal",
                                         silent = TRUE)
    bridge_warp3_dots <- bridge_sampler(samples = x, log_posterior = log_densityRcpp_mu,
                                        mu, data = NULL, lb = lb, ub = ub,
                                        method = "warp3",
                                        silent = TRUE)

    expect_equal(bridge_normal_dots$logml, expected = log(2*pi), tolerance = 0.01)
    expect_equal(bridge_warp3_dots$logml, expected = log(2*pi), tolerance = 0.01)

  }

})
