
context('basic bridge sampling behavior normal')

test_that("bridge sampler matches anlytical value normal example", {

  # library(bridgesampling)
  library(mvtnorm)

  x <- rmvnorm(1e4, mean = rep(0, 2), sigma = diag(2))
  colnames(x) <- c("x1", "x2")
  log_density <- function(s, data) {
    -.5*t(s)%*%s
  }
  assign(x = "log_density", value = log_density, envir = .GlobalEnv)

  lb <- rep(-Inf, 2)
  ub <- rep(Inf, 2)
  names(lb) <- names(ub) <- colnames(x)

  # check repetitions > 1
  bridge_normal <- bridge_sampler(samples = x, log_posterior = log_density,
                                  data = NULL, lb = lb, ub = ub,
                                  method = "normal", silent = TRUE, repetitions = 2)
  bridge_warp3 <- bridge_sampler(samples = x, log_posterior = log_density,
                                 data = NULL, lb = lb, ub = ub,
                                 method = "warp3", silent = TRUE, repetitions = 2)
  bridge_normal_c <- bridge_sampler(samples = x, log_posterior = "log_density",
                                    data = NULL, lb = lb, ub = ub,
                                    method = "normal", silent = TRUE, repetitions = 2)
  bridge_warp3_c <- bridge_sampler(samples = x, log_posterior = "log_density",
                                   data = NULL, lb = lb, ub = ub,
                                   method = "warp3", silent = TRUE, repetitions = 2)


  expect_equal(bridge_normal$logml, expected = rep(log(2*pi), length(bridge_normal$logml)), tolerance = 0.01)
  expect_equal(bridge_warp3$logml, expected = rep(log(2*pi), length(bridge_warp3$logml)), tolerance = 0.01)
  expect_equal(bridge_normal_c$logml, expected = rep(log(2*pi), length(bridge_normal_c$logml)), tolerance = 0.01)
  expect_equal(bridge_warp3_c$logml, expected = rep(log(2*pi), length(bridge_warp3_c$logml)), tolerance = 0.01)

  expect_equal(bf(bridge_normal, bridge_warp3)$bf, expected = rep(1, 2), tolerance = 0.1)

  # check repetitions = 1
  bridge_normal <- bridge_sampler(samples = x, log_posterior = log_density,
                                  data = NULL, lb = lb, ub = ub,
                                  method = "normal", silent = TRUE)
  bridge_warp3 <- bridge_sampler(samples = x, log_posterior = log_density,
                                 data = NULL, lb = lb, ub = ub,
                                 method = "warp3", silent = TRUE)
  bridge_normal_c <- bridge_sampler(samples = x, log_posterior = "log_density",
                                    data = NULL, lb = lb, ub = ub,
                                    method = "normal", silent = TRUE)
  bridge_warp3_c <- bridge_sampler(samples = x, log_posterior = "log_density",
                                   data = NULL, lb = lb, ub = ub,
                                   method = "warp3", silent = TRUE)

  expect_equal(bridge_normal$logml, expected = log(2*pi), tolerance = 0.01)
  expect_equal(bridge_warp3$logml, expected = log(2*pi), tolerance = 0.01)
  expect_equal(bridge_normal_c$logml, expected = log(2*pi), tolerance = 0.01)
  expect_equal(bridge_warp3_c$logml, expected = log(2*pi), tolerance = 0.01)


  # check using dots repetitions > 1
  mu <- c(1, 2)
  x <- rmvnorm(1e4, mean = mu, sigma = diag(2))
  colnames(x) <- c("x1", "x2")
  log_density <- function(s, data, ...) {
    -.5*t(s - mu) %*% (s - mu)
  }
  assign(x = "log_density", value = log_density, envir = .GlobalEnv)
  lb <- rep(-Inf, 2)
  ub <- rep(Inf, 2)
  names(lb) <- names(ub) <- colnames(x)
  bridge_normal_dots <- bridge_sampler(samples = x, log_posterior = log_density, mu,
                                       data = NULL, lb = lb, ub = ub, method = "normal",
                                       silent = TRUE, repetitions = 2)
  bridge_warp3_dots <- bridge_sampler(samples = x, log_posterior = log_density, mu,
                                      data = NULL, lb = lb, ub = ub, method = "normal",
                                      silent = TRUE, repetitions = 2)
  bridge_normal_c_dots <- bridge_sampler(samples = x, log_posterior = "log_density",
                                         mu, data = NULL, lb = lb, ub = ub,
                                         method = "normal", silent = TRUE, repetitions = 2)
  bridge_warp3_c_dots <- bridge_sampler(samples = x, log_posterior = "log_density",
                                        mu, data = NULL, lb = lb, ub = ub,
                                        method = "warp3", silent = TRUE, repetitions = 2)

  expect_equal(bridge_normal_dots$logml, expected = rep(log(2*pi), length(bridge_normal_dots$logml)), tolerance = 0.01)
  expect_equal(bridge_warp3_dots$logml, expected = rep(log(2*pi), length(bridge_warp3_dots$logml)), tolerance = 0.01)
  expect_equal(bridge_normal_c_dots$logml, expected = rep(log(2*pi), length(bridge_normal_c_dots$logml)), tolerance = 0.01)
  expect_equal(bridge_warp3_c_dots$logml, expected = rep(log(2*pi), length(bridge_warp3_c_dots$logml)), tolerance = 0.01)

  # check using dots
  mu <- c(1, 2)
  x <- rmvnorm(1e4, mean = mu, sigma = diag(2))
  colnames(x) <- c("x1", "x2")
  log_density <- function(s, data, ...) {
    -.5*t(s - mu) %*% (s - mu)
  }
  assign(x = "log_density", value = log_density, envir = .GlobalEnv)
  lb <- rep(-Inf, 2)
  ub <- rep(Inf, 2)
  names(lb) <- names(ub) <- colnames(x)
  bridge_normal_dots <- bridge_sampler(samples = x, log_posterior = log_density, mu,
                                       data = NULL, lb = lb, ub = ub, method = "normal",
                                       silent = TRUE)
  bridge_warp3_dots <- bridge_sampler(samples = x, log_posterior = log_density, mu,
                                      data = NULL, lb = lb, ub = ub, method = "normal",
                                      silent = TRUE)
  bridge_normal_c_dots <- bridge_sampler(samples = x, log_posterior = "log_density",
                                         mu, data = NULL, lb = lb, ub = ub,
                                         method = "normal", silent = TRUE)
  bridge_warp3_c_dots <- bridge_sampler(samples = x, log_posterior = "log_density",
                                        mu, data = NULL, lb = lb, ub = ub,
                                        method = "warp3", silent = TRUE)

  expect_equal(bridge_normal_dots$logml, expected = log(2*pi), tolerance = 0.01)
  expect_equal(bridge_warp3_dots$logml, expected = log(2*pi), tolerance = 0.01)
  expect_equal(bridge_normal_c_dots$logml, expected = log(2*pi), tolerance = 0.01)
  expect_equal(bridge_warp3_c_dots$logml, expected = log(2*pi), tolerance = 0.01)


  # check error_measures
  err <- error_measures(bridge_normal)
  expect_equal(names(err), c("re2", "cv", "percentage"))
  expect_is(unlist(err), "character")

  expect_error(error_measures(bridge_warp3), "not implemented for warp3")

  ### these are meant to check the bf and post_prob functions and not as a meaningful comparisons
  bf <- bf(bridge_normal, bridge_warp3)
  expect_is(bf$bf, "numeric")

  # without prior_prob
  post1 <- post_prob(bridge_normal, bridge_warp3, bridge_normal_c, bridge_warp3_c)
  expect_equal(sum(post1), 1)

  # with prior_prob
  post2 <- post_prob(bridge_normal, bridge_warp3, bridge_normal_c,
                             bridge_warp3_c, prior_prob = c(0.2, 0.1, 0.25, 0.45))
  expect_equal(sum(post2), 1)

  # with incorrect prior_prob
  expect_error(post_prob(bridge_normal, bridge_warp3, bridge_normal_c,
                                 bridge_warp3_c, prior_prob = c(0.2, 0.1, 0.25, 0.55)),
               "do not sum to one")

})


context('non-standard parameter spaces')

test_that("bridge sampler functions for non-standard parameter spaces", {


  # Test with only simplex
  ru <- replicate(10, runif(10))
  theta <- (ru / rowSums(ru))[, -10]
  colnames(theta) <- paste0("sim", 1:9)
  theta_t <- .transform2Real(theta,
                                   lb = rep(0, 9), ub = rep(1, 9),
                                   theta_types = rep("simplex", 9))

  expect_equal(theta_t$transTypes[1], c(sim1 = "simplex"))

  theta_t_t <- .invTransform2Real(theta_t$theta_t,
                                              lb = rep(0, 9), ub = rep(1, 9),
                                              theta_types = rep("simplex", 9))
  expect_equal(theta, theta_t_t)


  # tranformations work for different input shapes
  nsimp <- 4
  n     <- 100
  sum_to_one <- function(x) x / sum(x)
  ru <- t(replicate(n, c(rnorm(2), # unbounded
                       sum_to_one(runif(nsimp)), # simplex
                       runif(3), # double-bounded
                       abs(rnorm(1)), # lower-bounded
                       rnorm(2) %% (2*pi)))) # circular
  theta_original <- ru[, -(nsimp + 2)]
  pt <- c(rep("real", 2),
          rep("simplex", nsimp - 1),
          rep("real", 4),
          rep("circular", 2))
  lb <- c(rep(-Inf, 2),
          rep(0, nsimp - 1),
          rep(0, 4),
          rep(0, 2))
  ub <- c(rep(Inf, 2),
          rep(1, nsimp - 1),
          rep(1, 3),
          rep(Inf, 1),
          rep(2*pi, 2))
  nm <- c(paste0("unbounded", 1:2),
          paste0("simplex", 1:(nsimp - 1)),
          paste0("doublebounded", 1:3),
          paste0("lower", 1),
          paste0("circular", 1:2))
  colnames(theta_original) <- names(lb) <- names(ub) <- names(pt) <- nm


  theta_t <- .transform2Real(theta_original, lb, ub, pt)
  theta_t_t <- .invTransform2Real(theta_t$theta_t, lb, ub, pt)

  # The modulus is to force the circular variables to be equal if they lie on
  # the same place on the circle. The modulus is also taken for the linear
  # variables, for simplicity of programming.
  expect_equal(theta_original %% (2*pi), theta_t_t %% (2*pi))

  # Works with one row
  theta <- theta_original[1, , drop = FALSE]
  theta_t <- .transform2Real(theta, lb, ub, pt)
  theta_t_t <- .invTransform2Real(theta_t$theta_t, lb, ub, pt)

  # The modulus is to force the circular variables to be equal if they lie on
  # the same place on the circle.
  expect_equal(theta %% (2*pi), theta_t_t %% (2*pi))


  # Test bridge sampler function with non-standard sample spaces
  bs_ns <- bridge_sampler.matrix(
    theta_original,
    data = rnorm(10),
    log_posterior = function(s, data) -.5*t(s) %*% s,
    lb = lb, ub = ub, silent = TRUE, verbose = FALSE)

  expect_true(class(bs_ns) == "bridge")



  ############ TEST JACOBIAN
  n <- 2

  theta_full <- t(c(.4, .6))
  theta <- theta_full[, -n, drop = FALSE]
  colnames(theta) <- paste0("sim", (1:(n - 1)))

  y <- bridgesampling:::.transform2Real(theta,
                                        lb = rep(0, n - 1),
                                        ub = rep(1, n - 1),
                                        theta_types = rep("simplex", n - 1))$theta_t
  tt <- rep("simplex", n - 1)
  colnames(y) <- paste0("trans_sim", (1:(n - 1)))
  names(tt) <- paste0("sim", (1:(n - 1)))
  jacob <- .logJacobian(y, tt,
                        lb = rep(0, n),
                        ub = rep(1, n))

  expect_true(is.numeric(jacob))


  skip_if_not_installed("MCMCpack")

  invsimplex <- function(y) {
    y <- as.matrix(y)
    n <- length(y)
    colnames(y) <- paste0("trans_sim", (1:n))
    out1 <- .invTransform2Real(y,
                                        lb = rep(0, n),
                                        ub = rep(1, n),
                                     theta_types = rep("simplex", n))
    c(out1, 1 - sum(out1))
  }
  invsimplex(100)
  p_y <- function(y) {
    y <- as.matrix(y)
    n <- length(y)
    tt <- rep("simplex", n)
    colnames(y) <- paste0("trans_sim", (1:n))
    names(tt) <- paste0("sim", (1:n))
    MCMCpack::ddirichlet(invsimplex(y), theta_full*10) *
      exp(.logJacobian(y,
                                        tt,
                                        lb = rep(0, n),
                                        ub = rep(1, n)))
  }

  # The jaobian corrects for the transformation
  expect_equal(integrate(Vectorize(p_y), -100, 100)$value, 1)

})



