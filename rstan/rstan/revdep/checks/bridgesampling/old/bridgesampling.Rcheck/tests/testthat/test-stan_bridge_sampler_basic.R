
context('bridge_sampler.stanfit works.')


### H0: mu = 0
mH0 <- function(y, sigma2 = 1, alpha = 2, beta = 3, rel.tol = 10^(-10)) {
  n <- length(y)
  mH0integrand <- function(tau2, y, sigma2, alpha, beta) {
    (sigma2 + tau2)^(-n/2) * exp(-(n*mean(y)^2 + (n - 1)*sd(y)^2)/(2*(sigma2 + tau2))) *
      tau2^(-alpha - 1) * exp(-beta/tau2)
  }
  (2*pi)^(-n/2) * beta^alpha/gamma(alpha) * integrate(mH0integrand, 0, Inf, rel.tol = rel.tol,
                                                      y = y, sigma2 = sigma2, alpha = alpha,
                                                      beta = beta)$value
}


test_that("stan_bridge_sampler", {
  if (require(rstan)) {
    set.seed(12345)

    mu <- 0
    tau2 <- 0.5
    sigma2 <- 1

    n <- 20
    theta <- rnorm(n, mu, sqrt(tau2))
    y <- rnorm(n, theta, sqrt(sigma2))

    ### set prior parameters ###
    mu0 <- 0
    tau20 <- 1
    alpha <- 1
    beta <- 1

    # models
    stancodeH0 <- 'data {
    int<lower=1> n; // number of observations
    vector[n] y; // observations
    real<lower=0> alpha;
    real<lower=0> beta;
    real<lower=0> sigma2;
    }
    parameters {
    real<lower=0> tau2; // group-level variance
    vector[n] theta; // participant effects
    }
    model {
    target += inv_gamma_lpdf(tau2 | alpha, beta);
    target += normal_lpdf(theta | 0, sqrt(tau2));
    target += normal_lpdf(y | theta, sqrt(sigma2));
    }
    '

    # compile models
    tmp <- capture.output(suppressMessages(
    stanmodelH0 <- stan_model(model_code = stancodeH0, model_name="stanmodel")
    ))

    # fit models
    tmp <- capture.output(
    stanobjectH0 <- sampling(stanmodelH0, data = list(y = y, n = n,
                                                      alpha = alpha,
                                                      beta = beta),
                             iter = 3500, warmup = 500, chains = 4, show_messages = FALSE))
    expect_is(
    H0_bridge_norm <- bridge_sampler(stanobjectH0, method = "normal", silent = TRUE)
    , "bridge")

    expect_is(
      H0_bridge_norm_rep <-bridge_sampler(stanobjectH0, method = "normal", repetitions = 2, silent = TRUE)
      , "bridge_list")

    expect_is(
    H0_bridge_warp3 <- bridge_sampler(stanobjectH0, method = "warp3", silent = TRUE)
    , "bridge")

    expect_is(
      H0_bridge_warp3_rep <- bridge_sampler(stanobjectH0, method = "warp3", repetitions = 2, silent = TRUE)
      , "bridge_list")

    expect_equal(H0_bridge_norm$logml, log(mH0(y = y, sigma2 = sigma2, alpha = alpha, beta = beta)), tolerance = 0.1)
    expect_equal(H0_bridge_warp3$logml, log(mH0(y = y, sigma2 = sigma2, alpha = alpha, beta = beta)), tolerance = 0.1)

    expect_equal(H0_bridge_norm_rep$logml, rep(log(mH0(y = y, sigma2 = sigma2, alpha = alpha, beta = beta)), 2),
                 tolerance = 0.1)
    expect_equal(H0_bridge_warp3_rep$logml, rep(log(mH0(y = y, sigma2 = sigma2, alpha = alpha, beta = beta)), 2),
                 tolerance = 0.1)

  }
})


test_that("stan_bridge_sampler in multicore", {
  testthat::skip_on_cran()
  testthat::skip_on_travis()
  #testthat::skip_on_os("windows")
  if (require(rstan)) {
    set.seed(12345)

    mu <- 0
    tau2 <- 0.5
    sigma2 <- 1

    n <- 20
    theta <- rnorm(n, mu, sqrt(tau2))
    y <- rnorm(n, theta, sqrt(sigma2))

    ### set prior parameters ###
    mu0 <- 0
    tau20 <- 1
    alpha <- 1
    beta <- 1

    # models
    stancodeH0 <- 'data {
    int<lower=1> n; // number of observations
    vector[n] y; // observations
    real<lower=0> alpha;
    real<lower=0> beta;
    real<lower=0> sigma2;
    }
    parameters {
    real<lower=0> tau2; // group-level variance
    vector[n] theta; // participant effects
    }
    model {
    target += inv_gamma_lpdf(tau2 | alpha, beta);
    target += normal_lpdf(theta | 0, sqrt(tau2));
    target += normal_lpdf(y | theta, sqrt(sigma2));
    }
    '

    # compile models
    tmp <- capture.output(suppressMessages(
    stanmodelH0 <- stan_model(model_code = stancodeH0, model_name="stanmodel")
    ))

    # fit models
    tmp <- capture.output(
    stanobjectH0 <- sampling(stanmodelH0, data = list(y = y, n = n,
                                                      alpha = alpha,
                                                      beta = beta),
                             iter = 2500, warmup = 500, chains = 4, show_messages = FALSE))
    expect_is(
    H0_bridge_norm <- bridge_sampler(stanobjectH0, method = "normal", silent = TRUE, cores = 2)
    , "bridge")

    expect_is(
    H0_bridge_warp3 <- bridge_sampler(stanobjectH0, method = "warp3", silent = TRUE, cores = 2)
    , "bridge")

    expect_equal(H0_bridge_norm$logml, log(mH0(y = y, sigma2 = sigma2, alpha = alpha, beta = beta)), tolerance = 0.1)
    expect_equal(H0_bridge_warp3$logml, log(mH0(y = y, sigma2 = sigma2, alpha = alpha, beta = beta)), tolerance = 0.1)

  }
})
