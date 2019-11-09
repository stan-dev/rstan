
context('test vignette bridgesampling_example_stan.Rmd')

test_that("bridge sampler yields correct results", {

  testthat::skip_on_cran()
  testthat::skip_on_travis()

  # library(bridgesampling)
  if (require(rstan)) {

    ### generate data ###
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
    stancodeH1 <- 'data {
    int<lower=1> n; // number of observations
    vector[n] y; // observations
    real mu0;
    real<lower=0> tau20;
    real<lower=0> alpha;
    real<lower=0> beta;
    real<lower=0> sigma2;
    }
    parameters {
    real mu;
    real<lower=0> tau2; // group-level variance
    vector[n] theta; // participant effects
    }
    model {
    target += normal_lpdf(mu | mu0, sqrt(tau20));
    target += inv_gamma_lpdf(tau2 | alpha, beta);
    target += normal_lpdf(theta | mu, sqrt(tau2));
    target += normal_lpdf(y | theta, sqrt(sigma2));
    }
    '
    # compile models
    stanmodelH0 <- stan_model(model_code = stancodeH0, model_name="stanmodel")
    stanmodelH1 <- stan_model(model_code = stancodeH1, model_name="stanmodel")

    # fit models
    stanfitH0 <- sampling(stanmodelH0, data = list(y = y, n = n,
                                                   alpha = alpha,
                                                   beta = beta,
                                                   sigma2 = sigma2),
                          iter = 50000, warmup = 1000, chains = 3, cores = 1)
    stanfitH1 <- sampling(stanmodelH1, data = list(y = y, n = n,
                                                   mu0 = mu0,
                                                   tau20 = tau20,
                                                   alpha = alpha,
                                                   beta = beta,
                                                   sigma2 = sigma2),
                          iter = 50000, warmup = 1000, chains = 3, cores = 1)

    # compute log marginal likelihood via bridge sampling for H0
    H0.bridge <- bridge_sampler(stanfitH0, silent = TRUE)

    # compute log marginal likelihood via bridge sampling for H1
    H1.bridge <- bridge_sampler(stanfitH1, silent = TRUE)

    # compute percentage errors
    H0.error <- error_measures(H0.bridge)$percentage
    H1.error <- error_measures(H1.bridge)$percentage

    # compute Bayes factor
    BF01 <- bf(H0.bridge, H1.bridge)

    # compute posterior model probabilities (assuming equal prior model probabilities)
    post1 <- post_prob(H0.bridge, H1.bridge)

    # compute posterior model probabilities (using user-specified prior model probabilities)
    post2 <- post_prob(H0.bridge, H1.bridge, prior_prob = c(.6, .4))

    # "exact" ml H1
    mH1 <- function(data, rel.tol = 1e-10) {

      y <- data$y
      n <- data$n
      mu0 <- data$mu0
      tau20 <- data$tau20
      alpha <- data$alpha
      beta <- data$beta
      sigma2 <- data$sigma2

      mH1integrand <- function(tau2, y, sigma2, mu0, tau20, alpha, beta) {

        (sigma2 + tau2)^(-n/2) *
          exp(-1/2 * ((n*mean(y)^2 + (n - 1)*sd(y)^2)/(sigma2 + tau2) +
                        mu0^2/tau20 - ((n*mean(y))/(sigma2 + tau2) +
                                         mu0/tau20)^2 /
                        (n/(sigma2 + tau2) + 1/tau20))) *
          (n/(sigma2 + tau2) + 1/tau20)^(-1/2) * tau2^(-alpha - 1) *
          exp(-beta/tau2)

      }

      (2*pi)^(-n/2) * (tau20)^(-1/2) * beta^alpha/gamma(alpha) *
        integrate(mH1integrand, 0, Inf,
                  rel.tol = rel.tol, y = y,
                  sigma2 = sigma2, mu0 = mu0,
                  tau20 = tau20, alpha = alpha,
                  beta = beta)$value

    }

    exact_logmlH1 <- log(mH1(list(y = y, n = n,
                                  mu0 = mu0,
                                  tau20 = tau20,
                                  alpha = alpha,
                                  beta = beta,
                                  sigma2 = sigma2)))

    # "exact" ml H1
    mH0 <- function(data, rel.tol = 1e-10) {

      y <- data$y
      n <- data$n
      alpha <- data$alpha
      beta <- data$beta
      sigma2 <- data$sigma2

      mH0integrand <- function(tau2, y, sigma2, alpha, beta) {

        n <- length(y)

        (sigma2 + tau2)^(-n/2) * exp(-(n*mean(y)^2 + (n - 1)*sd(y)^2)/
                                       (2*(sigma2 + tau2))) *
          tau2^(-alpha - 1) * exp(-beta/tau2)

      }

      (2*pi)^(-n/2) * beta^alpha/gamma(alpha) *
        integrate(mH0integrand, 0, Inf, rel.tol = rel.tol,
                  y = y, sigma2 = sigma2, alpha = alpha,
                  beta = beta)$value

    }

    exact_logmlH0 <- log(mH0(list(y = y, n = n,
                                  alpha = alpha,
                                  beta = beta,
                                  sigma2 = sigma2)))

    exact_BF01 <- exp(exact_logmlH0 - exact_logmlH1)

    H0.bridge.curr <- H0.bridge
    H1.bridge.curr <- H1.bridge
    BF01.curr <- BF01
    post1.curr <- post1
    post2.curr <- post2

    load(system.file("extdata/", "vignette_example_stan.RData",
                     package = "bridgesampling"))

    expect_equal(
      H0.bridge.curr$logml,
      expected = exact_logmlH0,
      tolerance = 0.01
    )
    expect_equal(
      H1.bridge.curr$logml,
      expected = exact_logmlH1,
      tolerance = 0.01
    )
    expect_equal(
      BF01.curr$bf,
      expected = exact_BF01,
      tolerance = 0.01
    )
    expect_equal(
      H0.bridge.curr$logml,
      expected = H0.bridge$logml,
      tolerance = 0.01
    )
    expect_equal(
      H1.bridge.curr$logml,
      expected = H1.bridge$logml,
      tolerance = 0.01
    )
    expect_equal(
      BF01.curr$bf,
      expected = BF01$bf,
      tolerance = 0.01
    )
    expect_equal(
      post1.curr,
      expected = post1,
      tolerance = 0.01
    )
    expect_equal(
      post2.curr,
      expected = post2,
      tolerance = 0.01
    )

  }

})
