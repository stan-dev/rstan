
context('test vignette bridgesampling_example_jags.Rmd')

test_that("bridge sampler yields correct results", {

  testthat::skip_on_cran()
  testthat::skip_on_travis()

  # library(bridgesampling)
  if (require(R2jags)) {

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

    ### functions to get posterior samples ###

    # H0: mu = 0
    getSamplesModelH0 <- function(data,
                                  niter = 52000,
                                  nburnin = 2000,
                                  nchains = 3) {

      model <- "
      model {
      for (i in 1:n) {
      theta[i] ~ dnorm(0, invTau2)
      y[i] ~ dnorm(theta[i], 1/sigma2)
      }
      invTau2 ~ dgamma(alpha, beta)
      tau2 <- 1/invTau2
      }"

      s <- jags(data, parameters.to.save = c("theta", "invTau2"),
                model.file = textConnection(model),
                n.chains = nchains, n.iter = niter,
                n.burnin = nburnin, n.thin = 1)

      return(s)

    }

    # H1: mu != 0
    getSamplesModelH1 <- function(data,
                                  niter = 52000,
                                  nburnin = 2000,
                                  nchains = 3) {

      model <- "
      model {
      for (i in 1:n) {
      theta[i] ~ dnorm(mu, invTau2)
      y[i] ~ dnorm(theta[i], 1/sigma2)
      }
      mu ~ dnorm(mu0, 1/tau20)
      invTau2 ~ dgamma(alpha, beta)
      tau2 <- 1/invTau2
      }"

      s <- jags(data, parameters.to.save = c("theta", "mu", "invTau2"),
                model.file = textConnection(model),
                n.chains = nchains, n.iter = niter,
                n.burnin = nburnin, n.thin = 1)

      return(s)

    }

    ### get posterior samples ###

    # create data lists for JAGS
    data_H0 <- list(y = y, n = length(y), alpha = alpha,
                    beta = beta, sigma2 = sigma2)
    data_H1 <- list(y = y, n = length(y), mu0 = mu0,
                    tau20 = tau20, alpha = alpha,
                    beta = beta, sigma2 = sigma2)

    # fit models
    samples_H0 <- getSamplesModelH0(data_H0)
    samples_H1 <- getSamplesModelH1(data_H1)

    ### functions for evaluating the unnormalized posteriors on log scale ###

    log_posterior_H0 <- function(samples.row, data) {

      mu <- 0
      invTau2 <- samples.row[[ "invTau2" ]]
      theta <- samples.row[ paste0("theta[", seq_along(data$y), "]") ]

      sum(dnorm(data$y, theta, data$sigma2, log = TRUE)) +
        sum(dnorm(theta, mu, 1/sqrt(invTau2), log = TRUE)) +
        dgamma(invTau2, data$alpha, data$beta, log = TRUE)

    }

    log_posterior_H1 <- function(samples.row, data) {

      mu <- samples.row[[ "mu" ]]
      invTau2 <- samples.row[[ "invTau2" ]]
      theta <- samples.row[ paste0("theta[", seq_along(data$y), "]") ]

      sum(dnorm(data$y, theta, data$sigma2, log = TRUE)) +
        sum(dnorm(theta, mu, 1/sqrt(invTau2), log = TRUE)) +
        dnorm(mu, data$mu0, sqrt(data$tau20), log = TRUE) +
        dgamma(invTau2, data$alpha, data$beta, log = TRUE)

    }

    # specify parameter bounds H0
    cn <- colnames(samples_H0$BUGSoutput$sims.matrix)
    cn <- cn[cn != "deviance"]
    lb_H0 <- rep(-Inf, length(cn))
    ub_H0 <- rep(Inf, length(cn))
    names(lb_H0) <- names(ub_H0) <- cn
    lb_H0[[ "invTau2" ]] <- 0

    # specify parameter bounds H1
    cn <- colnames(samples_H1$BUGSoutput$sims.matrix)
    cn <- cn[cn != "deviance"]
    lb_H1 <- rep(-Inf, length(cn))
    ub_H1 <- rep(Inf, length(cn))
    names(lb_H1) <- names(ub_H1) <- cn
    lb_H1[[ "invTau2" ]] <- 0

    # compute log marginal likelihood via bridge sampling for H0
    H0.bridge <- bridge_sampler(samples = samples_H0, data = data_H0,
                                log_posterior = log_posterior_H0, lb = lb_H0,
                                ub = ub_H0, silent = TRUE)

    # compute log marginal likelihood via bridge sampling for H1
    H1.bridge <- bridge_sampler(samples = samples_H1, data = data_H1,
                                log_posterior = log_posterior_H1, lb = lb_H1,
                                ub = ub_H1, silent = TRUE)

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

    exact_logmlH1 <- log(mH1(data_H1))

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

    exact_logmlH0 <- log(mH0(data_H0))

    exact_BF01 <- exp(exact_logmlH0 - exact_logmlH1)

    H0.bridge.curr <- H0.bridge
    H1.bridge.curr <- H1.bridge
    BF01.curr <- BF01
    post1.curr <- post1
    post2.curr <- post2

    load(system.file("extdata/", "vignette_example_jags.RData",
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
