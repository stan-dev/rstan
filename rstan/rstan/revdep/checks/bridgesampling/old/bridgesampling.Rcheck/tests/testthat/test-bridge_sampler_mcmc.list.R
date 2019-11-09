
context('test bridge_sampler mcmc.list method')

test_that("bridge sampler matches analytical value", {

  testthat::skip_on_cran()
  testthat::skip_on_travis()

  # library(bridgesampling)
  if (require(R2jags) && require(runjags)) {

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

    ### function to get posterior samples ###

    # H1: mu != 0
    getSamplesModelH1 <- function(data, niter = 12000, nburnin = 2000,
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
                n.burnin = nburnin, n.thin = 1, progress.bar = "none")
      return(s)

    }

    getSamplesModelH1_runjags <- function(data, niter = 12000, nburnin = 2000,
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

      s <- suppressWarnings(runjags::run.jags(model = model, data = data,
                                     monitor = c("theta", "mu", "invTau2"),
                                     n.chains = 3, burnin = 2000,
                                     sample = 10000, silent.jags = TRUE))
      return(s)

    }

    ### get posterior samples ###

    # create data list for Jags
    data_H1 <- list(y = y, n = length(y), mu0 = mu0, tau20 = tau20, alpha = alpha,
                    beta = beta, sigma2 = sigma2)

    # fit model
    samples_H1 <- getSamplesModelH1(data_H1)
    samples_runjags <- getSamplesModelH1_runjags(data_H1)

    ### function for evaluating the unnormalized posterior on log scale ###
    log_posterior_H1 <- function(samples.row, data) {

      mu <- samples.row[[ "mu" ]]
      invTau2 <- samples.row[[ "invTau2" ]]
      theta <- samples.row[ paste0("theta[", seq_along(data$y), "]") ]

      sum(dnorm(data$y, theta, data$sigma2, log = TRUE)) +
        sum(dnorm(theta, mu, 1/sqrt(invTau2), log = TRUE)) +
        dnorm(mu, data$mu0, sqrt(data$tau20), log = TRUE) +
        dgamma(invTau2, data$alpha, data$beta, log = TRUE)

    }

    # specify parameter bounds
    cn <- colnames(samples_H1$BUGSoutput$sims.matrix)
    lb_H1 <- rep(-Inf, length(cn) - 1)
    ub_H1 <- rep(Inf, length(cn) - 1)
    names(lb_H1) <- names(ub_H1) <- cn[cn != "deviance"]
    lb_H1[[ "invTau2" ]] <- 0
    samples1 <- coda::as.mcmc(samples_H1)
    samples1 <- samples1[,cn != "deviance"]

    # mcmc.list
    bridge_normal <- bridge_sampler(samples = samples1, log_posterior = log_posterior_H1,
                                    data = data_H1, lb = lb_H1, ub = ub_H1,
                                    method = "normal", silent = TRUE, repetitions = 2)
    bridge_warp3 <- bridge_sampler(samples = samples1, log_posterior = log_posterior_H1,
                                   data = data_H1, lb = lb_H1, ub = ub_H1,
                                   method = "warp3", silent = TRUE, repetitions = 2)
    bridge_normal_m <- bridge_sampler(samples = samples1, log_posterior = log_posterior_H1,
                                      data = data_H1, lb = lb_H1, ub = ub_H1,
                                      method = "normal", silent = TRUE, repetitions = 2,
                                      cores = 2)
    bridge_warp3_m <- bridge_sampler(samples = samples1, log_posterior = log_posterior_H1,
                                     data = data_H1, lb = lb_H1, ub = ub_H1,
                                     method = "warp3", silent = TRUE, repetitions = 2,
                                     cores = 2)

    # mcmc
    bridge_normal_s <- bridge_sampler(samples = samples1[[1]], log_posterior = log_posterior_H1,
                                      data = data_H1, lb = lb_H1, ub = ub_H1,
                                      method = "normal", silent = TRUE, repetitions = 2)
    bridge_warp3_s <- bridge_sampler(samples = samples1[[1]], log_posterior = log_posterior_H1,
                                   data = data_H1, lb = lb_H1, ub = ub_H1,
                                   method = "warp3", silent = TRUE, repetitions = 2)
    bridge_normal_m_s <- bridge_sampler(samples = samples1[[1]], log_posterior = log_posterior_H1,
                                      data = data_H1, lb = lb_H1, ub = ub_H1,
                                      method = "normal", silent = TRUE, repetitions = 2,
                                      cores = 2)
    bridge_warp3_m_s <- bridge_sampler(samples = samples1[[1]], log_posterior = log_posterior_H1,
                                     data = data_H1, lb = lb_H1, ub = ub_H1,
                                     method = "warp3", silent = TRUE, repetitions = 2,
                                     cores = 2)

    # rjags
    bridge_normal_j <- bridge_sampler(samples = samples_H1, log_posterior = log_posterior_H1,
                                    data = data_H1, lb = lb_H1, ub = ub_H1,
                                    method = "normal", silent = TRUE, repetitions = 2)
    bridge_warp3_j <- bridge_sampler(samples = samples_H1, log_posterior = log_posterior_H1,
                                   data = data_H1, lb = lb_H1, ub = ub_H1,
                                   method = "warp3", silent = TRUE, repetitions = 2)
    bridge_normal_jm <- bridge_sampler(samples = samples_H1, log_posterior = log_posterior_H1,
                                      data = data_H1, lb = lb_H1, ub = ub_H1,
                                      method = "normal", silent = TRUE, repetitions = 2,
                                      cores = 2)
    bridge_warp3_jm <- bridge_sampler(samples = samples_H1, log_posterior = log_posterior_H1,
                                     data = data_H1, lb = lb_H1, ub = ub_H1,
                                     method = "warp3", silent = TRUE, repetitions = 2,
                                     cores = 2)

    # runjags
    bridge_normal_r <- bridge_sampler(samples = samples_runjags, log_posterior = log_posterior_H1,
                                      data = data_H1, lb = lb_H1, ub = ub_H1,
                                      method = "normal", silent = TRUE, repetitions = 2)
    bridge_warp3_r <- bridge_sampler(samples = samples_runjags, log_posterior = log_posterior_H1,
                                     data = data_H1, lb = lb_H1, ub = ub_H1,
                                     method = "warp3", silent = TRUE, repetitions = 2)
    bridge_normal_rm <- bridge_sampler(samples = samples_runjags, log_posterior = log_posterior_H1,
                                       data = data_H1, lb = lb_H1, ub = ub_H1,
                                       method = "normal", silent = TRUE, repetitions = 2,
                                       cores = 2)
    bridge_warp3_rm <- bridge_sampler(samples = samples_runjags, log_posterior = log_posterior_H1,
                                      data = data_H1, lb = lb_H1, ub = ub_H1,
                                      method = "warp3", silent = TRUE, repetitions = 2,
                                      cores = 2)

    # "exact" ml
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

    exact_logml <- log(mH1(data_H1))

    expect_equal(class(samples1), expected = "mcmc.list")
    expect_equal(
      bridge_normal$logml,
      expected = rep(exact_logml, length(bridge_normal$logml)),
      tolerance = 0.01
    )
    expect_equal(
      bridge_warp3$logml,
      expected = rep(exact_logml, length(bridge_warp3$logml)),
      tolerance = 0.01
    )
    expect_equal(
      bridge_normal_m$logml,
      expected = rep(exact_logml, length(bridge_normal_m$logml)),
      tolerance = 0.01
    )
    expect_equal(
      bridge_warp3_m$logml,
      expected = rep(exact_logml, length(bridge_warp3_m$logml)),
      tolerance = 0.01
    )

    expect_equal(class(samples1[[1]]), expected = "mcmc")
    expect_equal(
      bridge_normal_s$logml,
      expected = rep(exact_logml, length(bridge_normal_s$logml)),
      tolerance = 0.01
    )
    expect_equal(
      bridge_warp3_s$logml,
      expected = rep(exact_logml, length(bridge_warp3_s$logml)),
      tolerance = 0.01
    )
    expect_equal(
      bridge_normal_m_s$logml,
      expected = rep(exact_logml, length(bridge_normal_m_s$logml)),
      tolerance = 0.01
    )
    expect_equal(
      bridge_warp3_m_s$logml,
      expected = rep(exact_logml, length(bridge_warp3_m_s$logml)),
      tolerance = 0.01
    )

    expect_equal(class(samples_H1), expected = "rjags")
    expect_equal(
      bridge_normal_j$logml,
      expected = rep(exact_logml, length(bridge_normal_j$logml)),
      tolerance = 0.01
    )
    expect_equal(
      bridge_warp3_j$logml,
      expected = rep(exact_logml, length(bridge_warp3_j$logml)),
      tolerance = 0.01
    )
    expect_equal(
      bridge_normal_jm$logml,
      expected = rep(exact_logml, length(bridge_normal_jm$logml)),
      tolerance = 0.01
    )
    expect_equal(
      bridge_warp3_jm$logml,
      expected = rep(exact_logml, length(bridge_warp3_jm$logml)),
      tolerance = 0.01
    )

    expect_equal(class(samples_runjags), expected = "runjags")
    expect_equal(
      bridge_normal_r$logml,
      expected = rep(exact_logml, length(bridge_normal_r$logml)),
      tolerance = 0.01
    )
    expect_equal(
      bridge_warp3_r$logml,
      expected = rep(exact_logml, length(bridge_warp3_r$logml)),
      tolerance = 0.01
    )
    expect_equal(
      bridge_normal_rm$logml,
      expected = rep(exact_logml, length(bridge_normal_rm$logml)),
      tolerance = 0.01
    )
    expect_equal(
      bridge_warp3_rm$logml,
      expected = rep(exact_logml, length(bridge_warp3_rm$logml)),
      tolerance = 0.01
    )

    ### check that wrong lb and ub produce errors:
    ub_H0 <- ub_H1[-2]
    lb_H0 <- lb_H1[-1]
    expect_error(
      bridge_sampler(
        samples = samples_runjags,
        log_posterior = log_posterior_H1,
        data = data_H1,
        lb = lb_H1,
        ub = ub_H0
      ),
      "ub does not contain all parameters"
    )
    expect_error(
      bridge_sampler(
        samples = samples_runjags,
        log_posterior = log_posterior_H1,
        data = data_H1,
        lb = lb_H0,
        ub = ub_H1
      ),
      "lb does not contain all parameters"
    )
    expect_error(
      bridge_sampler(
        samples = samples1,
        log_posterior = log_posterior_H1,
        data = data_H1,
        lb = lb_H1,
        ub = ub_H0
      ),
      "ub does not contain all parameters"
    )
    expect_error(
      bridge_sampler(
        samples = samples1,
        log_posterior = log_posterior_H1,
        data = data_H1,
        lb = lb_H0,
        ub = ub_H1
      ),
      "lb does not contain all parameters"
    )
    expect_error(
      bridge_sampler(
        samples = samples_H1,
        log_posterior = log_posterior_H1,
        data = data_H1,
        lb = lb_H0,
        ub = ub_H1
      ),
      "lb does not contain all parameters"
    )
    expect_error(
      bridge_sampler(
        samples = samples_H1,
        log_posterior = log_posterior_H1,
        data = data_H1,
        lb = lb_H1,
        ub = ub_H0
      ),
      "ub does not contain all parameters"
    )

  }

})
