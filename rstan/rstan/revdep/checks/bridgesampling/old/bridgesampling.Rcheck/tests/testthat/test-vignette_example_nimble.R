
context('test vignette bridgesampling_example_nimble.Rmd')

test_that("bridge sampler yields correct results", {

  testthat::skip_on_cran()
  testthat::skip_on_travis()

  # library(bridgesampling)
  if (require(nimble)) {

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
    codeH0 <- nimbleCode({
      invTau2 ~ dgamma(1, 1)
      tau2 <- 1/invTau2
      for (i in 1:20) {
        theta[i] ~ dnorm(0, sd = sqrt(tau2))
        y[i] ~ dnorm(theta[i], sd = 1)
      }
    })
    codeH1 <- nimbleCode({
      mu ~ dnorm(0, sd = 1)
      invTau2 ~ dgamma(1, 1)
      tau2 <- 1/invTau2
      for (i in 1:20) {
        theta[i] ~ dnorm(mu, sd = sqrt(tau2))
        y[i] ~ dnorm(theta[i], sd = 1)
      }
    })

    ## steps for H0:
    modelH0 <- nimbleModel(codeH0)
    modelH0$setData(y = y) # set data
    cmodelH0 <- compileNimble(modelH0) # make compiled version from generated C++

    ## steps for H1:
    modelH1 <- nimbleModel(codeH1)
    modelH1$setData(y = y) # set data
    cmodelH1 <- compileNimble(modelH1) # make compiled version from generated C++

    # build MCMC functions, skipping customization of the configuration.
    mcmcH0 <- buildMCMC(modelH0,
                        monitors = modelH0$getNodeNames(stochOnly = TRUE,
                                                        includeData = FALSE))
    mcmcH1 <- buildMCMC(modelH1,
                        monitors = modelH1$getNodeNames(stochOnly = TRUE,
                                                        includeData = FALSE))
    # compile the MCMC function via generated C++
    cmcmcH0 <- compileNimble(mcmcH0, project = modelH0)
    cmcmcH1 <- compileNimble(mcmcH1, project = modelH1)

    # run the MCMC.  This is a wrapper for cmcmc$run() and extraction of samples.
    # the object samplesH1 is actually not needed as the samples are also in cmcmcH1
    samplesH0 <- runMCMC(cmcmcH0, niter = 1e5, nburnin = 1000, nchains = 2,
                         progressBar = FALSE)
    samplesH1 <- runMCMC(cmcmcH1, niter = 1e5, nburnin = 1000, nchains = 2,
                         progressBar = FALSE)

    # compute log marginal likelihood via bridge sampling for H0
    H0.bridge <- bridge_sampler(cmcmcH0, silent = TRUE)

    # compute log marginal likelihood via bridge sampling for H1
    H1.bridge <- bridge_sampler(cmcmcH1, silent = TRUE)

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

    # load(system.file("extdata/", "vignette_example_nimble.RData",
    #                  package = "bridgesampling"))

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
