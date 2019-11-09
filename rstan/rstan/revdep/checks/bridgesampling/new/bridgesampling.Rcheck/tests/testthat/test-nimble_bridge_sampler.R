
context('bridge_sampler.nimble works.')

test_that("nimble support works", {
  testthat::skip_on_cran()
  testthat::skip_on_travis()

  testthat::skip_if_not_installed("nimble")
  if (require(nimble)) {
    set.seed(12345)
    mu <- 0
    tau2 <- 0.5
    sigma2 <- 1
    n <- 20
    theta <- rnorm(n, mu, sqrt(tau2))
    y <- rnorm(n, theta, sqrt(sigma2))
    # create model
    codeH1 <- nimbleCode({
      mu ~ dnorm(0, sd = 1)
      invTau2 ~ dgamma(1, 1)
      tau2 <- 1/invTau2
      for (i in 1:20) {
        theta[i] ~ dnorm(mu, sd = sqrt(tau2))
        y[i] ~ dnorm(theta[i], sd = 1)
      }
    })

    modelH1 <- nimbleModel(codeH1)
    modelH1$setData(y = y) # set data

    # make compiled version from generated C++
    cmodelH1 <- compileNimble(modelH1)

    # build an MCMC, skipping customization of the configuration.
    mcmcH1 <- buildMCMC(modelH1,
                        monitors = modelH1$getNodeNames(stochOnly = TRUE,
                                                        includeData = FALSE))
    # compile the MCMC via generated C++
    cmcmcH1 <- compileNimble(mcmcH1, project = modelH1)

    # run the MCMC.  This is a wrapper for cmcmc$run() and extraction of samples.
    # the object samplesH1 is actually not needed as the samples are also in cmcmcH1
    samplesH1 <- runMCMC(cmcmcH1, niter = 1e5, nburnin = 1000, nchains = 2,
                         progressBar = FALSE)

    # bridge sampling
    bridge_H1 <- bridge_sampler(samples = cmcmcH1,
                                cores = 1,
                                method = "warp3",
                                repetitions = 2)

    expect_equal(bridge_H1$logml, rep(-37.7983064265064, 2), tolerance = 0.01)

  }


})
