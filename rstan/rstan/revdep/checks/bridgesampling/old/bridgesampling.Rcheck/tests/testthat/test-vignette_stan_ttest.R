
context('test vignette bridgesampling_stan_ttest.Rmd')

test_that("bridge sampler yields correct results", {

  testthat::skip_on_cran()
  testthat::skip_on_travis()

  # library(bridgesampling)
  if (require(rstan) && require(BayesFactor)) {

    set.seed(12345)

    # Sleep data from t.test example
    data(sleep)

    # compute difference scores
    y <- sleep$extra[sleep$group == 2] - sleep$extra[sleep$group == 1]
    n <- length(y)

    # models
    stancodeH0 <- '
    data {
    int<lower=1> n; // number of observations
    vector[n] y; // observations
    }
    parameters {
    real<lower=0> sigma2; // variance parameter
    }
    model {
    target += log(1/sigma2); // Jeffreys prior on sigma2
    target += normal_lpdf(y | 0, sqrt(sigma2)); // likelihood
    }
    '
    stancodeH1 <- '
    data {
    int<lower=1> n; // number of observations
    vector[n] y; // observations
    real<lower=0> r; // Cauchy prior scale
    }
    parameters {
    real delta;
    real<lower=0> sigma2;// variance parameter
    }
    model {
    target += cauchy_lpdf(delta | 0, r); // Cauchy prior on delta
    target += log(1/sigma2); // Jeffreys prior on sigma2
    target += normal_lpdf(y | delta*sqrt(sigma2), sqrt(sigma2));  // likelihood
    }
    '
    # compile models
    stanmodelH0 <- stan_model(model_code = stancodeH0, model_name="stanmodel")
    stanmodelH1 <- stan_model(model_code = stancodeH1, model_name="stanmodel")

    # fit models
    stanfitH0 <- sampling(stanmodelH0, data = list(y = y, n = n),
                          iter = 20000, warmup = 1000, chains = 4, cores = 1,
                          control = list(adapt_delta = .99))
    stanfitH1 <- sampling(stanmodelH1, data = list(y = y, n = n, r = 1/sqrt(2)),
                          iter = 20000, warmup = 1000, chains = 4, cores = 1,
                          control = list(adapt_delta = .99))

    set.seed(12345)
    suppressWarnings(H0 <- bridge_sampler(stanfitH0, silent = TRUE))
    H1 <- bridge_sampler(stanfitH1, silent = TRUE)

    # compute percentage errors
    H0.error <- error_measures(H0)$percentage
    H1.error <- error_measures(H1)$percentage

    # compute Bayes factor
    BF10 <- bf(H1, H0)

    # BayesFactor result
    BF10.BayesFactor <- extractBF(ttestBF(y), onlybf = TRUE, logbf = FALSE)

    # one-sided test
    stancodeHplus <- '
    data {
    int<lower=1> n; // number of observations
    vector[n] y; // observations
    real<lower=0> r; // Cauchy prior scale
    }
    parameters {
    real<lower=0> delta; // constrained to be positive
    real<lower=0> sigma2;// variance parameter
    }
    model {
    target += cauchy_lpdf(delta | 0, r) - cauchy_lccdf(0 | 0, r); // Cauchy prior on delta
    target += log(1/sigma2); // Jeffreys prior on sigma2
    target += normal_lpdf(y | delta*sqrt(sigma2), sqrt(sigma2));  // likelihood
    }
    '
    # compile and fit model
    stanmodelHplus <- stan_model(model_code = stancodeHplus, model_name="stanmodel")
    stanfitHplus <- sampling(stanmodelHplus, data = list(y = y, n = n, r = 1/sqrt(2)),
                             iter = 30000, warmup = 1000, chains = 4,
                             control = list(adapt_delta = .99))

    Hplus <- bridge_sampler(stanfitHplus, silent = TRUE)
    Hplus.error <- error_measures(Hplus)$percentage

    # compute Bayes factor
    BFplus0 <- bf(Hplus, H0)

    BFplus0.BayesFactor <- extractBF(ttestBF(y, nullInterval = c(0, Inf)),
                                     onlybf = TRUE, logbf = FALSE)[1]

    H0.curr <- H0
    H1.curr <- H1
    Hplus.curr <- Hplus
    BF10.curr <- BF10
    BFplus0.curr <- BFplus0

    load(system.file("extdata/", "vignette_stan_ttest.RData",
                     package = "bridgesampling"))

    expect_equal(
      H0.curr$logml,
      expected = H0$logml,
      tolerance = 0.01
    )
    expect_equal(
      H1.curr$logml,
      expected = H1$logml,
      tolerance = 0.01
    )
    expect_equal(
      BF10.curr$bf,
      expected = BF10$bf,
      tolerance = 0.01
    )
    expect_equal(
      BF10.curr$bf,
      expected = BF10.BayesFactor,
      tolerance = 0.03
    )
    expect_equal(
      BFplus0.curr$bf,
      expected = BFplus0$bf,
      tolerance = 0.01
    )
    expect_equal(
      BFplus0.curr$bf,
      expected = BFplus0.BayesFactor,
      tolerance = 0.03
    )
  }

})
