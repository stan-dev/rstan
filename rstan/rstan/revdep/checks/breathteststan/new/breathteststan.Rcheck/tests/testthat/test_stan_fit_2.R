context("Bayesian fit with Stan")

test_that("Data that cannot be fitted with nls_list/nlme work with stan_fit", {
  # with this seed, cf[10] does not fit with nls_list
  library(breathtestcore)

#  library(rstan)
#  library(dplyr)
#  library(rstan)
#  library(stringr)
#  library(testthat)
#  library(breathteststan)

  chains = 1
  student_t_df = 10
  dose = 100
  iter = 100
  sample_minutes = 15
  data = cleanup_data(simulate_breathtest_data(seed = 100)$data)
  comment(data) = "comment"
  fit = stan_fit(data, dose = dose, student_t_df = student_t_df,
                 chains = chains, iter = iter  )
  expect_is(fit, "breathtestfit")
  expect_is(fit, "breathteststanfit")
  expect_is(fit$stan_fit, "stanfit" )
  expect_identical(names(fit), c("coef", "data", "stan_fit", "coef_chain"))
  expect_equal(names(fit$data), names(data))
  expect_gt(sigma(fit), 0.9)
  expect_identical(comment(fit), "comment")

  cf = fit$coef
  expect_identical(unique(cf$group), "A")
  expect_identical(unique(cf$parameter), c("beta", "k", "m", "t50", "tlag"))
  expect_identical(unique(cf$stat), c("estimate", "q_0275", "q_25", "q_75", "q_975"))
  expect_equal(nrow(cf), 395)
  expect_equal(ncol(cf), 6)

  cf = coef(fit) # This is the "mean" group only
  expect_identical(unique(cf$group), "A")
  expect_identical(unique(cf$parameter), c("beta", "k", "m", "t50", "tlag"))
  expect_equal(nrow(cf), 79)
  expect_equal(ncol(cf), 5)
})

