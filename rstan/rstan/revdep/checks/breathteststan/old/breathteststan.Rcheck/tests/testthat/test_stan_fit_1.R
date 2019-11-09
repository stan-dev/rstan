context("Bayesian fit with Stan")

test_that("stanmodels exist", {
  expect_is(breathteststan:::stanmodels,"list")
  expect_is(breathteststan:::stanmodels$breath_test_1, "stanmodel")
  expect_gte(length(breathteststan:::stanmodels), 1L)
})

