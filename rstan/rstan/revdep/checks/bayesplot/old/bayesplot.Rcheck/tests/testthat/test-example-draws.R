library(bayesplot)
context("Example draws")

test_that("example_mcmc_draws throws correct errors", {
  expect_error(example_mcmc_draws(chains = 5), "chains <= 4")
  expect_error(example_mcmc_draws(chains = 0), "chains >= 1")
  expect_error(example_mcmc_draws(params = 7), "params <= 6")
  expect_error(example_mcmc_draws(params = 0), "params >= 1")
})
test_that("example_mcmc_draws returns correct structure", {
  expect_identical(dim(example_mcmc_draws()),
                   c(250L, 4L, 4L))
  expect_identical(dim(example_mcmc_draws(chains = 1, params = 6)),
                   c(250L, 6L))
  expect_identical(dim(example_mcmc_draws(params = 1)),
                   c(250L, 4L, 1L))
  expect_identical(dimnames(example_mcmc_draws(4, 6))[[3]],
                   c("alpha", "sigma", paste0("beta[", 1:4,"]")))
})


test_that("example ppc data works", {
  y <- example_y_data()
  expect_type(y, "integer")
  expect_true(is_vector_or_1Darray(y))

  yrep <- example_yrep_draws()
  expect_type(yrep, "double")
  expect_is(yrep, "matrix")
  expect_equal(ncol(yrep), length(y))

  group <- example_group_data()
  expect_s3_class(group, "factor")
  expect_equal(length(group), length(y))

  x <- example_x_data()
  expect_type(x, "double")
  expect_true(is_vector_or_1Darray(x))
  expect_equal(length(x), length(y))
})
