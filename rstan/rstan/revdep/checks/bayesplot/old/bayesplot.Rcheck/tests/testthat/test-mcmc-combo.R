library(bayesplot)
context("MCMC: combo")

source(test_path("data-for-mcmc-tests.R"))

test_that("mcmc_combo returns a gtable object", {
  expect_gtable(mcmc_combo(arr, regex_pars = "beta"))
  expect_gtable(mcmc_combo(mat, regex_pars = "beta",
                           binwidth = 1/20, combo = c("dens", "hist"),
                           facet_args = list(nrow = 2)))

  expect_gtable(mcmc_combo(dframe, regex_pars = "Intercept"))
  expect_gtable(mcmc_combo(dframe_multiple_chains, regex_pars = "Intercept",
                           combo = c("trace_highlight", "dens_overlay")))
  expect_gtable(mcmc_combo(arr1chain, regex_pars = "Intercept",
                           combo = c("trace", "hist")))

  expect_gtable(mcmc_combo(arr1, pars = "(Intercept)"))
  expect_gtable(mcmc_combo(mat1))
  expect_gtable(mcmc_combo(dframe1))
})

# functions that require multiple chains ----------------------------------
test_that("mcmc_combo throws error if 1 chain but multiple chains required", {
  expect_error(mcmc_combo(arr1chain, regex_pars = "beta",
               combo = c("trace_highlight", "dens")),
               "requires multiple chains")
  expect_error(mcmc_combo(mat, regex_pars = "beta",
               combo = c("trace_highlight", "hist")),
               "requires multiple chains")
  expect_error(mcmc_combo(dframe, regex_pars = "beta",
               combo = c("dens_overlay", "trace")),
               "requires multiple chains")
})

# other errors ------------------------------------------------------------
test_that("mcmc_combo throws errors", {
  expect_error(mcmc_combo(arr, combo = c("trace_highlight")),
               "'combo' should have at least two elements")
  expect_error(mcmc_combo(arr, regex_pars = "beta",
                          combo = c("animal", "hist", "tornado")),
               "The following functions were not found: mcmc_animal, mcmc_tornado")
})
