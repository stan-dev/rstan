library(bayesplot)
context("available_mcmc and available_ppc")


test_that("available_mcmc works", {
  a <- available_mcmc()
  expect_s3_class(a, "bayesplot_function_list")
  expect_s3_class(a, "character")
  expect_identical(
    as.character(a),
    sort(grep("^mcmc_", getNamespaceExports("bayesplot"), value = TRUE))
  )

  b <- available_mcmc("trace|dens")
  expect_s3_class(b, "bayesplot_function_list")
  expect_identical(
    as.character(b),
    sort(grep("^mcmc_dens|^mcmc_trace", getNamespaceExports("bayesplot"), value = TRUE))
  )

  expect_length(available_mcmc(pattern = "99999"), 0)
})

test_that("available_ppc works", {
  a <- available_ppc()
  expect_s3_class(a, "bayesplot_function_list")
  expect_s3_class(a, "character")
  expect_identical(
    as.character(a),
    sort(grep("^ppc_", getNamespaceExports("bayesplot"), value = TRUE))
  )

  b <- available_ppc("grouped")
  expect_s3_class(b, "bayesplot_function_list")
  expect_identical(
    as.character(b),
    sort(grep("_grouped$", getNamespaceExports("bayesplot"), value = TRUE))
  )

  c <- available_ppc("grouped", invert = TRUE)
  expect_false(any(grepl("grouped", c)))

  expect_length(available_ppc(pattern = "99999"), 0)
})

test_that("print.bayesplot_function_list works", {
  expect_output(print(available_ppc()), "bayesplot PPC module:")
  expect_output(print(available_mcmc()), "bayesplot MCMC module:")

  expect_output(print(available_ppc("ribbon")), "(matching pattern 'ribbon')")
  expect_output(print(available_mcmc("trace")), "trace_highlight")

  expect_output(print(available_ppc("grouped", invert = TRUE)),
                "excluding pattern 'grouped'")
})
