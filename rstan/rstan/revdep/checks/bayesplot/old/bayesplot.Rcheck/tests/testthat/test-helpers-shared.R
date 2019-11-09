library(bayesplot)
context("Shared: misc. functions")

# suggested packages ------------------------------------------------------
test_that("suggested_package throws correct errors", {
  expect_error(suggested_package("NOPACKAGE"),
               "Please install the NOPACKAGE package")
  expect_error(suggested_package(c("testthat", "gridExtra")), "length")
  expect_silent(suggested_package("testthat"))
  expect_silent(suggested_package("testthat", min_version = "0.0.1"))
  expect_error(suggested_package("testthat", min_version = "100000.0.0"))
})


# check_ignored_arguments -------------------------------------------------
test_that("check_ignored_arguments throws correct warnings", {
  expect_warning(check_ignored_arguments(a = 1, b = "2"),
               "The following arguments were unrecognized and ignored: a, b")
  expect_warning(check_ignored_arguments(a = 1, b = "2", ok_args = c("a", "c")),
                 "The following arguments were unrecognized and ignored: b")
  expect_silent(check_ignored_arguments(a = 1, b = "2", ok_args = c("a", "b", "c")))
})


# parameter selection -----------------------------------------------------
all_pars <- c("param_1", "param_2",
              "param[1]", "param[2]",
              "param[1,3,5]", "param[2,4,5]",
              "alpha", "beta")

test_that("select_parameters throws errors if 'explicit' not found", {
  expect_error(select_parameters(explicit = c("alpha", "ALPHA"),
                                 complete = all_pars),
               "don't match parameter names: ALPHA")
  expect_error(select_parameters(c("BETA", "ALPHA"), complete = all_pars),
               "don't match parameter names: BETA, ALPHA")
})
test_that("select_parameters throws errors if no regex matches", {
  expect_error(select_parameters(explicit = c("alpha", "beta"),
                                 patterns = "tomato|apple",
                                 complete = all_pars),
               "No matches for 'regex_pars'")
})

test_that("select_parameters works with regex", {
  expect_identical(select_parameters(patterns = "param", complete = all_pars),
                   all_pars[-c(7:8)])

  expect_identical(select_parameters(patterns = c("param", "tomato"), complete = all_pars),
                   all_pars[-c(7:8)])

  expect_identical(select_parameters(patterns = c("param\\[", "tomato"), complete = all_pars),
                   all_pars[3:6])

  expect_identical(select_parameters(patterns = c("param\\_"), complete = all_pars),
                   all_pars[1:2])
})

test_that("select_parameters works without regex", {
  expect_identical(select_parameters(explicit = "alpha", complete = all_pars),
                   "alpha")
  expect_identical(select_parameters(c("alpha", "param[1,3,5]"), complete = all_pars),
                   c("alpha", "param[1,3,5]"))
})

test_that("select_parameters works with both explicit and regex", {
  expect_identical(select_parameters(explicit = "alpha",
                                     patterns = "param",
                                     complete = all_pars),
                   c("alpha", all_pars[-c(7:8)]))
  expect_identical(select_parameters(explicit = "alpha",
                                     patterns = "alpha",
                                     complete = all_pars),
                   "alpha")
  expect_identical(select_parameters(explicit = c("alpha", "beta"),
                                     patterns = "param\\[|param\\_",
                                     complete = all_pars),
                   c(all_pars[7:8], all_pars[-c(7:8)]))
})

