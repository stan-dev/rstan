
test_that("Getting/setting 0 or 1 options works", {
  rstan:::rstan_options(a = 22, b = 23)
  on.exit(rstan::rstan_options(a = NULL, b = NULL))

  o <- rstan:::rstan_options()
  expect_null(o)

  rstan:::rstan_options(testname = 22)
  expect_equal(rstan:::rstan_options("testname"), 22)
})

test_that("Getting/setting multiple options works", {
  rstan:::rstan_options(a = 22, b = 23)
  on.exit(rstan::rstan_options(a = NULL, b = NULL, c = NULL, d = NULL))

  o <- rstan:::rstan_options("a", "b")
  expect_equal(o$a, 22)
  expect_equal(o$b, 23)
  ov <- rstan:::rstan_options(a = 34)
  expect_equal(ov, 22)
  o <- rstan:::rstan_options("a", "b")
  expect_equal(o$a, 34)
  # Warning that "rstan option 'c' not found"
  # FIXME: Should the warning be tested?
  o <- suppressWarnings(rstan:::rstan_options("a", "b", "c"))
  expect_equal(o$c, NA)
  # Warning that "rstan option 'c' not found"
  # FIXME: Should the warning be tested?
  o <- suppressWarnings(rstan:::rstan_options("a", "b", "c", d = 38))
  expect_equal(o$d, NA)
  expect_equal(rstan:::rstan_options("d"), 38)
})

test_that("Option plot_rhat_breaks works", {
  rstan_options(plot_rhat_breaks = c(2, 1.2))
  o <- rstan_options("plot_rhat_breaks")
  expect_equal(o, c(1.2, 2))

  rstan_options(plot_rhat_breaks = c(1.5, 2, 1.2))
  o <- rstan_options("plot_rhat_breaks")
  expect_equal(o, c(1.2, 1.5, 2))

  on.exit(rstan_options(plot_rhat_break = NULL))
})

# All options used in rstan.
# Note: In the RUnit tests this test case had to checks, so I decided to skip
# it when converting to testthat.
test_that("All rstan options works", {
  skip("Need to add expectations")

  rhat_nan_col <- rstan_options("plot_rhat_nan_col")
  rhat_large_col <- rstan_options("plot_rhat_large_col")
  rhat_breaks <- rstan_options("plot_rhat_breaks")
  rhat_colors <- rstan_options("plot_rhat_cols")
  rhat_breaks <- rstan_options("plot_rhat_breaks")
  rhat_colors <- rstan_options("plot_rhat_cols")
  rhat_legend_cols <- c(
    rhat_colors, rstan_options("plot_rhat_large_col"),
    rstan_options("plot_rhat_nan_col")
  )
  alert_col <- rstan_options("rstan_alert_col")
  chain_cols <- rstan_options("rstan_chain_cols")
  standard_width <- rstan_options("plot_standard_npar")
  max_width <- rstan_options("plot_max_npar")
  rstan_options("eigen_lib")
  rstan_options("boost_lib")
  rstan_options("rstan_chain_cols")
  warmup_col <- rstan_options("rstan_warmup_bg_col")
})
