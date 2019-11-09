library(bayesplot)
context("MCMC: intervals")

source(test_path("data-for-mcmc-tests.R"))


test_that("mcmc_intervals_data computes quantiles", {
  xs <- melt_mcmc(merge_chains(prepare_mcmc_array(arr, pars = "beta[1]")))
  d <- mcmc_intervals_data(arr, pars = "beta[1]",
                           prob = .3, prob_outer = .5)

  qs <- unlist(d[, c("ll", "l", "m", "h", "hh")])
  by_hand <- quantile(xs$Value, c(.25, .35, .5, .65, .75))
  expect_equivalent(qs, by_hand)

  expect_equal(d$parameter, factor("beta[1]"))
  expect_equal(d$outer_width, .5)
  expect_equal(d$inner_width, .3)
  expect_equal(d$point_est, "median")

  d2 <- mcmc_areas_data(arr, pars = "beta[1]", prob = .3, prob_outer = .5)
  sets <- split(d2, d2$interval)

  expect_equal(range(sets$inner$x), c(d$l, d$h))
  expect_equal(range(sets$outer$x), c(d$ll, d$hh))
})

test_that("mcmc_intervals_data computes point estimates", {
  xs <- melt_mcmc(merge_chains(prepare_mcmc_array(arr, pars = "beta[2]")))
  d <- mcmc_intervals_data(arr, pars = "beta[2]",
                           prob = .3, prob_outer = .5, point_est = "mean")

  expect_equivalent(d$m, mean(xs$Value))
  expect_equal(d$parameter, factor("beta[2]"))
  expect_equal(d$point_est, "mean")

  d <- mcmc_intervals_data(arr, pars = "(Intercept)",
                           prob = .3, prob_outer = .5,
                           point_est = "none")
  expect_true(!("m" %in% names(d)))
  expect_equal(d$point_est, "none")
})

test_that("mcmc_intervals returns a ggplot object", {
  expect_gg(mcmc_intervals(arr, pars = "beta[1]", regex_pars = "x\\:"))
  expect_gg(mcmc_intervals(arr1chain, pars = "beta[1]", regex_pars = "Intercept"))
  expect_gg(mcmc_intervals(mat, regex_pars = "beta"))
  expect_gg(mcmc_intervals(dframe))
  expect_gg(mcmc_intervals(dframe_multiple_chains))

  expect_gg(mcmc_intervals(arr1))
  expect_gg(mcmc_intervals(mat1))
  expect_gg(mcmc_intervals(dframe1))
})

test_that("mcmc_areas returns a ggplot object", {
  expect_gg(mcmc_areas(arr, pars = "beta[2]", regex_pars = "x\\:"))
  expect_gg(mcmc_areas(arr1chain, regex_pars = c("beta", "x\\:")))
  expect_gg(mcmc_areas(mat))
  expect_gg(mcmc_areas(dframe))
  expect_gg(mcmc_areas(dframe_multiple_chains))

  expect_gg(mcmc_areas(arr1))
  expect_gg(mcmc_areas(mat1))
  expect_gg(mcmc_areas(dframe1))
})

test_that("mcmc_areas_ridges returns a ggplot object", {
  expect_gg(mcmc_areas_ridges(arr, pars = "beta[2]", regex_pars = "x\\:"))
  expect_gg(mcmc_areas_ridges(arr1chain, regex_pars = c("beta", "x\\:")))
  expect_gg(mcmc_areas_ridges(mat))
  expect_gg(mcmc_areas_ridges(dframe))
  expect_gg(mcmc_areas_ridges(dframe_multiple_chains))

  expect_gg(mcmc_areas_ridges(arr1))
  expect_gg(mcmc_areas_ridges(mat1))
  expect_gg(mcmc_areas_ridges(dframe1))
})

test_that("mcmc_intervals/areas with rhat", {
  r <- runif(ncol(mat), 0.9, 1.3)
  rbad <- c(NA, r[-1])

  expect_error(mcmc_intervals(arr, rhat = r[-1]), "'rhat' has length")
  expect_error(expect_warning(mcmc_intervals(arr, rhat = rbad)))

  expect_gg(g <- mcmc_intervals(arr, rhat = r))
  if (utils::packageVersion("ggplot2") >= "3.0.0") {
    rhat_map <- g$layers[[3]][["mapping"]]
    expect_identical(as.character(rhat_map[["colour"]]), c("~", "rhat_rating"))
  }

  # areas with rhat.

  # layer 1 is maybe line vertical line. [skip]

  # layer 2 is inner interval.
  expect_gg(g2 <- mcmc_areas(arr, rhat = r))
  if (utils::packageVersion("ggplot2") >= "3.0.0") {
    rhat_map2 <- g2$layers[[2]][["mapping"]]
    expect_identical(as.character(rhat_map2$fill), c("~", "rhat_rating"))
    expect_identical(as.character(rhat_map2$colour), c("~", "rhat_rating"))
    # layer 3 is point estimate. manually colored. [skip]

    # layer 4 is outer interval.
    rhat_map4 <- g2$layers[[4]][["mapping"]]
    expect_identical(as.character(rhat_map4$colour), c("~", "rhat_rating"))

    # layer 5 is bottom line.

    rhat_map5 <- g2$layers[[5]][["mapping"]]
    expect_identical(as.character(rhat_map5$colour), c("~", "rhat_rating"))
  }
})

test_that("mcmc_areas_data computes density", {
  areas_data <- mcmc_areas_data(arr, point_est = "none")
  areas_data <- areas_data[areas_data$interval_width == 1, ]
  by_parameter <- split(areas_data, areas_data$parameter)

  # Manually compute the same
  raw_values <- melt_mcmc(merge_chains(prepare_mcmc_array(arr)))
  raw_values <- split(raw_values, interaction(raw_values$Parameter))

  do_dens <- function(df, interval_width, n) {
    x <- df$Value
    tail_width <- (1 - interval_width) / 2
    qs <- quantile(x, probs = c(tail_width, 1 - tail_width))
    dens <- density(x = x, from = min(qs), to = max(qs), n = n)
    data.frame(Parameter = unique(df$Parameter), x = dens$x, y = dens$y)
  }

  densities <- lapply(raw_values, do_dens, 1, 1024)

  for (name in names(by_parameter)) {
    expect_equivalent(by_parameter[[name]][["density"]],
                      densities[[name]][["y"]])
  }
})

test_that("compute_column_density can use density options (#118)", {
  # n_dens affects the number of rows in the return data-frame
  areas_data <- mcmc_areas_data(arr, point_est = "none", n_dens = 100)
  pars <- length(unique(areas_data$parameter))
  intervals <- length(unique(areas_data$interval))
  expect_equal(nrow(areas_data), 100 * intervals * pars)

  # If these raise errors, they are being evaluated
  expect_error(mcmc_areas_data(arr, bw = stop()))
  expect_error(mcmc_areas_data(arr, adjust = stop()))
  expect_error(mcmc_areas_data(arr, kernel = stop()))
})




# Visual tests -----------------------------------------------------------------

test_that("mcmc_intervals renders correctly", {
  testthat::skip_on_cran()

  p_base <- mcmc_intervals(vdiff_dframe)
  vdiffr::expect_doppelganger("mcmc intervals (default)", p_base)

  p_outer <- mcmc_intervals(vdiff_dframe, prob_outer = .8)
  vdiffr::expect_doppelganger("mcmc intervals (outer)", p_outer)

  p_inner <- mcmc_intervals(vdiff_dframe, prob = .1)
  vdiffr::expect_doppelganger("mcmc intervals (inner)", p_inner)

  rhats <- seq(from = 1, to = 1.15, length = 5)
  p_rhats <- mcmc_intervals(vdiff_dframe, rhat = rhats)
  vdiffr::expect_doppelganger("mcmc intervals (rhats)", p_rhats)

  p_no_points <- mcmc_intervals(vdiff_dframe, point_est = "none")
  vdiffr::expect_doppelganger("mcmc intervals (no points)", p_no_points)

  p_mean_points <- mcmc_intervals(vdiff_dframe, point_est = "mean")
  vdiffr::expect_doppelganger("mcmc intervals (means)", p_mean_points)
})

test_that("mcmc_areas renders correctly", {
  testthat::skip_on_cran()

  p_base <- mcmc_areas(vdiff_dframe)
  vdiffr::expect_doppelganger("mcmc areas (default)", p_base)

  p_outer <- mcmc_areas(vdiff_dframe, prob_outer = .8)
  vdiffr::expect_doppelganger("mcmc areas (outer)", p_outer)

  p_inner <- mcmc_areas(vdiff_dframe, prob = .1)
  vdiffr::expect_doppelganger("mcmc areas (inner)", p_inner)

  rhats <- seq(from = 1, to = 1.15, length = 5)
  p_rhats <- mcmc_areas(vdiff_dframe, rhat = rhats)
  vdiffr::expect_doppelganger("mcmc areas (rhats)", p_rhats)

  p_no_points <- mcmc_areas(vdiff_dframe, point_est = "none")
  vdiffr::expect_doppelganger("mcmc areas (no points)", p_no_points)

  p_mean_points <- mcmc_areas(vdiff_dframe, point_est = "mean")
  vdiffr::expect_doppelganger("mcmc areas (means)", p_mean_points)
})

test_that("mcmc_areas_ridges renders correctly", {
  testthat::skip_on_cran()

  p_base <- mcmc_areas_ridges(vdiff_dframe)
  vdiffr::expect_doppelganger("mcmc areas ridges (default)", p_base)

  p_inner <- mcmc_areas_ridges(vdiff_dframe, prob = .5)
  vdiffr::expect_doppelganger("mcmc areas ridges (inner)", p_inner)

  p_outer <- mcmc_areas_ridges(vdiff_dframe, prob = .5, prob_outer = .95)
  vdiffr::expect_doppelganger("mcmc areas ridges (outer)", p_outer)
})
