library(bayesplot)
context("PPC: intervals & ribbon")

source(test_path("data-for-ppc-tests.R"))

test_that("ppc_intervals returns ggplot object", {
  expect_gg(ppc_intervals(y, yrep))
  expect_gg(ppc_intervals(y, yrep, size = 2, fatten = 1))
  expect_gg(ppc_intervals(y, yrep, x = seq(1, 2 * length(y), by = 2)))
  expect_gg(ppc_intervals(y2, yrep2))
})

test_that("ppc_ribbon returns ggplot object", {
  expect_gg(ppc_ribbon(y, yrep, prob = 0.5))
  expect_gg(ppc_ribbon(y, yrep, alpha = 0, size = .5))
  expect_gg(ppc_ribbon(y2, yrep2, x = rnorm(length(y2)), prob = 0.5))
})


y <- rnorm(50)
yrep <- matrix(rnorm(500, 0, 2), ncol = 50)
x <- rep(1:10, each = 5)
group <- gl(5, 1, length = 50, labels = LETTERS[1:5])

test_that("ppc_intervals_grouped returns ggplot object", {
  expect_gg(ppc_intervals_grouped(y, yrep, x, group))
})

test_that("ppc_ribbon_grouped returns ggplot object", {
  expect_gg(
    ppc_ribbon_grouped(y, yrep, x, group, facet_args = list(scales = "fixed")))
})

test_that("ppc_intervals_data returns correct structure", {
  d <- ppc_intervals_data(y, yrep, x = 1:length(y), prob = .9)
  d_group <- ppc_intervals_data(y, yrep, x, group)
  expect_named(d, c("y_id", "y_obs", "x",
                    "outer_width", "inner_width",
                    "ll", "l", "m", "h", "hh"))
  expect_named(d_group, c("y_id", "y_obs", "group", "x",
                          "outer_width", "inner_width",
                          "ll", "l", "m", "h", "hh"))

  expect_error(
    ppc_intervals_data(y, yrep, x = 1:length(y), prob = 0), "prob")

  expect_error(
    ppc_intervals_data(y, yrep, x = 1:length(y), prob = 1.01), "prob")

  expect_error(
    ppc_intervals_data(y, yrep, x = 1:length(y), prob_outer = 0), "prob_outer")

  expect_error(
    ppc_intervals_data(y, yrep, x = 1:length(y), prob_outer = 1.01), "prob_outer")
})

test_that("ppc_intervals_data does math correctly", {
  d <- ppc_intervals_data(y, yrep, prob = .4, prob_outer = .8)
  qs <- unname(quantile(yrep[, 1], c(.1, .3, .5, .7, .9)))
  expect_equal(d$ll[1], qs[1])
  expect_equal(d$l[1],  qs[2])
  expect_equal(d$m[1],  qs[3])
  expect_equal(d$h[1],  qs[4])
  expect_equal(d$hh[1], qs[5])

  # Testing groups and known quantiles
  y <- rep(10, 4)
  group <- c("a", "a", "b", "b")
  yrep_g1 <- matrix(rep((0:20), 2), ncol = 2)
  yrep_g2 <- yrep_g1 - 10
  yrep <- cbind(yrep_g1, yrep_g2)

  d <- ppc_intervals_data(y, yrep, group = group, prob = .5, prob_outer = .9)
  expect_equal(unique(d$inner_width), .5)
  expect_equal(unique(d$outer_width), .9)
  expect_equal(d$ll, c( 1,  1, -9, -9))
  expect_equal(d$l,  c( 5,  5, -5, -5))
  expect_equal(d$m,  c(10, 10,  0,  0))
  expect_equal(d$h,  c(15, 15,  5,  5))
  expect_equal(d$hh, c(19, 19,  9,  9))
})




# Visual tests -----------------------------------------------------------------

test_that("ppc_intervals renders correctly", {
  testthat::skip_on_cran()

  p_base <- ppc_intervals(vdiff_y, vdiff_yrep)
  vdiffr::expect_doppelganger("ppc intervals (default)", p_base)

  p_x <- ppc_intervals(vdiff_y, vdiff_yrep, x = vdiff_y)
  vdiffr::expect_doppelganger("ppc intervals (x values)", p_x)

  p_50 <- ppc_intervals(vdiff_y, vdiff_yrep, prob = .50)
  vdiffr::expect_doppelganger("ppc intervals (interval width)", p_50)
})

test_that("ppc_intervals_grouped renders correctly", {
  testthat::skip_on_cran()

  p_base <- ppc_intervals_grouped(vdiff_y, vdiff_yrep, group = vdiff_group)
  vdiffr::expect_doppelganger("ppc intervals grouped (default)", p_base)

  p_x <- ppc_intervals_grouped(
    y = vdiff_y,
    yrep = vdiff_yrep,
    x = vdiff_y,
    group = vdiff_group)
  vdiffr::expect_doppelganger("ppc intervals grouped (x values)", p_x)
})

test_that("ppc_ribbon renders correctly", {
  testthat::skip_on_cran()

  p_base <- ppc_ribbon(vdiff_y, vdiff_yrep)
  vdiffr::expect_doppelganger("ppc ribbon (default)", p_base)

  p_x <- ppc_ribbon(vdiff_y, vdiff_yrep, x = vdiff_y)
  vdiffr::expect_doppelganger("ppc ribbon (x values)", p_x)

  p_50 <- ppc_ribbon(vdiff_y, vdiff_yrep, prob = 0.5)
  vdiffr::expect_doppelganger("ppc ribbon (interval width)", p_50)
})
