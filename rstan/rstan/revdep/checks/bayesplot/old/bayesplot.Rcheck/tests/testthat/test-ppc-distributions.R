library(bayesplot)
context("PPC: distributions")

source(test_path("data-for-ppc-tests.R"))

test_that("ppc_dens_overlay returns a ggplot object", {
  expect_gg(ppc_dens_overlay(y, yrep))
  expect_gg(ppc_dens_overlay(y2, yrep2, size = 0.5, alpha = 0.2))
})

test_that("ppc_ecdf_overlay returns a ggplot object", {
  expect_gg(ppc_ecdf_overlay(y, yrep, size = 0.5, alpha = 0.2))
  expect_gg(ppc_ecdf_overlay(y2, yrep2))
})

test_that("ppc_dens,pp_hist,ppc_freqpoly,ppc_boxplot return ggplot objects", {
  expect_gg(ppc_hist(y, yrep[1,, drop = FALSE]))
  expect_gg(ppc_hist(y, yrep[1:8, ]))
  expect_gg(ppc_hist(y2, yrep2))

  expect_gg(ppc_boxplot(y, yrep[1,, drop = FALSE]))
  expect_gg(ppc_boxplot(y, yrep[1:8, ]))
  expect_gg(ppc_boxplot(y2, yrep2, notch = FALSE))

  expect_gg(ppc_dens(y, yrep[1:8, ]))
  expect_gg(ppc_dens(y2, yrep2))

  expect_gg(ppc_freqpoly(y, yrep[1:8, ], binwidth = 2, size = 2, alpha = 0.1))
  expect_gg(ppc_freqpoly(y2, yrep2))

  expect_gg(p <- ppc_hist(y, yrep[1:8, ], binwidth = 3))
  if (utils::packageVersion("ggplot2") >= "3.0.0") {
    facet_var <- "~rep_label"
    expect_equal(as.character(p$facet$params$facets[1]), facet_var)
  }
})

test_that("ppc_freqpoly_grouped returns a ggplot object", {
  expect_gg(ppc_freqpoly_grouped(y, yrep[1:4, ], group))
  expect_gg(ppc_freqpoly_grouped(y, yrep[1:4, ], group,
                                 freq = TRUE, alpha = 0.5))

  expect_error(ppc_freqpoly_grouped(y2, yrep2, group2),
               "'group' must have more than one unique value")
})

test_that("ppc_violin_grouped returns a ggplot object", {
  expect_gg(ppc_violin_grouped(y, yrep, group))
  expect_gg(ppc_violin_grouped(y, yrep, as.numeric(group)))
  expect_gg(ppc_violin_grouped(y, yrep, as.integer(group)))
  expect_gg(ppc_violin_grouped(y, yrep, group, y_draw = "both", y_jitter = 0.3))

  expect_error(ppc_violin_grouped(y2, yrep2, group2),
               "'group' must have more than one unique value")
})




# Visual tests -----------------------------------------------------------------

test_that("ppc_hist renders correctly", {
  testthat::skip_on_cran()

  p_base <- ppc_hist(vdiff_y, vdiff_yrep[1:8, ])
  vdiffr::expect_doppelganger("ppc hist (default)", p_base)

  p_binwidth <- ppc_hist(vdiff_y, vdiff_yrep[1:8, ], binwidth = 3)
  vdiffr::expect_doppelganger("ppc hist (binwidth)", p_binwidth)
})

test_that("ppc_freqpoly renders correctly", {
  testthat::skip_on_cran()

  p_base <- ppc_freqpoly(vdiff_y, vdiff_yrep[1:8, ])
  vdiffr::expect_doppelganger("ppc freqpoly (default)", p_base)

  p_custom <- ppc_freqpoly(
    y = vdiff_y,
    yrep = vdiff_yrep[1:8, ],
    binwidth = 2,
    size = 2,
    alpha = 0.1)

  vdiffr::expect_doppelganger(
    title = "ppc freqpoly (alpha, binwidth, size)",
    fig = p_custom)
})

test_that("ppc_freqpoly_grouped renders correctly", {
  testthat::skip_on_cran()

  p_base <- ppc_freqpoly_grouped(vdiff_y, vdiff_yrep[1:3, ], vdiff_group)
  vdiffr::expect_doppelganger("ppc freq poly grouped (default)", p_base)
})

test_that("ppc_boxplot renders correctly", {
  testthat::skip_on_cran()

  p_base <- ppc_boxplot(vdiff_y, vdiff_yrep[1:8, ])
  vdiffr::expect_doppelganger("ppc boxplot (default)", p_base)

  p_no_notch <- ppc_boxplot(vdiff_y, vdiff_yrep[1:8, ], notch = FALSE)
  vdiffr::expect_doppelganger("ppc boxplot (no notch)", p_no_notch)

  p_custom <- ppc_boxplot(vdiff_y, vdiff_yrep[1:8, ], size = 1.5, alpha = .5)
  vdiffr::expect_doppelganger("ppc boxplot (alpha, size)", p_custom)
})

test_that("ppc_ecdf_overlay renders correctly", {
  testthat::skip_on_cran()

  p_base <- ppc_ecdf_overlay(vdiff_y2, vdiff_yrep2)
  vdiffr::expect_doppelganger("ppc ecdf overlay (default)", p_base)

  p_custom <- ppc_ecdf_overlay(vdiff_y2, vdiff_yrep2, size = 2, alpha = .2)
  vdiffr::expect_doppelganger("ppc ecdf overlay (alpha, size)", p_custom)
})

test_that("ppc_dens renders correctly", {
  testthat::skip_on_cran()

  p_base <- ppc_dens(vdiff_y, vdiff_yrep[1:8, ])
  vdiffr::expect_doppelganger("ppc density (default)", p_base)
})

test_that("ppc_dens_overlay renders correctly", {
  testthat::skip_on_cran()

  p_base <- ppc_dens_overlay(vdiff_y, vdiff_yrep)
  vdiffr::expect_doppelganger("ppc density overlay (default)", p_base)

  p_custom <- ppc_dens_overlay(vdiff_y, vdiff_yrep, size = 1, alpha = 0.2)
  vdiffr::expect_doppelganger("ppc density overlay (alpha, size)", p_custom)
})

test_that("ppc_violin_grouped renders correctly", {
  testthat::skip_on_cran()

  p_base <- ppc_violin_grouped(vdiff_y, vdiff_yrep, vdiff_group)
  vdiffr::expect_doppelganger("ppc violin grouped (default)", p_base)

  # lock in jitter
  set.seed(100)
  p_dots <- ppc_violin_grouped(
    y = vdiff_y,
    yrep = vdiff_yrep,
    group = vdiff_group,
    y_draw = "both")

  vdiffr::expect_doppelganger("ppc violin grouped (with points)", p_dots)

  p_dots_jitter <- ppc_violin_grouped(
    y = vdiff_y,
    yrep = vdiff_yrep,
    group = vdiff_group,
    y_draw = "points",
    y_jitter = 0.01)

  vdiffr::expect_doppelganger(
    "ppc violin grouped (points, low jitter)",
    p_dots_jitter)

  set.seed(seed = NULL)
})
