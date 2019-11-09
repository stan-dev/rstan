library(bayesplot)
suppressPackageStartupMessages(library(rstanarm))
context("PPC: discrete")


# bar plots ---------------------------------------------------------------
data("esoph", package = "datasets")
capture.output(
  fit <- stan_polr(tobgp ~ agegp, data = esoph, method = "probit",
                   prior = R2(0.2, "mean"), init_r = 0.1, seed = 12345,
                   algorithm = "fullrank") # for speed only
)
y <- as.integer(fit$y)
yrep_char <- posterior_predict(fit, draws = 10)
yrep <- sapply(data.frame(yrep_char, stringsAsFactors = TRUE), as.integer)

test_that("ppc_bars & ppc_bars_grouped return a ggplot object", {
  expect_gg(ppc_bars(y, yrep))
  expect_gg(ppc_bars(y, yrep, prob = 0))
  expect_gg(ppc_bars_grouped(y, yrep, group = esoph$agegp))
})

test_that("freq argument to ppc_bars works", {
  p_freq <- ggplot2::ggplot_build(ppc_bars(y, yrep, freq = TRUE))
  p_prop <- ggplot2::ggplot_build(ppc_bars(y, yrep, freq = FALSE))

  y_freq <- p_freq$data[[1]]$y
  y_prop <- p_prop$data[[1]]$y
  expect_equal(y_freq, as.integer(y_freq))
  expect_true(all(y_prop < 1) && all(y_prop > 0))
})

test_that("ppc_bars errors if y/yrep not natural numbers", {
  expect_error(ppc_bars(y + 0.5, yrep),
               "ppc_bars expects only non-negative integers in 'y'")
  expect_error(ppc_bars(y, yrep + 0.5),
               "ppc_bars expects only non-negative integers in 'yrep'")
})


# rootograms -----------------------------------------------------------
rm(list = ls())
source(test_path("data-for-ppc-tests.R"))

yrep3 <- matrix(yrep2, nrow = 5, ncol = ncol(yrep2), byrow = TRUE)

test_that("ppc_rootogram returns a ggplot object", {
  expect_gg(ppc_rootogram(y2, yrep2))
  expect_gg(ppc_rootogram(y2, yrep3, style = "hanging", prob = 0.5))
  expect_gg(ppc_rootogram(y2, yrep3, style = "suspended"))
})

test_that("ppc_rootogram errors if y/yrep not counts", {
  expect_error(ppc_rootogram(y, yrep),
               "ppc_rootogram expects counts as inputs to 'y'")
  expect_error(ppc_rootogram(y2, yrep[1:5, seq_along(y2)]),
               "ppc_rootogram expects counts as inputs to 'yrep'")
  expect_error(ppc_rootogram(y, yrep3),
               "ncol(yrep) must be equal to length(y)", fixed = TRUE)
})

