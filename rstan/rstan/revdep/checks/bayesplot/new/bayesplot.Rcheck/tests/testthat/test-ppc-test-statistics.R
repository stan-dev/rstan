library(bayesplot)
context("PPC: test-statistics")

source(test_path("data-for-ppc-tests.R"))

q25 <- function(x) quantile(x, 0.25)
prop0 <- function(x) mean(x == 0)

test_that("ppc_stat throws errors if function not found", {
  expect_error(ppc_stat(y, yrep, stat = "9999"), "not found")
  expect_error(ppc_stat_grouped(y, yrep, group, stat = "9999"), "not found")
  expect_error(ppc_stat_freqpoly_grouped(y, yrep, group, stat = "9999"), "not found")
})

test_that("ppc_stat throws errors if 'stat' wrong length", {
  expect_error(ppc_stat(y, yrep, stat = c("mean", "sd")), "not a function")
  expect_error(ppc_stat_grouped(y, yrep, group, stat = c("mean", "sd")), "not a function")
  expect_error(ppc_stat_freqpoly_grouped(y, yrep, group, stat = c(mean, sd)), "not a function")
})

test_that("ppc_stat returns ggplot object", {
  expect_gg(ppc_stat(y, yrep))
  expect_gg(ppc_stat(y, yrep, stat = "sd"))
  expect_gg(ppc_stat(y, yrep, stat = sd))
  expect_gg(ppc_stat(y, yrep, stat = "q25"))
  expect_gg(ppc_stat(y, yrep, stat = q25))
  expect_gg(ppc_stat(y, yrep, stat = function(x) median(x)))

  expect_gg(ppc_stat(y2, yrep2))
  expect_gg(ppc_stat(y2, yrep2, stat = "prop0"))
})

test_that("ppc_stat_2d returns ggplot object", {
  expect_gg(ppc_stat_2d(y, yrep))
  expect_gg(ppc_stat_2d(y, yrep, stat = c("q25", "median")))
  expect_gg(ppc_stat_2d(y, yrep, stat = c("q25", median)))
  expect_gg(ppc_stat_2d(y, yrep, stat = c(function(x) mean(x), function(y) sd(y))))

  expect_gg(ppc_stat_2d(y2, yrep2))
})

test_that("ppc_stat_grouped returns ggplot object", {
  expect_gg(ppc_stat_grouped(y, yrep, group))
  expect_gg(ppc_stat_grouped(y, yrep, as.numeric(group), stat = function(z) var(z)))
  expect_gg(ppc_stat_grouped(y, yrep, as.integer(group), stat = "sd"))

  expect_error(ppc_stat_grouped(y2, yrep2, group2),
               "'group' must have more than one unique value")
})
test_that("ppc_stat_freqpoly_grouped returns ggplot object", {
  expect_gg(ppc_stat_freqpoly_grouped(y, yrep, group, stat = "sd", freq = FALSE))
  expect_gg(ppc_stat_freqpoly_grouped(y, yrep, group, stat = function(x) sd(x), freq = TRUE))
})
