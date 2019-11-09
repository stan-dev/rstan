library(bayesplot)
suppressPackageStartupMessages(library(rstanarm))
suppressPackageStartupMessages(library(loo))
context("PPC: loo")

options(useFancyQuotes = FALSE)

ITER <- 1000
CHAINS <- 3
capture.output(
  fit <- stan_glm(mpg ~ wt + am, data = mtcars,
                  iter = ITER, chains = CHAINS, refresh = 0)
)
y <- fit$y
yrep <- posterior_predict(fit)
suppressWarnings(
  psis1 <- psis(-log_lik(fit), cores = 2)
)
lw <- weights(psis1)
suppressWarnings(
  pits <- rstantools::loo_pit(yrep, y, lw)
)


test_that("ppc_loo_pit gives deprecation warning but still works", {
  expect_warning(p1 <- ppc_loo_pit(y, yrep, lw), "deprecated")
  expect_gg(p1)
})

test_that("ppc_loo_pit_overlay returns ggplot object", {
  expect_gg(p1 <- ppc_loo_pit_overlay(y, yrep, lw, samples = 25))
  expect_gg(p2 <- ppc_loo_pit_qq(y, yrep, lw, compare = "normal"))
  expect_equal(p2$labels$x, "Normal")
})

test_that("ppc_loo_pit_qq returns ggplot object", {
  expect_gg(p1 <- ppc_loo_pit_qq(y, yrep, lw))
  expect_equal(p1$labels$x, "Uniform")
})

test_that("ppc_loo_pit functions work when pit specified instead of y,yrep,lw", {
  expect_gg(ppc_loo_pit_qq(pit = pits))
  expect_message(
    ppc_loo_pit_qq(y = y, yrep = yrep, lw = lw, pit = pits),
    "'pit' specified so ignoring 'y','yrep','lw' if specified"
  )

  expect_gg(ppc_loo_pit_overlay(pit = pits))
  expect_message(
    ppc_loo_pit_overlay(y = y, yrep = yrep, lw = lw, pit = pits),
    "'pit' specified so ignoring 'y','yrep','lw' if specified"
  )
})



test_that("ppc_loo_intervals returns ggplot object", {
  expect_gg(ppc_loo_intervals(y, yrep, psis_object = psis1))
  expect_gg(g <- ppc_loo_intervals(y, yrep, psis_object = psis1, order = "median"))
  expect_s3_class(g$data$x, "factor")
  expect_equal(nlevels(g$data$x), length(g$data$x))

  # subset argument
  expect_gg(g <- ppc_loo_intervals(y, yrep, psis_object = psis1, subset = 1:25))
  expect_equal(nrow(g$data), 25)
})

test_that("ppc_loo_ribbon returns ggplot object", {
  expect_gg(ppc_loo_ribbon(y, yrep, psis_object = psis1, prob = 0.7, alpha = 0.1))
  expect_gg(g <- ppc_loo_ribbon(y, yrep, psis_object = psis1, subset = 1:25))
  expect_equal(nrow(g$data), 25)
})

test_that("ppc_loo_intervals/ribbon work when 'intervals' specified", {
  intervals <- t(apply(yrep, 2, quantile, probs = c(0.1, 0.25, 0.5, 0.75, 0.9)))
  expect_gg(ppc_loo_intervals(y, intervals = intervals))
  expect_gg(ppc_loo_ribbon(y, intervals = intervals))
  expect_message(ppc_loo_ribbon(y, intervals = intervals),
                 "'intervals' specified so ignoring 'yrep', 'psis_object', 'subset', if specified")
  expect_message(ppc_loo_intervals(y, yrep, psis_object = psis1, intervals = intervals),
                 "'intervals' specified so ignoring 'yrep', 'psis_object', 'subset', if specified")
})

test_that("ppc_loo_intervals/ribbon work when 'intervals' has 3 columns", {
  intervals <- t(apply(yrep, 2, quantile, probs = c(0.1, 0.5, 0.9)))
  expect_gg(ppc_loo_intervals(y, intervals = intervals))
  expect_gg(ppc_loo_ribbon(y, intervals = intervals))
})

test_that("errors if dimensions of yrep and lw don't match", {
  expect_error(
    ppc_loo_pit_overlay(y, yrep, lw[, 1:5]),
    "identical(dim(yrep), dim(lw)) is not TRUE",
    fixed = TRUE
  )
})

test_that("error if subset is bigger than num obs", {
  expect_error(.psis_subset(psis1, 1:1000), "too many elements")
  expect_error(
    ppc_loo_intervals(y, yrep, psis_object = psis1, subset = 1:1000),
    "length(y) >= length(subset) is not TRUE",
    fixed = TRUE
  )
})


