library(bayesplot)
suppressPackageStartupMessages(library(rstanarm))
context("MCMC: recover")

alpha <- 1; beta <- c(-.5, .5); sigma <- 2
X <- matrix(rnorm(200), 100, 2)
y <- rnorm(100, mean = c(alpha + X %*% beta), sd = sigma)
capture.output(
  fit <- stan_glm(y ~ ., data = data.frame(y, X))
)
draws <- as.matrix(fit)
true <- c(alpha, beta, sigma)

test_that("mcmc_recover_intervals throws correct errors", {
  expect_error(
    mcmc_recover_intervals(draws, letters[1:ncol(draws)]),
    "is.numeric(true) is not TRUE",
    fixed = TRUE
  )
  expect_error(
    mcmc_recover_intervals(draws, true[-1]),
    "ncol(x) == length(true) is not TRUE",
    fixed = TRUE
  )
  expect_error(
    mcmc_recover_intervals(draws, true, batch = 1:3),
    "length(batch) == length(true) is not TRUE",
    fixed = TRUE
  )
  expect_error(
    mcmc_recover_intervals(draws, true, prob = 0.8, prob_outer = 0.5),
    "prob_outer >= prob is not TRUE",
    fixed = TRUE
  )
  expect_error(
    mcmc_recover_intervals(draws, true, prob = 0, prob_outer = 0.5),
    "prob > 0 is not TRUE",
    fixed = TRUE
  )
  expect_error(
    mcmc_recover_intervals(draws, true, prob = .5, prob_outer = 1.1),
    "prob_outer <= 1 is not TRUE",
    fixed = TRUE
  )
})

test_that("mcmc_recover_intervals returns a ggplot object", {
  expect_gg(mcmc_recover_intervals(draws, true))
  expect_gg(mcmc_recover_intervals(draws, true, batch = c(1, 2, 2, 1),
                                   point_est = "mean"))
  expect_gg(mcmc_recover_intervals(draws, true, batch = grepl("X", colnames(draws))))
  expect_gg(mcmc_recover_intervals(draws, true, batch = grepl("X", colnames(draws)),
                                   facet_args = list(ncol = 1)))
})

test_that("mcmc_recover_intervals works when point_est = 'none'", {
  a <- mcmc_recover_intervals(draws, true, batch = 1:4, point_est = "none")
  expect_gg(a)
  expect_equal(a$data$Point, rep(NA, ncol(draws)))
})


test_that("mcmc_recover_scatter returns a ggplot object", {
  expect_gg(mcmc_recover_scatter(draws, true))
  expect_gg(mcmc_recover_scatter(draws, true, batch = 1:4,
                                 point_est = "mean"))
  expect_gg(mcmc_recover_scatter(draws, true, batch = c(1, 2, 2, 1),
                                 point_est = "mean"))
  expect_gg(mcmc_recover_scatter(draws, true, batch = grepl("X", colnames(draws))))
  expect_gg(mcmc_recover_scatter(draws, true, batch = grepl("X", colnames(draws)),
                                 facet_args = list(ncol = 1)))
})


test_that("mcmc_recover_hist returns a ggplot object", {
  expect_gg(mcmc_recover_hist(draws, true))
  expect_gg(mcmc_recover_hist(draws, true, binwidth = .1,
                              facet_args = list(nrow = 1)))
})
