library(bayesplot)
suppressPackageStartupMessages(library(rstanarm))
context("Extractors")

ITER <- 1000
CHAINS <- 3
capture.output(
  fit <- stan_glm(mpg ~ wt + am, data = mtcars,
                  iter = ITER, chains = CHAINS, refresh = 0)
)

x <- list(cbind(a = 1:3, b = rnorm(3)), cbind(a = 1:3, b = rnorm(3)))

# nuts_params and log_posterior methods -----------------------------------
test_that("nuts_params.list throws errors", {
  x[[3]] <- c(a = 1:3, b = rnorm(3))
  expect_error(nuts_params.list(x), "list elements should be matrices")

  x[[3]] <- cbind(a = 1:3, d = rnorm(3))
  expect_error(nuts_params.list(x), "same column names")

  x[[3]] <- cbind(a = 1:4, b = rnorm(4))
  expect_error(nuts_params.list(x), "same dimensions")
})

test_that("nuts_params.list parameter selection ok", {
  expect_error(nuts_params.list(x, pars = "apple"), "subscript out of bounds")

  np <- nuts_params.list(x, pars = "b")
  expect_true(all(np$Parameter == "b"))
})

test_that("all nuts_params methods identical", {
  expect_identical(
    nuts_params(fit),
    nuts_params(fit$stanfit)
  )
  expect_identical(
    nuts_params(fit),
    nuts_params(rstan::get_sampler_params(fit$stanfit, inc_warmup = FALSE))
  )
})

test_that("nuts_params.stanreg returns correct structure", {
  np <- nuts_params(fit)
  expect_identical(colnames(np), c("Iteration", "Parameter", "Value", "Chain"))

  np_names <- paste0(c("accept_stat", "stepsize", "treedepth", "n_leapfrog",
                       "divergent", "energy"), "__")
  expect_identical(levels(np$Parameter), np_names)

  expect_equal(length(unique(np$Iteration)), floor(ITER / 2))
  expect_equal(length(unique(np$Chain)), CHAINS)
})

test_that("log_posterior.stanreg returns correct structure", {
  lp <- log_posterior(fit)
  expect_identical(colnames(lp), c("Iteration", "Value", "Chain"))
  expect_equal(length(unique(lp$Iteration)), floor(ITER / 2))
  expect_equal(length(unique(lp$Chain)), CHAINS)
})

test_that("rhat.stanreg returns correct structure", {
  r <- rhat(fit)
  expect_named(r)
  expect_equal(r, summary(fit)[1:length(r), "Rhat"])

  expect_identical(names(rhat(fit, regex_pars = c("wt", "am"))),
                   c("wt", "am"))
})

test_that("neff_ratio.stanreg returns correct structure", {
  ratio <- neff_ratio(fit)
  expect_named(ratio)
  ans <- summary(fit)[1:length(ratio), "n_eff"] / (floor(ITER / 2) * CHAINS)
  expect_equal(ratio, ans, tol = 0.001)
})

test_that("rhat.stanfit returns correct structure", {
  r <- rhat(fit$stanfit)
  expect_named(r)
  expect_equal(r, summary(fit)[, "Rhat"])

  r2 <- rhat(fit$stanfit, pars = c("wt", "sigma"))
  expect_named(r2)
  expect_equal(r2, summary(fit, pars = c("wt", "sigma"))[, "Rhat"])
})

test_that("neff_ratio.stanreg returns correct structure", {
  denom <- floor(ITER / 2) * CHAINS

  ratio <- neff_ratio(fit$stanfit)
  expect_named(ratio)
  ans <- summary(fit)[, "n_eff"] / denom
  expect_equal(ratio, ans, tol = 0.001)

  ratio2 <- neff_ratio(fit$stanfit, pars = c("wt", "sigma"))
  expect_named(ratio2)
  ans2 <- summary(fit, pars = c("wt", "sigma"))[, "n_eff"] / denom
  expect_equal(ratio2, ans2, tol = 0.001)
})
