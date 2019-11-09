library(bayesplot)
suppressPackageStartupMessages(library(rstanarm))
context("MCMC: nuts")

ITER <- 1000
CHAINS <- 3
capture.output(
  fit <- stan_glm(mpg ~ wt + am, data = mtcars,
                  iter = ITER, chains = CHAINS, refresh = 0)
)
np <- nuts_params(fit)
lp <- log_posterior(fit)

test_that("all mcmc_nuts_* (except energy) return gtable objects", {
  expect_gtable(mcmc_nuts_acceptance(np, lp))
  expect_gtable(mcmc_nuts_acceptance(np, lp, chain = CHAINS))

  expect_gtable(mcmc_nuts_treedepth(np, lp))
  expect_gtable(mcmc_nuts_treedepth(np, lp, chain = CHAINS))

  expect_gtable(mcmc_nuts_stepsize(np, lp))
  expect_gtable(mcmc_nuts_stepsize(np, lp, chain = CHAINS))

  expect_gtable(mcmc_nuts_divergence(np, lp))
  expect_gtable(mcmc_nuts_divergence(np, lp, chain = CHAINS))
})
test_that("all mcmc_nuts_* (except energy) error if chain argument is bad", {
  funs <- c("acceptance", "divergence", "treedepth", "stepsize")
  for (f in paste0("mcmc_nuts_", funs)) {
    expect_error(do.call(f, list(x=np, lp=lp, chain = CHAINS + 1)),
                 regexp = paste("only", CHAINS, "chains found"),
                 info = f)
    expect_error(do.call(f, list(x=np, lp=lp, chain = 0)),
                 regexp = "chain >= 1",
                 info = f)
  }
})

test_that("mcmc_nuts_energy returns a ggplot object", {
  p <- mcmc_nuts_energy(np)
  expect_gg(p)
  expect_s3_class(p$facet, "FacetWrap")
  expect_equal(names(p$facet$params$facets), "Chain")

  p <- mcmc_nuts_energy(np, merge_chains = TRUE)
  expect_gg(p)
  expect_s3_class(p$facet, "FacetNull")
})
test_that("mcmc_nuts_energy throws correct warnings", {
  expect_warning(mcmc_nuts_energy(np, chain = 1), "ignored: chain")
})


test_that("validate_nuts_data_frame throws errors", {
  expect_error(
    validate_nuts_data_frame(list(Iteration = 1, Chain = 1)),
    "NUTS parameters should be in a data frame"
  )
  expect_error(
    validate_nuts_data_frame(data.frame(Iteration = 1, apple = 2)),
    "NUTS parameter data frame must have columns: Iteration, Parameter, Value, Chain"
  )
  expect_error(
    validate_nuts_data_frame(np, as.matrix(lp)),
    "lp should be in a data frame"
  )

  lp2 <- lp
  colnames(lp2)[3] <- "Chains"
  expect_error(
    validate_nuts_data_frame(np, lp2),
    "lp data frame must have columns: Iteration, Value, Chain"
  )

  lp2 <- subset(lp, Chain %in% 1:2)
  expect_error(
    validate_nuts_data_frame(np, lp2),
    "Number of chains"
  )
})
