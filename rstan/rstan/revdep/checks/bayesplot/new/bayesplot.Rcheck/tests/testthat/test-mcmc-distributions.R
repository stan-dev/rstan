library(bayesplot)
context("MCMC: distributions")

source(test_path("data-for-mcmc-tests.R"))

get_palette <- function(ggplot, n) {
  scale <- ggplot$scales$get_scales("colour")
  scale$palette(n)
}

test_that("mcmc_hist returns a ggplot object", {
  expect_gg(mcmc_hist(arr, pars = "beta[1]", regex_pars = "x\\:"))
  expect_gg(mcmc_hist(arr1chain, regex_pars = "beta"))
  expect_gg(mcmc_hist(mat))
  expect_gg(mcmc_hist(dframe))
  expect_gg(mcmc_hist(dframe_multiple_chains))

  expect_gg(mcmc_hist(arr1))
  expect_gg(mcmc_hist(mat1))
  expect_gg(mcmc_hist(dframe1))
})

test_that("mcmc_dens returns a ggplot object", {
  expect_gg(mcmc_dens(arr, pars = "beta[2]", regex_pars = "x\\:"))
  expect_gg(mcmc_dens(arr1chain, regex_pars = "beta"))
  expect_gg(mcmc_dens(mat))

  expect_gg(mcmc_dens(dframe, transformations = list(sigma = function(x) x^2)))
  expect_gg(mcmc_dens(
    dframe_multiple_chains,
    transformations =
      list(sigma = function(x) x ^ 2, 'beta[1]' = "exp")
  ))

  expect_gg(mcmc_dens(arr1))
  expect_gg(mcmc_dens(mat1))
  expect_gg(mcmc_dens(dframe1))
})


# functions that require multiple chains ----------------------------------
test_that("mcmc_hist_by_chain returns a ggplot object", {
  expect_gg(mcmc_hist_by_chain(arr, pars = "beta[1]", regex_pars = "x\\:"))
  expect_gg(mcmc_hist_by_chain(dframe_multiple_chains,
                               regex_pars = c("(Intercept)", "beta")))
})

test_that("mcmc_dens_overlay returns a ggplot object", {
  expect_gg(mcmc_dens_overlay(arr, pars = "beta[1]", regex_pars = "x\\:"))
  expect_gg(mcmc_dens_overlay(dframe_multiple_chains,
                              pars = c("(Intercept)", "beta[2]")))
})

test_that("mcmc_dens_chains returns a ggplot object", {
  p <- mcmc_dens_chains(arr, pars = "beta[1]", regex_pars = "x\\:",
                        color_chains = FALSE)
  expect_gg(p)

  p2 <- mcmc_dens_overlay(dframe_multiple_chains,
                          pars = c("(Intercept)", "beta[2]"),
                          color_chains = TRUE)
  expect_gg(p2)
})

test_that("mcmc_dens_chains/mcmc_dens_overlay color chains", {
  p1 <- mcmc_dens_chains(arr, pars = "beta[1]", regex_pars = "x\\:",
                         color_chains = FALSE)
  p2 <- mcmc_dens_overlay(arr, pars = "beta[1]", regex_pars = "x\\:",
                          color_chains = FALSE)

  # Only one color when set not to color chains
  expect_equal(length(unique(get_palette(p1, 4))), 1)
  expect_equal(length(unique(get_palette(p2, 4))), 1)

  p3 <- mcmc_dens_chains(arr, pars = "beta[1]", regex_pars = "x\\:",
                         color_chains = TRUE)
  p4 <- mcmc_dens_overlay(arr, pars = "beta[1]", regex_pars = "x\\:",
                          color_chains = TRUE)

  # Chain coloring works
  expect_equal(get_palette(p3, 4), chain_colors(4))
  expect_equal(get_palette(p4, 4), chain_colors(4))
})

test_that("mcmc_violin returns a ggplot object", {
  expect_gg(mcmc_violin(arr, pars = "beta[2]", regex_pars = "x\\:"))
  expect_gg(mcmc_violin(dframe_multiple_chains,
                        regex_pars = c("\\(Intercept\\)$", "beta")))
})

test_that("mcmc_* throws error if 1 chain but multiple chains required", {
  expect_error(mcmc_hist_by_chain(mat), "requires multiple chains")
  expect_error(mcmc_hist_by_chain(dframe), "requires multiple chains")
  expect_error(mcmc_hist_by_chain(arr1chain), "requires multiple chains")

  expect_error(mcmc_dens_overlay(mat), "requires multiple chains")
  expect_error(mcmc_dens_overlay(dframe), "requires multiple chains")
  expect_error(mcmc_dens_overlay(arr1chain), "requires multiple chains")

  expect_error(mcmc_dens_chains(mat), "requires multiple chains")
  expect_error(mcmc_dens_chains(dframe), "requires multiple chains")
  expect_error(mcmc_dens_chains(arr1chain), "requires multiple chains")

  expect_error(mcmc_violin(mat), "requires multiple chains")
  expect_error(mcmc_violin(dframe), "requires multiple chains")
  expect_error(mcmc_violin(arr1chain), "requires multiple chains")
})
