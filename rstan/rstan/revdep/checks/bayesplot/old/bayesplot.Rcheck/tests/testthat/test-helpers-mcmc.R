library(bayesplot)
context("MCMC: misc. functions")

source(test_path("data-for-mcmc-tests.R"))


# validate_mcmc_x ----------------------------------------------------------
test_that("validate_mcmc_x works", {
  expect_identical(validate_mcmc_x(mat), mat)
  expect_identical(validate_mcmc_x(mat1), mat1)
  expect_identical(validate_mcmc_x(arr), arr)
  expect_identical(validate_mcmc_x(arr1), arr1)
  expect_identical(validate_mcmc_x(arr1chain), arr1chain)

  # error if df_with_chain
  expect_error(validate_mcmc_x(dframe_multiple_chains), "is_df_with_chain")

  # converts regular df to matrix
  expect_identical(validate_mcmc_x(dframe), as.matrix(dframe))

  # NAs
  mat[1, 2] <- NA
  arr[1, 2, 3] <- NA
  expect_error(validate_mcmc_x(mat), "NAs not allowed")
  expect_error(validate_mcmc_x(arr), "NAs not allowed")
})


# 3-D array helpers --------------------------------------------------------
test_that("is_mcmc_array works", {
  expect_false(is_mcmc_array(mat))
  expect_false(is_mcmc_array(dframe))
  expect_false(is_mcmc_array(dframe_multiple_chains))
  expect_false(is_mcmc_array(arr))
  arr2 <- set_mcmc_dimnames(arr, parnames = dimnames(arr)[[3]])
  expect_mcmc_array(arr2)
})

test_that("parameter_names works", {
  x <- example_mcmc_draws()
  expect_identical(parameter_names(x), dimnames(x)[[3]])

  dimnames(x) <- list(a = NULL, b = NULL, c = letters[1:dim(x)[3]])
  expect_identical(parameter_names(x), dimnames(x)[[3]])

  dimnames(x) <- NULL
  expect_error(parameter_names(x), "No parameter names found")
  expect_error(parameter_names(x[, 1, ]), "No parameter names found")
})

test_that("has_multiple_chains works", {
  expect_error(has_multiple_chains(mat), "is_3d_array")
  expect_error(has_multiple_chains(dframe_multiple_chains), "is_3d_array")
  expect_error(has_multiple_chains(chainlist), "is_3d_array")

  expect_true(has_multiple_chains(arr))

  arr2 <- set_mcmc_dimnames(arr, parnames = dimnames(arr)[[3]])
  expect_true(has_multiple_chains(arr2))

  arr1chain2 <- set_mcmc_dimnames(arr1chain, parnames = dimnames(arr1chain)[[3]])
  expect_false(has_multiple_chains(arr1chain2))
})

test_that("has_multiple_params works", {
  expect_error(has_multiple_params(mat), "is_3d_array")
  expect_error(has_multiple_params(dframe_multiple_chains), "is_3d_array")

  expect_true(has_multiple_params(arr), "is_3d_array")

  arr2 <- set_mcmc_dimnames(arr, parnames = dimnames(arr)[[3]])
  expect_true(has_multiple_params(arr2))

  arr2 <- arr2[, , 3, drop = FALSE]
  expect_false(has_multiple_params(arr2))
})

# data frame with ‘chain’ variable ----------------------------------------
test_that("is_df_with_chain works", {
  expect_false(is_df_with_chain(arr))
  expect_false(is_df_with_chain(mat))
  expect_false(is_df_with_chain(chainlist))
  expect_false(is_df_with_chain(dframe))
  expect_true(is_df_with_chain(dframe_multiple_chains))

  mat2 <- cbind(mat, chain = dframe_multiple_chains$chain)
  expect_false(is_df_with_chain(mat2))

  dframe_multiple_chains2 <-
    cbind(dframe_multiple_chains, Chain = dframe_multiple_chains$chain)
  dframe_multiple_chains2$chain <- NULL
  expect_true(is_df_with_chain(dframe_multiple_chains2))
})

test_that("validate_df_with_chain works", {
  expect_error(validate_df_with_chain(mat), "is_df_with_chain")

  dframe_multiple_chains2 <-
    cbind(dframe_multiple_chains, Chain = dframe_multiple_chains$chain)
  dframe_multiple_chains2$chain <- NULL

  expect_identical(validate_df_with_chain(dframe_multiple_chains),
                   dframe_multiple_chains2)

  dframe_multiple_chains2$Chain <-
    factor(dframe_multiple_chains2$Chain, labels = letters[1:4])
  a <- validate_df_with_chain(dframe_multiple_chains2)
  expect_type(a$Chain, "integer")
})

test_that("df_with_chain2array works", {
  a <- df_with_chain2array(dframe_multiple_chains)
  expect_mcmc_array(a)

  expect_error(df_with_chain2array(dframe), "is_df_with_chain")
})



# list of chains ----------------------------------------------------------
test_that("is_chain_list works", {
  expect_false(is_chain_list(arr))
  expect_false(is_chain_list(mat))
  expect_false(is_chain_list(dframe))
  expect_false(is_chain_list(dframe_multiple_chains))
  expect_true(is_chain_list(chainlist))
  expect_true(is_chain_list(chainlist1))
  expect_true(is_chain_list(chainlist1chain))
})

test_that("validate_chain_list works", {
  expect_error(validate_chain_list(mat), "is_chain_list")
  expect_identical(validate_chain_list(chainlist), chainlist)
  expect_identical(validate_chain_list(chainlist1), chainlist1)
  expect_identical(validate_chain_list(chainlist1chain), chainlist1chain)

  chainlist2 <- chainlist
  colnames(chainlist2[[1]]) <- colnames(chainlist[[1]])
  colnames(chainlist2[[1]])[1] <- "AAA"
  expect_error(validate_chain_list(chainlist2), "parameters for each chain")

  chainlist[[1]] <- chainlist[[1]][-1, ]
  expect_error(validate_chain_list(chainlist),
               "Each chain should have the same number of iterations")
})

test_that("chain_list2array works", {
  expect_mcmc_array(chain_list2array(chainlist))
  expect_mcmc_array(chain_list2array(chainlist1))
  expect_mcmc_array(chain_list2array(chainlist1chain))

  expect_error(chain_list2array(dframe), "is_chain_list")
})



# transformations ---------------------------------------------------------
test_that("validate_transformations throws correct works", {
  trans <- list(x = "exp", 'beta[1]' = function(x) x^2, sigma = log)
  expect_silent(
    validate_transformations(trans, pars = c("x", "beta[1]", "sigma"))
  )

  trans2 <- trans
  trans2[[1]] <- match.fun(trans[[1]])
  expect_equal(
    validate_transformations(trans, pars = c("x", "beta[1]", "sigma")),
    trans2
  )
})
test_that("validate_transformations throws correct errors", {
  expect_error(
    validate_transformations(list("log", exp)),
    "must be a _named_ list"
  )
  expect_error(
    validate_transformations(list(x = "log", function(x) x^2)),
    "Each element of 'transformations' must have a name"
  )
  expect_error(
    validate_transformations(list(x = "log", 'beta[2]' = exp),
                             pars = c("x", "beta[1]")),
    regexp = "don't match parameter names: beta[2]", fixed = TRUE
  )
})

test_that("apply_transformations works", {
  trans <- list('beta[1]' = "exp", sigma = function(x) x^2)

  arr_trans <- apply_transformations(arr, trans)
  expect_equal(arr_trans[,, "sigma"], arr[,, "sigma"]^2)
  expect_equal(arr_trans[,, "beta[1]"], exp(arr[,, "beta[1]"]))

  mat_trans <- apply_transformations(mat, trans)
  expect_equal(mat_trans[, "sigma"], mat[, "sigma"]^2)
  expect_equal(mat_trans[, "beta[1]"], exp(mat[, "beta[1]"]))
})

test_that("transformations recycled properly if not a named list", {
  # if transformations is a single string naming a function
  x <- prepare_mcmc_array(arr, regex_pars = "beta", transformations = "exp")
  expect_identical(parameter_names(x), c("exp(beta[1])", "exp(beta[2])"))

  # if transformations is a single function
  x <- prepare_mcmc_array(arr, pars = c("beta[1]", "sigma"), transformations = exp)
  expect_identical(parameter_names(x), c("t(beta[1])", "t(sigma)"))
})

# rhat and neff helpers ---------------------------------------------------
test_that("diagnostic_factor.rhat works", {
  rhats <- new_rhat(c(low = 0.99, low = 1, low = 1.01,
                      ok = 1.06, ok = 1.09, ok = 1.1,
                      high = 1.2, high = 1.7))

  r <- diagnostic_factor(unname(rhats))
  expect_equivalent(r, as.factor(names(rhats)))
  expect_identical(levels(r), c("low", "ok", "high"))
})
test_that("diagnostic_factor.neff_ratio works", {
  ratios <- new_neff_ratio(c(low = 0.05, low = 0.01,
                             ok = 0.2, ok = 0.49,
                             high = 0.51, high = 0.99, high = 1))

  r <- diagnostic_factor(unname(ratios))
  expect_equivalent(r, as.factor(names(ratios)))
  expect_identical(levels(r), c("low", "ok", "high"))
})

