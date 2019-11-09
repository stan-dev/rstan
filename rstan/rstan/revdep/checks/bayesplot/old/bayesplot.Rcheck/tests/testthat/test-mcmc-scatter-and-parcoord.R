library(bayesplot)
suppressPackageStartupMessages(library(rstanarm))
context("MCMC: scatter and parallel coordinates plots")

source(test_path("data-for-mcmc-tests.R"))

# also fit an rstanarm model to use with mcmc_pairs
capture.output(
  fit <- stan_glm(mpg ~ wt + am, data = mtcars, iter = 1000, chains = 2, refresh = 0)
)
post <- as.array(fit)
lp <- log_posterior(fit)
np <- nuts_params(fit)
divs <- sample(c(0,1), size = 1000, prob = c(0.25, 0.75), replace = TRUE)
np$Value[np$Parameter=="divergent__"] <- divs # fake divergences


test_that("mcmc_scatter returns a ggplot object", {
  expect_gg(mcmc_scatter(arr, pars = c("beta[1]", "beta[2]")))
  expect_gg(mcmc_scatter(arr1chain, regex_pars = "beta", size = 3, alpha = 0.5))
  expect_gg(mcmc_scatter(mat, pars = c("sigma", "(Intercept)")))
  expect_gg(mcmc_scatter(dframe, regex_pars = "x:[2,4]"))
  expect_gg(mcmc_scatter(dframe_multiple_chains,
                         pars = c("sigma", "(Intercept)")))
})

test_that("mcmc_hex returns a ggplot object", {
  expect_gg(mcmc_hex(arr, pars = c("beta[1]", "beta[2]")))
  expect_gg(mcmc_hex(arr1chain, regex_pars = "beta", binwidth = c(.5,.5)))
})

test_that("mcmc_scatter & mcmc_hex throw error if only 1 parameter", {
  expect_error(mcmc_scatter(arr, pars = "sigma"), "exactly 2 parameters")
  expect_error(mcmc_scatter(arr1), "exactly 2 parameters")
  expect_error(mcmc_scatter(mat1), "exactly 2 parameters")
  expect_error(mcmc_hex(dframe1), "exactly 2 parameters")
  expect_error(mcmc_hex(chainlist1), "exactly 2 parameters")
})

test_that("mcmc_scatter accepts NUTS info", {
  expect_gg(mcmc_scatter(post, pars = c("wt", "sigma"), np = np))

  div_style <- scatter_style_np(div_color = "orange", div_size = 2,
                                div_shape = 3, div_alpha = 0.5)
  g <- mcmc_scatter(post, pars = c("wt", "sigma"), np = np, np_style = div_style)
  expect_gg(g)
  expect_named(g$data, c("x", "y", "Divergent"))
})


# mcmc_pairs  -------------------------------------------------------------
test_that("mcmc_pairs returns a bayesplot_grid object", {
  expect_bayesplot_grid(mcmc_pairs(arr, pars = c("(Intercept)", "sigma")))
  expect_bayesplot_grid(mcmc_pairs(arr, pars = "sigma", regex_pars = "beta"))
  expect_bayesplot_grid(mcmc_pairs(arr, regex_pars = "x:[1-3]",
                                   transformations = "exp",
                                   diag_fun = "dens", off_diag_fun = "hex",
                                   diag_args = list(trim = FALSE),
                                   off_diag_args = list(binwidth = c(0.5, 0.5))))

  expect_bayesplot_grid(suppressWarnings(mcmc_pairs(arr1chain, regex_pars = "beta")))
  expect_bayesplot_grid(suppressWarnings(mcmc_pairs(mat, pars = c("(Intercept)", "sigma"))))
  expect_bayesplot_grid(suppressWarnings(mcmc_pairs(dframe, pars = c("(Intercept)", "sigma"))))
  expect_bayesplot_grid(mcmc_pairs(dframe_multiple_chains, regex_pars = "beta"))
})

test_that("no mcmc_pairs non-NUTS 'condition's fail", {
  expect_bayesplot_grid(
    mcmc_pairs(arr, pars = "sigma", regex_pars = "beta",
               condition = pairs_condition(chains = list(1, 2:4)))
    )
  expect_bayesplot_grid(
    mcmc_pairs(arr, pars = "sigma", regex_pars = "beta",
               condition = pairs_condition(draws = rep(c(T,F), length.out = prod(dim(arr)[1:2]))))
    )
  expect_bayesplot_grid(
    mcmc_pairs(arr, pars = "sigma", regex_pars = "beta",
               condition = pairs_condition(draws = 1/3))
  )
  expect_bayesplot_grid(
    mcmc_pairs(arr, pars = "sigma", regex_pars = "beta",
               condition = pairs_condition(chains = c(1,3)))
  )
})

test_that("mcmc_pairs works with NUTS info", {
  expect_bayesplot_grid(mcmc_pairs(post, pars = c("wt", "am", "sigma"), np = np))
  expect_bayesplot_grid(mcmc_pairs(post, pars = c("wt", "am"),
                                   condition = pairs_condition(nuts="energy__"), np = np))
  expect_bayesplot_grid(mcmc_pairs(post, pars = c("wt", "am"),
                                   condition = pairs_condition(nuts="divergent__"), np = np))
  expect_bayesplot_grid(mcmc_pairs(post, pars = c("wt", "am"),
                                   condition = pairs_condition(nuts = "lp__"), lp=lp, np = np,
                                   max_treedepth = 2))

  p <- mcmc_pairs(
    post,
    pars = c("wt", "am"),
    off_diag_fun = "hex",
    condition = pairs_condition(nuts = "lp__"),
    lp = lp,
    np = np,
    np_style = pairs_style_np(div_color = "firebrick", td_color = "dodgerblue", div_size = 2, td_size = 2),
    max_treedepth = with(np, max(Value[Parameter == "treedepth__"]) - 1)
  )
  expect_bayesplot_grid(p)
})


test_that("mcmc_pairs throws correct warnings and errors", {
  expect_warning(mcmc_pairs(arr1chain, regex_pars = "beta"),
                 "This plot is more useful with multiple chains")
  expect_error(mcmc_pairs(arr, pars = "sigma"),
               "requires at least two parameters")

  expect_error(
    mcmc_pairs(arr, condition = pairs_condition(draws = c(T, F))),
    "length(condition) == (n_iter * n_chain) is not TRUE",
    fixed = TRUE
  )
  expect_error(
    mcmc_pairs(arr, condition = pairs_condition(nuts = "accept_stat__")),
    "the 'np' argument to 'mcmc_pairs' must also be specified"
  )
  expect_error(
    mcmc_pairs(arr, condition = pairs_condition(nuts = "lp__")),
    "the 'lp' argument to 'mcmc_pairs' must also be specified"
  )
  expect_error(
    mcmc_pairs(arr, condition = "lp__"),
    'inherits(condition, "pairs_condition") is not TRUE',
    fixed = TRUE
  )

  expect_error(
    mcmc_pairs(post, pars = c("wt", "am"), max_treedepth = 2, np = np,
               np_style = list(color = "green")),
    'inherits(np_style, "nuts_style") is not TRUE',
    fixed = TRUE
  )

  post2 <- post
  post2[,1:2,"wt"] <- 0
  expect_warning(
    mcmc_pairs(post2, pars = c("wt", "am", "sigma")),
    "parameters were dropped because they are constant: wt"
  )

  post[,, "sigma"] <- post[,, "am"]
  expect_warning(
    mcmc_pairs(post, pars = c("wt", "sigma", "am")),
    "parameters were dropped because they are duplicative: am"
  )
})


# pairs_style_np -------------------------------------------------------
test_that("pairs_style_np returns correct structure", {
  style <- pairs_style_np(div_size = 3, td_color = "gray", td_shape = 1)
  expect_s3_class(style, "nuts_style")
  expect_named(style, c("color", "shape", "size", "alpha"), ignore.order = TRUE)
  expect_named(style$color, c("div", "td"))
  expect_named(style$size, c("div", "td"))
  expect_named(style$shape, c("div", "td"))
  expect_named(style$alpha, c("div", "td"))
})

test_that("pairs_style_np throws correct errors", {
  expect_error(
    pairs_style_np(div_size = "3"),
    "is.numeric(div_size) is not TRUE",
    fixed = TRUE
  )
  expect_error(
    pairs_style_np(td_color = 1),
    "is.character(td_color) is not TRUE",
    fixed = TRUE
  )
})


# pairs_condition ---------------------------------------------------------
test_that("pairs_condition returns correct structure", {
  # default
  cond0 <- pairs_condition()
  expect_s3_class(cond0, "pairs_condition")
  expect_equivalent(unclass(cond0), list())
  expect_equal(attr(cond0, "type"), "default")

  # chains
  cond1 <- pairs_condition(chains = 1:4)
  expect_s3_class(cond1, "integer")
  expect_s3_class(cond1, "pairs_condition")
  expect_equivalent(unclass(cond1), 1:4)
  expect_equal(attr(cond1, "type"), "chain_vector")

  cond2 <- pairs_condition(chains = list(1:4, 5:6))
  expect_s3_class(cond2, "list")
  expect_s3_class(cond2, "pairs_condition")
  expect_equivalent(unclass(cond2), list(upper=1:4, lower=5:6))
  expect_equal(attr(cond2, "type"), "chain_list")

  # draws
  cond3 <- pairs_condition(draws = 0.7)
  expect_s3_class(cond3, "numeric")
  expect_s3_class(cond3, "pairs_condition")
  expect_equivalent(unclass(cond3), 0.7)
  expect_equal(attr(cond3, "type"), "draws_proportion")

  cond4 <- pairs_condition(draws = c(T, F, T))
  expect_s3_class(cond4, "logical")
  expect_s3_class(cond4, "pairs_condition")
  expect_equivalent(unclass(cond4), c(T, F, T))
  expect_equal(attr(cond4, "type"), "draws_selection")

  # nuts
  cond5 <- pairs_condition(nuts = "lp__")
  expect_s3_class(cond5, "character")
  expect_s3_class(cond5, "pairs_condition")
  expect_equivalent(unclass(cond5), "lp__")
  expect_equal(attr(cond5, "type"), "nuts")
})

test_that("pairs_condition throws correct errors", {
  # chain
  expect_error(
    pairs_condition(chains = "abc"),
    "must be an integer vector or a list of two integer vectors"
  )
  expect_error(
    pairs_condition(chains = list(1:2, 3:4, 5:6)),
    "length(chains) == 2 is not TRUE",
    fixed = TRUE
  )
  expect_error(
    pairs_condition(chains = list(1:2, 2:3)),
    "Each chain can only be specified once"
  )
  expect_error(
    pairs_condition(chains = c(1:3, 2)),
    "Each chain can only be specified once"
  )

  # draws
  expect_error(
    pairs_condition(draws = "abc"),
    "must be a single proportion or a logical vector"
  )
  expect_error(
    pairs_condition(draws = 2),
    "draws > 0 && draws < 1 is not TRUE",
    fixed = TRUE
  )

  # nuts
  expect_error(
    pairs_condition(nuts = 2),
    "must be a single string"
  )
  expect_error(
    pairs_condition(nuts = c("lp__", "energy__")),
    "must be a single string"
  )
  expect_error(
    pairs_condition(nuts = "step_size__"),
    "stepsize__"
  )
})

test_that("pairs_condition message if multiple args specified", {
  options(useFancyQuotes = FALSE)
  expect_message(
    pairs_condition(chains = 2, draws = 0.5, nuts = "lp__"),
    "because they are superseded by 'chains': 'draws', 'nuts'",
    fixed = TRUE
  )
  expect_message(
    pairs_condition(chains = 2, nuts = "lp__"),
    "because they are superseded by 'chains': 'nuts'",
    fixed = TRUE
  )
  expect_message(
    pairs_condition(draws = 0.5, nuts = "lp__"),
    "because they are superseded by 'draws': 'nuts'",
    fixed = TRUE
  )
})



# mcmc_parcoord -----------------------------------------------------------
test_that("mcmc_parcoord returns a ggplot object", {
  expect_gg(mcmc_parcoord(arr, pars = c("(Intercept)", "sigma")))
  expect_gg(mcmc_parcoord(arr, pars = "sigma", regex_pars = "beta"))

  # with nuts info
  expect_gg(mcmc_parcoord(post, pars = c("wt", "am", "sigma"), np = np))
})

test_that("mcmc_parcoord throws correct warnings and errors", {
  expect_error(mcmc_parcoord(arr, pars = "sigma"),
               "requires at least two parameters")

  expect_error(
    mcmc_parcoord(post, np = np[, -1]),
    "NUTS parameter data frame must have columns: Iteration, Parameter, Value, Chain",
    fixed = TRUE
  )

  expect_error(
    mcmc_parcoord(post, np = np, np_style = list(div_color = "green")),
    'inherits(np_style, "nuts_style") is not TRUE',
    fixed = TRUE
  )
})


# parcoord_style_np -------------------------------------------------------
test_that("parcoord_style_np returns correct structure", {
  style <- parcoord_style_np()
  expect_s3_class(style, "nuts_style")
  expect_named(style, c("color", "alpha", "size"), ignore.order = TRUE)
  expect_named(style$color, c("div"))
  expect_named(style$size, c("div"))
  expect_named(style$alpha, c("div"))
})

test_that("parcoord_style_np throws correct errors", {
  expect_error(
    parcoord_style_np(div_size = "3"),
    "is.numeric(div_size) is not TRUE",
    fixed = TRUE
  )
  expect_error(
    parcoord_style_np(td_color = 1),
    "unused argument (td_color = 1)",
    fixed = TRUE
  )
})
