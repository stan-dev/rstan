## test tidy and glance methods from rstanarm_tidiers.R
stopifnot(require("testthat"), require("broom.mixed"), require("broom"))

if (suppressPackageStartupMessages(require(rstanarm, quietly = TRUE))) {
  fit <<- readRDS(system.file("extdata", "rstanarm_example.rds", package = "broom.mixed"))
  ## fit <- stan_glmer(mpg ~ wt + (1|cyl) + (1+wt|gear), data = mtcars,
  ## iter = 200, chains = 2)

  context("rstanarm tidiers")
  test_that("tidy works on rstanarm fits", {
    td1 <- tidy(fit)
    td2 <- tidy(fit, effects = "ran_vals")
    td3 <- tidy(fit, effects = "ran_pars")
    td4 <- tidy(fit, effects = "auxiliary")
    expect_equal(colnames(td1), c("term", "estimate", "std.error"))
  })

  test_that("tidy with multiple 'effects' selections works on rstanarm fits", {
    td1 <- tidy(fit, effects = c("ran_vals", "auxiliary"))
    expect_true(all(c("sigma", "mean_PPD") %in% td1$term))
    expect_equal(colnames(td1), c("term", "estimate", "std.error", "level", "group"))
  })

  test_that("conf.int works on rstanarm fits", {
    td1 <- tidy(fit, conf.int = TRUE, prob = 0.8)
    td2 <- tidy(fit, effects = "ran_vals", conf.int = TRUE, prob = 0.5)
    nms <- c("level", "group", "term", "estimate", "std.error", "conf.low", "conf.high")
    expect_equal(colnames(td2), nms)
  })

  test_that("glance works on rstanarm fits", {
    g1 <- glance(fit)
    g2 <- glance(fit, looic = TRUE, cores = 1, k_threshold = 0.7)
    expect_equal(colnames(g1), c("algorithm", "pss", "nobs", "sigma"))
    expect_equal(colnames(g2), c(colnames(g1), "looic", "elpd_loo", "p_loo"))
  })
} ## rstanarm available
