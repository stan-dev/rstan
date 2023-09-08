test_that("read_stan_csv works", {
  exfit <- read_stan_csv(dir(system.file("misc", package = "rstan"),
    pattern = "rstan_doc_ex_[[:digit:]].csv",
    full.names = TRUE
  ))
  expect_equal(dim(as.array(exfit)), c(100, 4, 10))
  expect_equal(exfit@model_pars, c("mu", "sigma", "z", "alpha", "lp__"))
  t <- matrix(c(0.005308, 0.003964, 0.004776, 0.003744, 0.010337, 0.002867, 0.004711, 0.004117),
    byrow = TRUE, ncol = 2
  )
  colnames(t) <- c("warmup", "sample")
  rownames(t) <- paste0("chain:", 1:nrow(t))
  expect_equal(get_elapsed_time(exfit), t)
})

test_that("read_stan_csv works with an incomplete csv (0 samples)", {
  skip("Backwards compatibility")

  exfit <- read_stan_csv(dir(system.file("misc", package = "rstan"),
    pattern = "rstan_doc_ex_incomplete_1.csv",
    full.names = TRUE
  ))
  expect_equal(exfit@model_pars, c("mu", "sigma", "z", "alpha", "lp__"))
  expect_true(is.null(get_elapsed_time(exfit)))
})

test_that("read_stan_csv works with an incomplete csv (few samples)", {
  skip("Backwards compatibility")

  exfit <- read_stan_csv(dir(system.file("misc", package = "rstan"),
    pattern = "rstan_doc_ex_incomplete_1.csv",
    full.names = TRUE
  ))
  expect_equal(exfit@model_pars, c("mu", "sigma", "z", "alpha", "lp__"))
})
