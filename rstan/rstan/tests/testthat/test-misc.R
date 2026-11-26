test_that("mklist works", {
  skip("Backwards compatibility")

  c <- 4
  z <- matrix(0, ncol = 2, nrow = 3)
  x <- list()
  x[[1]] <- 1:2
  x[[2]] <- 3:4
  y <- 3
  fun1 <- function(n) {
    a <- 3
    mklist(n)[[n]]
  }

  expect_equal(fun1("a"), 3)
  expect_equal(fun1("c"), 4)

  # No list.
  L <- mklist(c("y", "z"))
  expect_equal(L$y, y)
  expect_equal(L$z, z)

  L <- mklist(c("x", "y", "z"))
  expect_equal(L$x, x)
  expect_equal(L$y, y)
  expect_equal(L$z, z)

  # Only a list.
  L <- mklist("x")
  expect_equal(L$x, x)
})

test_that("get_time_from_csv works", {
  skip("Backwards compatibility")

  expect_true(all(is.na(get_time_from_csv(""))))

  t_lines <- c("aa", "aa")
  expect_true(all(is.na(get_time_from_csv(t_lines))))

  t_lines <- c(
    "# Elapsed Time: 0.005308 seconds (Warm-up)",
    "#               0.003964 seconds (Sampling)"
  )
  t <- rstan:::get_time_from_csv(t_lines)

  expect_equal(unname(t), c(0.005308, 0.003964))
  expect_named(t, c("warmup", "sample"))
})

test_that("parse_data ignores transformed data variables", {
  code <- "
    data { int<lower=1> N; vector[N] y; }
    transformed data { row_vector[N] x = y'; }
    parameters { real mu; }
    model { y ~ normal(mu, 1); }
  "
  cppcode <- stanc(model_code = code)$cppcode
  data <- list(N = 3L, y = c(1, 2, 3))
  # neither x in a calling frame should be picked up for transformed data x.
  # f's argument x is a promise under evaluation while parse_data() runs,
  # like withVisible(x) when sampling() is called inside source()
  x <- 99
  f <- function(x) x
  expect_equal(f(with(data, rstan:::parse_data(cppcode))), data)
})

test_that("parse_data skips promises under evaluation", {
  code <- "
    data { real parse_data_test_x; }
    parameters { real mu; }
    model { mu ~ normal(parse_data_test_x, 1); }
  "
  cppcode <- stanc(model_code = code)$cppcode
  f <- function(parse_data_test_x) parse_data_test_x
  expect_equal(f(with(list(), rstan:::parse_data(cppcode))), list())
})

test_that("parse_data does not clobber data named like its own variables", {
  code <- "
    data { real stuff; }
    parameters { real mu; }
    model { mu ~ normal(stuff, 1); }
  "
  cppcode <- stanc(model_code = code)$cppcode
  expect_equal(with(list(stuff = 2), rstan:::parse_data(cppcode)),
               list(stuff = 2))
})

test_that("parse_data falls back to the global environment", {
  code <- "
    data { real parse_data_test_global; }
    parameters { real mu; }
    model { mu ~ normal(parse_data_test_global, 1); }
  "
  cppcode <- stanc(model_code = code)$cppcode
  assign("parse_data_test_global", 3, envir = globalenv())
  on.exit(rm("parse_data_test_global", envir = globalenv()))
  expect_equal(with(list(), rstan:::parse_data(cppcode)),
               list(parse_data_test_global = 3))
})
