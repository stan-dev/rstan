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
