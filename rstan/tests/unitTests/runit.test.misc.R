test_mklist_fun <- function() {
  c <- 4
  z <- matrix(0, ncol = 2, nrow = 3)
  x <- list()
  x[[1]] <- 1:2
  x[[2]] <- 3:4
  y <- 3
  fun1 <- function(n) {
    a <- 3
    rstan:::mklist(n)[[n]]
  }
  checkEquals(fun1("a"), 3, checkNames = FALSE)
  checkEquals(fun1("c"), 4, checkNames = FALSE)

  L <- rstan:::mklist(c("y", "z")) # no list
  checkEquals(L$y, y)
  checkEquals(L$z, z)

  L <- rstan:::mklist(c("x", "y", "z"))
  checkEquals(L$x, x)
  checkEquals(L$y, y)
  checkEquals(L$z, z)

  L <- rstan:::mklist("x") # only list
  checkEquals(L$x, x)
}

test_get_time_from_csv <- function() {
  tlines <- character(2)
  tlines[1] <- "aa"
  tlines[2] <- "aa"
  checkTrue(all(is.na(rstan:::get_time_from_csv(""))))
  checkTrue(all(is.na(rstan:::get_time_from_csv(tlines))))
  tlines[1] <- "# Elapsed Time: 0.005308 seconds (Warm-up)"
  tlines[2] <- "#               0.003964 seconds (Sampling)"
  t <- rstan:::get_time_from_csv(tlines)
  checkEquals(t, c(0.005308, 0.003964), checkNames = FALSE)
  checkEquals(names(t), c("warmup", "sample"))
}
