test_mklist_fun <- function() {
  c <- 4
  fun1 <- function(n) {
    a <- 3
    rstan:::mklist(n)[[n]]
  }
  checkEquals(fun1("a"), 3, checkNames = FALSE)
  checkEquals(fun1("c"), 4, checkNames = FALSE)
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
