test_mklist_fun <- function() {
  c <- 4
  fun1 <- function(n) {
    a <- 3
    rstan:::mklist(n)[[n]]
  }
  checkEquals(fun1("a"), 3, checkNames = FALSE)
  checkEquals(fun1("c"), 4, checkNames = FALSE)
}
