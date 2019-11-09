LAG <- function(x, k=1, na.pad=TRUE) {
  z <- lag(as.zoo(x), -k, na.pad)
  dimnames(z) <- NULL
  as.xts(z)
}

# POSIXct index
test.lag_integer_POSIXt <- function() {
  x <- .xts(1:5, 1:5 + 0.0)
  checkIdentical(lag(x), LAG(x))
}
test.lag_numeric_POSIXt <- function() {
  x <- .xts(1:5 + 1.0, 1:5 + 0.0)
  checkIdentical(lag(x), LAG(x))
}
test.lag_logical_POSIXt <- function() {
  x <- .xts(1:5 > 2, 1:5 + 0.0)
  checkIdentical(lag(x), LAG(x))
}

# Date index
test.lag_integer_Date <- function() {
  x <- xts(1:5, as.Date("2016-01-01") - 5:1)
  checkIdentical(lag(x), LAG(x))
}
test.lag_numeric_Date <- function() {
  x <- xts(1:5 + 1.0, as.Date("2016-01-01") - 5:1)
  checkIdentical(lag(x), LAG(x))
}
test.lag_logical_Date <- function() {
  x <- xts(1:5 > 2, as.Date("2016-01-01") - 5:1)
  checkIdentical(lag(x), LAG(x))
}

# Type-check failure errors
test.lag_k_NA <- function() {
  x <- .xts(1:5, 1:5)
  checkException(lag(x, "a"), "'k' must be integer", TRUE)
}
test.lag_k_zero_length <- function() {
  x <- .xts(1:5, 1:5)
  checkException(lag(x, 1L, "a"), "'na.pad' must be logical", TRUE)
}

