
# POSIXct index
test.diff_integer_POSIXt <- function() {
  x <- .xts(1:5, 1:5 + 0.0)
  dx <- xts(rbind(NA_integer_, diff(coredata(x))), index(x))
  checkIdentical(diff(x), dx)
}
test.diff_numeric_POSIXt <- function() {
  x <- .xts(1:5 + 1.0, 1:5 + 0.0)
  dx <- xts(rbind(NA_real_, diff(coredata(x))), index(x))
  checkIdentical(diff(x), dx)
}
test.diff_logical_POSIXt <- function() {
  x <- .xts(1:5 > 2, 1:5 + 0.0)
  dx <- xts(rbind(NA, diff(coredata(x))), index(x))
  checkIdentical(diff(x), dx)
}

# Date index
test.diff_integer_Date <- function() {
  x <- xts(1:5, as.Date("2016-01-01") - 5:1)
  dx <- xts(rbind(NA_integer_, diff(coredata(x))), index(x))
  checkIdentical(diff(x), dx)
}
test.diff_numeric_Date <- function() {
  x <- xts(1:5 + 1.0, as.Date("2016-01-01") - 5:1)
  dx <- xts(rbind(NA_real_, diff(coredata(x))), index(x))
  checkIdentical(diff(x), dx)
}
test.diff_logical_Date <- function() {
  x <- xts(1:5 > 2, as.Date("2016-01-01") - 5:1)
  dx <- xts(rbind(NA, diff(coredata(x))), index(x))
  checkIdentical(diff(x), dx)
}

# Type-check failure errors
test.diff_differences_NA <- function() {
  x <- .xts(1:5, 1:5)
  checkException(diff(x, 1L, "a"), "'differences' must be integer")
}
test.diff_lag_NA <- function() {
  x <- .xts(1:5, 1:5)
  checkException(diff(x, "a", 1L), "'lag' must be integer")
}
test.diff_differences_LT1 <- function() {
  x <- .xts(1:5, 1:5)
  checkException(diff(x, 1L, -1L), "'diff.xts' defined only for positive lag and differences arguments")
}
test.diff_lag_LT1 <- function() {
  x <- .xts(1:5, 1:5)
  checkException(diff(x, -1L, 1L), "'diff.xts' defined only for positive lag and differences arguments")
}
