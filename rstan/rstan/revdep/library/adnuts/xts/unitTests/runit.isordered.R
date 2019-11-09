# Tests for isOrdered()
#

# Utility functions for tests {{{
run.isOrdered <- function(x) {
  c(isOrdered(x,  TRUE,  TRUE),
    isOrdered(x,  TRUE, FALSE),
    isOrdered(x, FALSE, FALSE),
    isOrdered(x, FALSE,  TRUE))
}
check.isOrdered <- function(x, v = rep(TRUE, 4)) {
  xc <- paste(capture.output(dput(x)), collapse = " ")
  checkIdentical(v[1], isOrdered(x,  TRUE,  TRUE), paste(xc, v[1], "increasing, strictly"))
  checkIdentical(v[2], isOrdered(x,  TRUE, FALSE), paste(xc, v[2], "increasing"))
  checkIdentical(v[3], isOrdered(x, FALSE, FALSE), paste(xc, v[3], "decreasing"))
  checkIdentical(v[4], isOrdered(x, FALSE,  TRUE), paste(xc, v[4], "decreasing, strictly"))
}
# }}}

TTTT <- rep(TRUE, 4)
FFFF <- !TTTT
TTFF <- c(TRUE, TRUE, FALSE, FALSE)
FFTT <- !TTFF

# Increasing {{{
test.isOrdered_incr <- function() {
  check.isOrdered(1:3, TTFF)
  check.isOrdered(-1:1, TTFF)
  check.isOrdered(c(1, 2, 3), TTFF)
  check.isOrdered(c(-1, 0, 1), TTFF)
}
### NA, NaN, Inf
# beg
test.isOrdered_incr_begNA <- function() {
  check.isOrdered(c(NA_integer_, 1L, 2L), FFFF)
  check.isOrdered(c(NA_real_, 1, 2), TTFF)
  check.isOrdered(c(NaN, 1, 2), TTFF)
  check.isOrdered(c(Inf, 1, 2), FFFF)
  check.isOrdered(c(-Inf, 1, 2), TTFF)
}
# mid
test.isOrdered_incr_midNA <- function() {
  check.isOrdered(c(1L, NA_integer_, 2L), FFFF)
  check.isOrdered(c(1, NA_real_, 2), TTTT)
  check.isOrdered(c(1, NaN, 2), TTTT)
  check.isOrdered(c(1, Inf, 2), FFFF)
  check.isOrdered(c(1, -Inf, 2), FFFF)
}
# end
test.isOrdered_incr_endNA <- function() {
  check.isOrdered(c(1L, 2L, NA_integer_), TTFF)
  check.isOrdered(c(1, 2, NA_real_), TTFF)
  check.isOrdered(c(1, 2, NaN), TTFF)
  check.isOrdered(c(1, 2, Inf), TTFF)
  check.isOrdered(c(1, 2, -Inf), FFFF)
}
###
# }}}

# Decreasing {{{
test.isOrdered_decr <- function() {
  check.isOrdered(1:-1, FFTT)
  check.isOrdered(3:1, FFTT)
  check.isOrdered(c(3, 2, 1), FFTT)
  check.isOrdered(c(1, 0, -1), FFTT)
}
### NA, NaN, Inf
# beg
test.isOrdered_decr_begNA <- function() {
  check.isOrdered(c(NA_integer_, 2L, 1L), FFTT)
  check.isOrdered(c(NA_real_, 2, 1), FFTT)
  check.isOrdered(c(NaN, 2, 1), FFTT)
  check.isOrdered(c(Inf, 2, 1), FFTT)
  check.isOrdered(c(-Inf, 2, 1), FFFF)
}
# mid
test.isOrdered_decr_midNA <- function() {
  check.isOrdered(c(2L, NA_integer_, 1L), FFFF)
  check.isOrdered(c(2, NA_real_, 1), TTTT)
  check.isOrdered(c(2, NaN, 1), TTTT)
  check.isOrdered(c(2, Inf, 1), FFFF)
  check.isOrdered(c(2, -Inf, 1), FFFF)
}
# end
test.isOrdered_decr_endNA <- function() {
  check.isOrdered(c(2L, 1L, NA_integer_), FFFF)
  check.isOrdered(c(2, 1, NA_real_), FFTT)
  check.isOrdered(c(2, 1, NaN), FFTT)
  check.isOrdered(c(2, 1, Inf), FFFF)
  check.isOrdered(c(2, 1, -Inf), FFTT)
}
###
# }}}
