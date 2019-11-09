# make.index.unique

test.make.index.unique_1us_default_eps <- function() {
  x <- .xts(1:5, rep(1e-6, 5))
  y <- make.index.unique(x)
  checkEqualsNumeric(.index(y), cumsum(rep(1e-6, 5)))
}

test.make.index.unique_returns_sorted_index <- function() {
  x <- .xts(1:5, c(rep(1e-6, 4), 3e-6))
  y <- make.index.unique(x, eps = 1e-6)
  checkEqualsNumeric(.index(y), cumsum(rep(1e-6, 5)))
}

test.make.index.unique_adds_eps_to_duplicates <- function() {
  epsilon <- c(1e-6, 1e-7, 1e-8)
  for (eps in epsilon) {
    x <- .xts(1:5, rep(eps, 5))
    y <- make.index.unique(x, eps = eps)
    checkEqualsNumeric(.index(y), cumsum(rep(eps, 5)))
  }
}

test.make.index.unique_no_warn_if_unique_timestamps_unchanged <- function() {
  x <- .xts(1:10, c(rep(1e-6, 9), 1e-5))
  y <- make.index.unique(x, eps = 1e-6)
  checkEqualsNumeric(.index(y), cumsum(rep(1e-6, 10)))
}

test.make.index.unique_warns_if_unique_timestamp_changes <- function() {
  # There should be a warning if the cumulative epsilon for a set of duplicate
  # index values is larger than the first unique index value that follows.
  # When this happens, we will overwrite that non-duplicate index value with
  # the prior index value + eps.
  eps <- 1e-6
  x <- .xts(1:5, c(rep(0, 4), 2*eps))
  orig <- options(warn = 2)
  on.exit(options(warn = orig$warn))
  checkException(y <- make.index.unique(x, eps = eps))
}

test.make.index.unique_warns_ONCE_if_unique_timestamp_changes <- function() {
  # There should be a warning if the cumulative epsilon for a set of duplicate
  # index values is larger than the first unique index value that follows.
  # When this happens, we will overwrite that non-duplicate index value with
  # the prior index value + eps.
  eps <- 1e-6
  x <- .xts(1:5, c(rep(0, 3), 2, 3) * eps)
  count <- 0L
  withCallingHandlers(make.index.unique(x, eps = eps),
                      warning = function(w) { count <<- count + 1L })
  checkEquals(count, 1L)
}


test.make.index.unique_converts_date_index_to_POSIXct <- function() {
  # It doesn't make sense to add a small epsilon to a date index. The C code
  # converts the integer index to a double, but it keeps the same index class.
  # The index class should be converted to POSIXct.
}
