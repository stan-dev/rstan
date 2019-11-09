#
# RUnit tests for the following 'xts' methods:
# rbind
# cbind
#
test.rbind_zero_length_non_zero_length_POSIXct_errors <- function() {
  xpz <- xts( , as.POSIXct("2017-01-01"))
  xp1 <- xts(1, as.POSIXct("2017-01-02"))
  zpz <- as.zoo(xpz)
  zp1 <- as.zoo(xp1)
  zpe <- tryCatch(rbind(zpz, zp1), error = identity)
  xpe <- tryCatch(rbind(xpz, xp1), error = identity)
  checkIdentical(zpe$message, xpe$message)
}
test.rbind_zero_length_non_zero_length_Date_errors <- function() {
  xpz <- xts( , as.Date("2017-01-01"))
  xp1 <- xts(1, as.Date("2017-01-02"))
  zpz <- as.zoo(xpz)
  zp1 <- as.zoo(xp1)
  zpe <- tryCatch(rbind(zpz, zp1), error = identity)
  xpe <- tryCatch(rbind(xpz, xp1), error = identity)
  checkIdentical(zpe$message, xpe$message)
}

# Test that as.Date.numeric() works at the top level (via zoo::as.Date()),
# and for functions defined in the xts namespace even if xts::as.Date.numeric()
# is not formally registered as an S3 method.
test.as.Date.numeric <- function() {
  # Define function that calls as.Date.numeric() ...
  f <- function(d) {
    as.Date(d)
  }
  # ... in xts' namespace
  environment(f) <- as.environment("package:xts")

  dd <- as.Date("2017-12-13")
  dn <- unclass(dd)
  checkIdentical(dd, as.Date(dn))  # via zoo::as.Date()
  checkIdentical(dd, f(dn))
}

# .subset.xts
# window.xts
# .toPOSIXct (indirectly)

test.window <- function() {
  # window function for xts series, use basic logic for testing & debugging
  # null start and end not supported
  window_dbg <- function(x, index. = index(x), start, end)
  {
    start <- xts:::.toPOSIXct(start, indexTZ(x))
    end <- xts:::.toPOSIXct(end, indexTZ(x))
    index. <- as.POSIXct(index., tz=indexTZ(x))
    all.indexes <- .index(x)
    in.index <- all.indexes %in% index.
    matches <- (in.index & all.indexes >= start & all.indexes <= end)
    x[matches,]
  }

  DAY = 24*3600
  base <- as.POSIXct("2000-12-31")
  dts <- base + c(1:10, 12:15, 17:20)*DAY
  x <- xts(1:length(dts), dts)

  # Range over gap
  start <- base + 11*DAY
  end <- base + 16*DAY
  bin <- window(x, start = start, end = end)
  reg <- window_dbg(x, start = start, end = end)
  checkIdentical(bin, reg, "Range over gap")

  # Range over one day
  start <- base + 12*DAY
  end <- base + 12*DAY
  bin <- window(x, start = start, end = end)
  reg <- window_dbg(x, start = start, end = end)
  checkIdentical(bin, reg, "Range over one day")

  # Empty Range over one day
  start <- base + 11*DAY
  end <- base + 11*DAY
  bin <- window(x, start = start, end = end)
  reg <- window_dbg(x, start = start, end = end)
  checkIdentical(bin, reg, "Empty Range over one day")

  # Range containing all dates
  start <- base
  end <- base + 21*DAY
  bin <- window(x, start = start, end = end)
  reg <- window_dbg(x, start = start, end = end)
  checkIdentical(bin, reg, "Range containing all dates")

  # Range past end
  start <- base + 16*DAY
  end <- base + 30*DAY
  bin <- window(x, start = start, end = end)
  reg <- window_dbg(x, start = start, end = end)
  checkIdentical(bin, reg, "Range past end")

  # Range before begin
  start <- base
  end <- base + 3*DAY
  bin <- window(x, start = start, end = end)
  reg <- window_dbg(x, start = start, end = end)
  checkIdentical(bin, reg, "Range before begin")

  # Test just start, end = NULL
  start <- base + 13 * DAY
  end <- base + 30*DAY
  bin <- window(x, start = start)
  reg <- window_dbg(x, start = start, end = end)
  checkIdentical(bin, reg, "Test just start, end = NULL")

  # Test just start, end = NULL, empty range
  start <- base + 25 * DAY
  end <- base + 30*DAY
  bin <- window(x, start = start)
  reg <- window_dbg(x, start = start, end = end)
  checkIdentical(bin, reg, "Test just start, end = NULL, empty range")

  # Test just end, start = NULL
  end <- base + 13 * DAY
  start <- base
  bin <- window(x, end = end)
  reg <- window_dbg(x, start = start, end = end)
  checkIdentical(bin, reg, "Test just end, start = NULL")

  # Test just end, start = NULL, empty range
  end <- base
  start <- base
  bin <- window(x, end = end)
  reg <- window_dbg(x, start = start, end = end)
  checkIdentical(bin, reg, "Test just end, start = NULL, empty range")

  # Test end = NULL, start = NULL
  start <- base
  end <- base + 30*DAY
  bin <- window(x)
  reg <- window_dbg(x, start = start, end = end)
  checkIdentical(bin, reg, "Test end = NULL, start = NULL")

  #######################################
  # Test for index. parameter
  start <- base
  end <- base + 30*DAY
  idx = index(x)[c(2,4,6)]
  bin <- window(x, index. = idx)
  reg <- window_dbg(x, index. = idx, start = start, end = end)
  checkIdentical(bin, reg, "Test for index. parameter")

  # Test index. outside range of dates in xts series
  start <- base
  end <- base + 30*DAY
  idx = c(start, index(x)[c(2,4,6)], end)
  bin <- window(x, index. = idx)
  reg <- window_dbg(x, index. = idx, start = start, end = end)
  checkIdentical(bin, reg, "Test index. outside range of dates in xts series")

  # Test NA in index
  start <- base
  end <- base + 30*DAY
  idx = c(start, index(x)[c(2,4,6)], end, NA)
  bin <- window(x, index. = idx)
  reg <- window_dbg(x, index. = idx, start = start, end = end)
  checkIdentical(bin, reg, "Test NA in index ")

  # Next 3 adapted from window.zoo example
  # Test basic window.zoo example
  x.date <- as.Date(paste(2003, rep(1:4, 4:1), seq(1,19,2), sep = "-"))
  x <- xts(matrix(1:20, ncol = 2), x.date)
  bin <- window(x, start = as.Date("2003-02-01"), end = as.Date("2003-03-01"))
  reg <- window_dbg(x, start = as.Date("2003-02-01"), end = as.Date("2003-03-01"))
  checkIdentical(bin, reg, "Test basic window.zoo example")

  # Test index + start
  bin <- window(x, index = x.date[1:6], start = as.Date("2003-02-01"))
  reg <- window_dbg(x, index = x.date[1:6], start = as.Date("2003-02-01"), end = as.Date("2004-01-01"))
  checkIdentical(bin, reg, "Test index + start")

  # Test just index
  bin <- window(x, index = x.date[c(4, 8, 10)])
  reg <- window_dbg(x, index = x.date[c(4, 8, 10)], start = as.Date("2003-01-01"), end = as.Date("2004-01-01"))
  checkIdentical(bin, reg, "Test just index")

  # Test decreasing index
  bin <- window(x, index = x.date[c(10, 8, 4)])
  reg <- window_dbg(x, index = x.date[c(10, 8, 4)], start = as.Date("2003-01-01"), end = as.Date("2004-01-01"))
  checkIdentical(bin, reg, "Test decreasing index")

  # Test index parameter with repeated dates in xts series
  idx <- sort(rep(1:5, 5))
  x <- xts(1:length(idx), as.Date("1999-12-31")+idx)
  bin <- window(x, index = as.Date("1999-12-31")+c(1,3,5))
  reg <- window_dbg(x, index = as.Date("1999-12-31")+c(1,3,5), start = as.Date("2000-01-01"), end = as.Date("2000-01-05"))
  checkIdentical(bin, reg, "Test index parameter with repeated dates in xts series")
  checkTrue(nrow(bin) == 3*5, "Test index parameter with repeated dates in xts series")

  # Test performance difference
  DAY = 24*3600
  base <- as.POSIXct("2000-12-31")
  dts <- base + c(1:10, 12:15, 17:20)*DAY
  x <- xts(1:length(dts), dts)
  start <- base + 14*DAY
  end <- base + 14*DAY
  cat("\n")
  print("performance:")
  print("binary search")
  print(system.time(replicate(1000, window(x, start = start, end = end)))) # Binary search is about 2x faster than regular
  print("regular search")
  print(system.time(replicate(1000, window_dbg(x, start = start, end = end))))
}

# test subset.xts for date subsetting by row
test.subset_i_datetime_or_character <- function() {
  base <- as.POSIXct("2000-12-31")
  dts <- base + c(1:10, 12:15, 17:20) * 24L * 3600L
  x <- xts(seq_along(dts), dts)

  # Note that "2001-01-11" is not in the series. Skipped by convention.
  d <- c("2001-01-10", "2001-01-11", "2001-01-12", "2001-01-13")

  for (type in c("double", "integer")) {
    storage.mode(.index(x)) <- type

    # Test scalar
    msg <- paste("scalar,", type, "index")
    bin <- window(x, start = d[1], end = d[1])
    checkIdentical(bin, x[d[1], ], paste("character", msg))
    checkIdentical(bin, x[I(d[1]), ], paste("as-is character", msg))
    checkIdentical(bin, x[as.POSIXct(d[1]), ], paste("POSIXct", msg))
    checkIdentical(bin, x[as.Date(d[1]), ], paste("Date", msg))

    # Test vector
    msg <- paste("vector,", type, "index")
    bin <- window(x, start = d[1], end = d[length(d)])
    checkIdentical(bin, x[d, ], paste("character", msg))
    checkIdentical(bin, x[I(d), ], paste("as-is character", msg))
    checkIdentical(bin, x[as.POSIXct(d), ], paste("POSIXct", msg))
    checkIdentical(bin, x[as.Date(d), ], paste("Date", msg))

    # Test character dates, and single column selection
    y <- xts(rep(2, length(dts)), dts)
    z <- xts(rep(3, length(dts)), dts)
    x2 <- cbind(y, x, z)
    sub <- x2[d, 2]  # Note that "2001-01-11" is not in the series. Skipped by convention.
    bin <- window(x, start = d[1], end = d[length(d)])
    checkTrue(nrow(sub) == nrow(bin), "Test character dates, and single column selection")
    checkTrue(all(sub == bin), "Test character dates, and single column selection")
  }
}

test.subset_i_ISO8601 <- function() {
  x <- xts(1:1000, as.Date("2000-01-01")+1:1000)

  for (type in c("double", "integer")) {
    storage.mode(.index(x)) <- type

    fmt <- paste("Test date range, %s;", type, "index")
    # Test Date Ranges
    sub <- x['200001'] # January 2000
    bin <- window(x, start = "2000-01-01", end = "2000-01-31")
    checkIdentical(bin, sub, sprintf(fmt, "2000-01"))

    # Test Date Ranges 2
    sub <- x['1999/2000'] # All of 2000 (note there is no need to use the exact start)
    bin <- window(x, start = "2000-01-01", end = "2000-12-31")
    checkIdentical(bin, sub, sprintf(fmt, "1999/2000"))

    # Test Date Ranges 3
    sub <- x['1999/200001'] # January 2000
    bin <- window(x, start = "2000-01-01", end = "2000-01-31")
    checkIdentical(bin, sub, sprintf(fmt, "1999/2000-01"))
  }
}
