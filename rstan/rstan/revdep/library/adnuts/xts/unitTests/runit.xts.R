# Tests for xts constructors
#

### NA in order.by {{{
# .xts()
test..xts_order.by_NA_integer <- function() {
  checkException(.xts(1:3, c(1L, 2L, NA)))
  checkException(.xts(1:3, c(NA, 2L, 3L)))
  checkException(.xts(1:3, c(1L, NA, 3L)))
}
test..xts_order.by_NA_double <- function() {
  checkException(.xts(1:3, c(1, 2, NA)))
  checkException(.xts(1:3, c(NA, 2, 3)))
  checkException(.xts(1:3, c(1, NA, 3)))
}
test..xts_order.by_NaN_double <- function() {
  checkException(.xts(1:3, c(1, 2, NaN)))
  checkException(.xts(1:3, c(NaN, 2, 3)))
  checkException(.xts(1:3, c(1, NaN, 3)))
}
test..xts_order.by_Inf_double <- function() {
  checkException(.xts(1:3, c(1, 2,  Inf)))
  checkException(.xts(1:3, c(-Inf, 2, 3)))
}
# xts()
test.xts_order.by_NA_integer <- function() {
  checkException(xts(1:3, as.Date(c(1L, 2L, NA), origin = "1970-01-01")))
  checkException(xts(1:3, as.Date(c(NA, 2L, 3L), origin = "1970-01-01")))
  checkException(xts(1:3, as.Date(c(1L, NA, 3L), origin = "1970-01-01")))
}
test.xts_order.by_NA_double <- function() {
  checkException(xts(1:3, .POSIXct(c(1, 2, NA))))
  checkException(xts(1:3, .POSIXct(c(NA, 2, 3))))
  checkException(xts(1:3, .POSIXct(c(1, NA, 3))))
}
test.xts_order.by_NaN_double <- function() {
  checkException(xts(1:3, .POSIXct(c(1, 2, NaN))))
  checkException(xts(1:3, .POSIXct(c(NaN, 2, 3))))
  checkException(xts(1:3, .POSIXct(c(1, NaN, 3))))
}
test.xts_order.by_Inf_double <- function() {
  checkException(xts(1:3, .POSIXct(c(1, 2,  Inf))))
  checkException(xts(1:3, .POSIXct(c(-Inf, 2, 3))))
}
### }}}

# Test that only first tzone element is stored
test.xts_only_use_first_tzone_element <- function() {
  tz <- "America/Chicago"
  i <- as.POSIXlt("2018-01-01", tz = tz)
  y <- xts(1, i)
  checkIdentical(tz, tzone(y))
}

# .xts()
test..xts_dimnames_in_dots <- function() {
  x <- .xts(1:5, 1:5, dimnames = list(NULL, "x"))
  y <- xts(1:5, index(x), dimnames = list(NULL, "x"))
  checkEquals(x, y)
}

checkXtsClass <- function(xts, class) {
  checkEquals(tclass(xts), class)
  checkEquals(indexClass(xts), class)
  checkEquals(attr(attr(xts, "index"), "tclass"), class)
}

### Check that .indexCLASS takes precedence over tclass when both specified
test..xts_class <- function() {
  checkXtsClass(.xts(1, 1), c("POSIXct", "POSIXt"))
  checkXtsClass(.xts(1, 1, tclass="timeDate"), "timeDate")
  checkXtsClass(.xts(1, 1, .indexCLASS="Date"), "Date")
  checkXtsClass(.xts(1, 1, tclass="timeDate", .indexCLASS="Date"), "Date")

  ## also check that tclass is ignored if specified as part of index
  checkXtsClass(.xts(1, structure(1, tzone="",tclass="yearmon")), c("POSIXct", "POSIXt"))
  checkXtsClass(.xts(1, structure(1, tzone="",tclass="yearmon"), tclass="timeDate"), "timeDate")
  checkXtsClass(.xts(1, structure(1, tzone="",tclass="yearmon"), .indexCLASS="Date"), "Date")
  checkXtsClass(.xts(1, structure(1, tzone="",tclass="yearmon"), tclass="timeDate", .indexCLASS="Date"), "Date")
}

checkXtsTz <- function(xts, tzone) {
  checkEquals(tzone(xts), tzone)
  checkEquals(indexTZ(xts), tzone)
  checkEquals(attr(attr(xts, "index"), "tzone"), tzone)
}

### Check that tzone is honoured and .indexTZ ignored
test..xts_tzone <- function() {
  sysTZ <- Sys.getenv("TZ")
  Sys.setenv(TZ = "UTC")
  on.exit(Sys.setenv(TZ = sysTZ), add = TRUE)

  checkXtsTz(.xts(1, 1), "UTC")
  checkXtsTz(.xts(1, 1, tzone="Europe/London"), "Europe/London")
  ## this case passes in 0.10-2 but looks wrong
  checkXtsTz(.xts(1, 1, .indexTZ="America/New_York"), "UTC")
  checkXtsTz(.xts(1, 1, tzone="Europe/London", .indexTZ="America/New_York"), "Europe/London")

  ## Cases where tzone is specified in the index
  checkXtsTz(.xts(1, structure(1, tzone="Asia/Tokyo",tclass="yearmon")), "Asia/Tokyo")
  checkXtsTz(.xts(1, structure(1, tzone="Asia/Tokyo",tclass="yearmon"), tzone="Europe/London"), "Europe/London")
  checkXtsTz(.xts(1, structure(1, tzone="Asia/Tokyo",tclass="yearmon"), .indexTZ="America/New_York"), "Asia/Tokyo")
  checkXtsTz(.xts(1, structure(1, tzone="Asia/Tokyo",tclass="yearmon"), tzone="Europe/London", .indexTZ="America/New_York"), "Europe/London")
}
