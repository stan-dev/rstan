# timeBasedSeq test
# 1999 to 2008 by year, Date
test.tbs_1999_to_2008_by_year_Date <- function() {
  tbs <- timeBasedSeq('1999/2008')
  bench <- seq(as.Date("1999-01-01"),as.Date("2008-01-01"),by='year')
  checkEqualsNumeric(tbs, bench)
}
# 1999 to 2008 by year, retclass='Date'
test.tbs_1999_to_2008_by_year_retclassDate <- function() {
  tbs <- timeBasedSeq('1999/2008', retclass='Date')
  bench <- seq(as.Date("1999-01-01"),as.Date("2008-01-01"),by='year')
  checkEqualsNumeric(tbs, bench)
}
# 1999 to 2008 by year, retclass="POSIXct"
test.tbs_1999_to_2008_by_year <- function() {
  tbs <- timeBasedSeq('1999/2008',retclass='POSIXct')
  bench <- seq(as.POSIXct("1999-01-01"),as.POSIXct("2008-01-01"),by='year')
  checkEqualsNumeric(tbs, bench)
}

# MONTHLY sequences
# defaults to yearmon from the zoo package
# NB: these differ by ~4.16e-5 on Solaris and rhub's windows-x86_64-devel
test.tbs_199901_to_200801_by_month <- function() {
  tbs <- timeBasedSeq('199901/200801')
  bench <- as.yearmon(seq(as.Date("1999-01-01"),as.Date("2008-01-01"),by='month'))
  checkEqualsNumeric(tbs, bench, tolerance = 1e-4)
}
test.tbs_199901_to_2008_by_month <- function() {
  tbs <- timeBasedSeq('199901/2008')
  bench <- as.yearmon(seq(as.Date("1999-01-01"),as.Date("2008-12-01"),by='month'))
  checkEqualsNumeric(tbs, bench, tolerance = 1e-4)
}
test.tbs_1999_to_200801_by_month <- function() {
  tbs <- timeBasedSeq('1999/200801')
  bench <- as.yearmon(seq(as.Date("1999-01-01"),as.Date("2008-01-01"),by='month'))
  checkEqualsNumeric(tbs, bench, tolerance = 1e-4)
}
# retclass=Date
test.tbs_199901_to_200801_by_month_Date <- function() {
  tbs <- timeBasedSeq('199901/200801', retclass='Date')
  bench <- seq(as.Date("1999-01-01"),as.Date("2008-01-01"),by='month')
  checkEqualsNumeric(tbs, bench)
}
test.tbs_199901_to_2008_by_month_Date <- function() {
  tbs <- timeBasedSeq('199901/2008', retclass='Date')
  bench <- seq(as.Date("1999-01-01"),as.Date("2008-12-01"),by='month')
  checkEqualsNumeric(tbs, bench)
}
test.tbs_1999_to_200801_by_month_Date <- function() {
  tbs <- timeBasedSeq('1999/200801', retclass='Date')
  bench <- as.Date(seq(as.Date("1999-01-01"),as.Date("2008-01-01"),by='month'))
  checkEqualsNumeric(tbs, bench)
}
# retclass=POSIXct
test.tbs_199901_to_200801_by_month_POSIXct <- function() {
  tbs <- timeBasedSeq('199901/200801', retclass='POSIXct')
  bench <- seq(as.POSIXct("1999-01-01"),as.POSIXct("2008-01-01"),by='month')
  checkEqualsNumeric(tbs, bench)
}
test.tbs_199901_to_2008_by_month_POSIXct <- function() {
  tbs <- timeBasedSeq('199901/2008', retclass='POSIXct')
  bench <- as.POSIXct(seq(as.POSIXct("1999-01-01"),as.POSIXct("2008-12-01"),by='month'),tzone='GMT')
  checkEqualsNumeric(tbs, bench)
}
test.tbs_1999_to_200801_by_month_POSIXct <- function() {
  tbs <- timeBasedSeq('1999/200801', retclass='POSIXct')
  bench <- seq(as.POSIXct("1999-01-01"),as.POSIXct("2008-01-01"),by='month')
  checkEqualsNumeric(tbs, bench)
}
