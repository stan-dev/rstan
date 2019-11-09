# Tests for plotting functions
data(sample_matrix)
x <- as.xts(sample_matrix, dateFormat = "Date")

# axTicksByTime
test.format_xts_yearqtr <- function() {
  xq <- to.quarterly(x)
  xtbt <- axTicksByTime(xq)
  checkIdentical(names(xtbt), c("2007-Q1", "2007-Q2"))
}

test.format_zoo_yearqtr <- function() {
  xq <- to.quarterly(x)
  xtbt <- axTicksByTime(as.zoo(xq))
  checkIdentical(names(xtbt), c("2007-Q1", "2007-Q2"))
}

test.axTicksByTime_ticks.on_quarter <- function() {
  tick_marks <- setNames(c(1, 4, 7, 10, 13, 16, 19, 22, 25, 25),
    c("\nJan\n2016", "\nApr\n2016", "\nJul\n2016", "\nOct\n2016",
      "\nJan\n2017", "\nApr\n2017", "\nJul\n2017", "\nOct\n2017",
      "\nJan\n2018", "\nJan\n2018"))
  if (.Platform$OS.type != "unix") {
    names(tick_marks) <- gsub("\n(.*)\n", "\\1 ", names(tick_marks))
  }
  ym <- as.yearmon("2018-01") - 24:0 / 12
  x <- xts(seq_along(ym), ym)

  xtbt <- axTicksByTime(x, ticks.on = "quarters")
  checkIdentical(xtbt, tick_marks)
}
