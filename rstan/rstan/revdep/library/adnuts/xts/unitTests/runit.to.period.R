# ensure first group is included in output
test.to.frequency_includes_first_group <- function() {
  data(sample_matrix)
  x <- as.xts(sample_matrix)
  x$Volume <- 1
  
  tf <- xts:::to.frequency(x, x$Volume, 90, name=NULL)
  tp <- .Call("toPeriod", x, c(0L, 90L, 180L), TRUE, 5L, FALSE, FALSE,
              c("Open", "High", "Low", "Close", "Volume") , PACKAGE="xts")

  checkIdentical(tf, tp)
}

