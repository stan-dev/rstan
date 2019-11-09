data(sample_matrix)
sample.data.frame <- data.frame(sample_matrix)
sample.xts <- as.xts(sample.data.frame)

test.convert_data.frame_to_xts <- function() {
  checkIdentical(sample.xts,as.xts(sample.data.frame))
}

test.convert_data.frame_to_xts_j1 <- function() {
  checkIdentical(sample.xts[,1],as.xts(sample.data.frame)[,1])
}
test.convert_data.frame_to_xts_i1 <- function() {
  checkIdentical(sample.xts[1,],as.xts(sample.data.frame)[1,])
}
test.convert_data.frame_to_xts_i1j1 <- function() {
  checkIdentical(sample.xts[1,1],as.xts(sample.data.frame)[1,1])
}
test.data.frame_reclass <- function() {
  checkIdentical(sample.data.frame,reclass(try.xts(sample.data.frame)))
}
test.data.frame_reclass_subset_reclass_j1 <- function() {
  checkIdentical(sample.data.frame[,1],reclass(try.xts(sample.data.frame))[,1])
}

# subsetting to 1 col converts to simple numeric - can't successfully handle
test.data.frame_reclass_subset_as.xts_j1 <- function() {
  checkIdentical(sample.data.frame[,1,drop=FALSE],reclass(try.xts(sample.data.frame)[,1]))
}
test.data.frame_reclass_subset_data.frame_j1 <- function() {
  # subsetting results in a vector, so can't be converted to xts
  checkException(try.xts(sample.data.frame[,1]))
}

# check for as.xts.data.frame when order.by is specified
test.convert_data.frame_to_xts_order.by_POSIXlt <- function() {
  orderby = as.POSIXlt(rownames(sample.data.frame))
  x <- as.xts(sample.data.frame, order.by = orderby)
  # tz = "" by default for as.POSIXlt.POSIXct
  y <- xts(coredata(sample.xts), as.POSIXlt(index(sample.xts)))
  checkIdentical(y, x)
}
test.convert_data.frame_to_xts_order.by_POSIXct <- function() {
  orderby = as.POSIXct(rownames(sample.data.frame))
  x <- as.xts(sample.data.frame, order.by = orderby)
  checkIdentical(sample.xts, x)
}
test.convert_data.frame_to_xts_order.by_Date <- function() {
  # tz = "UTC" by default for as.Date.POSIXct (y), but
  # tz = "" by default for as.Date.character (orderby)
  orderby = as.Date(rownames(sample.data.frame))
  x <- as.xts(sample.data.frame, order.by = orderby)
  y <- xts(coredata(sample.xts), as.Date(index(sample.xts), tz = ""))
  checkIdentical(y, x)
}
