if (requireNamespace("tseries", quietly = TRUE)) {
#
data(sample_matrix)

sample.irts <- tseries::irts(as.POSIXct(rownames(sample_matrix)),sample_matrix)
sample.irts.xts <- as.xts(sample.irts)

test.convert_irts_to_xts <- function() {
  checkIdentical(sample.irts.xts,as.xts(sample.irts))
}
test.convert_irts_to_xts_j1 <- function() {
  checkIdentical(sample.irts.xts[,1],as.xts(sample.irts)[,1])
}
test.convert_irts_to_xts_i1 <- function() {
  checkIdentical(sample.irts.xts[1,],as.xts(sample.irts)[1,])
}
test.convert_irts_to_xts_i1j1 <- function() {
  checkIdentical(sample.irts.xts[1,1],as.xts(sample.irts)[1,1])
}
test.irts_reclass <- function() {
  DEACTIVATED("irts forces rownames, xts disallows rownames. Unable to test")
  checkIdentical(sample.irts,reclass(try.xts(sample.irts)))
}
test.irts_reclass_subset_reclass_j1 <- function() {
  DEACTIVATED("irts forces rownames, xts disallows rownames. Unable to test")
  checkIdentical(sample.irts[,1],reclass(try.xts(sample.irts))[,1])
}
test.irts_reclass_subset_as.xts_j1 <- function() {
  DEACTIVATED("irts forces rownames, xts disallows rownames. Unable to test")
  checkIdentical(sample.irts[,1],reclass(try.xts(sample.irts)[,1]))
}
test.irts_reclass_subset_irts_j1 <- function() {
  DEACTIVATED("irts forces rownames, xts disallows rownames. Unable to test")
  checkIdentical(sample.irts[,1],reclass(try.xts(sample.irts[,1])))
}

}  # requireNamespace
