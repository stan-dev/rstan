data(sample_matrix)
sample.ts1 <- ts(sample_matrix,start=as.numeric(as.Date(rownames(sample_matrix)[1])))
sample.xts.ts1 <- as.xts(sample.ts1)


test.convert_ts_to_xts <- function() {
  checkIdentical(sample.xts.ts1,as.xts(sample.ts1))
}
test.convert_ts_to_xts_j1 <- function() {
  checkIdentical(sample.xts.ts1[,1],as.xts(sample.ts1)[,1])
}
test.convert_ts_to_xts_i1 <- function() {
  checkIdentical(sample.xts.ts1[1,],as.xts(sample.ts1)[1,])
}
test.convert_ts_to_xts_i1j1 <- function() {
  checkIdentical(sample.xts.ts1[1,1],as.xts(sample.ts1)[1,1])
}


test.ts_reclass <- function() {
  checkIdentical(sample.ts1,reclass(try.xts(sample.ts1)))
}
test.ts_reclass_subset_reclass_j1 <- function() {
  checkIdentical(sample.ts1[,1],reclass(try.xts(sample.ts1))[,1])
}
test.ts_reclass_subset_as.xts_j1 <- function() {
  checkIdentical(sample.ts1[,1],reclass(try.xts(sample.ts1)[,1]))
}
test.ts_reclass_subset_ts_j1 <- function() {
  checkIdentical(sample.ts1[,1],reclass(try.xts(sample.ts1[,1])))
}

# quarterly series
sample.ts4 <- ts(sample_matrix,start=1960,frequency=4)
sample.xts.ts4 <- as.xts(sample.ts4)


test.convert_ts4_to_xts <- function() {
  checkIdentical(sample.xts.ts4,as.xts(sample.ts4))
}
test.convert_ts4_to_xts_j1 <- function() {
  checkIdentical(sample.xts.ts4[,1],as.xts(sample.ts4)[,1])
}
test.convert_ts4_to_xts_i1 <- function() {
  checkIdentical(sample.xts.ts4[1,],as.xts(sample.ts4)[1,])
}
test.convert_ts4_to_xts_i1j1 <- function() {
  checkIdentical(sample.xts.ts4[1,1],as.xts(sample.ts4)[1,1])
}


test.ts4_reclass <- function() {
  checkIdentical(sample.ts4,reclass(try.xts(sample.ts4)))
}
test.ts4_reclass_subset_reclass_j1 <- function() {
  checkIdentical(sample.ts4[,1],reclass(try.xts(sample.ts4))[,1])
}
test.ts4_reclass_subset_as.xts_j1 <- function() {
  checkIdentical(sample.ts4[,1],reclass(try.xts(sample.ts4)[,1]))
}
test.ts4_reclass_subset_ts_j1 <- function() {
  checkIdentical(sample.ts4[,1],reclass(try.xts(sample.ts4[,1])))
}

# monthly series
sample.ts12 <- ts(sample_matrix,start=1990,frequency=12)
sample.xts.ts12 <- as.xts(sample.ts12)


test.convert_ts12_to_xts <- function() {
  checkIdentical(sample.xts.ts12,as.xts(sample.ts12))
}
test.convert_ts12_to_xts_j1 <- function() {
  checkIdentical(sample.xts.ts12[,1],as.xts(sample.ts12)[,1])
}
test.convert_ts12_to_xts_i1 <- function() {
  checkIdentical(sample.xts.ts12[1,],as.xts(sample.ts12)[1,])
}
test.convert_ts12_to_xts_i1j1 <- function() {
  checkIdentical(sample.xts.ts12[1,1],as.xts(sample.ts12)[1,1])
}


test.ts12_reclass <- function() {
  checkIdentical(sample.ts12,reclass(try.xts(sample.ts12)))
}
test.ts12_reclass_subset_reclass_j1 <- function() {
  checkIdentical(sample.ts12[,1],reclass(try.xts(sample.ts12))[,1])
}
test.ts12_reclass_subset_as.xts_j1 <- function() {
  checkIdentical(sample.ts12[,1],reclass(try.xts(sample.ts12)[,1]))
}
test.ts12_reclass_subset_ts_j1 <- function() {
  checkIdentical(sample.ts12[,1],reclass(try.xts(sample.ts12[,1])))
}
