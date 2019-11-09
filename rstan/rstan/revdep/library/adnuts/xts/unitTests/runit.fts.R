if (requireNamespace("fts", quietly = TRUE)) {

data(sample_matrix)
sample.fts1 <- ts(sample_matrix,start=as.Date(rownames(sample_matrix)[1]))
sample.fts1 <- fts::fts(index(sample.fts1), sample_matrix)
sample.xts.fts1 <- as.xts(sample.fts1)


test.convert_fts_to_xts <- function() {
  checkIdentical(sample.xts.fts1,as.xts(sample.fts1))
}
test.convert_fts_to_xts_j1 <- function() {
  checkIdentical(sample.xts.fts1[,1],as.xts(sample.fts1)[,1])
}
test.convert_fts_to_xts_i1 <- function() {
  checkIdentical(sample.xts.fts1[1,],as.xts(sample.fts1)[1,])
}
test.convert_fts_to_xts_i1j1 <- function() {
  checkIdentical(sample.xts.fts1[1,1],as.xts(sample.fts1)[1,1])
}


test.fts_reclass <- function() {
  checkIdentical(sample.fts1,reclass(try.xts(sample.fts1)))
}
test.fts_reclass_subset_reclass_j1 <- function() {
  checkIdentical(sample.fts1[,1],reclass(try.xts(sample.fts1))[,1])
}
test.fts_reclass_subset_as.xts_j1 <- function() {
  checkIdentical(sample.fts1[,1],reclass(try.xts(sample.fts1)[,1]))
}
test.fts_reclass_subset_fts_j1 <- function() {
  checkIdentical(sample.fts1[,1],reclass(try.xts(sample.fts1[,1])))
}

# quarterly series
sample.fts4 <- ts(sample_matrix,start=1960,frequency=4)
sample.fts4 <- fts::fts(index(sample.fts4), sample_matrix)
sample.xts.fts4 <- as.xts(sample.fts4)


test.convert_fts4_to_xts <- function() {
  checkIdentical(sample.xts.fts4,as.xts(sample.fts4))
}
test.convert_fts4_to_xts_j1 <- function() {
  checkIdentical(sample.xts.fts4[,1],as.xts(sample.fts4)[,1])
}
test.convert_fts4_to_xts_i1 <- function() {
  checkIdentical(sample.xts.fts4[1,],as.xts(sample.fts4)[1,])
}
test.convert_fts4_to_xts_i1j1 <- function() {
  checkIdentical(sample.xts.fts4[1,1],as.xts(sample.fts4)[1,1])
}


test.fts4_reclass <- function() {
  checkIdentical(sample.fts4,reclass(try.xts(sample.fts4)))
}
test.fts4_reclass_subset_reclass_j1 <- function() {
  checkIdentical(sample.fts4[,1],reclass(try.xts(sample.fts4))[,1])
}
test.fts4_reclass_subset_as.xts_j1 <- function() {
  checkIdentical(sample.fts4[,1],reclass(try.xts(sample.fts4)[,1]))
}
test.fts4_reclass_subset_fts_j1 <- function() {
  checkIdentical(sample.fts4[,1],reclass(try.xts(sample.fts4[,1])))
}

# monthly series
sample.fts12 <- ts(sample_matrix,start=1990,frequency=12)
sample.fts12 <- fts::fts(index(sample.fts12), sample_matrix)
sample.xts.fts12 <- as.xts(sample.fts12)


test.convert_fts12_to_xts <- function() {
  checkIdentical(sample.xts.fts12,as.xts(sample.fts12))
}
test.convert_fts12_to_xts_j1 <- function() {
  checkIdentical(sample.xts.fts12[,1],as.xts(sample.fts12)[,1])
}
test.convert_fts12_to_xts_i1 <- function() {
  checkIdentical(sample.xts.fts12[1,],as.xts(sample.fts12)[1,])
}
test.convert_fts12_to_xts_i1j1 <- function() {
  checkIdentical(sample.xts.fts12[1,1],as.xts(sample.fts12)[1,1])
}


test.fts12_reclass <- function() {
  checkIdentical(sample.fts12,reclass(try.xts(sample.fts12)))
}
test.fts12_reclass_subset_reclass_j1 <- function() {
  checkIdentical(sample.fts12[,1],reclass(try.xts(sample.fts12))[,1])
}
test.fts12_reclass_subset_as.xfts_j1 <- function() {
  checkIdentical(sample.fts12[,1],reclass(try.xts(sample.fts12)[,1]))
}
test.fts12_reclass_subset_fts_j1 <- function() {
  checkIdentical(sample.fts12[,1],reclass(try.xts(sample.fts12[,1])))
}

}  # requireNamespace
