#sample.zoo <-
#structure(c(43.46, 43.3, 43.95, 43.89, 44.01, 43.96, 44.71, 45.02, 
#45.35, 45.09), .Names = c("264", "263", "262", "261", "260", 
#"259", "258", "257", "256", "255"), index = structure(c(13516, 
#13517, 13518, 13521, 13522, 13523, 13524, 13525, 13529, 13530
#), class = "Date"), class = "zoo")
#
#sample.xts <-
#structure(c(43.46, 43.3, 43.95, 43.89, 44.01, 43.96, 44.71, 45.02, 
#45.35, 45.09), index = structure(c(13516, 13517, 13518, 13521, 
#13522, 13523, 13524, 13525, 13529, 13530), class = "Date"), class = c("xts", 
#"zoo"), .CLASS = "zoo", .Dim = c(10L, 1L), .Dimnames = list(c("264", 
#"263", "262", "261", "260", "259", "258", "257", "256", "255"
#), NULL), .ROWNAMES = c("264","263","262", "261", "260", "259",
#"258", "257", "256", "255"))
#
data(sample_matrix)
sample.zoo <- zoo(sample_matrix,as.Date(rownames(sample_matrix)))
sample.xts <- as.xts(sample.zoo)

test.convert_zoo_to_xts <- function() {
  checkIdentical(sample.xts,as.xts(sample.zoo))
}
test.convert_zoo_to_xts_j1 <- function() {
  checkIdentical(sample.xts[,1],as.xts(sample.zoo)[,1])
}
test.convert_zoo_to_xts_i1 <- function() {
  checkIdentical(sample.xts[1,],as.xts(sample.zoo)[1,])
}
test.convert_zoo_to_xts_i1j1 <- function() {
  checkIdentical(sample.xts[1,1],as.xts(sample.zoo)[1,1])
}
test.zoo_reclass <- function() {
  DEACTIVATED("rownames are not kept yet in current xts-dev")
  checkIdentical(sample.zoo,reclass(try.xts(sample.zoo)))
}
test.zoo_reclass_subset_reclass_j1 <- function() {
  DEACTIVATED("rownames are not kept yet in current xts-dev")
  checkIdentical(sample.zoo[,1],reclass(try.xts(sample.zoo))[,1])
}
test.zoo_reclass_subset_as.xts_j1 <- function() {
  checkIdentical(sample.zoo[,1],reclass(try.xts(sample.zoo)[,1]))
}
test.zoo_reclass_subset_zoo_j1 <- function() {
  checkIdentical(sample.zoo[,1],reclass(try.xts(sample.zoo[,1])))
}
