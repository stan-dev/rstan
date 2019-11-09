test.merge_empty_xts_with_2_scalars <- function() {
  m1 <- merge(xts(), 1, 1)
  m2 <- merge(merge(xts(), 1), 1)
  checkIdentical(m1, m2)
}

test.merge_more_than_2_zero_width_objects <- function() {
  zw1 <- xts()
  zw2 <- xts()
  zw3 <- xts()
  m1 <- merge(zw1, zw2, zw3)
  checkIdentical(m1, zw1)
}

# Tests for NA in index. Construct xts object using structure() because
# xts constructors should not allow users to create objects with NA in
# the index
indexHasNA_dbl <-
structure(1:5, .Dim = c(5L, 1L),
          index = structure(c(1, 2, 3, 4, NA), tzone = "",
                            tclass = c("POSIXct", "POSIXt")),
          .indexCLASS = c("POSIXct", "POSIXt"),
          .indexTZ = "", tclass = c("POSIXct", "POSIXt"),
          tzone = "", class = c("xts", "zoo"))

indexHasNA_int <-
structure(1:5, .Dim = c(5L, 1L),
          index = structure(c(1L, 2L, 3L, 4L, NA), tzone = "",
                            tclass = c("POSIXct", "POSIXt")),
          .indexCLASS = c("POSIXct", "POSIXt"),
          .indexTZ = "", tclass = c("POSIXct", "POSIXt"),
          tzone = "", class = c("xts", "zoo"))

test.merge_index_contains_NA_integer <- function() {
  checkException(merge(indexHasNA_int, indexHasNA_int), silent = TRUE)
}

test.merge_index_contains_NA_double <- function() {
  checkException(merge(indexHasNA_dbl, indexHasNA_dbl), silent = TRUE)
}

test.merge_index_contains_NaN <- function() {
  x <- indexHasNA_dbl
  idx <- attr(x, "index")
  idx[length(idx)] <- NaN
  attr(x, "index") <- idx
  checkException(merge(x, x), silent = TRUE)
}

test.merge_index_contains_Inf <- function() {
  x <- indexHasNA_dbl
  idx <- attr(x, "index")
  idx[length(idx)] <- Inf
  attr(x, "index") <- idx
  checkException(merge(x, x), silent = TRUE)

  idx <- rev(idx)
  idx[1L] <- -Inf
  attr(x, "index") <- idx
  checkException(merge(x, x), silent = TRUE)
}
# /end Tests for NA in index

# zero-length fill argument
test.merge_fill_NULL <- function() {
  x1 <- .xts(1, 1)
  x2 <- .xts(2, 2)
  x <- merge(x1, x2, fill = NULL)

  out <- .xts(matrix(c(1, NA, NA, 2), 2), c(1,2))
  colnames(out) <- c("x1", "x2")
  checkIdentical(x, out)
}

test.merge_fill_zero_length <- function() {
  x1 <- .xts(1, 1)
  x2 <- .xts(2, 2)
  x <- merge(x1, x2, fill = numeric())

  out <- .xts(matrix(c(1, NA, NA, 2), 2), c(1,2))
  colnames(out) <- c("x1", "x2")
  checkIdentical(x, out)
}

test.merge_with_zero_width_returns_original_type <- function() {
  M1 <- .xts(1:3, 1:3, dimnames = list(NULL, "m1"))
  for (m in c("double", "integer", "logical", "character")) {
    m1 <- M1
    storage.mode(m1) <- m
    e1 <- .xts(,1:3)
    m2 <- merge(m1, e1)
    checkIdentical(m1, m2)
  }
}

test.n_way_merge_on_all_types <- function() {
  D1 <- as.Date("2018-01-03")-2:0
  M1 <- xts(1:3, D1, dimnames = list(NULL, "m"))
  M3 <- xts(cbind(1:3, 1:3, 1:3), D1,
            dimnames = list(NULL, c("m", "m.1", "m.2")))

  types <- c("double", "integer", "logical", "character", "complex")
  for (type in types) {
    m1 <- M1
    m3 <- M3
    storage.mode(m1) <- storage.mode(m3) <- type
    m <- merge(m1, m1, m1)
    checkIdentical(m, m3)
  }
}

test.shorter_colnames_for_unnamed_args <- function() {
  X <- .xts(rnorm(10, 10), 1:10)

  types <- c("double", "integer", "logical", "character", "complex")
  for (type in types) {
    x <- X
    storage.mode(x) <- type
    mx <- do.call(merge, list(x, x))
    checkTrue(all(nchar(colnames(mx)) < 200), type)
  }
}
