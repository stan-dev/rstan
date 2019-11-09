na <- NA_integer_

# vector with even length, odd length

# no/yes result (potential infinite loop)
# https://www.topcoder.com/community/data-science/data-science-tutorials/binary-search/
test.integer_predicate_no_yes_stops <- function() {
  ans <- 2L
  ivec <- 3:4
  ikey <- ivec[ans]
  checkIdentical(ans, xts:::binsearch(ikey, ivec, TRUE))
  checkIdentical(ans, xts:::binsearch(ikey, ivec, FALSE))
}

# small steps between vector elements (test that we actually stop)
test.double_with_small_delta_stops <- function() {
  ans <- 10L
  dvec <- 1 + (-10:10 / 1e8)
  dkey <- dvec[ans]
  checkIdentical(ans, xts:::binsearch(dkey, dvec, TRUE))
  checkIdentical(ans, xts:::binsearch(dkey, dvec, FALSE))
}

test.find_first_zero_even_length <- function() {
  ivec <- sort(c(0L, -3:5L))
  dvec <- ivec * 1.0
  checkIdentical(4L, xts:::binsearch(0L,  ivec, TRUE))
  checkIdentical(4L, xts:::binsearch(0.0, dvec, TRUE))
}

test.find_last_zero_even_length <- function() {
  ivec <- sort(c(0L, -3:5L))
  dvec <- ivec * 1.0
  checkIdentical(5L, xts:::binsearch(0L,  ivec, FALSE))
  checkIdentical(5L, xts:::binsearch(0.0, dvec, FALSE))
}

test.find_first_zero_odd_length <- function() {
  ivec <- sort(c(0L, -3:5L))
  dvec <- ivec * 1.0
  checkIdentical(4L, xts:::binsearch(0L,  ivec, TRUE))
  checkIdentical(4L, xts:::binsearch(0.0, dvec, TRUE))
}

test.find_last_zero_odd_length <- function() {
  ivec <- sort(c(0L, -3:5L))
  dvec <- ivec * 1.0
  checkIdentical(5L, xts:::binsearch(0L,  ivec, FALSE))
  checkIdentical(5L, xts:::binsearch(0.0, dvec, FALSE))
}

# key is outside of vector
test.key_less_than_min <- function() {
  ivec <- 1:6
  checkIdentical(1L, xts:::binsearch(-9L, ivec, TRUE))
  checkIdentical(na, xts:::binsearch(-9L, ivec, FALSE))
  dvec <- ivec * 1.0
  checkIdentical(1L, xts:::binsearch(-9,  dvec, TRUE))
  checkIdentical(na, xts:::binsearch(-9,  dvec, FALSE))
}

test.key_greater_than_max <- function() {
  ivec <- 1:6
  checkIdentical(na, xts:::binsearch( 9L, ivec, TRUE))
  checkIdentical(6L, xts:::binsearch( 9L, ivec, FALSE))
  dvec <- ivec * 1.0
  checkIdentical(na, xts:::binsearch( 9,  dvec, TRUE))
  checkIdentical(6L, xts:::binsearch( 9,  dvec, FALSE))
}

# key is NA
test.key_is_NA <- function() {
  ivec <- 1:6
  ikey <- NA_integer_
  checkIdentical(na, xts:::binsearch(ikey, ivec, TRUE))
  checkIdentical(na, xts:::binsearch(ikey, ivec, FALSE))

  dvec <- ivec * 1.0
  dkey <- NA_real_
  checkIdentical(na, xts:::binsearch(dkey, dvec, TRUE))
  checkIdentical(na, xts:::binsearch(dkey, dvec, FALSE))
}

# key is zero-length
test.key_is_zero_length <- function() {
  # have empty key return NA
  ivec <- 1:6
  ikey <- integer()
  checkIdentical(na, xts:::binsearch(ikey, ivec, TRUE))
  checkIdentical(na, xts:::binsearch(ikey, ivec, FALSE))

  dvec <- ivec * 1.0
  dkey <- double()
  checkIdentical(na, xts:::binsearch(dkey, dvec, TRUE))
  checkIdentical(na, xts:::binsearch(dkey, dvec, FALSE))
}

# vec is zero-length
test.vec_is_zero_length <- function() {
  # have empty vector return NA
  ivec <- integer()
  ikey <- 0L
  checkIdentical(na, xts:::binsearch(ikey, ivec, TRUE))
  checkIdentical(na, xts:::binsearch(ikey, ivec, FALSE))

  dvec <- double()
  dkey <- 0.0
  checkIdentical(na, xts:::binsearch(dkey, dvec, TRUE))
  checkIdentical(na, xts:::binsearch(dkey, dvec, FALSE))
}

