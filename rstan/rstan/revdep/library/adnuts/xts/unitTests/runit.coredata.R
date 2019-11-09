test.coredata_vector <- function() {
  x <- xts(1, as.Date("2018-03-02"))
  z <- as.zoo(x)

  checkIdentical(coredata(x), coredata(z))
}

test.coredata_named_vector <- function() {
  x <- xts(c(hello = 1), as.Date("2018-03-02"))
  z <- as.zoo(x)

  checkIdentical(coredata(x), coredata(z))
}

test.coredata_matrix <- function() {
  x <- xts(cbind(1, 9), as.Date("2018-03-02"))
  z <- as.zoo(x)

  checkIdentical(coredata(x), coredata(z))
}

test.coredata_named_matrix <- function() {
  x <- xts(cbind(hello = 1, world = 9), as.Date("2018-03-02"))
  z <- as.zoo(x)

  checkIdentical(coredata(x), coredata(z))
}

test.coredata_data.frame <- function() {
  x <- xts(data.frame(hello = 1, world = 9), as.Date("2018-03-02"))
  z <- as.zoo(x)

  checkIdentical(coredata(x), coredata(z))
}

test.coredata_ts <- function() {
  x <- xts(ts(1), as.Date("2018-03-02"))
  z <- as.zoo(x)

  checkIdentical(coredata(x), coredata(z))
}

# empty objects
test.coredata_empty <- function() {
  x <- xts(, as.Date("2018-03-02"))
  z <- as.zoo(x)

  checkIdentical(coredata(x), coredata(z))
}

test.coredata_empty_dim <- function() {
  x <- xts(cbind(1, 9), as.Date("2018-03-02"))
  z <- as.zoo(x)
  x0 <- x[0,]
  z0 <- z[0,]

  checkIdentical(coredata(x0), coredata(z0))
}

test.coredata_empty_dim_dimnames <- function() {
  x <- xts(cbind(hello = 1, world = 9), as.Date("2018-03-02"))
  z <- as.zoo(x)
  x0 <- x[0,]
  z0 <- z[0,]

  checkIdentical(coredata(x0), coredata(z0))
}

