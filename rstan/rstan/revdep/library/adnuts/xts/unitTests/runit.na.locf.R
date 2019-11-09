XDAT <- .xts(c(1, NA, 3, 4, 5, 6), c(0, 4, 10, 19, 24, 29))
XIDX <- .xts(rep(0, 5), c(5, 10, 20, 25, 28))
MODES <- c("double", "integer", "character", "logical")

# na.locf.xts() on a univariate xts object
test.nalocf <- function() {
  for (m in MODES) {
    xdat <- XDAT
    storage.mode(xdat) <- m
    zdat <- as.zoo(xdat)

    x <- na.locf(xdat)
    z <- na.locf(zdat)
    #checkIdentical(x, as.xts(z))  # FALSE (attribute order differs)
    checkEquals(x, as.xts(z), check.attributes = TRUE)
  }
}

test.nalocf_leading_NA <- function() {
  for (m in MODES) {
    xdat <- XDAT
    storage.mode(xdat) <- m
    zdat <- as.zoo(xdat)

    xdat[1] <- NA
    zdat[1] <- NA

    x <- na.locf(xdat, na.rm = TRUE)
    z <- na.locf(zdat, na.rm = TRUE)
    checkEquals(x, as.xts(z), check.attributes = TRUE)

    x <- na.locf(xdat, na.rm = FALSE)
    z <- na.locf(zdat, na.rm = FALSE)
    checkEquals(x, as.xts(z), check.attributes = TRUE)
  }
}

test.nalocf_fromLast <- function() {
  for (m in MODES) {
    xdat <- XDAT
    storage.mode(xdat) <- m
    zdat <- as.zoo(xdat)

    x <- na.locf(xdat, fromLast = TRUE)
    z <- na.locf(zdat, fromLast = TRUE)
    checkEquals(x, as.xts(z), check.attributes = TRUE)
  }
}

test.nalocf_x <- function() {
  for (m in MODES) {
    xdat <- XDAT
    xidx <- XIDX
    storage.mode(xdat) <- storage.mode(xidx) <- m
    zdat <- as.zoo(xdat)
    zidx <- as.zoo(xidx)

    xidx <- rbind(xidx, .xts(0, 30))
    zidx <- as.zoo(xidx)

    x <- na.locf(xdat, x = index(xidx))
    z <- na.locf(zdat, x = index(zidx))
    checkEquals(x, as.xts(z), check.attributes = TRUE)
  }
}

test.nalocf_xout <- function() {
  for (m in MODES) {
    xdat <- XDAT
    xidx <- XIDX
    storage.mode(xdat) <- storage.mode(xidx) <- m
    zdat <- as.zoo(xdat)
    zidx <- as.zoo(xidx)

    x <- na.locf(xdat, xout = index(xidx))
    z <- na.locf(zdat, xout = index(zidx))
    checkEquals(x, as.xts(z), check.attributes = TRUE)
  }
}

# na.locf.xts() on a multivariate xts object
XDAT2 <- merge(one = XDAT, two = XDAT)

test.nalocf_by_column <- function() {
  for (m in MODES) {
    xdat <- XDAT2
    storage.mode(xdat) <- m
    zdat <- as.zoo(xdat)

    x <- na.locf(xdat)
    z <- na.locf(zdat)
    checkEquals(x, as.xts(z), check.attributes = TRUE)
  }
}

test.nalocf_by_column_leading_NA <- function() {
  for (m in MODES) {
    xdat <- XDAT2
    storage.mode(xdat) <- m
    zdat <- as.zoo(xdat)

    xdat[1] <- NA
    zdat[1] <- NA

if (FALSE) {
# bug w/zoo causes this to fail
# zoo:::na.locf.default() does not remove the first row
    x <- na.locf(xdat, na.rm = TRUE)
    z <- na.locf(zdat, na.rm = TRUE)
    checkEquals(x, as.xts(z), check.attributes = TRUE)
}

    x <- na.locf(xdat, na.rm = FALSE)
    z <- na.locf(zdat, na.rm = FALSE)
    checkEquals(x, as.xts(z), check.attributes = TRUE)
  }
}

test.nalocf_by_column_fromLast <- function() {
  for (m in MODES) {
    xdat <- XDAT2
    storage.mode(xdat) <- m
    zdat <- as.zoo(xdat)

    x <- na.locf(xdat, fromLast = TRUE)
    z <- na.locf(zdat, fromLast = TRUE)
    checkEquals(x, as.xts(z), check.attributes = TRUE)
  }
}

test.nalocf_by_column_x <- function() {
  for (m in MODES) {
    xdat <- XDAT2
    xidx <- XIDX
    storage.mode(xdat) <- storage.mode(xidx) <- m
    zdat <- as.zoo(xdat)
    zidx <- as.zoo(xidx)

    xidx <- rbind(xidx, .xts(0, 30))
    zidx <- as.zoo(xidx)

    x <- na.locf(xdat, x = index(xidx))
    z <- na.locf(zdat, x = index(zidx))
    checkEquals(x, as.xts(z), check.attributes = TRUE)
  }
}

test.nalocf_by_column_xout <- function() {
  for (m in MODES) {
    xdat <- XDAT2
    xidx <- XIDX
    storage.mode(xdat) <- storage.mode(xidx) <- m
    zdat <- as.zoo(xdat)
    zidx <- as.zoo(xidx)

    x <- na.locf(xdat, xout = index(xidx))
    z <- na.locf(zdat, xout = index(zidx))
    checkEquals(x, as.xts(z), check.attributes = TRUE)
  }
}

test.nalocf_by_column_1NA <- function() {
  narow <- 1L
  for (m in MODES) {
    xdrow <- XDAT2[narow,]
    xdat <- XDAT2 * NA
    xdat[narow,] <- xdrow
    storage.mode(xdat) <- m
    zdat <- as.zoo(xdat)

    x <- na.locf(xdat)
    z <- na.locf(zdat)
    checkEquals(x, as.xts(z), check.attributes = TRUE)
  }
}

test.nalocf_by_column_1NA_fromLast <- function() {
  narow <- nrow(XDAT2)
  for (m in MODES) {
    xdrow <- XDAT2[narow,]
    xdat <- XDAT2 * NA
    xdat[narow,] <- xdrow
    storage.mode(xdat) <- m
    zdat <- as.zoo(xdat)

    x <- na.locf(xdat, fromLast = TRUE)
    z <- na.locf(zdat, fromLast = TRUE)
    checkEquals(x, as.xts(z), check.attributes = TRUE)
  }
}

test.nalocf_first_column_all_NA <- function() {
  nacol <- 1L
  for (m in MODES) {
    xdat <- XDAT2
    xdat[,nacol] <- xdat[,nacol] * NA
    storage.mode(xdat) <- m
    zdat <- as.zoo(xdat)

    x <- na.locf(xdat)
    z <- na.locf(zdat)
    checkEquals(x, as.xts(z), check.attributes = TRUE)
  }
}

test.nalocf_last_column_all_NA <- function() {
  nacol <- NCOL(XDAT2)
  for (m in MODES) {
    xdat <- XDAT2
    xdat[,nacol] <- xdat[,nacol] * NA
    storage.mode(xdat) <- m
    zdat <- as.zoo(xdat)

    x <- na.locf(xdat)
    z <- na.locf(zdat)
    checkEquals(x, as.xts(z), check.attributes = TRUE)
  }
}
