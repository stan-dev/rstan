XDAT <- .xts(c(1, NA, 3, 4, 5, 6), c(0, 4, 10, 19, 24, 29))
XIDX <- .xts(rep(0, 5), c(5, 10, 20, 25, 28))
MODES <- c("double", "integer", "character", "logical")

test.naomit <- function() {
  for (m in MODES) {
    xdat <- XDAT
    xidx <- XIDX
    storage.mode(xdat) <- storage.mode(xidx) <- m
    zdat <- as.zoo(xdat)
    zidx <- as.zoo(xidx)

    x <- na.omit(xdat)
    z <- na.omit(zdat)
    # na.omit.xts adds "index" attribute to the "na.action" attribute
    attr(attr(x, "na.action"), "index") <- NULL
    #checkIdentical(x, as.xts(z))  # FALSE (attribute order differs)
    checkEquals(x, as.xts(z), check.attributes = TRUE)
  }
}

test.naomit_by_column <- function() {
  for (m in MODES) {
    xdat <- XDAT
    xidx <- XIDX
    storage.mode(xdat) <- storage.mode(xidx) <- m
    zdat <- as.zoo(xdat)
    zidx <- as.zoo(xidx)

    x <- na.omit(merge(one = xdat, two = xdat))
    z <- na.omit(merge(one = zdat, two = zdat))
    # na.omit.xts adds "index" attribute to the "na.action" attribute
    attr(attr(x, "na.action"), "index") <- NULL
    checkEquals(x, as.xts(z), check.attributes = TRUE)
  }
}
