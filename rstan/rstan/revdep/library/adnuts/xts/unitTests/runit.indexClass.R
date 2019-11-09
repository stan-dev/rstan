withPackage <- function(pkg, expr) {
  if (require(pkg, quietly = TRUE, character.only = TRUE)) {
    pkgpkg <- paste0("package:", pkg)
    on.exit(detach(pkgpkg, character.only = TRUE), add = TRUE)
    eval(expr)
  }
}

data(sample_matrix)

convert_xts <- as.xts(sample_matrix) # indexClass defaults to POSIXct

checkUTCindexTZ <- function(x) {
  if(any(indexClass(x) %in% xts:::.classesWithoutTZ)) {
    checkIdentical("UTC", attr(.index(x),'tzone'))
  }
}

# convert from 'POSIXct'
test.convert_POSIXct2Date <- function() {
  x <- convert_xts
  indexClass(x) <- 'Date'
  checkTrue(inherits(index(x),'Date'))
  checkUTCindexTZ(x)
}
test.convert_POSIXct2chron <- function() {
  withPackage("chron", {
  x <- convert_xts
  indexClass(x) <- 'chron'
  checkTrue(inherits(index(x),'dates'))
  checkUTCindexTZ(x)
  })
}
test.convert_POSIXct2yearmon <- function() {
  x <- convert_xts
  indexClass(x) <- 'yearmon'
  checkTrue(inherits(index(x),'yearmon'))
  checkUTCindexTZ(x)
}
test.convert_POSIXct2yearqtr <- function() {
  x <- convert_xts
  indexClass(x) <- 'yearqtr'
  checkTrue(inherits(index(x),'yearqtr'))
  checkUTCindexTZ(x)
}
test.convert_POSIXct2timeDate <- function() {
  withPackage("timeDate", {
  x <- convert_xts
  indexClass(x) <- 'timeDate'
  checkTrue(inherits(index(x),'timeDate'))
  checkUTCindexTZ(x)
  })
}
test.convert_POSIXct2POSIXct <- function() {
  x <- convert_xts
  indexClass(x) <- 'POSIXct'
  checkTrue(inherits(index(x),'POSIXct'))
  checkUTCindexTZ(x)
}

# Convert from 'Date'
indexClass(convert_xts) <- 'Date'

test.convert_Date2Date <- function() {
  x <- convert_xts
  indexClass(x) <- 'Date'
  checkTrue(inherits(index(x),'Date'))
  checkUTCindexTZ(x)
}
test.convert_Date2chron <- function() {
  withPackage("chron", {
  x <- convert_xts
  indexClass(x) <- 'chron'
  checkTrue(inherits(index(x),'dates'))
  checkUTCindexTZ(x)
  })
}
test.convert_Date2yearmon <- function() {
  x <- convert_xts
  indexClass(x) <- 'yearmon'
  checkTrue(inherits(index(x),'yearmon'))
  checkUTCindexTZ(x)
}
test.convert_Date2yearqtr <- function() {
  x <- convert_xts
  indexClass(x) <- 'yearqtr'
  checkTrue(inherits(index(x),'yearqtr'))
  checkUTCindexTZ(x)
}
test.convert_Date2timeDate <- function() {
  withPackage("timeDate", {
  x <- convert_xts
  indexClass(x) <- 'timeDate'
  checkTrue(inherits(index(x),'timeDate'))
  checkUTCindexTZ(x)
  })
}
test.convert_Date2POSIXct <- function() {
  x <- convert_xts
  indexClass(x) <- 'POSIXct'
  checkTrue(inherits(index(x),'POSIXct'))
  checkUTCindexTZ(x)
}

# Convert from 'chron'
if (requireNamespace("chron", quietly = TRUE)) {
indexClass(convert_xts) <- 'chron'

test.convert_chron2Date <- function() {
  x <- convert_xts
  indexClass(x) <- 'Date'
  checkTrue(inherits(index(x),'Date'))
  checkUTCindexTZ(x)
}
test.convert_chron2chron <- function() {
  withPackage("chron", {
  x <- convert_xts
  indexClass(x) <- 'chron'
  checkTrue(inherits(index(x),'dates'))
  checkUTCindexTZ(x)
  })
}
test.convert_chron2yearmon <- function() {
  x <- convert_xts
  indexClass(x) <- 'yearmon'
  checkTrue(inherits(index(x),'yearmon'))
  checkUTCindexTZ(x)
}
test.convert_chron2yearqtr <- function() {
  x <- convert_xts
  indexClass(x) <- 'yearqtr'
  checkTrue(inherits(index(x),'yearqtr'))
  checkUTCindexTZ(x)
}
test.convert_chron2timeDate <- function() {
  withPackage("timeDate", {
  x <- convert_xts
  indexClass(x) <- 'timeDate'
  checkTrue(inherits(index(x),'timeDate'))
  checkUTCindexTZ(x)
  })
}
test.convert_chron2POSIXct <- function() {
  x <- convert_xts
  indexClass(x) <- 'POSIXct'
  checkTrue(inherits(index(x),'POSIXct'))
  checkUTCindexTZ(x)
}
}

# Convert from 'yearmon'
indexClass(convert_xts) <- 'yearmon'

test.convert_yearmon2Date <- function() {
  x <- convert_xts
  indexClass(x) <- 'Date'
  checkTrue(inherits(index(x),'Date'))
  checkUTCindexTZ(x)
}
test.convert_yearmon2chron <- function() {
  withPackage("chron", {
  x <- convert_xts
  indexClass(x) <- 'chron'
  checkTrue(inherits(index(x),'dates'))
  checkUTCindexTZ(x)
  })
}
test.convert_yearmon2yearmon <- function() {
  x <- convert_xts
  indexClass(x) <- 'yearmon'
  checkTrue(inherits(index(x),'yearmon'))
  checkUTCindexTZ(x)
}
test.convert_yearmon2yearqtr <- function() {
  x <- convert_xts
  indexClass(x) <- 'yearqtr'
  checkTrue(inherits(index(x),'yearqtr'))
  checkUTCindexTZ(x)
}
test.convert_yearmon2timeDate <- function() {
  withPackage("timeDate", {
  x <- convert_xts
  indexClass(x) <- 'timeDate'
  checkTrue(inherits(index(x),'timeDate'))
  checkUTCindexTZ(x)
  })
}
test.convert_yearmon2POSIXct <- function() {
  x <- convert_xts
  indexClass(x) <- 'POSIXct'
  checkTrue(inherits(index(x),'POSIXct'))
  checkUTCindexTZ(x)
}

# Convert from 'yearqtr'
indexClass(convert_xts) <- 'yearqtr'

test.convert_yearqtr2Date <- function() {
  x <- convert_xts
  indexClass(x) <- 'Date'
  checkTrue(inherits(index(x),'Date'))
  checkUTCindexTZ(x)
}
test.convert_yearqtr2chron <- function() {
  withPackage("chron", {
  x <- convert_xts
  indexClass(x) <- 'chron'
  checkTrue(inherits(index(x),'dates'))
  checkUTCindexTZ(x)
  })
}
test.convert_yearqtr2yearmon <- function() {
  x <- convert_xts
  indexClass(x) <- 'yearmon'
  checkTrue(inherits(index(x),'yearmon'))
  checkUTCindexTZ(x)
}
test.convert_yearqtr2yearqtr <- function() {
  x <- convert_xts
  indexClass(x) <- 'yearqtr'
  checkTrue(inherits(index(x),'yearqtr'))
  checkUTCindexTZ(x)
}
test.convert_yearqtr2timeDate <- function() {
  withPackage("timeDate", {
  x <- convert_xts
  indexClass(x) <- 'timeDate'
  checkTrue(inherits(index(x),'timeDate'))
  checkUTCindexTZ(x)
  })
}
test.convert_yearqtr2POSIXct <- function() {
  x <- convert_xts
  indexClass(x) <- 'POSIXct'
  checkTrue(inherits(index(x),'POSIXct'))
  checkUTCindexTZ(x)
}

# Convert from 'timeDate'
if (requireNamespace("timeDate", quietly = TRUE)) {
indexClass(convert_xts) <- 'timeDate'

test.convert_timeDate2Date <- function() {
  x <- convert_xts
  indexClass(x) <- 'Date'
  checkTrue(inherits(index(x),'Date'))
  checkUTCindexTZ(x)
}
test.convert_timeDate2chron <- function() {
  withPackage("chron", {
  x <- convert_xts
  indexClass(x) <- 'chron'
  checkTrue(inherits(index(x),'dates'))
  checkUTCindexTZ(x)
  })
}
test.convert_timeDate2yearmon <- function() {
  x <- convert_xts
  indexClass(x) <- 'yearmon'
  checkTrue(inherits(index(x),'yearmon'))
  checkUTCindexTZ(x)
}
test.convert_timeDate2yearqtr <- function() {
  x <- convert_xts
  indexClass(x) <- 'yearqtr'
  checkTrue(inherits(index(x),'yearqtr'))
  checkUTCindexTZ(x)
}
test.convert_timeDate2timeDate <- function() {
  withPackage("timeDate", {
  x <- convert_xts
  indexClass(x) <- 'timeDate'
  checkTrue(inherits(index(x),'timeDate'))
  checkUTCindexTZ(x)
  })
}
test.convert_timeDate2POSIXct <- function() {
  x <- convert_xts
  indexClass(x) <- 'POSIXct'
  checkTrue(inherits(index(x),'POSIXct'))
  checkUTCindexTZ(x)
}
}

# set index and ensure TZ = "UTC"
test.checkUTC_set_index2Date <- function() {
  x <- .xts(1:2, 1:2)
  d <- c("2007-01-02", "2007-01-03")
  index(x) <- as.Date(d)
  checkUTCindexTZ(x)
}
test.checkUTC_set_index2chron <- function() {
  withPackage("chron", {
  x <- .xts(1:2, 1:2)
  d <- c("2007-01-02", "2007-01-03")
  index(x) <- chron::dates(d, format="Y-m-d")
  checkUTCindexTZ(x)
  })
}
test.checkUTC_set_index2yearmon <- function() {
  x <- .xts(1:2, 1:2)
  d <- c("2007-01-02", "2007-01-03")
  index(x) <- as.yearmon(d)
  checkUTCindexTZ(x)
}
test.checkUTC_set_index2yearqtr <- function() {
  x <- .xts(1:2, 1:2)
  d <- c("2007-01-02", "2007-01-03")
  index(x) <- as.yearqtr(d)
  checkUTCindexTZ(x)
}

# error checking

test.indexClass_NULL <- function() {
  checkException((indexClass(convert_xts) <- NULL))
}
test.indexClass_empty_string <- function() {
  checkException((indexClass(convert_xts) <- ''))
}
test.indexClass_full_index <- function() {
  checkException((indexClass(convert_xts) <- index(convert_xts)))
}
test.indexClass_unquoted_symbol <- function() {
  withPackage("timeDate", {
  checkException((indexClass(convert_xts) <- timeDate))
  })
}
test.indexClass_missing_object <- function() {
  rm(Date)
  checkException((indexClass(convert_xts) <- Date))
}
