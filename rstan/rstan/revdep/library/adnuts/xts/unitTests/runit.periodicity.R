P <- structure(
  list(difftime = structure(0, units = "secs", class = "difftime"),
       frequency = 0,
       start = structure(.POSIXct(1, "UTC"), tclass = c("POSIXct", "POSIXt")),
       end = structure(.POSIXct(1, "UTC"), tclass = c("POSIXct", "POSIXt")),
       units = "secs",
       scale = "seconds",
       label = "second"),
  class = "periodicity")

test.periodicity_on_one_observation_warns <- function() {
  x <- xts(1, .POSIXct(1, "UTC"))
  p <- periodicity(x)
  checkIdentical(p, P)

  opt <- options(warn = 2)
  on.exit(options(warn = opt$warn))

  checkException(p <- periodicity(x))
}
test.periodicity_on_zero_observations_warns <- function() {
  x <- xts(, .POSIXct(numeric(0), "UTC"))
  p <- periodicity(x)
  P$start <- NA
  P$end <- NA
  checkIdentical(p, P)

  opt <- options(warn = 2)
  on.exit(options(warn = opt$warn))

  checkException(p <- periodicity(x))
}
