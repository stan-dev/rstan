get_accumulator <- function(start=0L) {
  if (start != 0L) {
    if (length(start) != 1) 
      stop("Start value must be a length-1 integer vector.")
  }
  return(.Call('get_accumulator_', start))
}

check_accumulator <- function(accumulator = ACC) {
  return(.Call('check_accumulator_', accumulator))
}



