get_accumulator <- function() {
  return(.Call('get_accumulator_'))
}

check_accumulator <- function(accumulator = ACC) {
  return(.Call('check_accumulator_', accumulator))
}



