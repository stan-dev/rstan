get_rng <- function(seed=0L) {
  if (!identical(seed, 0L)) {
    if (length(seed) != 1) 
      stop("Seed must be a length-1 integer vector.")
  }
  return(.Call('get_rng_', seed))
}



