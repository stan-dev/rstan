get_rng <- function(seed=0L) {
  warning("The 'get_rng' function is deprecated and does not do anything useful.",
          "Exposed Stan functions now use R's pseudo-random number generator.",
          "Call 'set.seed' to achieve reproducibility.")
  if (seed != 0L) {
    if (length(seed) != 1) 
      stop("Seed must be a length-1 integer vector.")
  }
  return(.Call('get_rng_', seed))
}



