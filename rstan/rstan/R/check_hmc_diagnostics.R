check_hmc_diagnostics <- function(object) {
  cat("\nDivergences:\n")
    check_divergences(object)
  cat("\nTree depth:\n")
    check_treedepth(object)
  cat("\nEnergy:\n")
    check_energy(object)
}


# Check transitions that ended with a divergence
#
# @param object A stanfit object.
# @return nothing, just prints the number (and percentage) of iterations that 
#   ended with a divergence and, if any, suggests increasing adapt_delta.
#
check_divergences <- function(object) {
  divergent <- sampler_param_vector(object, "divergent__")
  n <- sum(divergent)
  N <- length(divergent)
  
  if (n == 0) {
    message(
      sprintf("%s of %s iterations ended with a divergence.", n, N)
    )
  } else {
    message(
      sprintf("%s of %s iterations ended with a divergence (%s%%).",
              n, N, 100 * n / N), 
      "\nTry increasing 'adapt_delta' to remove the divergences."
    )
  }
}

# Check transitions that ended prematurely due to maximum tree depth limit
#
# @param object A stanfit object.
# @return Nothing, just prints the number (and percentage) of iterations that
#   saturated the max treedepth and, if any, suggests increasing max_treedepth.
#
check_treedepth <- function(object) {
  max_depth <- object@stan_args[[1]]$control$max_treedepth
  if (is.null(max_depth)) {
    max_depth <- 10
  }
  treedepths <- sampler_param_vector(object, "treedepth__")
  n <- sum(treedepths == max_depth)
  N <- length(treedepths)
  
  if (n == 0) {
    message(
      sprintf("%s of %s iterations saturated the maximum tree depth of %s.", 
              n, N, max_depth)
    )
  } else {
    message(
      sprintf('%s of %s iterations saturated the maximum tree depth of %s (%s%%).',
              n, N, max_depth, 100 * n / N), 
      "\nTry increasing 'max_treedepth' to avoid saturation."
    )
  }
}


# Check the energy Bayesian fraction of missing information (E-BFMI)
#
# @param object A stanfit object.
# @return Nothing, just prints E-BFMI for chains with low E-BFMI and suggests 
#   reparameterizing.
#
check_energy <- function(object) {
  energies_by_chain <- sampler_param_matrix(object, "energy__")
  EBFMIs <- apply(energies_by_chain, 2, function(x) {
    numer <- sum(diff(x) ^ 2) / length(x)
    denom <- var(x)
    numer / denom
  })
  bad_chains <- which(EBFMIs < 0.2)
  if (!length(bad_chains)) {
    message("E-BFMI indicated no pathological behavior.")
  } else {
    message(
      "E-BFMI indicated possible pathological behavior:\n", 
      sprintf("  Chain %s: E-BFMI = %.3f\n", bad_chains, EBFMIs[bad_chains]),
      "E-BFMI below 0.2 indicates you may need to reparameterize your model."
    )
  }
}


# internal ----------------------------------------------------------------

# Extract single sampler parameters in conventient form
#
# @param object A stanfit object.
# @param param The name of a single sampler parameter (e.g. "divergent__").
# @return sampler_param_vector returns a single vector of length chains * iters
#   (post-warmup). sampler_param_matrix returns an iters (post-warmup) by chains
#   matrix.
#
sampler_param_vector <- function(object, param) {
  stopifnot(length(param) == 1, is.character(param))
  sampler_params <- get_sampler_params(object, inc_warmup = FALSE)
  do.call(rbind, sampler_params)[, param]
}
sampler_param_matrix <- function(object, param) {
  stopifnot(length(param) == 1, is.character(param))
  sampler_params <- get_sampler_params(object, inc_warmup=FALSE)
  sapply(sampler_params, function(x) x[, param])
}
