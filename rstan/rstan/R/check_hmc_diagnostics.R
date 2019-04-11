# This file is part of RStan
# Copyright (C) 2018 Trustees of Columbia University
#
# RStan is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# RStan is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

throw_sampler_warnings <- function(object) {
  if (!is(object, "stanfit"))
    stop("'object' must be of class 'stanfit'")
  sp <- get_sampler_params(object, inc_warmup = FALSE)
  n_d <- sum(sapply(sp, FUN = function(x) {
    if ("divergent__" %in% colnames(x)) return(sum(x[,"divergent__"]))
    else return(0)
  }))
  if (n_d > 0) {
    ad <- object@stan_args[[1]]$control
    if (is.null(ad)) ad <- 0.8
    else {
      ad <- ad$adapt_delta
      if (is.null(ad)) ad <- 0.8
    }
    warning("There were ", n_d, " divergent transitions after warmup.",
            " Increasing adapt_delta above ", ad, " may help. See\n",
            "http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup", call. = FALSE)
  }
  max_td <- object@stan_args[[1]]$control
  if (is.null(max_td)) max_td <- 10
  else {
    max_td <- max_td$max_treedepth
    if (is.null(max_td)) max_td <- 10
  }
  n_m <- sum(sapply(sp, FUN = function(x) {
    if ("treedepth__" %in% colnames(x)) return(sum(x[,"treedepth__"] >= max_td))
    else return(0)
  }))
  if (n_m > 0)
    warning("There were ", n_m,
            " transitions after warmup that exceeded the maximum treedepth.",
            " Increase max_treedepth above ", max_td, ". See\n",
            "http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded", call. = FALSE)
  n_e <- 0L
  if (is_sfinstance_valid(object) && all(sapply(sp, function(x) "energy__" %in% colnames(x)))) {
    E <- as.matrix(sapply(sp, FUN = function(x) x[,"energy__"]))
    threshold <- 0.2
    if (nrow(E) > 1) {
      EBFMI <- get_num_upars(object) / apply(E, 2, var)
      n_e <- sum(EBFMI < threshold, na.rm = TRUE)
    }
    else n_e <- 0L
    if (n_e > 0)
      warning("There were ", n_e, 
              " chains where the estimated Bayesian Fraction of Missing Information",
              " was low. See\n", 
              "http://mc-stan.org/misc/warnings.html#bfmi-low", call. = FALSE)
  }
  if (n_d > 0 || n_m > 0 || n_e > 0) 
    warning("Examine the pairs() plot to diagnose sampling problems\n",
            call. = FALSE, noBreaks. = TRUE)
  
  rhat <- apply(sims, MARGIN = 3, FUN = rhat)
  if (any(rhat > 1.01))
      warning("The largest R-hat is ", round(max(rhat,2)),
            ", indicating chains have not mixed.\n",
            "Running the chains for more iterations may help. See\n",
            "http://mc-stan.org/misc/warnings.html#r-hat")
  bulk_ess <- apply(sims, MARGIN = 3, FUN = ess_bulk)
  if (any(bulk_ess < 100 * ncol(object)))
    warning("Bulk Effective Samples Size (ESS) is too low, ",
            "indicating posterior means and medians may be unreliable.\n",
            "Running the chains for more iterations may help. See\n",
            "http://mc-stan.org/misc/warnings.html#bulk-ess")
  tail_ess <- apply(sims, MARGIN = 3, FUN = ess_tail)
  if (any(tail_ess < 100 * ncol(object)))
    warning("Tail Effective Samples Size (ESS) is too low, indicating",
            "posterior variances and tail quantiles may be unreliable.\n",
            "Running the chains for more iterations may help. See\n",
            "http://mc-stan.org/misc/warnings.html#tail-ess")
  
  return(invisible(NULL))
}


# Check divergences, treedepth, and energy diagnostics
#
# @param object A stanfit object.
# @return Nothing, just prints the output from the functions it calls internally.
#
check_hmc_diagnostics <- function(object) {
  stopifnot(is.stanfit(object))
  cat("\nDivergences:\n")
  check_divergences(object)
  cat("\nTree depth:\n")
  check_treedepth(object)
  cat("\nEnergy:\n")
  check_energy(object)
}


# Get a logical vector indicating whict transitions ended with a divergence
#
# @param object A stanfit object.
# @return a logical vector
#
get_divergent_iterations <- function(object) {
  stopifnot(is.stanfit(object))
  sampler_param_vector(object, "divergent__") > 0
}

# Get the number of transitions that ended with a divergence
#
# @param object A stanfit object.
# @return the number of divergent transitions
#
get_num_divergent <- function(object) {
  stopifnot(is.stanfit(object))
  divergent <- get_divergent_iterations(object)
  sum(divergent)
}

# Check transitions that ended with a divergence
#
# @param object A stanfit object.
# @return nothing, just prints the number (and percentage) of iterations that 
#   ended with a divergence and, if any, suggests increasing adapt_delta.
#
check_divergences <- function(object) {
  stopifnot(is.stanfit(object))
  divergent <- get_divergent_iterations(object)
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

get_treedepth_threshold <- function(object) {
  stopifnot(is.stanfit(object))
  max_depth <- object@stan_args[[1]]$control$max_treedepth
  if (is.null(max_depth)) {
    max_depth <- 10
  }
  max_depth
}

# Get a logical vector indicating transitions that ended prematurely 
# due to maximum tree depth limit
#
# @param object A stanfit object.
# @return a logical vector
#
get_max_treedepth_iterations <- function(object) {
  stopifnot(is.stanfit(object))
  max_depth <- get_treedepth_threshold(object)
  treedepths <- sampler_param_vector(object, "treedepth__") >= max_depth
  
  treedepths
}

# Get the number of transitions that ended prematurely due to maximum tree depth limit
#
# @param object A stanfit object.
# @return the number of affected transitions
#
get_num_max_treedepth <- function(object) {
  stopifnot(is.stanfit(object))
  sum(get_max_treedepth_iterations(object))
}

# Check transitions that ended prematurely due to maximum tree depth limit
#
# @param object A stanfit object.
# @return Nothing, just prints the number (and percentage) of iterations that
#   saturated the max treedepth and, if any, suggests increasing max_treedepth.
#
check_treedepth <- function(object) {
  stopifnot(is.stanfit(object))
  max_depth <- get_treedepth_threshold(object)
  treedepths <- get_max_treedepth_iterations(object)
  n <- sum(treedepths)
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

# Get the Bayesian fraction of missing information (E-BFMI)
#
# @param object A stanfit object.
# @return vector, one element per chain containing the E-BFMI
#
get_bfmi <- function(object) {
  stopifnot(is.stanfit(object))
  energies_by_chain <- sampler_param_matrix(object, "energy__")
  EBFMIs <- apply(energies_by_chain, 2, function(x) {
    numer <- sum(diff(x) ^ 2) / length(x)
    denom <- var(x)
    numer / denom
  })
  
  EBFMIs
}

# Get the chains with low Bayesian fraction of missing information (E-BFMI)
#
# @param object A stanfit object.
# @return a vector of IDs of chains with low E-BFMI 
#
get_low_bfmi_chains <- function(object) {
  stopifnot(is.stanfit(object))
  which(get_bfmi(object) < 0.2)
}

# Check the energy Bayesian fraction of missing information (E-BFMI)
#
# @param object A stanfit object.
# @return Nothing, just prints E-BFMI for chains with low E-BFMI and suggests 
#   reparameterizing.
#
check_energy <- function(object) {
  stopifnot(is.stanfit(object))
  EBFMIs <- get_bfmi(object)
  
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


# Get the number of actual leapfrog evaluations for each iteration
#
# @param object A stanfit object.
# @return an integer vector with the number of evaluations for each iteration
#
get_num_leapfrog_per_iteration <- function(object) {
  stopifnot(is.stanfit(object))
  sampler_param_vector(object,"n_leapfrog__")
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

