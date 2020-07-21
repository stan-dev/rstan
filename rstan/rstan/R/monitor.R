# This file is part of RStan
# Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018 Trustees of Columbia University
# Copyright (C) 2018, 2019 Aki Vehtari, Paul Bürkner
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

fft_next_good_size <- function(N) {
  # Find the optimal next size for the FFT so that
  # a minimum number of zeros are padded.
  if (N <= 2)
    return(2)
  while (TRUE) {
    m <- N
    while ((m %% 2) == 0) m <- m / 2
    while ((m %% 3) == 0) m <- m / 3
    while ((m %% 5) == 0) m <- m / 5
    if (m <= 1)
      return(N)
    N <- N + 1
  }
}

#' Autocovariance estimates
#'
#' Compute autocovariance estimates for every lag for the specified
#' input sequence using a fast Fourier transform approach. Estimate
#' for lag t is scaled by N-t.
#'
#' @param y A numeric vector forming a sequence of values.
#'
#' @return A numeric vector of autocovariances at every lag (scaled by N-lag).
autocovariance <- function(y) {
  N <- length(y)
  M <- fft_next_good_size(N)
  Mt2 <- 2 * M
  yc <- y - mean(y)
  yc <- c(yc, rep.int(0, Mt2 - N))
  transform <- fft(yc)
  ac <- fft(Conj(transform) * transform, inverse = TRUE)
  # use "biased" estimate as recommended by Geyer (1992)
  ac <- Re(ac)[1:N] / (N^2 * 2)
  ac
}

#' Autocorrelation estimates
#'
#' Compute autocorrelation estimates for every lag for the specified
#' input sequence using a fast Fourier transform approach. Estimate
#' for lag t is scaled by N-t.
#'
#' @param y A numeric vector forming a sequence of values.
#'
#' @return A numeric vector of autocorrelations at every lag (scaled by N-lag).
autocorrelation <- function(y) {
  ac <- autocovariance(y)
  ac <- ac / ac[1]
}

#' Rank normalization
#'
#' Compute rank normalization for a numeric array. First replace each
#' value by its rank. Average rank for ties are used to conserve the
#' number of unique values of discrete quantities. Second, normalize
#' ranks via the inverse normal transformation.
#'
#' @param x A numeric array of values.
#'
#' @return A numeric array of rank normalized values with the same
#'     size as input.
z_scale <- function(x) {
  S <- length(x)
  r <- rank(x, ties.method = 'average')
  z <- qnorm((r - 1 / 2) / S)
  z[is.na(x)] <- NA
  if (!is.null(dim(x))) {
    # output should have the input dimension
    z <- array(z, dim = dim(x), dimnames = dimnames(x))
  }
  z
}

#' Rank uniformization
#'
#' Compute rank uniformization for a numeric array. First replace each
#' value by its rank. Average rank for ties are used to conserve the
#' number of unique values of discrete quantities. Second, uniformize
#' ranks to scale \code{[1/(2S), 1-1/(2S)]}, where \code{S} is the the number
#' of values.
#'
#' @param x A numeric array of values.
#'
#' @return A numeric array of rank uniformized values with the same
#'     size as input.
#'
u_scale <- function(x) {
  S <- length(x)
  r <- rank(x, ties.method = 'average')
  u <- (r - 1 / 2) / S
  u[is.na(x)] <- NA
  if (!is.null(dim(x))) {
    # output should have the input dimension
    u <- array(u, dim = dim(x), dimnames = dimnames(x))
  }
  u
}

#' Rank values
#'
#' Compute ranks for a numeric array. First replace each
#' value by its rank. Average rank for ties are used to conserve the
#' number of unique values of discrete quantities. Second, normalize
#' ranks via the inverse normal transformation.
#'
#' @param x A numeric array of values.
#'
#' @return A numeric array of ranked values with the same
#'     size as input.
#'
r_scale <- function(x) {
  S <- length(x)
  r <- rank(x, ties.method = 'average')
  r[is.na(x)] <- NA
  if (!is.null(dim(x))) {
    # output should have the input dimension
    r <- array(r, dim = dim(x), dimnames = dimnames(x))
  }
  r
}

split_chains <- function(sims) {
  # split Markov chains
  # Args:
  #   sims: a 2D array of samples (# iter * # chains)
  if (is.vector(sims)) {
    dim(sims) <- c(length(sims), 1)
  }
  niter <- dim(sims)[1]
  if (niter == 1L) return(sims)
  half <- niter / 2
  cbind(sims[1:floor(half), ], sims[ceiling(half + 1):niter, ])
}

is_constant <- function(x, tol = .Machine$double.eps) {
  abs(max(x) - min(x)) < tol
}

#' Traditional Rhat convergence diagnostic
#'
#' Compute the Rhat convergence diagnostic for a single parameter
#' For split-Rhat, call this with split chains.
#'
#' @param sims A 2D array _without_ warmup samples (# iter * # chains).
#'
#' @return A single numeric value for Rhat.
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2019). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
rhat_rfun <- function(sims) {
  if (anyNA(sims)) {
    return(NA)
  }
  if (any(!is.finite(sims))) {
    return(NaN)
  }
  if (is_constant(sims)) {
    return(NA)
  }
  if (is.vector(sims)) {
    dim(sims) <- c(length(sims), 1)
  }
  chains <- ncol(sims)
  n_samples <- nrow(sims)
  chain_mean <- numeric(chains)
  chain_var <- numeric(chains)
  for (i in seq_len(chains)) {
    chain_mean[i] <- mean(sims[, i])
    chain_var[i] <- var(sims[, i])
  }
  var_between <- n_samples * var(chain_mean)
  var_within <- mean(chain_var)
  sqrt((var_between / var_within + n_samples - 1) / n_samples)
}

#' Effective sample size
#'
#' Compute the effective sample size estimate for a sample of several chains
#' for one parameter. For split-ESS, call this with split chains.
#'
#' @param sims A 2D array _without_ warmup samples (# iter * # chains).
#'
#' @return A single numeric value for the effective sample size.
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2019). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
ess_rfun <- function(sims) {
  if (is.vector(sims)) {
    dim(sims) <- c(length(sims), 1)
  }
  chains <- ncol(sims)
  n_samples <- nrow(sims)
  if (n_samples < 3L || anyNA(sims)) {
    return(NA)
  }
  if (any(!is.finite(sims))) {
    return(NaN)
  }
  if (is_constant(sims)) {
    return(NA)
  }
  acov <- lapply(seq_len(chains), function(i) autocovariance(sims[, i]))
  acov <- do.call(cbind, acov)
  chain_mean <- apply(sims, 2, mean)
  mean_var <- mean(acov[1, ]) * n_samples / (n_samples - 1)
  var_plus <- mean_var * (n_samples - 1) / n_samples
  if (chains > 1)
    var_plus <- var_plus + var(chain_mean)

  # Geyer's initial positive sequence
  rho_hat_t <- rep.int(0, n_samples)
  t <- 0
  rho_hat_even <- 1
  rho_hat_t[t + 1] <- rho_hat_even
  rho_hat_odd <- 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
  rho_hat_t[t + 2] <- rho_hat_odd
  while (t < nrow(acov) - 5 && !is.nan(rho_hat_even + rho_hat_odd) &&
         (rho_hat_even + rho_hat_odd > 0)) {
    t <- t + 2
    rho_hat_even = 1 - (mean_var - mean(acov[t + 1, ])) / var_plus
    rho_hat_odd = 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
    if ((rho_hat_even + rho_hat_odd) >= 0) {
      rho_hat_t[t + 1] <- rho_hat_even
      rho_hat_t[t + 2] <- rho_hat_odd
    }
  }
  max_t <- t
  # this is used in the improved estimate
  if (rho_hat_even>0)
      rho_hat_t[max_t + 1] <- rho_hat_even

  # Geyer's initial monotone sequence
  t <- 0
  while (t <= max_t - 4) {
    t <- t + 2
    if (rho_hat_t[t + 1] + rho_hat_t[t + 2] >
        rho_hat_t[t - 1] + rho_hat_t[t]) {
      rho_hat_t[t + 1] = (rho_hat_t[t - 1] + rho_hat_t[t]) / 2;
      rho_hat_t[t + 2] = rho_hat_t[t + 1];
    }
  }
  ess <- chains * n_samples
  # Geyer's truncated estimate
  # tau_hat <- -1 + 2 * sum(rho_hat_t[1:max_t])
  # Improved estimate reduces variance in antithetic case
  tau_hat <- -1 + 2 * sum(rho_hat_t[1:max_t]) + rho_hat_t[max_t+1]
  # Safety check for negative values and with max ess equal to ess*log10(ess)
  tau_hat <- max(tau_hat, 1/log10(ess))
  ess <- ess / tau_hat
  ess
}

#' Rhat convergence diagnostic
#'
#' Compute Rhat convergence diagnostic as the maximum of rank normalized
#' split-Rhat and rank normalized folded-split-Rhat for one parameter.
#'
#' @param sims A 2D array _without_ warmup samples (# iter * # chains).
#'
#' @return A single numeric value for the effective sample size.
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2019). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
#'
#' @export
Rhat <- function(sims) {
  bulk_rhat <- rhat_rfun(z_scale(split_chains(sims)))
  sims_folded <- abs(sims - median(sims))
  tail_rhat <- rhat_rfun(z_scale(split_chains(sims_folded)))
  max(bulk_rhat, tail_rhat)
}

#' Bulk effective sample size (bulk-ESS)
#'
#' Compute bulk effective sample size estimate (bulk-ESS) for one parameter.
#' Bulk-ESS is useful as a generic diagnostic for the sampling
#' efficiency in the bulk of the posterior. It is defined as the
#' effective sample size for rank normalized values using split chains.
#'
#' @param sims A 2D array _without_ warmup samples (# iter * # chains).
#'
#' @return A single numeric value for the bulk effective sample size.
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2019). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
#'
#' @export
ess_bulk <- function(sims) {
  ess_rfun(z_scale(split_chains(sims)))
}

#' Tail effective sample size (tail-ESS)
#'
#' Compute tail effective sample size estimate (tail-ESS) for one parameter.
#' Tail-ESS is useful for generic diagnostic for the sampling
#' efficiency in the tails of the posterior. It is defined as
#' the minimum of the effective sample sizes for 5% and 95% quantiles.
#'
#' @param sims A 2D array _without_ warmup samples (# iter * # chains).
#'
#' @return A single numeric value for the tail effective sample size.
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2019). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
#'
#' @export
ess_tail <- function(sims) {
  I05 <- sims <= quantile(sims, 0.05, na.rm = TRUE)
  q05_ess <- ess_rfun(split_chains(I05))
  I95 <- sims <= quantile(sims, 0.95, na.rm = TRUE)
  q95_ess <- ess_rfun(split_chains(I95))
  min(q05_ess, q95_ess)
}

#' Quantile effective sample size
#'
#' Compute effective sample size estimate for a quantile estimate of
#' one parameter.
#'
#' @param sims A 2D array _without_ warmup samples (# iter * # chains).
#' @param prob A single numeric value of probability.
#'
#' @return A single numeric value for the effective sample size for a
#'     quantile estimate corresponding to the probability.
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2019). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
#'
#' @export
ess_quantile <- function(sims, prob) {
  I <- sims <= quantile(sims, prob, na.rm = TRUE)
  ess_rfun(split_chains(I))
}

#' Effective sample size
#'
#' Compute effective sample size estimate for a mean (expectation)
#' estimate of one parameter.
#'
#' @param sims A 2D array _without_ warmup samples (# iter * # chains).
#'
#' @return A single numeric value for the effective sample size
#'     estimate for mean estimate.
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2019). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
#'
#' @export
ess_mean <- function(sims) {
  ess_rfun(split_chains(sims))
}

#' Effective sample size
#'
#' Compute effective sample size estimate for standard deviation (s)
#' estimate of one parameter. This is defined as minimum of effective
#' sample size estimate for mean and mean of squared value.
#'
#' @param sims A 2D array _without_ warmup samples (# iter * # chains).
#'
#' @return A single numeric value for the effective sample size
#'     estimate for standard deviation estimate.
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2019). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
#'
#' @export
ess_sd <- function(sims) {
  min(ess_rfun(split_chains(sims)), ess_rfun(split_chains(sims^2)))
}

#' Monte Carlo diagnostics for a quantile
#'
#' Compute Monte Carlo standard error, 5%-quantile, 95%-quantile, and
#' effective sample size estimate for a quantile estimate of a single
#' parameter.
#'
#' @param sims A 2D array _without_ warmup samples (# iter * # chains).
#' @param prob A single numeric value of probability.
#'
#' @return A data frame with Monte Carlo standard error (mcse),
#'     5%-quantile (Q05), 95%-quantile (Q95), and effective sample
#'     size estimate (ess).
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2019). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
conv_quantile <- function(sims, prob) {
  if (is.vector(sims)) {
    dim(sims) <- c(length(sims), 1)
  }
  ess <- ess_quantile(sims, prob)
  p <- c(0.1586553, 0.8413447, 0.05, 0.95)
  a <- qbeta(p, ess * prob + 1, ess * (1 - prob) + 1)
  ssims <- sort(sims)
  S <- length(ssims)
  th1 <- ssims[max(round(a[1] * S), 1)]
  th2 <- ssims[min(round(a[2] * S), S)]
  mcse <- (th2 - th1) / 2
  th1 <- ssims[max(round(a[3] * S), 1)]
  th2 <- ssims[min(round(a[4] * S), S)]
  data.frame(mcse = mcse, Q05 = th1, Q95 = th2, ess = ess)
}

#' Monte Carlo standard error for a quantile
#'
#' Compute Monte Carlo standard error for a quantile estimate of a
#' single parameter.
#'
#' @param sims A 2D array _without_ warmup samples (# iter * # chains).
#' @param prob A single numeric value of probability.
#'
#' @return A single numeric value for Monte Carlo standard error for a
#'     quantile estimate corresponding to the probability.
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2019). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
#'
#' @export
mcse_quantile <- function(sims, prob) {
  conv_quantile(sims, prob)$mcse
}

#' Monte Carlo standard error for mean
#'
#' Compute Monte Carlo standard error for mean (expectation) of a
#' single parameter.
#'
#' @param sims A 2D array _without_ warmup samples (# iter * # chains).
#'
#' @return A single numeric value for Monte Carlo standard error
#'     for mean estimate.
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2019). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2019). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
#'
#' @export
mcse_mean <- function(sims) {
  sd(sims) / sqrt(ess_mean(sims))
}

#' Monte Carlo standard error for standard error
#'
#' Compute Monte Carlo standard error for standard deviation (sd) of a
#' single parameter using Stirling's approximation and assuming
#' approximate normality.
#'
#' @param sims A 2D array _without_ warmup samples (# iter * # chains).
#'
#' @return A single numeric value for Monte Carlo standard error
#'     for standard deviation estimate.
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2019). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
#'
#' @export
mcse_sd <- function(sims) {
  # assumes normality of sims and uses Stirling's approximation
  ess_sd <- ess_sd(sims)
  sd(sims) * sqrt(exp(1) * (1 - 1 / ess_sd)^(ess_sd - 1) - 1)
}

#' Summary of General Simulation Results
#'
#' Create a summary for general simulation results. Computed
#' quantities are specified quantiles, mean, standard deviation,
#' corresponding Monte Carlo standard errors, Rhat, Bulk-ESS and
#' Tail-ESS.
#'
#' @param sims A 3-dimensional array of simulation results. The first
#'     dimension is the number of iterations per chain, the second its
#'     the number of chains and the third is the number of
#'     parameters. Alternatively, \code{sims} can be a \code{stanfit}
#'     object from which the simulation results will be extracted.
#' @param warmup The number of iterations used for warmup. These will
#'     be removed before computing summary. Default is 0.
#' @param probs A vector of probabilities defining summarizing
#'     quantiles. Default is c(0.05, 0.50, 0.95).
#' @param se Logical, defaulting to FALSE. If TRUE print also MCSEs.
#' @param print Logical, defaulting to PRINT. If TRUE print summary table.
#' @param digits Positive integer, a number of digits printedm
#'     defaulting to 1.
#'
#' @return A summary matrix.
#'
#' @references
#' Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
#' Paul-Christian Bürkner (2019). Rank-normalization, folding, and
#' localization: An improved R-hat for assessing convergence of
#' MCMC. \emph{arXiv preprint} \code{arXiv:1903.08008}.
#'
#' @export
monitor <- function(sims, warmup = floor(dim(sims)[1] / 2),
                    probs = c(0.025, 0.25, 0.50, 0.75, 0.975),
                    digits_summary = 1, print = TRUE, ...) {
  if (inherits(sims, "stanfit")) {
    chains <- sims@sim$chains
    iter <- sims@sim$iter
    warmup <- 0L
    parnames <- names(sims)
    sims <- as.array(sims)
  } else {
    dim_sims <- dim(sims)
    if (is.null(dim_sims)) {
      dim(sims) <- c(length(sims), 1, 1)
    } else if (length(dim_sims) == 2) {
      dim(sims) <- c(dim_sims, 1)
    } else if (length(dim_sims) > 3) {
      stop("'sims' has more than 3 dimensions")
    }
    parnames <- dimnames(sims)[[3]]
    if (is.null(parnames)) {
      parnames <- paste0("V", seq_len(dim(sims)[3]))
    }
    iter <- dim(sims)[1]
    chains <- dim(sims)[2]
    if (warmup > dim(sims)[1]) {
      stop("warmup is larger than the total number of iterations")
    }
    if (warmup >= 1) {
      sims <- sims[-seq_len(warmup), , , drop = FALSE]
    }
  }

  out <- vector("list", length(parnames))
  out <- setNames(out, parnames)
  # loop over parameters
  for (i in seq_along(out)) {
    sims_i <- sims[, , i]
    valid <- all(is.finite(sims_i))
    quan <- unname(quantile(sims_i, probs = probs, na.rm = TRUE))
    quan2 <- quantile(sims_i, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)
    mean <- mean(sims_i)
    sd <- sd(sims_i)
    mcse_quan <- sapply(probs, mcse_quantile, sims = sims_i)
    mcse_mean <- mcse_mean(sims_i)
    mcse_sd <- mcse_sd(sims_i)
    rhat <- Rhat(sims_i)
    ess_bulk <- round(ess_bulk(sims_i))
    ess_tail <- round(ess_tail(sims_i))
    ess <- round(ess_rfun(sims_i))
    out[[i]] <- c(
      mean, mcse_mean, sd, quan, ess, rhat,
      valid, quan2, mcse_quan, mcse_sd, ess_bulk, ess_tail
    )
  }

  out <- as.data.frame(do.call(rbind, out))
  probs_str <- names(quantile(sims_i, probs = probs, na.rm = TRUE))
  str_quan <- paste0("Q", probs * 100)
  str_quan2 <- paste0("Q", c(0.05, 0.5, 0.95) * 100)
  str_mcse_quan <- paste0("MCSE_", str_quan)
  colnames(out) <- c("mean", "se_mean", "sd", probs_str, "n_eff", "Rhat",
                     "valid", str_quan2, str_mcse_quan, "MCSE_SD", "Bulk_ESS", "Tail_ESS")
  rownames(out) <- parnames

  # replace NAs with appropriate values if draws are valid
  S <- prod(dim(sims)[1:2])
  valid <- out[, "valid"]
  out[valid & !is.finite(out[, "Rhat"]), "Rhat"] <- 1
  out[valid & !is.finite(out[, "Bulk_ESS"]), "Bulk_ESS"] <- S
  out[valid & !is.finite(out[, "Tail_ESS"]), "Tail_ESS"] <- S
  SE_vars <- colnames(out)[grepl("^SE_", colnames(out), ignore.case = TRUE)]
  for (v in SE_vars) {
  	out[valid & !is.finite(out[, v]), v] <- 0
  }

  out <- structure(
    out,
    chains = chains,
    iter = iter,
    warmup = warmup,
    class = c("simsummary", "data.frame")
  )
  if (print) print.simsummary(out, digits = digits_summary, ...)
  return(invisible(out))
}

print.simsummary <- function(x, digits = 3, se = FALSE, ...) {
  atts <- attributes(x)
  px <- x
  class(px) <- "data.frame"
  quan2 <- grepl("^Q", colnames(px))
  if (se) {
    px <- cbind(px[ , quan2], Mean = px$mean, SD = px$sd,
                px[ , grepl("^MCSE", colnames(px))], " Rhat" = px$Rhat,
                Bulk_ESS = px$Bulk_ESS, Tail_ESS = px$Tail_ESS)
  } else {
    px <- cbind(px[ , quan2], Mean = px$mean, SD = px$sd, " Rhat" = px$Rhat,
                Bulk_ESS = px$Bulk_ESS, Tail_ESS = px$Tail_ESS)
  }

  decimal_places <- max(1, digits - 1)
  px$` Rhat` <- round(px$` Rhat`, digits = max(2, decimal_places))
  estimates <- setdiff(names(px), c(" Rhat", "Bulk_ESS", "Tail_ESS"))
  px[, estimates] <- round(px[, estimates], digits = decimal_places)
  cat(
    "Inference for the input samples (", atts$chains,
    " chains: each with iter = ", atts$iter,
    "; warmup = ", atts$warmup, "):\n\n", sep = ""
  )
  print(px, ...)
  if (!isTRUE(atts$extra)) {
  	cat(
  		"\nFor each parameter, Bulk_ESS and Tail_ESS are crude measures of \n",
  		"effective sample size for bulk and tail quantities respectively (an ESS > 100 \n",
  		"per chain is considered good), and Rhat is the potential scale reduction \n",
  		"factor on rank normalized split chains (at convergence, Rhat <= 1.05).\n", sep = ""
  	)
  }
  invisible(x)
}

# this is needed to make DeLorean's vignette build
`[.simsummary` <- function (x, i, j, drop = if (missing(i)) TRUE else length(j) == 1) {
  out <- `[.data.frame`(x, i, j, drop)
  nms <- rownames(x)[i]
  if (drop) names(out) <- nms
  else rownames(out) <- nms
  return(out)
}
