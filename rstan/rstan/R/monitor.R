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
  posterior::rhat(sims)
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
  posterior::ess_bulk(sims)
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
  posterior::ess_tail(sims)
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
  posterior::ess_quantile(sims, prob)
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
  posterior::ess_mean(sims)
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
  # TODO: are these two equivalent/ok to change to posterior's implementation?
  # min(ess_rfun(split_chains(sims)), ess_rfun(split_chains(sims^2)))
  posterior::ess_sd(sims)
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
  posterior::mcse_quantile(sims, prob)
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
  posterior::mcse_sd(sims)
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
    quan <- unname(posterior::quantile2(sims_i, probs = probs))
    #quan2 <- posterior::quantile2(sims_i, probs = c(0.05, 0.5, 0.95))
    mean <- mean(sims_i)
    sd <- sd(sims_i)
    mcse_quan <- sapply(probs, function(p) posterior::mcse_quantile(sims_i, probs = p))
    mcse_mean <- posterior::mcse_mean(sims_i)
    mcse_sd <- posterior::mcse_sd(sims_i)
    rhat <- posterior::rhat(sims_i)
    ess_bulk <- round(posterior::ess_bulk(sims_i))
    ess_tail <- round(posterior::ess_tail(sims_i))
    ess <- round(posterior::ess_bulk(sims_i))
    out[[i]] <- c(
      mean, mcse_mean, sd, quan, ess, rhat,
      valid, #quan2,
      mcse_quan, mcse_sd, ess_bulk, ess_tail
    )
  }
  
  out <- as.data.frame(do.call(rbind, out))
  #probs_str <- names(quantile(sims_i, probs = probs, na.rm = TRUE))
  str_quan <- paste0("Q", probs * 100)
  #str_quan2 <- paste0("Q", c(0.05, 0.5, 0.95) * 100)
  str_mcse_quan <- paste0("MCSE_", str_quan)
  colnames(out) <- c("mean", "se_mean", "sd", str_quan, "n_eff", "Rhat",
                     "valid", #str_quan2,
                     str_mcse_quan, "MCSE_SD", "Bulk_ESS", "Tail_ESS")
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

# should NA be returned by a convergence diagnostic?
should_return_NA <- function(x) {
  anyNA(x) || any(!is.finite(x)) || is_constant(x)
}
