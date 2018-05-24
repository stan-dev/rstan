# This file is part of RStan
# Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017 Trustees of Columbia University
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

# stan_prob_autocovariance <- function(v) { 
#   .Call("stan_prob_autocovariance", v)
# }

fft_next_good_size <- function(N) {
  # Find the optimal next size for the FFT so that
  # a minimum number of zeros are padded.
  if (N <= 2)
    return(2)
  while (TRUE) {
    m = N
    while ((m %% 2) == 0) m = m / 2
    while ((m %% 3) == 0) m = m / 3
    while ((m %% 5) == 0) m = m / 5
    if (m <= 1)
      return(N)
    N = N + 1
  }
}

autocovariance <- function(y) {
  # Compute autocovariance estimates for every lag for the specified
  # input sequence using a fast Fourier transform approach.
  N <- length(y)
  M <- fft_next_good_size(N)
  Mt2 <- 2 * M
  yc <- y - mean(y)
  yc <- c(yc, rep.int(0, Mt2-N))
  transform <- fft(yc)
  ac <- fft(Conj(transform) * transform, inverse = TRUE)
  ac <- Re(ac)[1:N]/(N*2*seq(N, 1, by = -1))
  ac
}
   
autocorrelation <- function(y) {
  # Compute autocorrelation estimates for every lag for the specified
  # input sequence using a fast Fourier transform approach.
  ac <- autocovariance(y)
  ac <- ac/ac[1]
}

ess_rfun <- function(sims) {
  # Compute the effective sample size for samples of several chains 
  # for one parameter; see the C++ code of function  
  # effective_sample_size in chains.cpp 
  # 
  # Args:
  #   sims: a 2-d array _without_ warmup samples (# iter * # chains) 
  # 
  if (is.vector(sims)) dim(sims) <- c(length(sims), 1)
  chains <- ncol(sims)
  n_samples <- nrow(sims)

  acov <- lapply(1:chains, 
                 FUN = function(i) autocovariance(sims[,i])) 
  acov <- do.call(cbind, acov)
  chain_mean <- apply(sims, 2, mean)
  mean_var <- mean(acov[1,]) * n_samples / (n_samples - 1) 
  var_plus <- mean_var * (n_samples - 1) / n_samples
  if (chains > 1) 
    var_plus <- var_plus + var(chain_mean)
  # Geyer's initial positive sequence
  rho_hat_t <- rep.int(0, n_samples)
  t <- 0
  rho_hat_even <- 1;
  rho_hat_t[t+1] <- rho_hat_even
  rho_hat_odd <- 1 - (mean_var - mean(acov[t+2, ])) / var_plus
  rho_hat_t[t+2] <- rho_hat_odd
  t <- 2  
  while (t < nrow(acov)-1 && !is.nan(rho_hat_even + rho_hat_odd) && (rho_hat_even + rho_hat_odd > 0)) {
    rho_hat_even = 1 - (mean_var - mean(acov[t+1, ])) / var_plus
    rho_hat_odd = 1 - (mean_var - mean(acov[t+2, ])) / var_plus
    if ((rho_hat_even + rho_hat_odd) >= 0) {
      rho_hat_t[t+1] <- rho_hat_even
      rho_hat_t[t+2] <- rho_hat_odd
    }
    t <- t + 2
  }
  max_t <- t
  # Geyer's initial monotone sequence
  t <- 2
  while (t <= max_t - 2) {  
    if (rho_hat_t[t + 1] + rho_hat_t[t + 2] >
  	rho_hat_t[t - 1] + rho_hat_t[t]) {
      rho_hat_t[t + 1] = (rho_hat_t[t - 1] + rho_hat_t[t]) / 2;
      rho_hat_t[t + 2] = rho_hat_t[t + 1];
    }
    t <- t + 2
  }
  ess <- chains * n_samples
  ess <- ess / (-1 + 2 * sum(rho_hat_t[1:max_t]))
  ess 
} 

split_rhat_rfun <- function(sims) {
  # Compute the split rhat for the diagnostics of converging; 
  # see the C++ code of split_potential_scale_reduction in chains.cpp.  
  # 
  # Args:
  #   sims: a 2-d array _without_ warmup samples (# iter * # chains) 
  # 
  # Note: 
  #   The R function wrapping the C++ implementation is defined 
  #   in chains.R with name rstan_splitrhat2_cpp 
  if (is.vector(sims)) dim(sims) <- c(length(sims), 1)
  chains <- ncol(sims)
  n_samples <- nrow(sims)
  half_n <- floor(n_samples / 2)
  # cat("n_samples=", n_samples, "\n"); cat("chains=", chains, "\n")
  # cat("half_n=", half_n, "\n")
  idx_2nd <- n_samples - half_n + 1
  
  split_chain_mean <- numeric(chains * 2)
  split_chain_var <- numeric(chains * 2)
  
  for (i in 1:chains) {
    split_chain_mean[i] <- mean(sims[1:half_n, i])
    split_chain_var[i] <- var(sims[1:half_n, i])
    split_chain_mean[chains + i] <- mean(sims[idx_2nd:n_samples, i])
    split_chain_var[chains + i] <- var(sims[idx_2nd:n_samples, i])
  } 
  var_between <- half_n * var(split_chain_mean)
  var_within <- mean(split_chain_var) 
  sqrt((var_between/var_within + half_n -1)/half_n)
} 

monitor <- function(sims, warmup = floor(dim(sims)[1] / 2), 
                    probs = c(0.025, 0.25, 0.50, 0.75, 0.975), 
                    digits_summary = 1, 
                    print = TRUE, ...) { 
  # print the summary for a general simulation results 
  # of 3-d array: # iter * # chains * # parameters 
  # Args:
  #   sims: a 3-d array described above 
  #   warmup: the number of iterations used for warmup 
  #   probs: probs of summarizing quantiles 
  #   print: print out the results
  # 
  # Return: 
  #   A summary array  
  dim_sims <- dim(sims)
  if (is.null(dim_sims))
      dim(sims) <- c(length(sims), 1, 1)
  if (length(dim_sims) == 2)
      dim(sims) <- c(dim_sims, 1)
  if (length(dim_sims) > 3) 
    stop("'sims' has more than 3 dimensions")
  dim_sims <- dim(sims)
  if (warmup > dim_sims[1])
    stop("warmup is larger than the total number of iterations")
  if (is(sims, "stanfit")) {
    warmup <- 0L
    sims <- as.array(sims)
  }
  dimnames_sims <- dimnames(sims)
  parnames <- dimnames_sims[[3]]
  num_par <- dim_sims[3]
  
  if (is.null(parnames)) parnames <- paste0("V", 1:num_par)
  sims_wow <- if (warmup >= 1) apply(sims, c(2, 3), FUN = function(x) x[-(1:warmup)]) else sims 
  m <- apply(sims_wow, 3, mean)
  sd <- sapply(1:num_par, FUN = function(i) sd(as.vector(sims_wow[,,i]))) 
  quan <- lapply(1:num_par, FUN = function(i) quantile(sims_wow[,,i], probs = probs))
  probs_str <- names(quan[[1]])
  quan <- do.call(rbind, quan)
  rhat <- sapply(1:num_par, FUN = function(i) split_rhat_rfun(sims_wow[,,i]))
  ess <- sapply(1:num_par, FUN = function(i) ess_rfun(sims_wow[,,i]))
  sem <- sd / sqrt(ess)
 
  summary <- cbind(m, sem, sd, quan, ess, rhat)
  colnames(summary) <- c("mean", "se_mean", "sd", probs_str, "n_eff", "Rhat")
  rownames(summary) <- parnames 
  if (print) {
    cat("Inference for the input samples (")
    cat(dim_sims[2], " chains: each with iter=", dim_sims[1], "; warmup=", warmup, "):\n\n", sep = "")
    # round n_eff to integers
    summary[, 'n_eff'] <- round(summary[, 'n_eff'], 0)
    print(round(summary, digits_summary), ...)
 
    cat("\nFor each parameter, n_eff is a crude measure of effective sample size,\n", 
        "and Rhat is the potential scale reduction factor on split chains (at \n",
        "convergence, Rhat=1).\n", sep = '')
  } 
  invisible(summary) 
} 

