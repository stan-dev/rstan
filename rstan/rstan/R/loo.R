# Leave-one-out cross-validation
# 
# The \code{loo} method for stanfit objects ---a wrapper around the
# \code{loo.array} method from the \pkg{loo} package--- computes approximate
# leave-one-out cross-validation using Pareto smoothed importance sampling
# (PSIS-LOO CV).
# 
# @param x stanfit object
# @param pars Name of parameter, transformed parameter, or generated quantity in
#   the Stan program corresponding to the pointwise log-likelihood. If not
#   specified the default behavior is to look for \code{"log_lik"}.
# @param ... Ignored.
# @param cores Number of cores to use for parallelization. Passed to
#   \code{\link[loo]{loo.array}}.
# @param save_psis Should the intermediate results from \code{\link[loo]{psis}}
#   be saved in the returned object? The default is \code{FALSE}. This can be
#   useful to avoid repeated computation when using other functions in the
#   \pkg{loo} and \pkg{bayesplot} packages.
# 
# @details Stan does not automatically compute and store the log-likelihood. It
#   is up to the user to incorporate it into the Stan program if it is to be
#   extracted after fitting the model. In a Stan model, the pointwise log
#   likelihood can be coded as a vector in the transformed parameters block
#   (and then summed up in the model block) or it can be coded entirely in the
#   generated quantities block. We recommend using the generated quantities
#   block so that the computations are carried out only once per iteration
#   rather than once per HMC leapfrog step.
#
#   For example, the following is the \code{generated quantities} block for
#   computing and saving the log-likelihood for a linear regression model with
#   \code{N} data points, outcome \code{y}, predictor matrix \code{X},
#   coefficients \code{beta}, and standard deviation \code{sigma}:
#
#  \code{vector[N] log_lik;}
#
#  \code{for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | X[n, ] * beta, sigma);}
#
loo.stanfit <-
  function(x,
           pars = "log_lik",
           ...,
           save_psis = FALSE,
           cores = getOption("mc.cores", 1)) {
    stopifnot(length(pars) == 1L)
    LLarray <- loo::extract_log_lik(stanfit = x,
                                    parameter_name = pars,
                                    merge_chains = FALSE)
    r_eff <- loo::relative_eff(x = exp(LLarray), cores = cores)
    loo::loo.array(LLarray,
                   r_eff = r_eff,
                   cores = cores,
                   save_psis = save_psis)
  }

setMethod("loo", "stanfit", loo.stanfit)
