\name{sbc}
\alias{sbc}
\alias{print.sbc}
\alias{plot.sbc}
\title{Simulation Based Calibration (sbc)}
\description{
\if{html}{\figure{stanlogo.png}{options: width="25px" alt="http://mc-stan.org/about/logo/"}}
Check whether a model is well-calibrated with respect to the
prior distribution and hence possibly amenable to obtaining
a posterior distribution conditional on observed data.
}
\usage{
  sbc(stanmodel, data, M, ..., save_progress)
  \method{plot}{sbc}(x, thin = 3, ...)
  \method{print}{sbc}(x, ...)

}
\arguments{
  \item{stanmodel}{
    An object of \code{\link{stanmodel-class}} that is first created by
    calling the \code{\link{stan_model}} function}
  \item{data}{A named \code{list} or \code{environment} providing the data for 
    the model, or a character vector for all the names of objects to use as data. 
    This is the same format as in \code{\link{stan}} or \code{\link{sampling}}.}
  \item{M}{
    The number of times to condition on draws from the prior predictive
    distribution}
  \item{\dots}{Additional arguments that are passed to \code{\link{sampling}},
    such as \code{refresh = 0} when calling \code{sbc}. For the \code{plot}
    and \code{print} methods, the additional arguments are not used.}
  \item{x}{An object produced by \code{sbc}}
  \item{thin}{An integer vector of length one indicating the thinning interval
    when plotting, which defaults to 3}
  \item{save_progress}{If a directory is provided, stanfit objects
    are saved to disk making it easy to resume a partial \code{sbc}
    run after interruption.}
}
\details{
  This function assumes adherence to the following conventions in the 
  underlying Stan program:
  \enumerate{
  \item Realizations of the unknown parameters are drawn in the \code{transformed data}
    block of the Stan program and are postfixed with an underscore, such as
    \code{theta_}. These are considered the \dQuote{true} parameters being estimated by
    the corresponding symbol declared in the \code{parameters} block, which
    should have the same name except for the trailing underscore, such as \code{theta}.
  \item The realizations of the unknown parameters are then conditioned on when drawing from
    the prior predictive distribution, also in the \code{transformed data} block.
    There is no restriction on the symbol name that holds the realizations from
    the prior predictive distribution but for clarity, it should not end with
    a trailing underscore.
  \item The realizations of the unknown parameters should be copied into a \code{vector}
    in the \code{generated quantities} block named \code{pars_}.
  \item The realizations from the prior predictive distribution should be copied 
    into an object (of the same type) in the \code{generated quantities} block
    named \code{y_}. Technically, this step is optional and could be omitted
    to conserve RAM, but inspecting the realizations from the prior predictive distribution
    is a good way to judge whether the priors are reasonable.
  \item The \code{generated quantities} block must contain an integer array named 
    \code{ranks_} whose only values are zero or one, depending on whether the realization of a 
    parameter from the posterior distribution exceeds the corresponding \dQuote{true} 
    realization, such as \code{theta > theta_;}. These are not actually "ranks"
    but can be used afterwards to reconstruct (thinned) ranks.
  \item The \code{generated quantities} block may contain a vector named \code{log_lik}
    whose values are the contribution to the log-likelihood by each observation. This
    is optional but facilitates calculating Pareto k shape parameters to judge whether
    the posterior distribution is sensitive to particular observations.
  }
  
  Although the user can pass additional arguments to \code{\link{sampling}} through the \dots,
  the following arguments are hard-coded and should not be passed through the \dots:
  \enumerate{
  \item \code{pars = "ranks_"} because nothing else needs to be stored for each posterior draw
  \item \code{include = TRUE} to ensure that \code{"ranks_"} is included rather than excluded
  \item \code{chains = 1} because only one chain is run for each integer less than \code{M}
  \item \code{seed} because a sequence of seeds is used across the \code{M} runs to preserve
    independence across runs
  \item \code{save_warmup = FALSE} because the warmup realizations are not relevant
  \item \code{thin = 1} because thinning can and should be done after the Markov Chain is
    finished, as is done by the \code{thin} argument to the \code{plot} method in order to
    make the histograms consist of approximately independent realizations
  }
  Other arguments will take the default values used by \code{\link{sampling}} unless 
  passed through the \dots. Specifying \code{refresh = 0} is recommended to avoid printing
  a lot of intermediate progress reports to the screen. It may be necessary to pass a 
  list to the \code{control} argument of \code{\link{sampling}} with elements \code{adapt_delta}
  and / or \code{max_treedepth} in order to obtain adequate results.
  
  Ideally, users would want to see the absence of divergent transitions (which is shown
  by the \code{print} method) and other warnings, plus an approximately uniform histogram 
  of the ranks for each parameter (which are shown by the \code{plot} method). See the 
  vignette for more details.
}
\value{
  The \code{sbc} function outputs a list of S3 class \code{"sbc"}, which contains the
  following elements:
  \enumerate{
  \item \code{ranks} A list of \code{M} matrices, each with number of 
    rows equal to the number of saved iterations and number of columns equal to 
    the number of unknown parameters. These matrices contain the realizations
    of the \code{ranks_} object from the \code{generated quantities} block of the
    Stan program.
  \item \code{Y} If present, a matrix of realizations from the prior predictive
    distribution whose rows are equal to the number of observations and whose columns
    are equal to \code{M}, which are taken from the \code{y_} object in the 
    \code{generated quantities} block of the Stan program.
  \item \code{pars} A matrix of realizations from the prior distribution whose rows
    are equal to the number of parameters and whose columns are equal to \code{M}, 
    which are taken from the \code{pars_} object in the \code{generated quantities}
    block of the Stan program.
  \item \code{pareto_k} A matrix of Pareto k shape parameter estimates or \code{NULL}
    if there is no \code{log_lik} symbol in the \code{generated quantities} block
    of the Stan program
  \item \code{sampler_params} A three-dimensional array that results from combining
    calls to \code{\link{get_sampler_params}} for each of
    the \code{M} runs. The resulting matrix has rows equal to the number of 
    post-warmup iterations, columns equal to six, and \code{M} floors. The columns
    are named \code{"accept_stat__"}, \code{"stepsize__"}, \code{"treedepth__"},
    \code{"n_leapfrog__"}, \code{"divergent__"}, and \code{"energy__"}. The most
    important of which is \code{"divergent__"}, which should be all zeros and
    perhaps \code{"treedepth__"}, which should only rarely get up to the value
    of \code{max_treedepth} passed as an element of the \code{control} list
    to \code{\link{sampling}} or otherwise defaults to \eqn{10}.
  }
  The \code{print} method outputs the number of divergent transitions and
  returns \code{NULL} invisibly.
  The \code{plot} method returns a \code{\link[ggplot2]{ggplot}} object
  with histograms whose appearance can be further customized.
}


\references{
  The Stan Development Team 
  \emph{Stan Modeling Language User's Guide and Reference Manual}. 
  \url{http://mc-stan.org}. 

  Talts, S., Betancourt, M., Simpson, D., Vehtari, A., and Gelman, A. (2018).
  Validating Bayesian Inference Algorithms with Simulation-Based Calibration.
  arXiv preprint arXiv:1804.06788. \url{https://arxiv.org/abs/1804.06788}
}
\seealso{
  \code{\link{stan_model}} and \code{\link{sampling}}
}
\examples{
scode <- "
data {
  int<lower = 1> N;
  real<lower = 0> a;
  real<lower = 0> b;
}
transformed data { // these adhere to the conventions above
  real pi_ = beta_rng(a, b);
  int y = binomial_rng(N, pi_);
}
parameters {
  real<lower = 0, upper = 1> pi;
}
model {
  target += beta_lpdf(pi | a, b);
  target += binomial_lpmf(y | N, pi);
}
generated quantities { // these adhere to the conventions above
  int y_ = y;
  vector[1] pars_;
  int ranks_[1] = {pi > pi_};
  vector[N] log_lik;
  pars_[1] = pi_;
  for (n in 1:y) log_lik[n] = bernoulli_lpmf(1 | pi);
  for (n in (y + 1):N) log_lik[n] = bernoulli_lpmf(0 | pi);
}
"
}
