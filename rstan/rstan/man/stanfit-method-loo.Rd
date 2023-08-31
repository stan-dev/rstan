\name{loo.stanfit}
\alias{loo.stanfit}
\alias{loo,stanfit-method}
\alias{loo}

\title{Approximate leave-one-out cross-validation}

\description{
A \code{\link[loo]{loo}} method that is customized for stanfit objects.
The \code{loo} method for stanfit objects ---a wrapper around the \code{array}
method for \code{\link[loo]{loo}} in the  \pkg{loo} package --- computes PSIS-LOO CV,
approximate leave-one-out cross-validation using Pareto smoothed importance
sampling (Vehtari, Gelman, and Gabry, 2017a,2017b).
}

\usage{
  \method{loo}{stanfit}(x,
    pars = "log_lik",
    save_psis = FALSE,
    cores = getOption("mc.cores", 1),
    moment_match = FALSE,
    k_threshold = 0.7,
    \dots)
}

\arguments{
  \item{x}{An object of S4 class \code{stanfit}.}
  \item{pars}{Name of transformed parameter or generated quantity in
  the Stan program corresponding to the pointwise log-likelihood. If not
  specified the default behavior is to look for \code{"log_lik"}.}
  \item{save_psis}{Should the intermediate results from \code{\link[loo]{psis}}
  be saved in the returned object? The default is \code{FALSE}. This can be
  useful to avoid repeated computation when using other functions in the
  \pkg{loo} and \pkg{bayesplot} packages.}
  \item{cores}{Number of cores to use for parallelization. The default is 1 unless
  \code{cores} is specified or the \code{mc.cores} \link[=options]{option}
  has been set.}
\item{moment_match}{Logical; Whether to use the moment matching algorithm for
  observations with high Pareto k values to improve accuracy. Note:
  because the moment matching algorithm relies on the \code{unconstrain_pars}
  method in RStan it is only available if run in the same R session as fitting the
  model.}
  \item{k_threshold}{Threshold value for Pareto k values above which
  the moment matching algorithm is used. If \code{moment_match} is \code{FALSE},
  this is ignored.}
  \item{\dots}{Ignored.}
}

\details{
  Stan does not automatically compute and store the log-likelihood. It
  is up to the user to incorporate it into the Stan program if it is to be
  extracted after fitting the model. In a Stan program, the pointwise log
  likelihood can be coded as a vector in the transformed parameters block
  (and then summed up in the model block) or it can be coded entirely in the
  generated quantities block. We recommend using the generated quantities
  block so that the computations are carried out only once per iteration
  rather than once per HMC leapfrog step.

  For example, the following is the \code{generated quantities} block for
  computing and saving the log-likelihood for a linear regression model with
  \code{N} data points, outcome \code{y}, predictor matrix \code{X} (including
  column of 1s for intercept), coefficients \code{beta},
  and standard deviation \code{sigma}:

 \code{vector[N] log_lik;}

 \code{for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | X[n, ] * beta, sigma);}

  This function automatically uses Pareto k diagnostics for assessing
  the accuracy of importance sampling for each observation. When the
  diagnostics indicate that importance sampling for certain observations
  is inaccurate, a moment matching algorithm can be used, which can
  improve the accuracy (Paananen et al., 2020).
}

\value{
  A list with class \code{c("psis_loo", "loo")}, as detailed in the
  \code{\link[loo]{loo}} documentation.
}

\references{
Vehtari, A., Gelman, A., and Gabry, J. (2017a).
Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC.
\emph{Statistics and Computing}. 27(5), 1413-1432.
\code{doi:10.1007/s11222-016-9696-4}.
\url{https://arxiv.org/abs/1507.04544},
\url{https://link.springer.com/article/10.1007/s11222-016-9696-4}

Vehtari, A., Gelman, A., and Gabry, J. (2017b).
Pareto smoothed importance sampling. arXiv preprint:
\url{https://arxiv.org/abs/1507.02646}

Yao, Y., Vehtari, A., Simpson, D., and Gelman, A. (2018).
Using stacking to average Bayesian predictive distributions.
Bayesian Analysis, advance publication, \code{doi:10.1214/17-BA1091}.
\url{https://projecteuclid.org/euclid.ba/1516093227}.

Paananen, T., Piironen, J., Buerkner, P.-C., Vehtari, A. (2020).
Implicitly Adaptive Importance Sampling.
arXiv preprint:
\url{https://arxiv.org/abs/1906.08850}.
}

\seealso{
 \itemize{
 \item The \pkg{\link[loo]{loo}} package documentation, including the vignettes for
 many examples (\url{https://mc-stan.org/loo/}).
 \item \code{\link[loo]{loo_moment_match}} for the moment matching algorithm.
 \item \code{\link[loo]{loo_model_weights}} for model averaging/weighting via
 stacking or pseudo-BMA weighting.
 }
}

\examples{
\dontrun{
# Generate a dataset from N(0,1)
N <- 100
y <- rnorm(N, 0, 1)

# Suppose we have three models for y:
#  1) y ~ N(-1, sigma)
#  2) y ~ N(0.5, sigma)
#  3) y ~ N(0.6,sigma)
#
stan_code <- "
data {
  int N;
  vector[N] y;
  real mu_fixed;
}
  parameters {
  real<lower=0> sigma;
}
model {
  sigma ~ exponential(1);
  y ~ normal(mu_fixed, sigma);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N) log_lik[n] = normal_lpdf(y[n]| mu_fixed, sigma);
}"

mod <- stan_model(model_code = stan_code)
fit1 <- sampling(mod, data=list(N=N, y=y, mu_fixed=-1))
fit2 <- sampling(mod, data=list(N=N, y=y, mu_fixed=0.5))
fit3 <- sampling(mod, data=list(N=N, y=y, mu_fixed=0.6))

# use the loo method for stanfit objects
loo1 <- loo(fit1, pars = "log_lik")
print(loo1)

# which is equivalent to
LL <- as.array(fit1, pars = "log_lik")
r_eff <- loo::relative_eff(exp(LL))
loo1b <- loo::loo.array(LL, r_eff = r_eff)
print(loo1b)

# compute loo for the other models
loo2 <- loo(fit2)
loo3 <- loo(fit3)

# stacking weights
wts <- loo::loo_model_weights(list(loo1, loo2, loo3), method = "stacking")
print(wts)

# use the moment matching for loo with a stanfit object
loo_mm <- loo(fit1, pars = "log_lik", moment_match = TRUE)
print(loo_mm)
}
}
