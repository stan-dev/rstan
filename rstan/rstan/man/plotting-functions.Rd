\name{rstan-plotting-functions}
\alias{rstan-plotting-functions}

\title{RStan Plotting Functions}
\description{
  List of RStan plotting functions that return ggplot objects
}
\section{RStan plotting functions}{
\describe{
  \item{Posterior intervals and point estimates}{\code{\link{stan_plot}}}
  \item{Traceplots}{\code{\link{stan_trace}}}
  \item{Histograms}{\code{\link{stan_hist}}}
  \item{Kernel density estimates}{\code{\link{stan_dens}}}
  \item{Scatterplots}{\code{\link{stan_scat}}}
  \item{Diagnostics for Hamiltonian Monte Carlo and the No-U-Turn Sampler}{\code{\link{stan_diag}}}
  \item{Rhat}{\code{\link{stan_rhat}}}
  \item{Ratio of effective sample size to total posterior sample size}{\code{\link{stan_ess}}}
  \item{Ratio of Monte Carlo standard error to posterior standard deviation}{\code{\link{stan_mcse}}}
  \item{Autocorrelation}{\code{\link[rstan]{stan_ac}}}
}
}
\seealso{
\code{\link[=rstan_gg_options]{Plot options}}
}
