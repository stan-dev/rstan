\name{plot-methods}
\alias{plot-methods} 
% \alias{plot,ANY,ANY-method}
\alias{plot,stanfit,missing-method}
\docType{methods}

\title{\code{plot}: parameter estimates} 
\description{Plot posterior intervals and point estimates for parameters and 
generated quantities.} 

\usage{
  \S4method{plot}{stanfit,missing}(x, ..., pars, include = TRUE, unconstrain = FALSE)
} 

\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "stanfit", y = "missing")}}{Plot posterior intervals for parameters and generated quantities. See \code{\link{stan_plot}} 
    (for which \code{plot} is a wrapper) for more details.}
  }
}

\arguments{
  \item{x}{An instance of class \code{\linkS4class{stanfit}}.}
  \item{...}{Optional arguments describing plot elements. See Details.}
  \item{pars}{A character vector of parameter names. Defaults to all parameters 
  or the first 10 parameters (if there are more than 10).}
  \item{include}{Should the parameters given by the \code{pars} argument be
    included (the default) or excluded from the plot? Only relevant if 
    \code{pars} is not missing.}
  \item{unconstrain}{Should parameters be plotted on the unconstrained space?
    Defaults to \code{FALSE}.}
}
\details{
The following arguments can be specied in \code{...}:
\describe{
  \item{point_est}{The point estimate to show. Either "mean" or "median".}
  \item{show_density}{Should kernel density estimates be plotted above 
  the intervals? Defaults to \code{FALSE}.}
  \item{ci_level}{The interval to highlight. Defaults to 0.5 (a 50 percent interval).}
  \item{outer_level}{The outer interval also shown but not highlighted.
  Defaults to 0.95.}
  \item{show_outer_line}{Should the \code{outer_level} interval be shown or hidden? 
  Defaults to = \code{TRUE} (to plot it).}
  \item{fill_color,outline_color,est_color}{Colors to override the defaults for 
  the highlighted interval, the outer interval (and density outline), and the 
  point estimate.}
  }
}
\value{A \code{\link[ggplot2]{ggplot}} object that can be further customized
  using the \pkg{ggplot2} package.} 
\seealso{
\code{\link[=plotting-functions]{List of RStan plotting functions}},
\code{\link[=rstan_gg_options]{Plot options}}
}
\examples{\dontrun{
library(rstan) 
fit <- stan(model_code = "parameters {real y;} model {y ~ normal(0,1);}") 
plot(fit)
plot(fit, show_density = TRUE)
}
} 
\keyword{methods}
