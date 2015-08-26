\name{stan_trace}
\alias{stan_trace}
\alias{stan_dens}
\alias{stan_hist}
\alias{stan_plot}
\alias{stan_ac}
\alias{stan_stepsize}
\alias{stan_sample}
\alias{stan_treedepth}
\alias{stan_divergence}
\alias{stan_param}
\alias{stan_rhat}
\alias{stan_ess}
\alias{stan_mcse}

\title{Visual Posterior Analysis}
\usage{
  stan_trace(object, pars, include = TRUE, unconstrain = FALSE,
             inc_warmup = FALSE, window = NULL, nrow = NULL, ncol = NULL,
             theme_elements = NULL, ...)
}
\arguments{
  \item{object}{A stanfit or stanreg object.}
  
  \item{pars}{A character vector of parameter names. If not specified, all
    parameters are used.}
  
  \item{include}{Should the parameters given by the \code{pars} argument be
    included (the default) or excluded from the plot?}
  
  \item{unconstrain}{Should parameters be plotted on the unconstrained space?
    Defaults to \code{FALSE}.}
  
  \item{inc_warmup}{Should warmup iterations be included? Defaults to
    \code{FALSE}.}
  
  \item{window}{The window of iterations to plot. The default of \code{NULL}
    means plotting all iterations. \code{window} needs to be an integer vector
    of at least length 1. The smaller value specifies the starting iteration
    and the larger value the last iteration. When there is only one value
    provided, it is taken to be the starting iteration. Examples:
      \code{window=100}, \code{window=c(1500,1600)}.}
  
  \item{nrow,ncol}{Passed to \code{\link[ggplot2]{facet_wrap}}.}
  
  \item{theme_elements}{Additional ggplot2 theme elements to use. This should
    be an object created by \code{\link[ggplot2]{theme}}. For example:
      \code{theme_elements = theme(panel.background = element_rect(fill =
                                                                     "gray"), legend.position = "top")}.}
  
  \item{...}{Additional named arguments to \code{\link[ggplot2]{geom_path}}
    (e.g. \code{linetype}, \code{size}, \code{alpha}, etc.).}
}
\value{
  A \code{\link[ggplot2]{ggplot}} object that can be further customized
  using the \pkg{ggplot2} package.
}
\description{
  ggplot2 for RStan
}