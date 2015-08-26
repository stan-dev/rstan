\name{Diagnostic plots}
\alias{stan_diag}
\alias{stan_par}
\alias{stan_rhat}
\alias{stan_ess}
\alias{stan_mcse}

\title{Diagnostic plots}
\usage{
  stan_diag(object, 
            information = c("sample","stepsize", "treedepth","divergence"), 
            chain = 0, ...)
  stan_par(object, par, chain = 0, ...)
  
  stan_rhat(object, pars, ...)
  stan_ess(object, pars, ...)
  stan_mcse(object, pars, ...)
}
\arguments{
  \item{object}{A stanfit or stanreg object.}
  
  \item{information}{The information to be contained in the diagnostic plot.}
  
  \item{par,pars}{The name of a single scalar parameter (\code{par}) or 
    one or more parameter names (\code{pars}).}
  
  \item{chain}{If \code{chain=0} (the default) all chains are combined.
    Otherwise the plot for \code{chain} is overlaid on the plot for all chains
    combined.}
  
  \item{...}{For \code{stan_diag} and \code{stan_par}, optional arguments to 
    \code{\link[gridExtra]{arrangeGrob}}. For \code{stan_rhat}, \code{stan_ess}, 
    and \code{stan_mcse}, optional arguments to \code{\link[ggplot2]{stat_bin}}.}
}
\value{
  For \code{stan_diag} and \code{stan_par}, a list containing the ggplot objects for 
  each of the displayed plots. For \code{stan_rhat}, \code{stan_ess}, 
  and \code{stan_mcse}, a single ggplot object.
}
\description{
  Diagnostic plots
}
\note{
  For \code{stan_par} red points are superimposed onto the
  scatterplots indicating which (if any) iterations encountered a divergent
  transition. Also for \code{stan_par}, yellow points indicate a transition
  that hit the maximum treedepth rather than terminated its evolution
  normally.
}
\examples{
fit <- stan_demo("eight_schools", iter = 200)
stan_diag(fit, info = 'sample', chain = 1)
samp_info <- stan_diag(fit, info = 'sample')
samp_info[[3]]
}
