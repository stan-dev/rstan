\name{Diagnostic plots}
\alias{stan_diag}
\alias{stan_par}
\alias{stan_rhat}
\alias{stan_ess}
\alias{stan_mcse}

\title{RStan Diagnostic plots}
\description{
  Diagnostic plots for HMC and NUTS using ggplot2. 
}
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
    and \code{stan_mcse}, optional arguments to \code{stat_bin} in the 
    \pkg{ggplot2} package.}
}
\value{
  For \code{stan_diag} and \code{stan_par}, a list containing the ggplot objects for 
  each of the displayed plots. For \code{stan_rhat}, \code{stan_ess}, 
  and \code{stan_mcse}, a single ggplot object.
}
\details{
 \describe{
  \item{\code{stan_rhat},\code{stan_ess},\code{stan_mcse}}{Respectively, 
  these plots show the 
  distribution of the Rhat statistic, the ratio of effective sample size to 
  total sample size, and the ratio of Monte Carlo standard error
  to posterior standard deviation for the estimated parameters. These 
  plots are not intended to identify individual parameters, but rather to allow
  for quickly identifying if the estimated values of these quantities are
  desireable for all parameters.}
  \item{\code{stan_par}}{Calling \code{stan_par} generates three plots: 
  (i) a scatterplot of \code{par} vs. the accumulated log-posterior (\code{lp__}), 
  (ii) a scatterplot of \code{par} vs. the average Metropolis acceptance rate 
  (\code{accept_stat}), and 
  (iii) a violin plot showing the distribution of \code{par} at each of the 
  sampled step sizes (one per chain).
  For the scatterplots, red points are superimposed to indicate which 
  (if any) iterations encountered a divergent transition. Yellow points indicate 
  a transition that hit the maximum treedepth rather than terminated its 
  evolution normally.}
  \item{\code{stan_diag}}{The \code{information} argument is used to specify which 
  plots \code{stan_diag} should generate:
  
  \itemize{
    \item \code{information='sample'} Histograms of \code{lp__} 
      and \code{accept_stat}, as well as a scatterplot showing their 
      joint distribution.
    \item \code{information='stepsize'} Violin plots showing the 
      distributions of \code{lp__} and \code{accept_stat} at each of the sampled
      step sizes (one per chain).
    \item \code{information='treedepth'} Histogram of \code{treedepth} and 
      violin plots showing the distributions of \code{lp__} and 
      \code{accept_stat} for each value of \code{treedepth}.
    \item \code{information='divergence'} Violin plots showing the 
      distributions of \code{lp__} and \code{accept_stat} for iterations that
      encountered divergent transitions (\code{divergent=1}) and those that
      did not (\code{divergent=0}).
    }
  }
}
}
\note{
  For details about the individual diagnostics and sampler parameters and their
  interpretations see the Stan Modeling Language User's Guide and Reference 
  Manual at \url{https://mc-stan.org/documentation/}.
}
\seealso{
\code{\link[=rstan-plotting-functions]{List of RStan plotting functions}},
\code{\link[=rstan_gg_options]{Plot options}}
}
\examples{
\dontrun{
fit <- stan_demo("eight_schools")

stan_diag(fit, info = 'sample') # shows three plots together
samp_info <- stan_diag(fit, info = 'sample') # saves the three plots in a list
samp_info[[3]] # access just the third plot

stan_diag(fit, info = 'sample', chain = 1) # overlay chain 1

stan_par(fit, par = "mu")
}
}
