
\name{pairs.stanfit}
\alias{pairs.stanfit}

\title{Create a matrix of output plots from a \code{stanfit} object} 

\description{A \code{\link{pairs}} method that is customized for MCMC output} 

\usage{
  \method{pairs}{stanfit}(x, labels = NULL, panel = NULL, \dots, lower.panel = NULL, 
    upper.panel = NULL, diag.panel = NULL, text.panel = NULL, 
    label.pos = 0.5 + has.diag/3, cex.labels = NULL, font.labels = 1, 
    row1attop = TRUE, gap = 1, pars = NULL, condition = "accept_stat__", 
    include = TRUE) 
} 

\arguments{
  \item{x}{An object of S4 class \code{stanfit}}
  \item{labels, panel, \dots, lower.panel, upper.panel, diag.panel}{Same as in 
    \code{\link{pairs}} syntactically but see the Details section for different
    default arguments}
  \item{text.panel, label.pos, cex.labels, font.labels, row1attop, gap}{Same as in 
    \code{\link{pairs}}}
  \item{pars}{If not \code{NULL}, a character vector indicating which quantities
    to include in the plots, which is passed to \code{\link{extract}}. Thus,
    by default, all unknown quantities are included, which may be far too many
    to visualize on a small computer screen. If \code{include = FALSE}, then
    the named parameters are excluded from the plot.}
  \item{condition}{If \code{NULL}, it will plot roughly half of the
    chains in the lower panel and the rest in the upper panel. An
    integer vector can be passed to select some subset of the chains, of which
    roughly half will be plotted in the lower panel and the rest in the upper
    panel. A list of two integer vectors can be passed, each specifying a subset
    of the chains to be plotted in the lower and upper panels respectively.
    A single number between zero and one exclusive can be passed, which is 
    interpreted as the proportion of realizations (among all chains) to plot
    in the lower panel starting with the first realization in each chain, 
    with the complement (from the end of each chain) plotted in the upper panel.
    A (possibly abbreviated) character vector of length one can be passed among 
    \code{"accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__"},
    \code{"n_divergent__"} or \code{"lp__"}, which are the variables produced by 
    \code{\link{get_sampler_params}} and \code{\link{get_logposterior}}. In that 
    case the lower panel will plot realizations that are below the median of the 
    indicated variable (or are zero in the case of \code{"n_divergent__"}) and 
    the upper panel will plot realizations that are greater than or equal to the 
    median of the indicated variable (or are one in the case of \code{"n_divergent__"}). 
    Finally, any logical vector whose length is equal to the product of the number 
    of iterations and the number of chains can be passed, in which case realizations
    corresponding to \code{FALSE} and \code{TRUE} will be plotted in the lower 
    and upper panel respectively. The default is \code{"accept_stat__"}.}
  \item{include}{Logical scalar indicating whether to include (the default) or
    excldue the parameters named in the \code{pars} argument from the plot.}
} 

\details{
  This method differs from the default \code{\link{pairs}} method in the following
  ways. If unspecified, the \code{\link{smoothScatter}} function is used for the
  off-diagonal plots, rather than \code{\link{points}}, since the former is more
  appropriate for visualizing thousands of draws from a posterior distribution.
  Also, if unspecified, histograms of the marginal distribution of each quantity
  are placed on the diagonal of the plot, after pooling all of the chains specified
  by the \code{chain\_id} argument.
  
  The draws from the warmup phase are always discarded before plotting.
  
  By default, the lower (upper) triangle of the plot contains draws with below
  (above) median acceptance probability. Also, if \code{condition} is not
  \code{"n_divergent__"}, red points will be superimposed onto the smoothed
  density plots indicating which (if any) iterations encountered a divergent 
  transition. Otherwise, yellow points indicate a transition that hit the
  maximum treedepth rather than terminated its evolution normally.
}
  
\seealso{
  S4 class \code{\linkS4class{stanfit}} and its method \code{\link{extract}} as
  well as the \code{\link{pairs}} generic function. Also, see
  \code{\link{get_sampler_params}} and \code{\link{get_logposterior}}.
} 

\examples{
example(read_stan_csv)
pairs(fit, pars = c("mu", "sigma", "alpha", "lp__"))
}
