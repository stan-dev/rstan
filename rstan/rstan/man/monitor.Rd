\name{monitor}
\alias{monitor}
\alias{print.simsummary}
\alias{[.simsummary}
\title{Compute summaries of MCMC draws and monitor convergence} 
\description{
  Similar to the \code{print} method for \code{stanfit} objects, but \code{monitor}
  takes an array of simulations as its argument rather than a \code{stanfit}
  object. For a 3-D array (iterations * chains * parameters) of MCMC draws, 
  \code{monitor} computes means, standard deviations, quantiles, Monte Carlo 
  standard errors, split Rhats, and effective sample sizes. By default, half of 
  the iterations are considered warmup and are excluded.
}

\usage{
monitor(sims, warmup = floor(dim(sims)[1]/2), 
        probs = c(0.025, 0.25, 0.5, 0.75, 0.975), 
        digits_summary = 1, print = TRUE, ...)
\method{print}{simsummary}(x, digits = 3, se = FALSE, ...)
\method{[}{simsummary}(x, i, j, drop = if (missing(i)) TRUE else length(j) == 1)
        
}
\arguments{
  \item{sims}{A 3-D array (iterations * chains * parameters) of MCMC 
    simulations from any MCMC algorithm.}
  \item{warmup}{The number of warmup iterations to be excluded 
    when computing the summaries. The default is half of the total number 
    of iterations. If \code{sims} doesn't include the warmup iterations then 
    \code{warmup} should be set to zero.}
  \item{probs}{A numeric vector specifying quantiles of interest. The 
    defaults is \code{c(0.025,0.25,0.5,0.75,0.975)}.}
  \item{digits_summary}{The number of significant digits to use when printing
    the summary, defaulting to 1. Applies to the quantities other than 
    the effective sample size, which is always rounded to the nearest integer.}
  \item{print}{Logical, indicating whether to print the summary after the 
    computations are performed.} 
  \item{\dots}{Additional arguments passed to the underlying \code{print} method.}
  \item{x}{An object of class \code{simsummary} created by \code{monitor}}
  \item{digits}{An integer scalar defaulting to 3 for the number of decimal
    places to print}
  \item{se}{A logical scalar defaulting to \code{FALSE} indicating whether
    to print the estimated standard errors of the estimates}
  \item{i}{A vector indicating which rows of the object created by \code{monitor}
    to select}
  \item{j}{A vector indicating which columns of the object crated by \code{monitor}
    to select}
  \item{drop}{A logical scalar indicating whether the resulting object 
    should return a vector where possible}
}

\value{
  A 2-D array with rows corresponding to parameters and columns to the 
  summary statistics that can be printed and subset.
}

\references{
  The Stan Development Team 
  \emph{Stan Modeling Language User's Guide and Reference Manual}. 
  \url{http://mc-stan.org}. 
}

\seealso{
  S4 class \code{\linkS4class{stanfit}} and particularly its
  \code{\link[=print.stanfit]{print}} method.
}

\examples{
csvfiles <- dir(system.file('misc', package = 'rstan'),
                pattern = 'rstan_doc_ex_[0-9].csv', full.names = TRUE)
fit <- read_stan_csv(csvfiles)
# The following is just for the purpose of giving an example
# since print can be used for a stanfit object.
monitor(extract(fit, permuted = FALSE, inc_warmup = TRUE))
} 
