\docType{methods}
\name{summary-methods}
\alias{summary,stanfit-method}

\title{Summary method for stanfit objects}
\description{
Summarize the distributions of estimated parameters and derived quantities 
using the posterior draws.
}

\usage{
  \S4method{summary}{stanfit}(object, pars, probs = c(0.025, 0.25, 0.50, 0.75, 0.975),
  use_cache = TRUE, ...) 
}

\arguments{
  \item{object}{An instance of class \code{\linkS4class{stanfit}}.}
  \item{pars}{A character vector of parameter names. Defaults to all parameters
  as well as the log-posterior (\code{lp__}).}
  \item{probs}{A numeric vector of \code{\link[stats]{quantile}}s of interest. 
    The default is \code{c(0.025,0.25,0.5,0.75,0.975)}.}
  \item{use_cache}{Logical, defaulting to \code{TRUE}. When \code{use_cache=TRUE}
    the summary quantities for all parameters are computed and cached for future
    use. Setting \code{use_cache=FALSE} can be used to avoid performing the 
    summary computations for all parameters if \code{pars} is given as some 
    specific parameters.}
  \item{\dots}{Currently unused.}
}
\value{
The \code{summary} method returns a named list with elements \code{summary} 
and \code{c_summary}, which contain summaries for for all chains merged and 
individual chains, respectively.
Included in the summaries are quantiles, means, standard deviations (sd),
effective sample sizes (n_eff), and split Rhats (the potential scale
reduction derived from all chains after splitting each chain in half and
treating the halves as chains). For the summary of all chains merged, 
Monte Carlo standard errors (se_mean) are also reported.
}

\seealso{
\itemize{
\item \code{\link{monitor}}, which computes similar summaries but accepts an 
array of MCMC draws as its input rather than a \code{stanfit} object.
\item The RStan vignettes for more example usage.
}
}

\examples{\dontrun{
ecode <- '
  parameters {
    array[2] real<lower=0> y;
  } 
  model {
    y ~ exponential(1);
  }
'
fit <- stan(model_code = ecode)
s <- summary(fit, probs = c(0.1, 0.9))
s$summary  # all chaines merged
s$c_summary  # individual chains
}
}
\keyword{methods}
