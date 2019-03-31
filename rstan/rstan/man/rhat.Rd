\name{rhat}
\alias{rhat}
\alias{ess_bulk}
\alias{ess_tail}
\title{
Convergence and efficiency diagnostics for Markov Chains
}
\description{
These functions are improved versions of the traditional
Rhat (for convergence) and Effective Sample Size (for
efficiency).
}
\usage{
rhat(sims)
ess_bulk(sims)
ess_tail(sims)
}
\arguments{
  \item{sims}{
  A two-dimensional array whose rows are equal to the
  number of iterations of the Markov Chain(s) and whose
  columns are equal to the number of Markov Chains
  (preferably more than one). The cells are the realized
  draws for a particular parameter or function of parameters.
}
}
\details{
Add details
}
\value{
The \code{rhat} function produces an estimated potential
scale factor reduction from running the chain for an infinite
number of iteraions and should be less than \eqn{1.015} at
convergence. The \code{ess_bulk} function produces an estimated
Effective Sample Size from near the center of the posterior
distribution, while the \code{ess_bulk} function produces an 
estimated Effective Sample Size from the tails of the posterior
distribution. Both Effective Sample Size estimates should be
at least \eqn{100} per Markov Chain in order to be reliable and
indicate that estimates of respective posterior quantiles are 
reliable.
}
\references{
Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and
Paul-Christian Bürkner (2019). Rank-normalization, folding, and
localization: An improved R-hat for assessing convergence of MCMC.
\emph{arXiv preprint} \code{arXiv:1903.08008}.
}
\author{
Paul-Christian Burkner and Aki Vehtari
}
\seealso{
\code{\link{monitor}}
}
\examples{
# pretend these draws came from five actual Markov Chins
sims <- matrix(rnorm(500), nrow = 100, ncol = 5)
rhat(sims)
ess_bulk(sims)
ess_tail(sims)
}
