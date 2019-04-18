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
\value{
The \code{rhat} function produces R-hat convergence diagnostic, which
compares the between- and within-chain estimates for model parameters
and other univariate quantities of interest. If chains have not mixed
well (ie, the between- and within-chain estimates don't agree), R-hat is
larger than 1. We recommend running at least four chains by default and
only using the sample if R-hat is less than 1.01. Stan reports R-hat
which is the maximum of rank normalized split-R-hat and rank normalized
folded-split-R-hat, which works for thick tailed distributions and is
sensitive also to differences in scale.

The \code{ess_bulk} function produces an estimated Bulk Effective Sample
Size (bulk-ESS) using rank normalized draws. Bulk-ESS is useful measure
for sampling efficiency in the bulk of the distribution (related e.g. to
efficiency of mean and median estimates), and is well defined even if
the chains do not have finite mean or variance.

The \code{ess_bulk} function produces an estimated Tail Effective Sample
Size (tail-ESS) by computing the minimum of effective sample sizes for
5\% and 95\% quantiles. Tail-ESS is useful measure for sampling
efficiency in the tails of the distribution (related e.g. to efficiency
of variance and tail quantile estimates).

Both bulk-ESS and tail-ESS should be at least \eqn{100} (approximately)
per Markov Chain in order to be reliable and indicate that estimates of
respective posterior quantiles are reliable.
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
