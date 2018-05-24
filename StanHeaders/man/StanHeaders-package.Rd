\name{StanHeaders-package}
\alias{StanHeaders-package}
\alias{StanHeaders}
\docType{package}
\title{\packageTitle{StanHeaders}}
\description{\packageDescription{StanHeaders}}
\author{\packageAuthor{StanHeaders}}
\section{Maintainer}{\packageMaintainer{StanHeaders}}
\details{
  There are no R functions, only C++ headers. To use the Stan Math Library,
  it is sufficient to include \code{LinkingTo: StanHeaders} in the 
  \code{DESCRIPTION} file of another package. The Stan Math Library is
  sufficient for many purposes. If, in addition, another package needs to
  utilize the MCMC, optimization, or parsing facilities of
  \pkg{StanHeaders}, then it is also necessary to include the \code{src}
  directory of \pkg{StanHeaders} in the other package's 
  \code{PKG_CPPFLAGS}. For an example of such a \code{Makevars}, see
  \url{https://raw.githubusercontent.com/stan-dev/rstan/develop/rstan/rstan/src/Makevars}.
}
\references{
  http://mc-stan.org
}
\keyword{ package }
\examples{
  # No examples because there are no R functions.
}
