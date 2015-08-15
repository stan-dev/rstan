\name{StanHeaders-package}
\alias{StanHeaders-package}
\alias{StanHeaders}
\docType{package}
\title{
  The Stan C++ Library Headers
}
\description{
  The Stan project implements Markov Chain Monte Carlo and optimization for 
  statistical models. This R package provides the Stan header files so that
  they can easily be used by other R packages.
}
\details{
\tabular{ll}{
Package: \tab StanHeaders\cr
Type: \tab Package\cr
Version: \tab 2.7.0\cr
Date: \tab 2015-07-16\cr
License: \tab BSD 3 clause\cr
}
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
\author{
  Ben Goodrich and Jiqiang Guo for the R package and the Stan Development 
  Team for Stan.
  Maintainer: Ben Goodrich <benjamin.goodrich@columbia.edu>
}
\references{
  http://mc-stan.org
}
\keyword{ package }
\examples{
  # No examples because there are no R functions.
}
