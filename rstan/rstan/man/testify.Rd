\name{testify}
\alias{testify} 
\title{
Export User-Defined Stan Functions to \R for Testing
}
\description{
  Stan's language allows users to define their own functions in a 
  functions block at the top of a Stan program. This utility uses
  \code{\link[Rcpp]{sourceCpp}} to export those user-defined functions
  to the \code{\link{.GlobalEnv}} for testing inside \R or for doing
  posterior predictive simulations outside of the 
  \code{generated quantities} block of a .stan program
}
\usage{
  testify(stanmodel) 
} 

\arguments{
  \item{stanmodel}{An object of \code{\link{stanmodel-class}},
    \code{\link{stanfit-class}}, a list produced by \code{\link{stanc}} 
    or a path to a .stan file. In any of these cases, the underlying 
    Stan program should contain a non-empty \code{functions} block}
}

\details{
  If a user-defined Stan function ends in \code{_rng}, then it can
  use the Boost pseudo-random number generator used by Stan. Note
  that a \code{seed} argument will be added to the \code{\link{formals}}
  on the R side that defaults to 0L but any non-negative integer
  can be passed as the \code{seed} the \strong{first} time \strong{any}
  user-defined function ending in \code{_rng} is called. In other words,
  the Boost pseudo-random number generator is initialized with the given
  \code{seed} but is declared with the \code{static} C++ keyword, meaning
  that it will not be reinitialized by subsequent calls to user-defined
  functions ending in \code{_rng}.
  
  If a user-defined Stan function ends in \code{_lp}, then it can
  modify the log-probability used by Stan to evaluate Metropolis
  proposals or as an objective function for optimization. Note that a
  \code{lp__} argument will be added to the \code{\link{formals}} 
  on the R side that defaults to zero. A \code{\link{double}} precision 
  scalar may be  passed to this argument when the function is called from R. 
  Such a user-defined Stan function can terminate with \code{return lp__;} or 
  can execute \code{print(lp__);} to verify that the calculation is correct.
}
\value{
  \code{NULL} invisibly but the side effect is that new functions
  are defined in \code{\link{.GlobalEnv}}
}
\seealso{
  \code{\link[Rcpp]{sourceCpp}}
}
\examples{
example(read_stan_csv)
testify(fit) # correct useage but throws warning because no C++ code
}
