\name{testify}
\alias{testify} 
\title{
Export User-Defined Stan Functions to \R for Testing
}
\description{
  Stan's language allows users to define their own functions in a 
  functions block at the top of a Stan program. This utility uses
  \code{\link[Rcpp]{sourceCpp}} to export those user-defined functions
  to the \code{\link{.GlobalEnv}} for testing inside \R.
}
\usage{
  testify(stanmodel) 
} 

\arguments{
  \item{stanmodel}{An object of \code{\link{stanmodel-class}} or
    \code{\link{stanfit-class}} whose associated Stan program
    contains a non-empty functions block}
}

\details{
  If a user-defined Stan function ends in \code{_rng}, then it can
  use the Boost pseudo-random number generator used by Stan. Note
  that a \code{seed} argument will be added to the \code{\link{formals}}
  on the R side that defaults to 0. Thus, by default, the same 
  pseudo-random number(s) will be returned each time the function is
  called from R. To change this behavior, an positive integer must be
  passed to the \code{seed} argument when the function is called from R.
  
  If a user-defined Stan function ends in \code{_lp}, then it can
  modify the log-probability used by Stan to evaluate Metropolis
  proposals or as an objective function for optimization. Note that a
  \code{lp__} argument will be added to the \code{\link{formals}} 
  on the R side that does not have a default argument. Thus, a 
  \code{\link{double}} precision scalar must be  passed to this 
  arguments when the function is called from R. The Stan function 
  can terminate with \code{return lp;} or can execute 
  \code{print(lp);} to verify that the calculation is correct.
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
