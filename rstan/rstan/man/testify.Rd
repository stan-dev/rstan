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
