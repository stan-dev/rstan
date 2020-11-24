\name{nlist}
\alias{nlist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Created named lists} 
\description{
Create a named list using specified names or, if names are omitted, using the 
names of the objects in the list. The code ]\code{list(a = a, b = b)} becomes 
\code{nlist(a,b)} and \code{list(a = a, b = 2)} becomes \code{nlist(a, b = 2)}, 
etc. This is convenient when creating the list of data to pass to Stan. 
}
\usage{
  nlist(...)
}
\arguments{
  \item{...}{The objects to include in the list.}
}
\value{
  A named list.
}
\examples{
# All variables already defined
x <- 1
y <- 2
nlist(x, y)

# Define some variables in the call and take the rest from the environment
nlist(x, y, z = 3)
}
\keyword{rstan} 
