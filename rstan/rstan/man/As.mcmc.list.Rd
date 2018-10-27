\name{As.mcmc.list}
\alias{As.mcmc.list}
\title{Create an mcmc.list from a stanfit object} 
\description{
  Create an \code{\link[coda]{mcmc.list}} (\pkg{coda}) 
  from a \code{stanfit} object.
}
\usage{
  As.mcmc.list(object, pars, include = TRUE, ...)
}
\arguments{
  \item{object}{object of class \code{"stanfit"}}
  \item{pars}{optional character vector of parameters to include}
  \item{include}{logical scalar indicating whether to include (the default) or 
    exclude the parameters named in \code{pars}}
  \item{\dots}{unused}
}
\value{
  An object of class \code{\link[coda]{mcmc.list}}.
}
