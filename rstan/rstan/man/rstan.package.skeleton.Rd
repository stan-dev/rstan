\name{rstan.package.skeleton}
\alias{rstan.package.skeleton}
\title{Create a Skeleton for a New Source Package with Stan Programs}
\description{This function is very similar to \code{\link[utils]{package.skeleton}}
  but is designed for source packages that want to include Stan Programs
  that can be built into binary versions.
}

\usage{
rstan.package.skeleton(name = "anRpackage", list = character(),
                       environment = .GlobalEnv, path = ".", 
                       force = FALSE, code_files = character(),
                       stan_files = character())
}

\arguments{
  \item{name, list, environment, path, force, code\_files}{Same is in 
    \code{\link[utils]{package.skeleton}}}
  \item{stan\_files}{A character vector with paths to .stan files
    to include in the package. Otherwise similar to \code{code_files}.}
}

\details{
  This function first calls \code{\link[utils]{package.skeleton}} and
  then adds the files listed in \code{stan_files} to an exec directory.
  Finally, it downloads several files from the rstanarm GitHub repository
  to facilitate building the resulting package. Note that rstanarm is
  licensed under the GPL >= 3, so package builders who do not want to be
  governed by that license should not use the downloaded files that contain
  R code. Otherwise, it may be worth considering whether it would be easier
  to include your .stan programs and supporting R code in the \pkg{rstanarm}
  package.
} 

\value{
  Used for its side-effects
} 

\seealso{
  \url{https://github.com/stan-dev/rstanarm/} for an example of how to
  create a package that contains pre-compilable .stan programs
} 
\examples{
  # See example(package.skeleton)
}
