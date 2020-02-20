\name{loo_moment_match.stanfit}
\alias{loo_moment_match.stanfit}
\alias{loo_moment_match,stanfit-method}
\alias{loo_moment_match}

\title{Moment matching for efficient approximate leave-one-out cross-validation (LOO)}

\description{
A \code{\link[loo_moment_match]{loo_moment_match}} method that
is customized for stanfit objects.
The \code{loo_moment_match} method for stanfit objects ---a wrapper around the
\code{\link[loo_moment_match]{loo_moment_match.default}} (\pkg{loo} package)--- updates
a loo object using moment matching (Paananen et al., 2019).
}


\usage{
\method{loo_moment_match}{stanfit}(x,
    loo = loo,
    \dots) 
}
\arguments{
\item{x}{An object of S4 class \code{stanfit}.}

\item{loo}{A loo object that is modified.}

\item{\dots}{Further arguments.}
}

\value{
The \code{loo_moment_match()} methods return an updated \code{loo} object.
}




\references{
Paananen, T., Piironen, J., Buerkner, P.-C., Vehtari, A. (2019).
Pushing the Limits of Importance Sampling through Iterative Moment Matching.
\href{https://arxiv.org/abs/1906.08850/}{preprint arXiv:1906.08850}
}
\seealso{
\code{\link[=loo]{loo()}}, \code{\link[=loo]{loo_moment_match()}}
}

