\name{as.array}
\alias{is.array.stanfit}
\alias{as.array.stanfit}
\alias{as.matrix.stanfit}
\alias{as.data.frame.stanfit}
\alias{dim.stanfit}
\alias{dimnames.stanfit}
\alias{names.stanfit}
\alias{names<-.stanfit}

\title{Create array, matrix, or data.frame objects from samples in a \code{stanfit} object} 


\description{The samples (without warmup) included in a \code{\linkS4class{stanfit}} 
  object can be coerced to an \code{array}, \code{matrix}, or \code{data.frame}.
  Methods are also provided for checking and setting names and dimnames.} 

\usage{
  \method{as.array}{stanfit}(x, \dots) 
  \method{as.matrix}{stanfit}(x, \dots)
  \method{as.data.frame}{stanfit}(x, \dots)
  \method{is.array}{stanfit}(x) 
  \method{dim}{stanfit}(x)
  \method{dimnames}{stanfit}(x)
  \method{names}{stanfit}(x)
  \method{names}{stanfit}(x) <- value
} 

\arguments{
  \item{x}{An object of S4 class \code{\linkS4class{stanfit}}.}
  \item{\dots}{Additional parameters that can be passed to \code{\link{extract}} 
    for extracting samples from \code{x}. For now, \code{pars}
    is the only additional parameter supported.
  }
  \item{value}{For the \code{names} replacement method, a character vector to 
    set/replace the parameter names in \code{x}.}
} 

\details{
  \code{as.array} and \code{as.matrix} can be applied to a \code{stanfit}
  object to coerce the samples without warmup to \code{array} or \code{matrix}.
  The \code{as.data.frame} method first calls \code{as.matrix} and then coerces
  this matrix to a \code{data.frame}.

  The array has three named dimensions: iterations, chains, parameters. 
  For \code{as.matrix}, all chains are combined, leaving a matrix of iterations
  by parameters.
}
\value{
  \code{as.array}, \code{as.matrix}, and \code{as.data.frame} return an array,
  matrix, and data.frame, respectively.

  \code{dim} and \code{dimnames} return the dim and dimnames of the 
  array object that could be created, while \code{names} returns the third
  element of the \code{dimnames}, which are the names of the margins of the
  posterior distribution. The \code{names} assignment method allows for 
  assigning more interpretable names to them.

  \code{is.array} returns \code{TRUE} for \code{stanfit} objects that include
  samples; otherwise \code{FALSE}. 

  When the \code{stanfit} object does not contain samples, empty objects
  are returned from \code{as.array}, \code{as.matrix}, \code{as.data.frame},
  \code{dim}, \code{dimnames}, and \code{names}. 
}
  
\seealso{
  S4 class \code{\linkS4class{stanfit}} and its method \code{\link{extract}}
} 

\examples{\dontrun{
ex_model_code <- '
  parameters {
    real alpha[2,3];
    real beta[2]; 
  } 
  model {
    for (i in 1:2) for (j in 1:3) 
      alpha[i, j] ~ normal(0, 1); 
    for (i in 1:2) 
      beta[i] ~ normal(0, 2); 
    # beta ~ normal(0, 2) // vectorized version
  } 
'

## fit the model 
fit <- stan(model_code = ex_model_code, chains = 4) 

dim(fit)
dimnames(fit)
is.array(fit) 
a <- as.array(fit)
m <- as.matrix(fit)
d <- as.data.frame(fit)
}}  
