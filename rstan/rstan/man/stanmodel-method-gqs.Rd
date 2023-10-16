\docType{methods}
\name{gqs}
\alias{gqs} 
\alias{gqs,stanmodel-method}

\title{Draw samples of generated quantities from a Stan model}
\description{
  Draw samples from the generated quantities block of a 
  \code{\linkS4class{stanmodel}}.
}

\usage{
  \S4method{gqs}{stanmodel}(object, data = list(), draws, 
    seed = sample.int(.Machine$integer.max, size = 1L))
} 

\section{Methods}{
  \describe{
    \item{\code{object}}{\code{signature(object = "stanmodel")} 
      Evaluate the generated quantities block of a Stan program
      by supplying \code{data} and the \code{draws} output from a
      previous Stan program.
    }
  }
}

\arguments{
  \item{object}{An object of class \code{\linkS4class{stanmodel}}.}

  \item{data}{A named \code{list} or \code{environment}
    providing the data for the model or a character vector 
    for all the names of objects used as data. 
    See the \strong{Passing data to Stan} section in \code{\link{stan}}.} 

  \item{draws}{A matrix of posterior draws, typically created by
    calling \code{as.matrix} on a \code{\linkS4class{stanfit}}.}

  \item{seed}{The seed for random number generation. The default is generated 
    from 1 to the maximum integer supported by \R on the machine.
    When a seed is specified by a number, \code{as.integer} will be applied to it. 
    If \code{as.integer} produces \code{NA}, the seed is generated randomly. 
    The seed can also be specified as a character string of digits, such as
    \code{"12345"}, which is converted to integer.}
}

\value{
   An object of S4 class \code{\linkS4class{stanmodel}} representing
   the fitted results. 
} 

\seealso{
  \code{\linkS4class{stanmodel}}, \code{\linkS4class{stanfit}}, \code{\link{stan}}
}

\examples{\dontrun{
m <- stan_model(model_code = 'parameters {real y;} model {y ~ normal(0,1);}')
f <- sampling(m, iter = 300)
mc <-
'
parameters {real y;}
generated quantities {real y_rep = normal_rng(y, 1);}
'
m2 <- stan_model(model_code = mc)
f2 <- gqs(m2, draws = as.matrix(f))
f2
}}
