\name{vb-methods}
\docType{methods}
\alias{vb} 
\alias{vb-methods}
\alias{vb,stanmodel-method}

\title{Run Stan's Variational Algorithm for Approximate Posterior Sampling}
\description{
  Approximately draw from a posterior distribution using variational inference.
  We recommend calling \code{\link{stan}} or \code{\link{sampling}} for 
  final inferences and only using \code{vb} to get a rough idea of the parameter
  distributions.
}

\usage{
  %% vb(object, \dots)  
  \S4method{vb}{stanmodel}(object, data = list(), pars = NA, include = TRUE,
    seed = sample.int(.Machine$integer.max, 1),
    check_data = TRUE, sample_file = tempfile(fileext = '.csv'),
    algorithm = c("meanfield", "fullrank"), \dots)
} 

\section{Methods}{
  \describe{
    \item{vb}{\code{signature(object = "stanmodel")}}{
      %%  ~~describe this method here~~
      Call Stan's variational Bayes methods
      for the model defined by S4 class \code{stanmodel} 
      given the data, initial values, etc. 
    }
  }
}

\arguments{
  \item{object}{An object of class \code{\linkS4class{stanmodel}}.}
  \item{data}{A named \code{list} or \code{environment}
    providing the data for the model or a character vector 
    for all the names of objects used as data. 
    See the notes in \code{\link{stan}}.}
  \item{pars}{If not \code{NA}, then a character vector naming parameters,
    which are included in the output if \code{include = TRUE} and excluded
    if \code{include = FALSE}. By default, all parameters are included.}
  \item{include}{Logical scalar indicating whether to include (the default)
    or exclude the parameters named in \code{pars}}
  \item{seed}{The seed for random number generation. The default is generated
    from 1 to the maximum integer supported by \R on the machine. 
    When a seed is specified by a number, \code{as.integer} will be applied to
    it. If \code{as.integer} produces \code{NA}, the seed is generated
    randomly.  We can also specify a seed using a character string of digits,
    such as \code{"12345"}, which will be converted to integer.}
  \item{check_data}{Logical: if \code{TRUE}, the data would be preprocessed; 
    otherwise not. If the data is not checked and preprocessed, it is safe 
    to leave it to be the default \code{TRUE}. See the notes in
    \code{\link{stan}}.}
  \item{sample_file}{A character string of file name for specifying where to 
    write samples for \emph{all} parameters and other saved quantities. 
    This defaults to a temporary file.}
  \item{algorithm}{One of "meanfield" and "fullrank" indicating which
    variational inference algorithm is used.  meanfield: mean-field 
    approximation; fullrank: full-rank covariance.  The default is 'meanfield'.}
  \item{\dots}{Other optional parameters, refer to the manuals for both
    CmdStan and Stan. 
    \enumerate{
      \item \code{iter} (positive \code{integer}), the maximum number of iterations,
      defaults to 10000.
      \item \code{grad_samples} (positive \code{integer}), the number of samples
      for Monte Carlo enumerate of gradients, defaults to 1.
      \item \code{elbo_samples} (positive \code{integer}), the number of samples
      for Monte Carlo estimate of ELBO (objective function), defaults to 100.
      (ELBO stands for "the evidence lower bound".)
      \item \code{eta} (\code{double}), positive stepsize weighting parameters
      for variational iteration but is ignored if adaptation is engaged, which
      is the case by default
      \item \code{adapt_engaged} (\code{logical}), flag indicating whether to
      automatically adapt the stepsize and defaults to \code{TRUE}
      \item \code{tol_rel_obj} (positive \code{double}), convergence tolerance on the 
      relative norm of the objective, defaults to 0.01. 
      \item \code{eval_elbo} (positive \code{integer}), evaluate ELBO every Nth iteration,
      defaults to 100.
      \item \code{output_samples} (positive \code{integer}), number of posterior
      samples to draw and save, defaults to 1000.
      \item \code{adapt_iter} (positive \code{integer}), number of iterations to
        adapt the stepsize if \code{adapt_engaged = TRUE} and ignored otherwise 
    }
  }
}
\value{
  An object of \code{\link{stanfit-class}}
} 

\seealso{
  \code{\linkS4class{stanmodel}}
 
  The manuals of CmdStan and Stan. 
}

\examples{\dontrun{
m <- stan_model(model_code = 'parameters {real y;} model {y ~ normal(0,1);}')
f <- vb(m)
}} 

