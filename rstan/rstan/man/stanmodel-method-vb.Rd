\docType{methods}
\name{vb}
\alias{vb} 
\alias{vb,stanmodel-method}

\title{Run Stan's variational algorithm for approximate posterior sampling}
\description{
  Approximately draw from a posterior distribution using variational inference.
  
  This is still considered an experimental feature.
  We recommend calling \code{\link{stan}} or \code{\link{sampling}} for 
  final inferences and only using \code{vb} to get a rough idea of the parameter
  distributions.
}

\usage{
  %% vb(object, \dots)  
  \S4method{vb}{stanmodel}(object, data = list(), pars = NA, include = TRUE,
    seed = sample.int(.Machine$integer.max, 1), 
    init = 'random', check_data = TRUE, 
    sample_file = tempfile(fileext = '.csv'),
    algorithm = c("meanfield", "fullrank"), 
    importance_resampling = FALSE, keep_every = 1,
    \dots)
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
    See the \strong{Passing data to Stan} section in \code{\link{stan}}.}
    
  \item{pars}{If not \code{NA}, then a character vector naming parameters,
    which are included in the output if \code{include = TRUE} and excluded
    if \code{include = FALSE}. By default, all parameters are included.}
    
  \item{include}{Logical scalar defaulting to \code{TRUE} indicating
    whether to include or exclude the parameters given by the 
    \code{pars} argument. If \code{FALSE}, only entire multidimensional
    parameters can be excluded, rather than particular elements of them.}
    
  \item{seed}{The seed for random number generation. The default is generated 
    from 1 to the maximum integer supported by \R on the machine. Even if 
    multiple chains are used, only one seed is needed, with other chains having 
    seeds derived from that of the first chain to avoid dependent samples.
    When a seed is specified by a number, \code{as.integer} will be applied to it. 
    If \code{as.integer} produces \code{NA}, the seed is generated randomly. 
    The seed can also be specified as a character string of digits, such as
    \code{"12345"}, which is converted to integer.}
    
  \item{init}{Initial values specification. See the detailed documentation for 
    the init argument in \code{\link{stan}}.}
      
  \item{check_data}{Logical, defaulting to \code{TRUE}. If \code{TRUE} 
    the data will be preprocessed; otherwise not.
    See the \strong{Passing data to Stan} section in \code{\link{stan}}.}
    
  \item{sample_file}{A character string of file name for specifying where to 
    write samples for \emph{all} parameters and other saved quantities. 
    This defaults to a temporary file.}
    
  \item{algorithm}{Either \code{"meanfield"} (the default) or \code{"fullrank"}, 
    indicating which variational inference algorithm is used. The \code{"meanfield"} 
    option uses a fully factorized Gaussian for the approximation whereas the 
    \code{fullrank} option uses a Gaussian with a full-rank covariance matrix 
    for the approximation. Details and additional references are available in 
    the Stan manual.}
    
  \item{importance_resampling}{Logical scalar (defaulting to \code{FALSE}) indicating
    whether to do importance resampling to adjust the draws at the optimum
    to be more like draws from the posterior distribution}
  \item{keep_every}{Integer scalar (defaulting to 1) indicating the interval
    by which to thin the draws when \code{imporance_resampling = TRUE}}
    
  \item{\dots}{Other optional parameters:
    \itemize{
      \item \code{iter} (positive \code{integer}), the maximum number 
      of iterations, defaulting to 10000.
      \item \code{grad_samples} (positive \code{integer}), the number of samples
      for Monte Carlo estimate of gradients, defaulting to 1.
      \item \code{elbo_samples} (positive \code{integer}), the number of samples
      for Monte Carlo estimate of ELBO (objective function), defaulting to 100.
      (ELBO stands for "the evidence lower bound".)
      \item \code{eta} (\code{double}), positive stepsize weighting parameter
      for variational inference but is ignored if adaptation is engaged, which
      is the case by default.
      \item \code{adapt_engaged} (\code{logical}), a flag indicating whether to
      automatically adapt the stepsize, defaulting to \code{TRUE}.
      \item \code{tol_rel_obj} (positive \code{double}), the convergence tolerance 
      on the relative norm of the objective, defaulting to 0.01. 
      \item \code{eval_elbo} (positive \code{integer}), evaluate ELBO every Nth 
      iteration, defaulting to 100.
      \item \code{output_samples} (positive \code{integer}), number of posterior
      samples to draw and save, defaults to 1000.
      \item \code{adapt_iter} (positive \code{integer}), the maximum number of 
      iterations to adapt the stepsize, defaulting to 50. Ignored if 
      \code{adapt_engaged = FALSE}.
    }
    Refer to the manuals for both CmdStan and Stan for more details. 
  }
}
\value{
  An object of \code{\link{stanfit-class}}.
} 

\seealso{
  \code{\linkS4class{stanmodel}}
 
  The manuals of CmdStan and Stan. 
}

\references{
  The Stan Development Team 
  \emph{Stan Modeling Language User's Guide and Reference Manual}. 
  \url{https://mc-stan.org}.

  The Stan Development Team
  \emph{CmdStan Interface User's Guide}.
  \url{https://mc-stan.org}.
}

\examples{\dontrun{
m <- stan_model(model_code = 'parameters {real y;} model {y ~ normal(0,1);}')
f <- vb(m)
}} 
