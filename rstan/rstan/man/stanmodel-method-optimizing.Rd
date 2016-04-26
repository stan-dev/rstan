\name{optimizing-methods}
\docType{methods}
\alias{optimizing} 
\alias{optimizing-methods}
\alias{optimizing,stanmodel-method}

\title{Obtain a point estimate by maximizing the joint posterior}
\description{
  Obtain a point estimate by maximizing the joint posterior 
  from the model defined by class \code{stanmodel}. 
  This method is a generic function of the S4 class \code{stanmodel}. 
}

\usage{
  %% optimizing(object, \dots)  
  \S4method{optimizing}{stanmodel}(object, data = list(), 
    seed = sample.int(.Machine$integer.max, 1), init = 'random', 
    check_data = TRUE, sample_file = NULL, 
    algorithm = c("LBFGS", "BFGS", "Newton"),
    verbose = FALSE, hessian = FALSE, as_vector = TRUE, 
    draws = 0, constrained = TRUE, \dots)   
} 

\section{Methods}{
  \describe{
    \item{optimizing}{\code{signature(object = "stanmodel")}}{
      %%  ~~describe this method here~~
      Call Stan's optimization methods to obtain a point estimate 
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
    See the Note section in \code{\link{stan}}.}
    
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
    See the Note section in \code{\link{stan}}.}
    
  \item{sample_file}{A character string of file name for specifying where to 
    write samples for \emph{all} parameters and other saved quantities. 
    If not provided, files are not created. When the folder specified 
    is not writable, \code{tempdir()} is used.}
    
  \item{algorithm}{One of \code{"Newton"}, \code{"BFGS"}, 
    and \code{"LBFGS"} (the default) indicating which optimization algorithm 
    to use.}
    
  \item{verbose}{\code{TRUE} or \code{FALSE} (the default): flag indicating 
    whether to print intermediate output from Stan on the console, which might
    be helpful for model debugging.}
     
  \item{hessian}{\code{TRUE} or \code{FALSE} (the default): flag indicating
     whether to calculate the Hessian (via numeric differentiation of the
     gradient function in the unconstrained parameter space).}
     
  \item{as_vector}{\code{TRUE} (the default) or \code{FALSE}: flag indicating
     whether a vector is used to store the point estimate found. A list can be 
     used instead by specifying it to be \code{FALSE}.}
     
  \item{draws}{A non-negative integer (that defaults to zero) indicating how
     many times to draw from a multivariate normal distribution whose parameters
     are the mean vector and the inverse negative Hessian in the unconstrained
     space.}
     
  \item{constrained}{A logical scalar indicating, if \code{draws > 0}, whether
     the draws should be transformed to the constrained space defined in the
     parameters block of the Stan program. Defaults to \code{TRUE}.}
     
  \item{\dots}{Other optional parameters:
    \itemize{
      \item \code{iter} (\code{integer}), the maximum number of iterations, 
      defaulting to 2000.
      \item \code{save_iterations} (logical), a flag indicating whether to save 
      the iterations, defaulting to \code{FALSE}.
      \item \code{refresh} (\code{integer}), the number of interations between 
      screen updates, defaulting to 100.
      \item \code{init_alpha} (\code{double}), for BFGS and LBFGS, 
      the line search step size for first iteration, defaulting to 0.001.
      \item \code{tol_obj} (\code{double}), for BFGS and LBFGS, 
      the convergence tolerance on changes in objective function value, 
      defaulting to 1e-12.
      \item \code{tol_rel_obj} (\code{double}), for BFGS and LBFGS, 
      the convergence tolerance on relative changes in objective function value,
      defaulting to 1e4.
      \item \code{tol_grad} (\code{double}), for BFGS and LBFGS, 
      the convergence tolerance on the norm of the gradient, defaulting to 1e-8.
      \item \code{tol_rel_grad} (\code{double}), for BFGS and LBFGS, 
      the convergence tolerance on the relative norm of the gradient, defaulting
      to 1e7.
      \item \code{tol_param} (\code{double}), for BFGS and LBFGS, 
      the convergence tolerance on changes in parameter value, defaulting to 1e-8.
      \item \code{history_size} (\code{integer}), for LBFGS, 
      the number of update vectors to use in Hessian approximations, 
      defaulting to 5.
    }
    Refer to the manuals for both CmdStan and Stan for more details.
  }
}
\value{
   If the optimization is done successfully, a list with named components:
   \item{par}{The point estimate found. Its form (vector or list) 
     is determined by the \code{as_vector} argument.}
   \item{value}{The value of the log-posterior (up to an additive constant, 
     the \code{"lp__"} in Stan) corresponding to \code{par}.}
   \item{hessian}{The Hessian matrix if \code{hessian} is \code{TRUE}}
   \item{theta_tilde}{If \code{draws > 0}, the matrix of parameter draws
    in the constrained or unconstrained space, depending on the value of 
    the \code{constrained} argument.}
   \item{log_p}{If \code{draws > 0}, a vector of length draws that contains
     the value of the log-posterior evaluated at each row of \code{theta_tilde}.}
   \item{log_g}{If \code{draws > 0}, a vector of length draws that contains
     the value of the logarithm of the multivariate normal density evaluated
     at each row of \code{theta_tilde}.}

  If the optimization is not completed for reasons such as feeding wrong data,
  it returns \code{NULL}. 
} 

\seealso{
  \code{\linkS4class{stanmodel}}
}

\examples{\dontrun{
m <- stan_model(model_code = 'parameters {real y;} model {y ~ normal(0,1);}')
f <- optimizing(m, hessian = TRUE)
}} 
