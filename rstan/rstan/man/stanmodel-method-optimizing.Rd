\name{optimizing-methods}
\docType{methods}
\alias{optimizing} 
\alias{optimizing-methods}
\alias{optimizing,stanmodel-method}

\title{\code{optimizing}: obtain a point estimate by maximizing the joint posterior}
\description{
  Obtain a point estimate by maximizing the joint posterior 
  from the model defined by class \code{stanmodel}. 
  This method is a generic function of the S4 class \code{stanmodel}. 
}

\usage{
  %% optimizing(object, \dots)  
  \S4method{optimizing}{stanmodel}(object, data = list(), 
    seed = sample.int(.Machine$integer.max, 1), init = 'random', 
    check_data = TRUE, sample_file, 
    algorithm = c("BFGS", "Nesterov", "Newton"),
    verbose = FALSE, hessian = FALSE, as_vector = TRUE, \dots)   
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
  \item{data}{An object of class \code{list}, \code{environment}
    providing the data for the model, or a vector character strings 
    for all the names of objects used as data in the working space. 
    See the notes in \code{\link{stan}}.} 
  \item{seed}{The seed for random number generation. The default is generated
    from 1 to the maximum integer supported by \R on the machine. 
    When a seed is specified by a number, \code{as.integer} will be applied to
    it. If \code{as.integer} produces \code{NA}, the seed is generated
    randomly.  We can also specify a seed using a character string of digits,
    such as \code{"12345"}, which will be converted to integer.}
  \item{init}{One of digit \code{0}, string \code{"0"} or \code{"random"}, 
    a function that returns a list, or a named list of initial parameter
    values.
    \code{"0"}: initialize all to be zero on the unconstrained support; 
    \code{"random"}: randomly generated;
    \code{list}: a list specifying the initial values of parameters by name.
    \code{function}: a function that returns a list for specifying 
      the initial values of parameters for a chain.}
  \item{check_data}{Logical: if \code{TRUE}, the data would be preprocessed; 
    otherwise not. If the data is not checked and preprocessed, it is safe 
    to leave it to be the default \code{TRUE}. See the notes in
    \code{\link{stan}}.}
  \item{sample_file}{A character string of file name for specifying where to 
    write samples for \emph{all} parameters and other saved quantities. 
    If not provided, files are not created. When the folder specified 
    is not writable, \code{tempdir()} is used.}
  \item{algorithm}{One of "Newton", "Nesterov", and "BFGS" indicating which
    optimization algorithm is used. The default is BFGS.}
  \item{verbose}{\code{TRUE} or \code{FALSE}: flag indicating whether 
     to print intermediate output from Stan on the console.}
  \item{hessian}{\code{TRUE} or \code{FALSE} (the default): flag indicating
     whether to calculate the Hessian (via numeric differentiation of the
     gradient function in the unconstrained parameter space)}   
  \item{as_vector}{\code{TRUE} (the default) or \code{FALSE}: flag indicating
     whether a vector is used to for the point estimate found. A list can be used
     instead by specifying it to be \code{FALSE}}
  \item{\dots}{Other optional parameters. 
    \enumerate{
      \item \code{iter} (\code{integer}), the maximum number of iterations
      \item \code{save_iterations} (logical), whether to save the iterations 
      \item \code{refresh} (\code{integer})
      \item \code{stepsize} (\code{double}), for Nesterov, see the manual of Stan
      \item \code{init_alpha} (\code{double}, default to 0.0001), for BFGS, see manual of Stan
      \item \code{tol_obj} (\code{double}, default to 1e-12), for BFGS, see the manual of Stan
      \item \code{tol_grad} (\code{double}, default to 1e-8), for BFGS, see the manual of Stan
      \item \code{tol_param} (\code{double}, default to 1e-8), for BFGS, see the manual of Stan
      \item \code{tol_rel_obj} (\code{double}, default to 1e4), for BFGS, see the manual of Stan
      \item \code{tol_rel_grad} (\code{double}, default to 1e3), for BFGS, see the manual of Stan
    }
  }
}
\value{
   A list with components if the optimization is done successfully:
   \item{par}{The point estimate found. Its form (vector or list) 
     is determined by argument \code{as_vector}.}
   \item{value}{The value of the log-posterior (up to an additive constant, 
   the \code{"lp__"} in Stan) corresponding to \code{par}.}

  If the optimization is not finished for reasons such as feeding wrong data,
  it returns \code{NULL}. 
} 

\seealso{
  \code{\linkS4class{stanmodel}}
}

\examples{\dontrun{
m <- stan_model(model_code = 'parameters {real y;} model {y ~ normal(0,1);}')
f <- optimizing(m, hessian = TRUE)
}} 

