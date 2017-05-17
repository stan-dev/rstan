\name{expose_stan_functions}
\alias{expose_stan_functions} 
\title{
Expose user-defined Stan functions to \R for testing and simulation
}
\description{
  The Stan modeling language allows users to define their own functions in a 
  \code{functions} block at the top of a Stan program. The 
  \code{expose_stan_functions} utility function uses
  \code{\link[Rcpp]{sourceCpp}} to export those user-defined functions
  to the specified environment for testing inside \R or for doing posterior
  predictive simulations in \R rather than in the \code{generated 
  quantities} block of a Stan program.
}
\usage{
  expose_stan_functions(stanmodel, env = globalenv())
}

\arguments{
  \item{stanmodel}{
    A \code{\linkS4class{stanmodel}} object, a \code{\linkS4class{stanfit}} 
    object, a list produced by \code{\link{stanc}} or the path to a Stan 
    program (\code{.stan} file). In any of these cases, the underlying Stan 
    program should contain a non-empty \code{functions} block.
  }
 \item{env}{
   An \code{\link[base]{environment}} passed to the \code{env} argument of 
   \code{\link[Rcpp]{sourceCpp}}. By default, the global environment
   is used.
  }
}

\details{
  There are a few special types of user-defined Stan functions for which some 
  additional details are relevant:
  
  \subsection{(P)RNG functions}{
  If a user-defined Stan function ends in \code{_rng}, then it can
  use the Boost pseudo-random number generator used by Stan. When exposing
  such functions to \R, a \code{seed} argument will be added to the 
  \code{\link{formals}}. This \code{seed} argument defaults to \code{0L}, 
  but any non-negative integer can be passed as the \code{seed} the 
  \strong{first} time \strong{any} user-defined function ending in \code{_rng} 
  is called. In other words, the Boost pseudo-random number generator is 
  initialized with the given \code{seed} but is declared with the \code{static} 
  C++ keyword, meaning that it will not be reinitialized by subsequent calls 
  to user-defined functions ending in \code{_rng}.
  }
  \subsection{LP functions}{
  If a user-defined Stan function ends in \code{_lp}, then it can
  modify the log-probability used by Stan to evaluate Metropolis
  proposals or as an objective function for optimization. When exposing
  such functions to \R, a \code{lp__} argument will be added to the 
  \code{\link{formals}}. This \code{lp__} argument defaults to zero, but a 
  \code{\link{double}} precision scalar may be  passed to this argument when the 
  function is called from \R. Such a user-defined Stan function can terminate 
  with \code{return target();} or can execute \code{print(target());} to verify that 
  the calculation is correct.
  }
}
\value{
  The names of the new functions in \code{env} are returned invisibly.
}
\seealso{
  \code{\link[Rcpp]{sourceCpp}}
}
\examples{
  mc <- 
  '
  functions {
    int fibonacci(int x);
    int fibonacci(int x) {
      if (x <= 0) reject("x must be positive");
      return x <= 2 ? 1 : fibonacci(x - 1) + fibonacci(x - 2);
    }
  }
  model {}
  '
  cppcode <- stanc(model_code = mc, model_name = "Fibonacci")
  \dontrun{
  expose_stan_functions(cppcode)
  }
}
