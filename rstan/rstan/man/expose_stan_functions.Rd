\name{expose_stan_functions}
\alias{expose_stan_functions}
\alias{get_stream}
\alias{get_rng}
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
  expose_stan_functions(stanmodel, includes = NULL, 
                        show_compiler_warnings = FALSE, ...)
  get_rng(seed = 0L)
  get_stream()
}

\arguments{
  \item{stanmodel}{
    A \code{\linkS4class{stanmodel}} object, a \code{\linkS4class{stanfit}} 
    object, a list produced by \code{\link{stanc}} or the path to a Stan 
    program (\code{.stan} file). In any of these cases, the underlying Stan 
    program should contain a non-empty \code{functions} block.
  }
  \item{includes}{If not \code{NULL} (the default), then a character vector of
    length one (possibly containing one or more \code{"\n"}) of the form 
    \code{'#include "/full/path/to/my_header.hpp"'}, which will be inserted
    into the C++ code in the model's namespace and can be used to provide definitions 
    of functions that are declared but not defined in \code{stanmodel}
  }
  \item{show_compiler_warnings}{Logical scalar defaulting to \code{FALSE} that
    controls whether compiler warnings, which can be numerous and have never
    been relevant, are shown}
  \item{seed}{An integer vector of length one indicating the state of Stan's
    pseudo-random number generator, but Stan's pseudo-random number generator
    is not used anymore in Stan functions that have been exposed to R}
  \item{\dots}{
    Further arguments passed to \code{\link[Rcpp]{sourceCpp}}.
  }
}

\details{
  The \code{expose_stan_functions} function requires as much compliance with
  the C++14 standard as is implemented in the RTools toolchain for Windows.

  Exposed user-defined Stan functions have an additional \code{pstream__}
  argument that is added to the \code{\link{formals}} that previously was
  used to manage the output stream if the function called \code{reject} or
  \code{print}. However, the output stream is not handled internally and the
  additional \code{pstream__} argument (which defaults to \code{NULL}) does
  nothing.

  There are a few special types of user-defined Stan functions for which some 
  additional details are relevant:
  
  \subsection{(P)RNG functions}{
  If a user-defined Stan function ends in \code{_rng}, then it now uses R's
  pseudo-random number stream rather than Boost's. However, for backward
  compatibility, a \code{base_rng__} argument is added to the \code{\link{formals}}
  of the R version of the function, which defaults to \code{NULL} and is not
  used. To reproduce the results, call \code{\link{set.seed}} as usual.
  }
  \subsection{LP functions}{
  If a user-defined Stan function ends in \code{_lp}, then it can
  modify the log-probability used by Stan to evaluate an objective function. When 
  exposing such functions to \R, a \code{lp__} argument will be added to the 
  \code{\link{formals}}. This \code{lp__} argument defaults to zero, but a 
  \code{\link{double}} precision scalar may be passed to this argument when the 
  function is called from \R. Such a user-defined Stan function can terminate 
  with \code{return target();} or can execute \code{print(target());} to verify that 
  the calculation is correct.
  }
}
\value{
  The names of the new functions in \code{env} are returned invisibly.
}
\seealso{
  \code{\link[Rcpp]{sourceCpp}} and the section in the Stan User Manual on
  user-defined functions
}
\examples{
\dontrun{
model_code <-
  '
  functions {
    real standard_normal_rng() {
      return normal_rng(0,1);
   }
  }
'
expose_stan_functions(stanc(model_code = model_code))
# For reproducibility, first call set.seed() with some integer argument
standard_normal_rng()
}
}
