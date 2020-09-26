\name{stanFunction}
\alias{stanFunction}
\title{Compile and Call a Stan Math Function
}
\description{
Call a function defined in the Stan Math Library from R using this wrapper around 
\code{\link[Rcpp]{cppFunction}}.
}
\usage{
  stanFunction(function_name, ..., env = parent.frame(), rebuild = FALSE,
               cacheDir = getOption("rcpp.cache.dir", tempdir()), 
               showOutput = verbose, verbose = getOption("verbose"))
}

\arguments{
  \item{function_name}{
    A \code{\link{character}} vector of length one that is the unscoped basename of a 
    C++ function under the \code{prim/} directory of the Stan Math Library that you 
    would like to evaluate. Functions (such as \code{integrate_1d}) of other functions
    are not permitted and neither are functions (such as \code{reject}) of characters.
  }
  \item{\dots}{
    Further arguments that are passed to \code{function_name} in \code{tag = value}
    form, which are passed to \code{function_name} by  \emph{position}. See the Details 
    and Examples sections.
  }
  \item{env,rebuild,cacheDir,showOutput,verbose}{
    The same as in \code{\link[Rcpp]{cppFunction}}
  }
}
\details{
  The \code{stanFunction} function essentially compiles and
  evaluates a C++ function of the form
  \preformatted{auto function_name(...) \{ return stan::math::function_name(...); \}}
  It is essential to pass all arguments to \code{function_name} through the \dots
  in order for the C++ wrapper to know what the argument types are. The mapping
  between R types and Stan types is
  \tabular{lr}{
    \bold{R type}  \tab \bold{Stan type} \cr
    \code{double}  \tab \code{real} \cr
    \code{integer} \tab \code{int}  \cr
    \code{vector}  \tab \code{vector} \cr
    \code{matrix(*, nrow = 1)} \tab \code{row_vector}\cr
    \code{matrix}  \tab \code{matrix}
  }
  and, in addition, lists of the aforementioned R types map
  to arrays of Stan types and thus must not be ragged if they are nested. The
  Stan version of the function is called with arguments specified by position,
  i.e. in the order that they appear in the \dots. However, the R wrapper
  function has arguments whose names are the same as the names passed through
  the \dots.
}
\value{
  The result of \code{function_name} evaluated at the arguments
  that are passed through the \dots, which could be of various
  R types. It also has the side effect of defining a function
  named \code{function_name} in the environment given by the
  \code{env} argument that can subsequently be called with
  inputs of the same type (but not necessarily the same value)
  that were passed through the \dots.
}
\examples{
  files <- dir(system.file("include", "stan", "math", "prim",
                           package = "StanHeaders"), 
               pattern = "hpp$", recursive = TRUE)
  functions <- sub("\\\\.hpp$", "", 
                   sort(unique(basename(files[dirname(files) != "."]))))
  length(functions) # you could call most of these Stan functions
  
  \dontrun{
    log(sum(exp(exp(1)), exp(pi))) # true value
    
    stanFunction("log_sum_exp", x = exp(1), y = pi)
    args(log_sum_exp) # now exists in .GlobalEnv
    log_sum_exp(x = pi, y = exp(1))
    
    # but log_sum_exp() was not defined for a vector or matrix
    x <- c(exp(1), pi)
    try(log_sum_exp(x))
    stanFunction("log_sum_exp", x = x) # now it is
    
    # log_sum_exp() is now also defined for a matrix
    log_sum_exp(as.matrix(x))
    log_sum_exp(t(as.matrix(x)))
    log_sum_exp(rbind(x, x))
    
    # but log_sum_exp() was not defined for a list
    try(log_sum_exp(as.list(x)))
    stanFunction("log_sum_exp", x = as.list(x)) # now it is

    # in rare cases, passing a nested list is needed
    stanFunction("dims", x = list(list(1:3)))
    
    # nullary functions work but are not that interesting
    stanFunction("negative_infinity")
    
    # PRNG functions work by adding a seed argument
    stanFunction("lkj_corr_rng", K = 3L, eta = 1)
    args(lkj_corr_rng) # has a seed argument
  }
}
