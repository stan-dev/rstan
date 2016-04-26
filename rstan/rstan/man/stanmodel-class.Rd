
\name{stanmodel-class} 
\Rdversion{1.1}

\docType{class}
\alias{stanmodel-class}
%% \alias{sampling,stanmodel-method}
\alias{get_cppcode}
\alias{get_cxxflags}
\alias{get_cppcode,stanmodel-method}
\alias{get_cxxflags,stanmodel-method}
\alias{show,stanmodel-method}

\title{Class representing model compiled from C++}

\description{
  A \code{stanmodel} object represents the model compiled from C++ code. 
  The \code{sampling} method defined in this class may be used to 
  draw samples from the model and \code{optimizing} method is for
  obtaining a point estimate by maximizing the log-posterior.
} 

\section{Objects from the Class}{
  Instances of \code{stanmodel} are usually created by calling 
  function \code{stan_model} or function \code{stan}.
}


\section{Slots}{
  \describe{
    \item{\code{model_name}:}{The model name, an object of type \code{character}.} 
    \item{\code{model_code}:}{The Stan model specification, an object of type \code{character}.} 
    \item{\code{model_cpp}:}{Object of type \code{list}
      that includes the C++ code for the model.} 
    \item{\code{mk_cppmodule}:}{A function to return a RCpp module. This function will be
      called in function \code{sampling} and \code{optimzing} with one
      argument (the instance of \code{stanmodel} itself).} 
    \item{\code{dso}:}{Object of S4 class \code{cxxdso}. The container for the dynamic 
      shared objects compiled from the C++ code of the model, returned from function
      \code{cxxfunction} in package \pkg{inline}.}
  }
}

\section{Methods}{
  \describe{
    \item{\code{show}}{\code{signature(object = "stanmodel")}: print the Stan model specification.}
    \item{\code{vb}}{\code{signature(object = "stanmodel")}: use the variational Bayes algorithms.}
    \item{\code{sampling}}{\code{signature(object = "stanmodel")}: draw samples for
      the model (see \code{\link{sampling}}).}
    \item{\code{optimizing}}{\code{signature(object = "stanmodel")}: obtain a point 
      estimate by maximizing the posterior (see \code{\link{optimizing}}).}
    \item{\code{get_cppcode}}{\code{signature(object = "stanmodel")}: return the C++ code for the model
      as a character string. This is part of the C++ code that is compiled to the dynamic 
      shared object for the model.}
    \item{\code{get_cxxflags}}{\code{signature(object = "stanmodel")}: return the \code{CXXFLAGS} 
      used for compiling the model. The returned string is like \code{CXXFLAGS = -O3}.}
  }
}

\note{
  Objects of class \code{stanmodel} can be saved for use across 
  \R sessions only if \code{save_dso = TRUE} is set during calling
  functions that create \code{stanmodel} objects (e.g., \code{stan} 
  and \code{stan_model}). 

  Even if \code{save_dso = TRUE}, the model cannot be loaded on
  a platform (operating system, 32 bits or 64 bits, etc.) that differs from
  the one on which it was compiled.
}

\seealso{
\code{\link{stan_model}}, \code{\link{stanc}}
\code{\link{sampling}}, \code{\link{optimizing}}, \code{\link{vb}}
}
\examples{
showClass("stanmodel")
}
\keyword{classes}
