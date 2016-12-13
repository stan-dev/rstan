\name{stan_model}
\alias{stan_model}
\title{Construct a Stan model} 
\description{
  Construct an instance of S4 class \code{stanmodel} from a model 
  specified in Stan's modeling language. A \code{stanmodel} object 
  can then be used to draw samples from the model. The Stan program
  (the model expressed in the Stan modeling language) is first translated to
  C++ code and then the C++ code for the model plus other auxiliary 
  code is compiled into a dynamic shared object (DSO) and then loaded. 
  The loaded DSO for the model can be executed to draw samples, allowing
  inference to be performed for the model and data.
}

\usage{
  stan_model(file, model_name = "anon_model", 
             model_code = "", stanc_ret = NULL, 
             boost_lib = NULL, eigen_lib = NULL, 
             save_dso = TRUE, verbose = FALSE, 
             auto_write = rstan_options("auto_write"), 
             obfuscate_model_name = TRUE, 
             allow_undefined = FALSE, includes = NULL)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{file}{A character string or a connection that \R supports 
    specifying the Stan model specification in Stan's modeling language. 
  }

  \item{model_name}{A character string naming the model; defaults 
    to \code{"anon_model"}. However, the model name will be derived from 
    \code{file} or \code{model_code} (if \code{model_code} is the name of a
    character string object) if \code{model_name} is not specified.}

  \item{model_code}{Either a character string containing the model 
    specification or the name of a character string object in the workspace.
    This is an alternative to specifying the model via the \code{file}
    or \code{stanc_ret} arguments. 
  }
  \item{stanc_ret}{A named list returned from a previous call to 
    the \code{\link{stanc}} function. The list can be used to specify the model  
    instead of using the \code{file} or \code{model_code} arguments. 
  }
  \item{boost_lib}{The path to a version of the Boost C++ library to
    use instead of the one in the \pkg{BH} package.}
  \item{eigen_lib}{The path to a version of the Eigen C++ library to
    use instead of the one in the \pkg{RcppEigen} package.} 
  \item{save_dso}{Logical, defaulting to \code{TRUE}, indicating 
    whether the  dynamic shared object (DSO) compiled from the C++ code for the 
    model will be saved or not. If \code{TRUE}, we can draw samples from
    the same model in another \R session using the saved DSO (i.e., 
    without compiling the C++ code again).}
  \item{verbose}{Logical, defaulting to \code{FALSE}, indicating whether
    to report additional intermediate output to the console, 
    which might be helpful for debugging.}
  \item{auto_write}{Logical, defaulting to the value of 
    \code{rstan_options("auto_write")}, indicating whether to write the
    object to the hard disk using \code{\link{saveRDS}}. Although this argument
    is \code{FALSE} by default, we recommend calling 
    \code{rstan_options("auto_write" = TRUE)} in order to avoid unnecessary 
    recompilations. If \code{file} is supplied and its \code{\link{dirname}} 
    is writable, then the object will be written to that same directory, 
    substituting a \code{.rda} extension for the \code{.stan} extension. 
    Otherwise, the object will be written to the \code{\link{tempdir}}.}
  \item{obfuscate_model_name}{A logical scalar that is \code{TRUE} by default and
    passed to \code{\link{stanc}}.}
  \item{allow_undefined}{A logical scalar that is \code{FALSE} by default and
    passed to \code{\link{stanc}}.}
  \item{includes}{If not \code{NULL} (the default), then a character vector of
    length one (possibly containing \code{"\n"}) of the form 
    \code{'#include "/full/path/to/my_header.hpp"'}, which will be inserted
    into the C++ code in the model's namespace and can be used to provide definitions 
    of functions that are declared but not defined in \code{file} or
    \code{model_code} when \code{allow_undefined = TRUE}}
}

\details{
  If a previously compiled \code{stanmodel} exists on the hard drive, its validity
  is checked and then returned without recompiling. The most common form of 
  invalidity seems to be Stan code that ends with a \code{\}} rather than a blank
  line, which causes the hash checker to think that the current model is different
  than the one saved on the hard drive. To avoid reading previously 
  compiled \code{stanmodel}s from the hard drive, supply the \code{stanc_ret}
  argument rather than the \code{file} or \code{model_code} arguments.
  
  There are three ways to specify the model's code for \code{stan_model}: 
  \enumerate{ 
  \item parameter \code{model_code}: a character string containing the 
        Stan model specification,
  \item parameter \code{file}: a file name (or a connection) from
       which to read the Stan model specification, or 
  \item parameter \code{stanc_ret}: a list returned by \code{stanc}
        to be reused.
  } 
}
\value{
  An instance of S4 class \code{\linkS4class{stanmodel}} that can be
  passed to the \code{\link{sampling}}, \code{\link{optimizing}}, and 
  \code{\link{vb}} functions.
}
\references{
  The Stan Development Team 
  \emph{Stan Modeling Language User's Guide and Reference Manual}. 
  \url{http://mc-stan.org/}. 
}
\seealso{
  \code{\linkS4class{stanmodel}} for details on the class. 
  
  \code{\link{sampling}}, \code{\link{optimizing}}, and \code{\link{vb}}, 
  which take a \code{stanmodel} object as input, for estimating the model 
  parameters.
  
  More details on Stan, including the full user's guide and reference 
  manual, can be found at \url{http://mc-stan.org/}. 
}
\examples{\dontrun{
stancode <- 'data {real y_mean;} parameters {real y;} model {y ~ normal(y_mean,1);}'
mod <- stan_model(model_code = stancode)
fit <- sampling(mod, data = list(y_mean = 0))
fit2 <- sampling(mod, data = list(y_mean = 5))
}}
\keyword{rstan}
