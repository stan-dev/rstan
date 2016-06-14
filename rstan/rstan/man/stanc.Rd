\name{stanc}
\alias{stanc}
\alias{stanc_builder}
\docType{package}
\title{
Translate Stan model specification to C++ code
}
\description{
  Translate a model specification in Stan code to C++ code, which can 
  then be compiled and loaded for sampling. 
}
\usage{
  stanc(file, model_code = '', model_name = "anon_model", verbose = FALSE,
        obfuscate_model_name = TRUE) 
  stanc_builder(file, isystem = dirname(file), 
                verbose = FALSE, obfuscate_model_name = FALSE)
} 

\arguments{
  \item{file}{A character string or a connection that \R supports 
    specifying the Stan model specification in Stan's modeling language.}
  \item{model_code}{Either a character string containing a Stan model 
    specification or the name of a character string object in the workspace. 
    This parameter is used only if parameter \code{file} is not specified, 
    so it defaults to the empty string.} 
  \item{model_name}{A character string naming the model. The 
    default is \code{"anon_model"}. However, the model name will be derived
    from \code{file} or \code{model_code} (if \code{model_code} is the name of a
    character string object) if \code{model_name} is not specified.} 
  \item{verbose}{Logical, defaulting to \code{FALSE}. If \code{TRUE} 
    more intermediate information is printed during the translation procedure.}
  \item{obfuscate_model_name}{Logical, defaulting to \code{TRUE}, 
    indicating whether to use a randomly-generated character string for the 
    name of the C++ class. This prevents name clashes when compiling multiple 
    models in the same R session.}
  \item{isystem}{A character vector of length one naming a path to look for 
    file paths in \code{file} that are to be included within the Stan program
    named by \code{file}. See the Details section below.}
} 

\details{
  The \code{stanc_builder} function supports the standard C++ convention of
  specifying something like \code{#include "my_includes.txt"} on an entire line
  within the file named by the \code{file} argument. In other words, 
  \code{stanc_builder} would look for \code{"my_includes.txt"} in (or under) the 
  directory named by the \code{isystem} argument and insert its contents verbatim 
  at that position  before calling \code{stanc} on the resulting \code{model_code}. 
  This mechanism reduces the need to copy common chunks of code across Stan 
  programs.
  
  Note that line numbers referred to in parser warnings or errors refer to the
  postprocessed Stan program rather than \code{file}. In the case of a parser
  error, the postprocessed Stan program will be printed after the error message.
  Line numbers referred to in messages while Stan is executing also refer to
  the postprocessed Stan program which can be obtained by calling
  \code{\link{get_stancode}}.
}

\value{
  A list with named entries:
  \enumerate{ 
    \item \code{model_name} Character string for the model name.
    \item \code{model_code} Character string for the model's Stan specification. 
    \item \code{cppcode}    Character string for the model's C++ code. 
    \item \code{status}     Logical indicating success/failure (\code{TRUE/FALSE}) of 
                            translating the Stan code. 
  } 
}

\note{
  Unlike \R, in which variable identifiers may contain dots (e.g. \code{a.1}),
  Stan prohibits dots from occurring in variable identifiers. Furthermore, 
  C++ reserved words and Stan reserved words may not be used for variable names; 
  see the Stan User's Guide for a complete list.  
} 

\references{
  The Stan Development Team 
  \emph{Stan Modeling Language User's Guide and Reference Manual}. 
  \url{http://mc-stan.org/}. 

  The Stan Development Team
  \emph{CmdStan Interface User's Guide}.
  \url{http://mc-stan.org}.
} 

\seealso{
  \code{\link{stan_model}} and \code{\link{stan}}
} 

\examples{
stanmodelcode <- "
data {
  int<lower=0> N;
  real y[N];
} 

parameters {
  real mu;
} 

model {
  mu ~ normal(0, 10);
  y ~ normal(mu, 1); 
} 
"

r <- stanc(model_code = stanmodelcode, model_name = "normal1") 
str(r)
}
