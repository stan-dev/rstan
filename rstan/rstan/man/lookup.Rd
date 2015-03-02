\name{lookup}
\alias{lookup} 
\title{
Look up the Stan function that corresponds to a \R function.
}
\description{
  This function helps to map between R functions and Stan functions.
}
\usage{
  lookup(FUN, ReturnType = character()) 
} 

\arguments{
  \item{FUN}{A character string naming a \R function or a \R function for
    which the (near) equivalent Stan function is sought. If \code{"~"},
    then all Stan sampling statements are returned.}
  \item{ReturnType}{A character string of positive length naming a valid
    return type for a Stan function: \code{int}, \code{int[]}, \code{matrix},
    \code{real}, \code{real[,]}, \code{real[]}, \code{row_vector}, 
    \code{T[]}, \code{vector}, or \code{void}. If \code{"ANY"} is passed,
    then the entire \code{\link{data.frame}} is returned and can be 
    inspected with the \code{\link{View}} function, for example.}
}

\value{
  Ordinarily, a data.frame with rows equal to the number of partial matches and four columns:
  \enumerate{
    \item \code{StanFunction} Character string for the Stan function's name.
    \item \code{Arguments} Character string indicating the arguments to that Stan function. 
    \item \code{ReturnType} Character string indicating the return type of that Stan function.
    \item \code{Page} Integer indicating the page of the Stan reference manual where  
                      that Stan function is defined. 
  }
  If \code{FUN = "~"}, a data.frame with rows equal to the number of sampling statements
  defined in Stan and four columns:
  \enumerate{
    \item \code{FirstArgument} Character string for the first argument for that sampling statement.
    \item \code{StanStatement} Character string indicating the name of the distribution in Stan. 
    \item \code{Arguments} Character string indicating the remaining arguments to that Stan function. 
    \item \code{Page} Integer indicating the page of the Stan reference manual where  
                      that Stan sampling statement is defined. 
  }
  If there are no matching Stan functions, a character string indicating so is returned.
}

\references{
  The Stan Development Team 
  \emph{Stan Modeling Language User's Guide and Reference Manual}. 
  \url{http://mc-stan.org/}. 

  The Stan Development Team
  \emph{CmdStan Interface User's Guide}.
  \url{http://mc-stan.org}.
} 

\examples{
lookup(dnorm)
nrow(lookup("~"))
lookup("foo")
}
