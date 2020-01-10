\name{CxxFlags}
\alias{CxxFlags}
\alias{LdFlags}
\title{Compilation flags for StanHeaders
}
\description{
Output the compiler or linker flags required to build with the \pkg{StanHeaders}
package
}
\usage{
  CxxFlags(as_character = FALSE)
  LdFlags(as_character = FALSE)
}

\arguments{
  \item{as_character}{
    A logical scalar that defaults to \code{\link{FALSE}} that indicates whether
    to return the compiler or linker flags as a \code{\link{character}} vector of
    length one. Otherwise, the compiler or linker flags are merely output to the
    screen, which is appropriate when called from a Makevars or Makevars.win file
  }
}
\details{
  These functions are currently not exported and are typically called from a Makevars 
  or a Makevars.win file of another package as follows:
  
  PKG_CXXFLAGS += $(shell "${R_HOME}/bin/Rscript" -e "StanHeaders:::CxxFlags()")
  PKG_LDLIBS += $(shell "${R_HOME}/bin/Rscript" -e "StanHeaders:::LdFlags()")
}
\value{
  If \code{as_character} is \code{\link{TRUE}}, then these functions return
  a character vector of length one. Otherwise, (which is the default) these
  functions return \code{\link{NULL}} invisibly after outputing the compiler
  or linker flags to the screen.
}
