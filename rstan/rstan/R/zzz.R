.onLoad <- function(libname, pkgname) { }

.onAttach <- function(...) {
  rstanLib <- dirname(system.file(package = "rstan"))
  pkgdesc <- packageDescription("rstan", lib.loc = rstanLib)
  builddate <- gsub(';.*$', '', pkgdesc$Packaged)
#  gitrev <- substring(git_head(), 0, 12)
  gitrev <- "CRAN"
  packageStartupMessage(paste("rstan (Version ", pkgdesc$Version, ", packaged: ", builddate, ", GitRev: ", gitrev, ")", sep = ""))
} 

