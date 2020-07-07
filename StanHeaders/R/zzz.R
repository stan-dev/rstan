.onAttach <- function(libname, pkgname) {
  old_pkg_cxxflags <- Sys.getenv("PKG_CXXFLAGS")
  plugin <- dir(system.file("include", "stan", "math", "prim", 
                            package = "StanHeaders", mustWork = TRUE),
                pattern = "Eigen.hpp$", recursive = TRUE, full.names = TRUE)[1]
  StanMath <- system.file("include", package = "StanHeaders", mustWork = FALSE)
  RcppEigen <- system.file("include", package = "RcppEigen", mustWork = FALSE)
  new_pkg_cxxflags <- paste(old_pkg_cxxflags, 
                            "-include", plugin,
                            "-I", shQuote(StanMath),
                            "-I", shQuote(RcppEigen))
  Sys.setenv(PKG_CXXFLAGS = new_pkg_cxxflags)
  return(invisible(NULL))  
}

.onDetach <- function(libpath) {
  new_pkg_cxxflags <- Sys.getenv("PKG_CXXFLAGS")
  plugin <- dir(system.file("include", "stan", "math", "prim", 
                            package = "StanHeaders", mustWork = TRUE),
                pattern = "Eigen.hpp$", recursive = TRUE, full.names = TRUE)[1]
  StanMath <- system.file("include", package = "StanHeaders", mustWork = FALSE)
  RcppEigen <- system.file("include", package = "RcppEigen", mustWork = FALSE)
  old_pkg_cxxflags <- sub(paste("-include", plugin, 
                                "-I", shQuote(StanMath),
                                "-I", shQuote(RcppEigen)), "", new_pkg_cxxflags)
  Sys.setenv(PKG_CPPFLAGS = old_pkg_cxxflags)
  return(invisible(NULL))
}
