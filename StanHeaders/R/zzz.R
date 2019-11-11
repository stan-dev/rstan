.onAttach <- function(libname, pkgname) {
  old_pkg_cxxflags <- Sys.getenv("PKG_CXXFLAGS")
  plugin <- system.file("include", "stan", "math", "prim", "mat", "fun", "Eigen.hpp",
                        package = "StanHeaders", mustWork = TRUE)
  StanMath <- system.file("include", package = "StanHeaders", mustWork = TRUE)
  RcppEigen <- system.file("include", package = "RcppEigen", mustWork = TRUE)
  new_pkg_cxxflags <- paste(old_pkg_cxxflags, 
                            "-include", plugin,
                            "-I", StanMath,
                            "-I", RcppEigen)
  Sys.setenv(PKG_CXXFLAGS = new_pkg_cxxflags)
  return(invisible(NULL))  
}

.onDetach <- function(libpath) {
  new_pkg_cxxflags <- Sys.getenv("PKG_CXXFLAGS")
  plugin <- system.file("include", "stan", "math", "prim", "mat", "fun", "Eigen.hpp",
                        package = "StanHeaders", mustWork = TRUE)
  StanMath <- system.file("include", package = "StanHeaders", mustWork = TRUE)
  RcppEigen <- system.file("include", package = "RcppEigen", mustWork = TRUE)
  old_pkg_cxxflags <- sub(paste("-include", plugin, 
                                "-I", StanMath,
                                "-I", RcppEigen), "", new_pkg_cxxflags)
  Sys.setenv(PKG_CPPFLAGS = old_pkg_cxxflags)
  return(invisible(NULL))
}
