.onAttach <- function(libname, pkgname) {
  pkg_cppflags <- Sys.getenv("PKG_CPPFLAGS")
  plugin <- system.file("include", "stan", "math", "prim", "mat", "fun", "Eigen.hpp",
                        package = "StanHeaders", mustWork = TRUE)
  Sys.setenv(PKG_CPPFLAGS = paste(pkg_cppflags, "-include", plugin))
  return(invisible(NULL))  
}

.onDetach <- function(libpath) {
  pkg_cppflags <- Sys.getenv("PKG_CPPFLAGS")
  plugin <- system.file("include", "stan", "math", "prim", "mat", "fun", "Eigen.hpp",
                        package = "StanHeaders", mustWork = TRUE)
  pkg_cppflags <- sub(paste("-include", plugin), "", pkg_cppflags, fixed = TRUE)
  Sys.setenv(PKG_CPPFLAGS = pkg_cppflags)
  return(invisible(NULL))
}
