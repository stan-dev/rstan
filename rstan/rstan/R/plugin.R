##
## Define rstan plugin for inline package.
## (original name: inline.R)
##

rstan_inc_path_fun <- function() { 
  system.file('include', package = 'rstan')
} 

# Using RcppEigen
eigen_path_fun <- function() {
  rstan_options("eigen_lib")
}

boost_path_fun <- function() {
  rstan_options("boost_lib")
}

boost_path_fun2 <- function() {
  rstan_options("boost_lib2")
}

# If included in RStan
# eigen_path_fun() <- paste0(rstan_inc_path_fun(), '/stanlib/eigen_3.1.0')

PKG_CPPFLAGS_env_fun <- function() {
   paste(' -I"', file.path(rstan_inc_path_fun(), '/stansrc" '),
         ' -isystem"', file.path(eigen_path_fun(), '" '),
         ' -isystem"', file.path(eigen_path_fun(), '/unsupported" '),
         ' -isystem"', boost_path_fun2(), '"', # boost_not_in_BH should come         
         ' -isystem"', boost_path_fun(), '"',  # before BH/include
         ' -I"', rstan_inc_path_fun(), '"', 
         ' -Wno-unused-function -Wno-uninitialized',
         ' -DBOOST_RESULT_OF_USE_TR1 -DBOOST_NO_DECLTYPE -DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG', sep = '')
}

legitimate_space_in_path <- function(path) {
  # For windows, use the short path name (8.3 format) 
  # 
  if (.Platform$OS.type == "windows") { 
    path <- normalizePath(path)
    if (grepl(" ", path, fixed = TRUE)) 
      path <- utils::shortPathName(path)
    # it is weird that the '\\' in the path name will be gone
    # when passed to cxxfunction, so change it to '/' 
    path <- gsub('\\\\', '/', path, perl = TRUE)
  }
  path 
} 

rstanplugin <- function() {
  Rcpp_plugin <- getPlugin("Rcpp")
  rcpp_pkg_libs <- Rcpp_plugin$env$PKG_LIBS
  rcpp_pkg_path <- system.file(package = 'Rcpp')
  rcpp_pkg_path2 <- legitimate_space_in_path(rcpp_pkg_path) 
 
  # In case  we have space (typicall on windows though not necessarily)
  # in the file path of Rcpp's library. 
  
  # If rcpp_PKG_LIBS contains space without preceding '\\', add `\\'; 
  # otherwise keept it intact
  if (grepl('[^\\\\]\\s', rcpp_pkg_libs, perl = TRUE))
    rcpp_pkg_libs <- gsub(rcpp_pkg_path, rcpp_pkg_path2, rcpp_pkg_libs, fixed = TRUE) 

  list(includes = '',
       body = function(x) x,
       LinkingTo = c("Rcpp"),
       ## FIXME see if we can use LinkingTo for RcppEigen's header files
       env = list(PKG_LIBS = paste(rcpp_pkg_libs),
                  PKG_CPPFLAGS = paste(Rcpp_plugin$env$PKG_CPPFLAGS,
                                        PKG_CPPFLAGS_env_fun(), collapse = " ")))
}


# inlineCxxPlugin would automatically get registered in inline's plugin list.
# Note that everytime rstan plugin is used, inlineCxxPlugin
# gets called so we can change some settings on the fly
# for example now by setting rstan_options(boost_lib=xxx)
inlineCxxPlugin <- function(...) {
  settings <- rstanplugin()
  settings
}

