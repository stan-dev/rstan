# This file is part of RStan
# Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017 Trustees of Columbia University
#
# RStan is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# RStan is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

##
## Define rstan plugin for inline package.
## (original name: inline.R)
##

inc_path_fun <- function(package) {
  system.file('include', package = package)
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

PKG_CPPFLAGS_env_fun <- function() {
   Eigen <- dir(system.file("include", "stan", "math", "prim", 
                            package = "StanHeaders", mustWork = TRUE), 
                pattern = "Eigen.hpp$", full.names = TRUE, recursive = TRUE)[1]
   paste(' -I"', file.path(inc_path_fun("Rcpp"), '" '),
         ' -I"', file.path(eigen_path_fun(), '" '),
         ' -I"', file.path(eigen_path_fun(), 'unsupported" '),
         # ' -I"', boost_path_fun2(), '"', # boost_not_in_BH should come
         ' -I"', boost_path_fun(), '"',  # before BH/include
         ' -I"', file.path(inc_path_fun("StanHeaders"), "src", '" '),
         ' -I"', file.path(inc_path_fun("StanHeaders"), '" '),
         ' -I"', file.path(inc_path_fun("RcppParallel"), '" '),
         ' -I"', inc_path_fun("rstan"), '"',
         ' -DEIGEN_NO_DEBUG ',
         ' -DBOOST_DISABLE_ASSERTS ',
         ' -DBOOST_PENDING_INTEGER_LOG2_HPP ',
         ' -DSTAN_THREADS ',
         ' -DBOOST_NO_AUTO_PTR ',
         ' -include ', shQuote(Eigen), ' ',
         ifelse (.Platform$OS.type == "windows", ' -std=c++1y',
                 ' -D_REENTRANT -DRCPP_PARALLEL_USE_TBB=1 '),
         sep = '')
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

  if (.Platform$OS.type == "windows") {
    StanHeaders_pkg_libs <- system.file("libs", .Platform$r_arch,
                                        package = "StanHeaders", mustWork = TRUE)
    RcppParallel_pkg_libs <- system.file("lib", .Platform$r_arch,
                                         package = "RcppParallel", mustWork = TRUE)
    rstan_StanServices <- system.file("lib", .Platform$r_arch, "libStanServices.a",
                                      package = "rstan", mustWork = TRUE)
  }
  else {
    StanHeaders_pkg_libs <- system.file("lib", .Platform$r_arch,
                                        package = "StanHeaders", mustWork = TRUE)
    RcppParallel_pkg_libs <- system.file("lib", .Platform$r_arch,
                                         package = "RcppParallel", mustWork = TRUE)
    rstan_StanServices <- system.file("lib", .Platform$r_arch, "libStanServices.a",
                                      package = "rstan", mustWork = TRUE)
  }

  # In case  we have space (typical on windows though not necessarily)
  # in the file path of Rcpp's library.

  # If rcpp_PKG_LIBS contains space without preceding '\\', add `\\';
  # otherwise keept it intact
  if (grepl('[^\\\\]\\s', rcpp_pkg_libs, perl = TRUE))
      rcpp_pkg_libs <- gsub(rcpp_pkg_path, rcpp_pkg_path2, rcpp_pkg_libs, fixed = TRUE)

  # cat("INFO: rcpp_pkg_libs = ", rcpp_pkg_libs, "\n")

  tbb_libs <- "-ltbb"
  if (!is.null(tbbmalloc_proxyDllInfo))
      tbb_libs <- paste(tbb_libs, "-ltbbmalloc -ltbbmalloc_proxy")

  PL <- paste(rcpp_pkg_libs,
              rstan_StanServices,
              paste0("-L", shQuote(StanHeaders_pkg_libs)),
              "-lStanHeaders",
              paste0("-L", shQuote(RcppParallel_pkg_libs)),
              tbb_libs)

  list(includes = '// [[Rcpp::plugins(cpp14)]]\n',
       body = function(x) x,
       env = list(PKG_LIBS = PL,
                  LOCAL_LIBS = if (.Platform$OS.type == "windows") PL,
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

