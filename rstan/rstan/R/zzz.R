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

rstan_load_time <- as.POSIXct("1970-01-01 00:00.00 UTC")
RNG <- 0
OUT <- 0

tbbmalloc_proxyDllInfo <- NULL

.onLoad <- function(libname, pkgname) {
  assignInMyNamespace("rstan_load_time", value = Sys.time())
  set_rstan_ggplot_defaults()
  assignInMyNamespace("RNG", value = get_rng(0))
  assignInMyNamespace("OUT", value = get_stream())
  Rcpp::loadModule("class_model_base", what = TRUE)
  Rcpp::loadModule("class_stan_fit", what = TRUE)
  ## the tbbmalloc_proxy is not loaded by RcppParallel which is linked
  ## in by default on macOS
  if(Sys.info()["sysname"] == "Darwin") {
      tbbmalloc_proxy  <- system.file("lib/libtbbmalloc_proxy.dylib", package="RcppParallel", mustWork=FALSE)
      tbbmalloc_proxyDllInfo <<- dyn.load(tbbmalloc_proxy, local = FALSE, now = TRUE)
  }
}

.onAttach <- function(...) {
  rstanLib <- dirname(system.file(package = "rstan"))
  pkgdesc <- packageDescription("rstan", lib.loc = rstanLib)
  gitrev <- substring(git_head(), 0, 12)
  packageStartupMessage(paste("rstan (Version ", pkgdesc$Version, ", GitRev: ", gitrev, ")", sep = ""))
  packageStartupMessage("For execution on a local, multicore CPU with excess RAM we recommend calling\n",
                        "options(mc.cores = parallel::detectCores()).\n",
                        "To avoid recompilation of unchanged Stan programs, we recommend calling\n",
                        "rstan_options(auto_write = TRUE)")
  if (.Platform$OS.type == "windows") {
    R_version <- as.integer(R.version$major)
    processor_msg <- ifelse(R_version >= 4, 
                            "Sys.setenv(LOCAL_CPPFLAGS = '-march=native -mtune=native')", 
                            "Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')")
    packageStartupMessage("For improved execution time, we recommend calling\n", processor_msg)
  }
}

.onUnload <- function(libpath) {
   # unload the package library
   library.dynam.unload("rstan", libpath)

   # unload tbbmalloc_proxy if we loaded it
   if (!is.null(tbbmalloc_proxyDllInfo)) {
       dyn.unload(tbbmalloc_proxyDllInfo[["path"]])
       tbbmalloc_proxyDllInfo <<- NULL
   }
}
