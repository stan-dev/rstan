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
stanc_ctx <- V8::v8()

tbbmalloc_proxyDllInfo <- NULL

.onLoad <- function(libname, pkgname) {
  assignInMyNamespace("stanc_ctx",  value= V8::v8())
  stanc_js <- system.file("stanc.js", package = "rstan")
  if (!file.exists(stanc_js)) {
    warning(paste0("Default stancjs compiler not found, ",
                   "downloading the current version from github."))
    stanc_js <- "https://github.com/stan-dev/stanc3/releases/download/v2.27.0/stanc.js"
  }
  rstan:::stanc_ctx$source(stanc_js)
  assignInMyNamespace("rstan_load_time", value = Sys.time())
  set_rstan_ggplot_defaults()
  assignInMyNamespace("RNG", value = get_rng(0))
  assignInMyNamespace("OUT", value = get_stream())
  Rcpp::loadModule("class_model_base", what = TRUE)
  Rcpp::loadModule("class_stan_fit", what = TRUE)
  ## the tbbmalloc_proxy is not loaded by RcppParallel which is linked
  ## in by default on macOS; unloading only works under R >= 4.0 so that
  ## this is only done for R >= 4.0
  if(FALSE && R.version$major >= 4 && Sys.info()["sysname"] == "Darwin") {
      tbbmalloc_proxy  <- system.file("lib/libtbbmalloc_proxy.dylib", package="RcppParallel", mustWork=FALSE)
      tbbmalloc_proxyDllInfo <<- dyn.load(tbbmalloc_proxy, local = FALSE, now = TRUE)
  }
}

.onAttach <- function(...) {
  packageStartupMessage("\nrstan version ",
                        utils::packageVersion("rstan"),
                        " (Stan version ",
                        stan_version(), ")\n")
  packageStartupMessage("For execution on a local, multicore CPU with excess RAM we recommend calling\n",
                        "options(mc.cores = parallel::detectCores()).\n",
                        "To avoid recompilation of unchanged Stan programs, we recommend calling\n",
                        "rstan_options(auto_write = TRUE)",
                        "\nFor within-chain threading using `reduce_sum()` or `map_rect()` Stan functions,\n",
                        "change `threads_per_chain` option:\n",
                        paste0("rstan_options(threads_per_chain = ",
                               rstan_options("threads_per_chain"), ")\n"))
  if (.Platform$OS.type == "windows") {
    packageStartupMessage("Do not specify '-march=native' in 'LOCAL_CPPFLAGS' or a Makevars file")
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
