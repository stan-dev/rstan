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

.onLoad <- function(libname, pkgname) {
  if (requireNamespace("V8", quietly = TRUE)) {
    assign("stanc_ctx", V8::v8(), envir = topenv())
  } else assign("stanc_ctx", QuickJSR::JSContext$new(stack_size = 4 * 1024 * 1024), envir = topenv())
  stanc_js <- system.file("stanc.js", package = "StanHeaders", mustWork = TRUE)
  test <- try(stanc_ctx$source(stanc_js), silent = TRUE)
  if (inherits(test, "try-error")) {
    stanc_js <- system.file("exec", "stanc.js", package = "rstan", mustWork = TRUE)
    stanc_ctx$source(stanc_js)
  }
  assignInMyNamespace("rstan_load_time", value = Sys.time())
  set_rstan_ggplot_defaults()
  assignInMyNamespace("RNG", value = get_rng(0))
  assignInMyNamespace("OUT", value = get_stream())
  Rcpp::loadModule("class_model_base", what = TRUE)
  Rcpp::loadModule("class_stan_fit", what = TRUE)
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
}
