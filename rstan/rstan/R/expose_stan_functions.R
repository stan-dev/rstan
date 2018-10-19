# This file is part of RStan
# Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018 Trustees of Columbia University
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

expose_stan_functions_hacks <- function(code, includes = NULL) {
  code <- paste("#include <exporter.h>\n#include <RcppEigen.h>", code, sep="\n")
  code <- gsub("// [[stan::function]]", 
               "// [[Rcpp::depends(rstan)]]\n// [[Rcpp::export]]", code, fixed = TRUE)
  code <- gsub("stan::math::accumulator<double>& lp_accum__, std::ostream* pstream__ = nullptr){", 
               "std::ostream* pstream__ = nullptr){\nstan::math::accumulator<double> lp_accum__;", 
               code, fixed = TRUE)
  if(is.null(includes)) return(code)
  code <- sub("\n\nstan::io::program_reader prog_reader__() {",
              paste0("\n", includes, "\nstan::io::program_reader prog_reader__() {"), 
              code, fixed = TRUE)
  return(code)
}

expose_stan_functions <- function(stanmodel, includes = NULL, ...) {
  mc <- NULL
  if(is(stanmodel, "stanfit")) {
    mc <- get_stancode(get_stanmodel(stanmodel))
  }
  else if(is.list(stanmodel)) {
    mc <- stanmodel$model_code
  }
  else if(is(stanmodel, "stanmodel")) {
    mc <- get_stancode(stanmodel)
  }
  else if(is.character(stanmodel)) {
    if(length(stanmodel) == 1) mc <- get_model_strcode(stanmodel, NULL)
    else mc <- get_model_strcode(model_code = stanmodel)
  }
  else stop("'stanmodel' is not a valid object")
  
  tf <- tempfile(fileext = ".stan")
  writeLines(mc, con = tf)
  md5 <- paste("user", tools::md5sum(tf), sep = "_")
  stopifnot(stanc(model_code = mc, model_name = "User-defined functions")$status)
  r <- .Call("stanfuncs", mc, md5, allow_undefined = TRUE)
  code <- expose_stan_functions_hacks(r$cppcode, includes)

  R_version <- with(R.version, paste(major, minor, sep = "."))
  if (.Platform$OS.type == "windows" && R_version < "3.6.0") {
    has_USE_CXX11 <- Sys.getenv("USE_CXX11") != ""
    Sys.setenv(USE_CXX11 = 1) # -std=c++1y gets added anyways
    if (!has_USE_CXX11) on.exit(Sys.unsetenv("USE_CXX11"))
  } else {
    has_USE_CXX14 <- Sys.getenv("USE_CXX14") != ""
    Sys.setenv(USE_CXX14 = 1)
    if (!has_USE_CXX14) on.exit(Sys.unsetenv("USE_CXX14"))
  }
  
  compiled <- suppressWarnings(pkgbuild::with_build_tools(
    Rcpp::sourceCpp(code = paste(code, collapse = "\n"), ...)) )
  DOTS <- list(...)
  ENV <- DOTS$env
  if (is.null(ENV)) ENV <- globalenv()
  for (x in compiled$functions) {
    FUN <- get(x, envir = ENV)
    args <- formals(FUN)
    args$pstream__ <- get_stream()
    if ("lp__" %in% names(args)) args$lp__ <- 0
    if ("base_rng__" %in% names(args)) args$base_rng__ <- get_rng()
    formals(FUN) <- args
    assign(x, FUN, envir = ENV)
  }
  return(invisible(compiled$functions))
}
