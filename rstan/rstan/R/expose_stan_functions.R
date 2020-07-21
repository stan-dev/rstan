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
  code <- paste("// [[Rcpp::depends(StanHeaders)]]",
                "// [[Rcpp::depends(rstan)]]",
                "// [[Rcpp::depends(RcppEigen)]]",
                "// [[Rcpp::depends(BH)]]",
                "#include <stan/math/prim/mat/fun/Eigen.hpp>", 
                "#include <boost/integer/integer_log2.hpp>",
                "#include <exporter.h>",
                "#include <RcppEigen.h>",
                code, sep = "\n")
  code <- gsub("// [[stan::function]]", 
               "// [[Rcpp::export]]", code, fixed = TRUE)
  code <- gsub("stan::math::accumulator<double>& lp_accum__, std::ostream* pstream__ = nullptr){", 
               "std::ostream* pstream__ = nullptr){\nstan::math::accumulator<double> lp_accum__;", 
               code, fixed = TRUE)
  code <- gsub("= nullptr", "= 0", code, fixed = TRUE)
  if(is.null(includes)) return(code)
  code <- sub("\n\nstan::io::program_reader prog_reader__() {",
              paste0("\n", includes, "\nstan::io::program_reader prog_reader__() {"), 
              code, fixed = TRUE)
  return(code)
}

expose_stan_functions <- function(stanmodel, includes = NULL, 
                                  show_compiler_warnings = FALSE, ...) {
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
  stopifnot(stanc(model_code = mc, model_name = "User-defined functions",
                  allow_undefined = TRUE)$status)
  r <- .Call("stanfuncs", mc, md5, allow_undefined = TRUE)
  code <- expose_stan_functions_hacks(r$cppcode, includes)

  WINDOWS <- .Platform$OS.type == "windows"
  R_version <- with(R.version, paste(major, minor, sep = "."))
  if (WINDOWS && R_version < "3.7.0") {
    has_USE_CXX11 <- Sys.getenv("USE_CXX11") != ""
    Sys.setenv(USE_CXX11 = 1) # -std=c++1y gets added anyways
    if (!has_USE_CXX11) on.exit(Sys.unsetenv("USE_CXX11"))
  } else {
    has_USE_CXX14 <- Sys.getenv("USE_CXX14") != ""
    Sys.setenv(USE_CXX14 = 1)
    if (!has_USE_CXX14) on.exit(Sys.unsetenv("USE_CXX14"))
  }

  if (rstan_options("required"))
    pkgbuild::has_build_tools(debug = FALSE) || 
    pkgbuild::has_build_tools(debug = TRUE)

  old_LOCAL_LIBS <- Sys.getenv("LOCAL_LIBS")
  if (WINDOWS) {
    TBB <- system.file("lib", .Platform$r_arch, package = "RcppParallel", mustWork = TRUE)
    SH  <- system.file("libs", .Platform$r_arch, package = "StanHeaders",  mustWork = TRUE)
    Sys.setenv(LOCAL_LIBS = paste0("-L", shQuote(TBB), " -tbb -tbbmalloc ", 
                                   "-L", shQuote(SH) , " -lStanHeaders "))
    on.exit(Sys.setenv(LOCAL_LIBS = old_LOCAL_LIBS))
    no_march_flags <- .remove_march_makevars()
  }

  has_LOCAL_CPPFLAGS <- WINDOWS && Sys.getenv("LOCAL_CPPFLAGS") != ""
  if (WINDOWS && !grepl("32", .Platform$r_arch) && !has_LOCAL_CPPFLAGS) {
    Sys.setenv(LOCAL_CPPFLAGS = "-march=core2")
    on.exit(Sys.unsetenv("LOCAL_CPPFLAGS"), add = TRUE)
  }
  
  if (!isTRUE(show_compiler_warnings)) {
    tf <- tempfile(fileext = ".warn")
    zz <- file(tf, open = "wt")
    sink(zz, type = "output")
    on.exit(close(zz), add = TRUE)
    on.exit(sink(type = "output"), add = TRUE)
  }
  compiled <- pkgbuild::with_build_tools(try(suppressWarnings(
    Rcpp::sourceCpp(code = paste(code, collapse = "\n"), ...)), silent = TRUE),
    required = rstan_options("required") &&
    # workaround for packages with src/install.libs.R
      identical(Sys.getenv("WINDOWS"), "TRUE") &&
      !identical(Sys.getenv("R_PACKAGE_SOURCE"), "") )
  if (!isTRUE(show_compiler_warnings)) {
    sink(type = "output")
    close(zz)
    try(file.remove(tf), silent = TRUE)
    on.exit(NULL)
    if (WINDOWS && R_version < "3.7.0") {
      if (!has_USE_CXX11) on.exit(Sys.unsetenv("USE_CXX11"), add = TRUE)
    } else {
      if (!has_USE_CXX14) on.exit(Sys.unsetenv("USE_CXX14"), add = TRUE)
    }
  }
  DOTS <- list(...)
  if (isTRUE(DOTS$dryRun)) return(code)
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
