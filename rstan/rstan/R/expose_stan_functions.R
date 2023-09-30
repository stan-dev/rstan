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
                "// [[Rcpp::plugins(rstan)]]",
                "// [[Rcpp::depends(RcppEigen)]]",
                "// [[Rcpp::depends(BH)]]",
                "#include <stan/math/prim/fun/Eigen.hpp>",
                "#include <stan/math/prim/meta.hpp>",
                "#include <boost/integer/integer_log2.hpp>",
                "#include <exporter.h>",
                "#include <RcppEigen.h>",
                includes,
                code, sep = "\n")
  code <- gsub("// [[stan::function]]",
               "// [[Rcpp::export]]", code, fixed = TRUE)
  code <- gsub("stan::math::accumulator<double>& lp_accum__,(\\n)?(\\s*)?std::ostream\\* pstream__ = (nullptr|0))(\\s*)?\\{",
               "std::ostream* pstream__ = nullptr){\nstan::math::accumulator<double> lp_accum__;",
               code)
  code <- gsub("pstream__(\\s*|)=(\\s*|)nullptr", "pstream__ = 0", code)
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
  r <- stanc(model_code = mc, model_name = "User-defined functions",
             allow_undefined = TRUE,
             standalone_functions = TRUE)
  code <- expose_stan_functions_hacks(r$cppcode, includes)

  WINDOWS <- .Platform$OS.type == "windows"
  if (WINDOWS && R.version$major < 4) {
    stop("expose_stan_functions requires R >= 4.0 on Windows to use C++17")
  } else {
    has_USE_CXX17 <- Sys.getenv("USE_CXX17") != ""
    Sys.setenv(USE_CXX17 = 1)
    if (!has_USE_CXX17) on.exit(Sys.unsetenv("USE_CXX17"))
  }

  if (rstan_options("required"))
    pkgbuild::has_build_tools(debug = FALSE) ||
    pkgbuild::has_build_tools(debug = TRUE)

  if (WINDOWS) {
    has_march = .warn_march_makevars()
    if (has_march) {
      user_makevar = Sys.getenv("R_MAKEVARS_USER")
      Sys.setenv(R_MAKEVARS_USER = NULL)
      on.exit(Sys.setenv(R_MAKEVARS_USER = user_makevar))
    }
  }

  if (!isTRUE(show_compiler_warnings)) {
    tf <- tempfile(fileext = ".warn")
    zz <- file(tf, open = "wt")
    sink(zz, type = "output")
    on.exit(close(zz), add = TRUE)
    on.exit(sink(type = "output"), add = TRUE)
  }
  Rcpp::registerPlugin("rstan", rstanplugin)
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
  }
  DOTS <- list(...)
  if (isTRUE(DOTS$dryRun)) return(code)
  if (inherits(compiled, "try-error")) stop("Compilation failed!")
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
