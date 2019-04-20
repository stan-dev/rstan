# This file is part of RStan
# Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019 Trustees of Columbia University
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
  lines <- scan(what = character(), sep = "\n", quiet = TRUE, text = code)
  lines <- append(lines, after = 0L,
                  values = c("// [[Rcpp::depends(rstan)]]", 
                             "#include <RcppEigen.h>", 
                             "#include <rstan/boost_random_R.hpp>"))
  lines <- gsub("// [[stan::function]]", 
                "// [[Rcpp::export]]", lines, fixed = TRUE)
  
  pat <- ", std::ostream* pstream__ = nullptr"
  mark <- grep(pat, lines, fixed = TRUE)
  lines[mark] <- gsub(pat, ", Rcpp::Nullable<Rcpp::IntegerVector> pstream__ = R_NilValue", 
                      lines[mark], fixed = TRUE)
  lines[mark + 2L] <- gsub("pstream__);", "&Rcpp::Rcout);", lines[mark + 2L], fixed = TRUE)
  pat <- "(std::ostream* pstream__ = nullptr)"
  mark <- grep(pat, lines, fixed = TRUE)
  if (length(mark) > 0L) {
    lines[mark] <- gsub(pat, "(Rcpp::Nullable<Rcpp::IntegerVector> pstream__ = R_NilValue)", 
                        lines[mark], fixed = TRUE)
    lines[mark + 2L] <- gsub("pstream__);", "&Rcpp::Rcout);", lines[mark + 2L], fixed = TRUE)
  }
  
  pat <- ", boost::ecuyer1988& base_rng__"
  mark <- grep(pat, lines, fixed = TRUE)
  if (length(mark) > 0L) {
    lines[mark] <- gsub(pat, ", Rcpp::Nullable<Rcpp::IntegerVector> base_rng__ = R_NilValue", 
                        lines[mark], fixed = TRUE)
    mark <- mark + 2L
    lines[mark] <- gsub("boost::ecuyer1988>(", "boost_random_R>(", lines[mark], fixed = TRUE)
    lines[mark] <- gsub("base_rng__", "brR", lines[mark], fixed = TRUE)
    mark <- mark - 2L
    for (j in rev(mark))
      lines <- append(lines, values = "boost_random_R brR;", after = j)
  }
  pat <- "(boost::ecuyer1988& base_rng__"
  mark <- grep(pat, lines, fixed = TRUE)
  if (length(mark) > 0L) {
    lines[mark] <- gsub(pat, "(Rcpp::Nullable<Rcpp::IntegerVector> base_rng__ = R_NilValue", 
                        lines[mark], fixed = TRUE)
    mark <- mark + 2L
    lines[mark] <- gsub("boost::ecuyer1988>(", "boost_random_R>(", lines[mark], fixed = TRUE)
    lines[mark] <- gsub("base_rng__", "brR", lines[mark], fixed = TRUE)
    mark <- mark - 2L
    for (j in rev(mark))
      lines <- append(lines, values = "boost_random_R brR;", after = j)
  }
  
  pat <- "double& lp__, stan::math::accumulator<double>& lp_accum__"
  mark <- grep(pat, lines, fixed = TRUE)
  if (length(mark) > 0L) {
    lines[mark] <- gsub(pat, "double lp__ = 0.0", lines[mark], fixed = TRUE)
    for (j in rev(mark))
      lines <- append(lines, "stan::math::accumulator<double> lp_accum__;", after = j)
  }
  # if there is lp_accum__, there is always and lp__ before it
  
  if(!is.null(includes)) {
    mark <- grep("stan::io::program_reader prog_reader__()", lines, fixed = TRUE)
    lines <- append(lines, values = includes, after = mark - 1L)
  }
  code <- paste(lines, collapse = "\n")
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
    pkgbuild::has_build_tools(debug = FALSE) || pkgbuild::has_build_tools(debug = TRUE)
  
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
  compiled <- pkgbuild::with_build_tools(suppressWarnings(
    Rcpp::sourceCpp(code = paste(code, collapse = "\n"), ...)),
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
  return(invisible(compiled$functions))
}
