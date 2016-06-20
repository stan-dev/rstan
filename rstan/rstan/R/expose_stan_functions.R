# This file is part of RStan
# Copyright (C) 2012, 2013, 2014, 2015, 2016 Jiqiang Guo and Benjamin Goodrich
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

expose_stan_functions <- function(stanmodel) {
  if(is(stanmodel, "stanfit")) {
    stanmodel <- get_stanmodel(stanmodel)
    stanmodel <- get_cppcode(stanmodel)
  }
  else if(is.list(stanmodel)) {
    stanmodel <- stanmodel$cppcode
  }
  else if(is(stanmodel, "stanmodel")) {
    stanmodel <- get_cppcode(stanmodel)
  }
  else if(is.character(stanmodel)) {
    if(length(stanmodel) == 1) stanmodel <- stanc(file = stanmodel)$cppcode
    else stanmodel <- stanc(model_code = stanmodel)$cppcode
  }
  else stop("'stanmodel' is not a valid object")
  
  if(is.null(stanmodel)) {
    warning("could not obtain C++ code for this 'stanmodel'")
    return(invisible(NULL))
  }
  
  lines <- scan(what = character(), sep = "\n", quiet = TRUE, text = stanmodel)
  end <- grep("^class", lines) - 1L # only care about things before Stan's class declaration
  
  # get rid of Stan's local namespace declaration
  lines <- grep("namespace \\{$", lines[1:end], value = TRUE, invert = TRUE)

  # deal with lack of PKG_CFLAGS (necessary stuff is #included in the last step before compilation)
  lines <- sub("#include <stan/model/model_header.hpp>", "", lines, fixed = TRUE)
  lines <- sub("#include <stan/services/command.hpp>", "", lines, fixed = TRUE)
    
  # get rid of templating and just use double because that is about all R can pass
  lines <- gsub("typename boost::math::tools::promote_args.*type ", "double ", lines)
  lines <- gsub("((std::vector<)+)typename boost::math::tools::promote_args.*(>::type)+", "\\1double", lines)
  lines <- gsub("vector<Eigen::Matrix<.*Eigen::Dynamic,1> >", "vector<vector_d>", lines)
  lines <- gsub("Eigen::Matrix<.*Eigen::Dynamic,1>", "vector_d", lines)
  lines <- gsub("vector<Eigen::Matrix<.*1,Eigen::Dynamic> >", "vector<row_vector_d>", lines)
  lines <- gsub("Eigen::Matrix<.*1,Eigen::Dynamic>", "row_vector_d", lines)
  lines <- gsub("vector<Eigen::Matrix<.*,Eigen::Dynamic> >", "vector<matrix_d>", lines)
  lines <- gsub("Eigen::Matrix<.*,Eigen::Dynamic>", "matrix_d", lines)
  
  # kill foo_lpdf<false> functions because of templating
  templated <- grep("_lp[dm]f<false>", lines)
  if(length(templated) > 0) for(i in rev(templated)) {
    end <- i + 1L
    while(!grepl("^}$", lines[end])) end <- end + 1L
    start <- i - 1L
    while(!grepl("^template", lines[start])) start <- start - 1L
    lines <- lines[-c(start:end)]
  }
  
  # stick using:: inside user-defined functions
  usings <- grep("^using", lines, value = TRUE)
  lines <- grep("^using", lines, value = TRUE, invert = TRUE)
  openings <- grep("std::ostream* pstream__) {", lines, fixed = TRUE)
  if(length(openings) == 0) {
    warning("no user-defined functions found")
    return(invisible(NULL))
  }
  for(i in rev(openings)) {
    # hard-code former arguments that cannot be passed from R
    if(grepl("_rng", lines[i], fixed = TRUE)) {
      lines[i] <- gsub("RNG\\&.*\\{$", "const int& seed = 0) {", lines[i])
      lines <- append(lines, c("static boost::ecuyer1988 base_rng__(seed);", 
                               "std::ostream* pstream__ = &Rcpp::Rcout;", 
                               "(void) pstream__;"), i)
    }
    else if(grepl("_lp", lines[i], fixed = TRUE)) {
      lines[i] <- gsub(", T_lp_accum__\\&.*\\{$", ") {", lines[i])
      lines <- append(lines, c("stan::math::accumulator<double> lp_accum__;",
                               "std::ostream* pstream__ = &Rcpp::Rcout;", 
                               "(void) pstream__;"), i)
    }
    else if(grepl("(std::ostream* pstream__) {", lines[i], fixed = TRUE)) {
      lines[i] <- gsub("(std::ostream* pstream__) {", "() {", lines[i], fixed = TRUE)
      lines <- append(lines, c("std::ostream* pstream__ = &Rcpp::Rcout;", 
                               "(void) pstream__;"), i)
    }
    else {
      lines[i] <- gsub(", std::ostream* pstream__) {", ") {", lines[i], fixed = TRUE)
      lines <- append(lines, c("std::ostream* pstream__ = &Rcpp::Rcout;", 
                               "(void) pstream__;"), i)
    }
    lines <- append(lines, usings, i) # make the usings:: local to the function
  }
  
  # yank unused using statements
  lines <- grep("^using stan::io::", lines, value = TRUE, invert = TRUE)
  lines <- grep("^using stan::model::prob_grad", lines, value = TRUE, invert = TRUE)
  
  # convert inline declarations to Rcpp export declarations
  lines <- gsub("^inline$", "// \\[\\[Rcpp::export\\]\\]", lines)

  # declare attributes for Rcpp for non-functor user-defined Stan functions
  templates <- grep("^template .*$", lines)
  for(i in rev(templates)) {
    if(!grepl("functor__", lines[i - 1L]) && lines[i + 1L] != "// [[Rcpp::export]]")
      lines <- append(lines, "// [[Rcpp::export]]", i - 1L)
  }
  
  # do not export function declarations created by the user
  declarations <- grep("std::ostream* pstream__);", lines, fixed = TRUE)
  if(length(declarations) > 0) for(i in rev(declarations)) { # walk back
    lines[i] <- gsub(", std::ostream* pstream__);", ");", lines[i], fixed = TRUE)
    lines[i] <- gsub("RNG& base_rng__", "const int& seed", lines[i], fixed = TRUE)
    j <- i - 1L
    while(lines[j] != "// [[Rcpp::export]]") j <- j - 1L
    lines <- lines[-j]
  }

  # special cases
  ODE_lines <- grep("integrate_ode(", lines, fixed = TRUE)
  ODE_statements <- grep("integrate_ode(", lines, fixed = TRUE, value = TRUE)

  print_lines <- grep("if (pstream__)", lines, fixed = TRUE)
  print_statements <- grep("if (pstream__)", lines, fixed = TRUE, value = TRUE)
  
  # handle more pstream__ arguments
  lines <- gsub(", pstream__)", ")", lines, fixed = TRUE)
  lines <- gsub("(pstream__)", "()", lines, fixed = TRUE)
  lines <- gsub(", std::ostream* pstream__) const {",
                ", std::ostream* pstream__ = &Rcpp::Rcout) const {",
                lines, fixed = TRUE)
  lines <- gsub("(std::ostream* pstream__)", 
                "(std::ostream* pstream__ = &Rcpp::Rcout)", 
                lines, fixed = TRUE)
  
  # put back pstream__ arguments
  if (length(ODE_lines) > 0) lines[ODE_lines] <- ODE_statements
  if (length(print_lines) > 0) lines[print_lines] <- print_statements
  
  lines <- gsub("typename boost::math::tools::promote_args.*(>::type)+", "double", lines)
  
  # remove more base_rng__ arguments
  lines <- gsub(", RNG& base_rng__", "", 
                lines, fixed = TRUE)
  lines <- gsub("(RNG& base_rng__,", "(",
                lines, fixed = TRUE)
  lines <- gsub("(RNG& base_rng__", "(",
                lines, fixed = TRUE)
  lines <- gsub("([[:space:]]+return .*_rng.*), base_rng__\\);",
                "\\1);", lines)
  lines <- gsub("([[:space:]]+return .*_rng)\\(base_rng__\\);",
                "\\1();", lines)
  lines <- gsub("_rng\\(base_rng__\\)", "_rng\\(seed, base_rng__\\)", lines)
  # lines <- gsub(", base_rng__\\)\\);", ", seed\\)\\);", lines)
  # lines <- gsub("stan::math::promote_scalar<fun_return_scalar_t__>\\((.*)_rng\\((.*), base_rng__",
  #               "stan::math::promote_scalar<fun_return_scalar_t__>\\(\\1_rng\\(\\2, seed", 
  #               lines)

  
  # remove line numbering things
  lines <- grep("current_statement_begin__", lines, 
                fixed = TRUE, value = TRUE, invert = TRUE)
                  
  # replace more templating with doubles
  lines <- gsub("const T[0-9]+__&", "const double&", lines)
  lines <- gsub("T_lp__& lp__", "double lp__ = 0.0", lines)
  lines <- gsub("^typename.*$", "double", lines)
  lines <- grep("^[[:space:]]*template", lines, invert = TRUE, value = TRUE)
  lines <- gsub("<T[0-9]+__>", "<double>", lines)
  
  # deal with accumulators
  lines <- gsub(", T_lp_accum__& lp_accum__", "", lines, fixed = TRUE)
  lines <- gsub(", lp_accum__", "", lines, fixed = TRUE)
  lines <- gsub("get_lp(lp__)", "get_lp(lp__, lp_accum__)", lines, fixed = TRUE)
  
  # make propto__ false to not skip anything that is double
  lines <- gsub("const static bool propto__ = true;",
                "const static bool propto__ = false;", lines, fixed = TRUE)
  
  # avoid catch messages that say to report a bug
  lines <- gsub('"*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"', "e.what()",
                lines, fixed = TRUE)

  # restore Stan's Eigen typedefs that were clobbered by the previous lines
  lines <- sub("typedef vector_d vector_d;", "using stan::math::vector_d;", lines)
  lines <- sub("typedef row_vector_d row_vector_d;", "using stan::math::row_vector_d;", lines)
  lines <- sub("typedef matrix_d matrix_d;", "using stan::math::matrix_d;", lines)
  
  # add dependencies
  extras <- dir(rstan_options("boost_lib2"), pattern = "hpp$", 
                full.names = TRUE, recursive = TRUE)
  has_model <- any(grepl("stan::model", lines, fixed = TRUE))
  lines <- c("// [[Rcpp::depends(rstan)]]",
             "#include <Rcpp.h>",
             "#include <RcppEigen.h>",
             if (length(extras) > 0) sapply(extras, FUN = function(x)
               paste0("#include<", x, ">")),             
             "#include <stan/math.hpp>",
             "#include <src/stan/lang/rethrow_located.hpp>",
             if (has_model) "#include <src/stan/model/indexing.hpp>",
             "#include <boost/exception/all.hpp>",
             "#include <boost/random/linear_congruential.hpp>",
             
             "#include <cmath>",
             "#include <cstddef>",
             "#include <fstream>",
             "#include <iostream>",
             "#include <sstream>",
             "#include <stdexcept>",
             "#include <utility>",
             "#include <vector>",

             "#include <boost/date_time/posix_time/posix_time_types.hpp>",
             "#include <boost/math/special_functions/fpclassify.hpp>",
             "#include <boost/random/additive_combine.hpp>",
             "#include <boost/random/uniform_real_distribution.hpp>",
             
             lines)
  
  # try to compile
  on.exit(message("Here is the C++ code that does not compile. Please report bug."))
  on.exit(print(lines), add = TRUE)
  compiled <- Rcpp::sourceCpp(code = paste(lines, collapse = "\n"))
  on.exit(NULL)
  return(invisible(compiled$functions))
}
