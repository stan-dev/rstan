# This file is part of RStan
# Copyright (C) 2012, 2013, 2014, 2015 Jiqiang Guo and Benjamin Goodrich
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

testify <- function(stanmodel) {
  if(is(stanmodel, "stanfit")) {
    stanmodel <- get_stanmodel(stanmodel)
    stanmodel <- get_cppcode(stanmodel)
  }
  else if(is.list(stanmodel)) {
    stanmodel <- stanmodel$cppcode
  }
  else {
    stopifnot(is(stanmodel, "stanmodel"))
    stanmodel <- get_cppcode(stanmodel)
  }  
  if(is.null(stanmodel)) {
    warning("could not obtain C++ code for this 'stanmodel'")
    return(invisible(NULL))
  }
  
  lines <- scan(what = character(), sep = "\n", quiet = TRUE, text = stanmodel)
  end <- grep("^class", lines) - 1L # only care about things before Stan's class declaration
  
  # get rid of Stan's local namespace declaration
  lines <- grep("namespace \\{$", lines[1:end], value = TRUE, invert = TRUE)
  
  # get rid of templating and just use double because that is about all R can pass
  lines <- gsub("typename boost::math::tools::promote_args.*type ", "double ", lines)  
  lines <- gsub("Eigen::Matrix<.*Eigen::Dynamic,1>", "vector_d", lines)
  lines <- gsub("Eigen::Matrix<.*1,Eigen::Dynamic>", "row_vector_d", lines)
  lines <- gsub("Eigen::Matrix<.*>", "matrix_d", lines)
  
  # restore Stan's Eigen typedefs that were clobbered by the previous lines
  lines <- gsub("typedef vector_d vector_d;", 
                "typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;", lines)
  lines <- gsub("typedef row_vector_d row_vector_d;", 
                "typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;", lines)
  lines <- gsub("typedef matrix_d matrix_d;",
                "typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;", lines)

  # kill foo_log<false> functions because of templating
  templated <- grep("_log<false>", lines, fixed = TRUE)
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
      lines <- append(lines, c("boost::ecuyer1988 base_rng__(seed);", 
                               "std::ostream* pstream__ = &Rcpp::Rcout;"), i)
    }
    else {
      lines[i] <- gsub(", std::ostream* pstream__) {", ") {", lines[i], fixed = TRUE)
      lines <- append(lines, "std::ostream* pstream__ = &Rcpp::Rcout;", i)
    }
    lines <- append(lines, usings, i) # make the usings:: local to the function
  }
  
  # get rid of inline declarations
  lines <- grep("^inline", lines, value = TRUE, invert = TRUE)
  
  # declare attributes for Rcpp for non-functor user-defined Stan functions
  templates <- grep("^template .*$", lines)
  for(i in rev(templates)) if(!grepl("functor__", lines[i - 1L])) {
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
  
  # remove more pstream__ arguments
  lines <- gsub(", std::ostream* pstream__) const {", ") const {", 
                lines, fixed = TRUE)
  lines <- gsub(", pstream__)", ")", lines, fixed = TRUE)
  
  # remove more base_rng__ arguments
  lines <- gsub(", RNG& base_rng__", "", 
                lines, fixed = TRUE)
  lines <- gsub("(RNG& base_rng__", "(",
                lines, fixed = TRUE)
  lines <- gsub("([[:space:]]+return .*_rng.*), base_rng__\\);",
                "\\1);", lines)
  lines <- gsub("([[:space:]]+return .*_rng)\\(base_rng__\\);",
                "\\1();", lines)
                  
  # replace more templating with doubles
  lines <- gsub("const T[0-9]+__&", "double&", lines)
  lines <- gsub("T_lp__&", "double&", lines)
  lines <- gsub("T_lp_accum__&", "double&", lines)
  lines <- gsub("^typename.*$", "double", lines)
  lines <- grep("^template", lines, invert = TRUE, value = TRUE)

  # deal with (lack of) accumulators
  lines <- gsub(", double& lp_accum__", "", lines, fixed = TRUE)
  lines <- gsub("lp_accum__.add", "lp__ += ", lines, fixed = TRUE)
  lines <- gsub(", lp_accum__", "", lines, fixed = TRUE)
  
  # add dependencies
  lines <- c("// [[Rcpp::depends(StanHeaders)]]", 
             "// [[Rcpp::depends(BH)]]",
             "// [[Rcpp::depends(RcppEigen)]]",
             "#include<Rcpp.h>",
             "#include<RcppEigen.h>",
             lines)
  
  # try to compile
  compiled <- sourceCpp(code = paste(lines, collapse = "\n"))
  return(invisible(NULL))
}
