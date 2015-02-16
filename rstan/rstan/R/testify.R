testify <- function(stanmodel) {
  if(is(stanmodel, "stanfit")) stanmodel <- get_stanmodel(stanmodel)
  else stopifnot(is(stanmodel, "stanmodel"))
  stanmodel <- get_cppcode(stanmodel)
  if(is.null(stanmodel)) {
    warning("could not obtain C++ code for this 'stanmodel'")
    return(invisible(NULL))
  }
  lines <- scan(what = character(), sep = "\n", quiet = TRUE, text = stanmodel)
  end <- grep("^class", lines) - 1L # only care about things before Stan's class declaration
  
  # get rid of local namespace
  lines <- grep("namespace \\{$", lines[1:end], value = TRUE, invert = TRUE)
  
  # get rid of templating and just use double because that is about all R can pass
  lines <- gsub("Eigen::Matrix<.*Eigen::Dynamic,1>", "vector_d", lines)
  lines <- gsub("Eigen::Matrix<.*1,Eigen::Dynamic>", "row_vector_d", lines)
  lines <- gsub("Eigen::Matrix<.*>", "matrix_d", lines)
  lines <- gsub("typename boost::math::tools::promote_args.*type ", "double ", lines)
  
  # restore Stan's Eigen typedefs that were clobbered by the previous lines
  lines <- gsub("typedef vector_d vector_d;", 
                "typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;", lines)
  lines <- gsub("typedef row_vector_d row_vector_d;", 
                "typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;", lines)
  lines <- gsub("typedef matrix_d matrix_d;",
                "typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;", lines)
                                  
  # stick using:: inside user-defined functions
  usings <- grep("^using", lines, value = TRUE)
  lines <- grep("^using", lines, value = TRUE, invert = TRUE)
  openings <- grep("std::ostream* pstream__) {", lines, fixed = TRUE)
  if(length(openings) == 0) {
    warning("no user-defined functions found")
    return(invisible(NULL))
  }
  for(i in rev(openings)) {
    lines[i] <- gsub(", std::ostream* pstream__) {", 
                     ") { std::ostream* pstream__ = &Rcpp::Rcout;", 
                     lines[i], fixed = TRUE)
    lines <- append(lines, usings, i)
  }
  
  # get rid of inline declarations
  lines <- grep("^inline", lines, value = TRUE, invert = TRUE)
  
  # declare attributes for Rcpp
  lines <- gsub("^template <typename.*$", "// [[Rcpp::export]]", lines)
  
  # do not export function declarations
  declarations <- grep("std::ostream* pstream__);", lines, fixed = TRUE)
  if(length(declarations) > 0) for(i in rev(declarations)) { # walk back
    lines[i] <- gsub(", std::ostream* pstream__);", ");", lines[i], fixed = TRUE)
    j <- i - 1L
    while(lines[j] != "// [[Rcpp::export]]") j <- j - 1L
    lines <- lines[-j]
  }
  
  # do not export functors
  functors <- grep("functor__ \\{$", lines)
  lines <- lines[-c(functors + 1L)]

  # remove more pstream__ arguments
  lines <- gsub(", std::ostream* pstream__) const {", ") const {", lines, fixed = TRUE)
  lines <- gsub(", pstream__)", ")", lines, fixed = TRUE)
  
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
