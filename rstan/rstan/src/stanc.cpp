// #include <Rcpp.h>
// #include <string>
// #include <iostream>
#include <stan/version.hpp>
#include <stan/gm/compiler.hpp>

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

#include <Rcpp.h>
#include <rstan/io/r_ostream.hpp> 

RcppExport SEXP CPP_stanc(SEXP model_stancode, SEXP model_name);
RcppExport SEXP CPP_stan_version(); 

SEXP CPP_stan_version() {
  BEGIN_RCPP
  std::string stan_version 
    = stan::MAJOR_VERSION + "." +
      stan::MINOR_VERSION + "." +
      stan::PATCH_VERSION;
  SEXP __sexp_result;
  PROTECT(__sexp_result = Rcpp::wrap(stan_version));
  UNPROTECT(1);
  return __sexp_result;
  END_RCPP
} 

SEXP CPP_stanc(SEXP model_stancode, SEXP model_name) { 
  BEGIN_RCPP;
  static const int SUCCESS_RC = 0;
  static const int EXCEPTION_RC = -1;
  static const int PARSE_FAIL_RC = -2;
  
  static const bool INCLUDE_MAIN = true; 
  /*
  std::string stan_version 
    = stan::MAJOR_VERSION + "." +
      stan::MINOR_VERSION + "." +
      stan::PATCH_VERSION;
  */

  std::string mcode_ = Rcpp::as<std::string>(model_stancode); 
  std::string mname_ = Rcpp::as<std::string>(model_name); 
   
  std::stringstream out;
  std::istringstream in(mcode_); 
  try {
    bool valid_model
      = stan::gm::compile(&rstan::io::rcerr,in,out,mname_,!INCLUDE_MAIN);
    if (!valid_model) {
      return Rcpp::List::create(Rcpp::Named("status") = PARSE_FAIL_RC); 

    }
  } catch(const std::exception& e) {
    // REprintf("\nERROR PARSING\n %s\n", e.what()); 
    return Rcpp::List::create(Rcpp::Named("status") = EXCEPTION_RC,
                              Rcpp::Named("msg") = Rcpp::wrap(e.what())); 
  }

  Rcpp::List lst = 
    Rcpp::List::create(Rcpp::Named("status") = SUCCESS_RC, 
                       Rcpp::Named("model_cppname") = mname_,
                       Rcpp::Named("cppcode") = out.str());
  SEXP __sexp_result;
  PROTECT(__sexp_result = lst);
  UNPROTECT(1);
  return __sexp_result;
  END_RCPP
}
