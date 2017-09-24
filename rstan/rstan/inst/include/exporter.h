#include <RcppCommon.h>
#include <boost/random/additive_combine.hpp>
#include <iostream>

namespace Rcpp {
  SEXP wrap(boost::ecuyer1988 RNG);
  SEXP wrap(boost::ecuyer1988& RNG);
  SEXP wrap(std::ostream stream);
  template <> boost::ecuyer1988 as(SEXP ptr_RNG);
  template <> boost::ecuyer1988& as(SEXP ptr_RNG);
  template <> std::ostream* as(SEXP ptr_stream);
  namespace traits {
    template <> struct input_parameter<boost::ecuyer1988&>;
  }
}


#include <Rcpp.h>

namespace Rcpp {
  SEXP wrap(boost::ecuyer1988 RNG){
    boost::ecuyer1988* ptr_RNG = &RNG;
    Rcpp::XPtr<boost::ecuyer1988> Xptr_RNG(ptr_RNG);
    return Xptr_RNG;
  }

  SEXP wrap(boost::ecuyer1988& RNG){
    boost::ecuyer1988* ptr_RNG = &RNG;
    Rcpp::XPtr<boost::ecuyer1988> Xptr_RNG(ptr_RNG);
    return Xptr_RNG;
  }

  SEXP wrap(std::ostream stream) {
    std::ostream* ptr_stream = &stream;
    Rcpp::XPtr<std::ostream> Xptr_stream(ptr_stream);
    return Xptr_stream;
  }

  template<> boost::ecuyer1988 as(SEXP ptr_RNG) {
    Rcpp::XPtr<boost::ecuyer1988> ptr(ptr_RNG);
    boost::ecuyer1988& RNG = *ptr; 
 		return RNG;
  }

  template<> boost::ecuyer1988& as(SEXP ptr_RNG) {
    Rcpp::XPtr<boost::ecuyer1988> ptr(ptr_RNG);
    boost::ecuyer1988& RNG = *ptr; 
 		return RNG;
  }

  template<> std::ostream* as(SEXP ptr_stream) {
    Rcpp::XPtr<std::ostream> ptr(ptr_stream);
    return ptr;
  }

  namespace traits {
    template <>
    struct input_parameter<boost::ecuyer1988&> {
      typedef typename Rcpp::ReferenceInputParameter<boost::ecuyer1988&> type ;
    };
  }

}

