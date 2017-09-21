#include <boost/random/additive_combine.hpp>
#include <RcppCommon.h>

namespace Rcpp {
  SEXP wrap(boost::ecuyer1988 RNG);
  template <> boost::ecuyer1988& as(SEXP ptr_RNG);
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

  template<> boost::ecuyer1988 as(SEXP ptr_RNG) {
    Rcpp::XPtr<boost::ecuyer1988> ptr(ptr_RNG);
    boost::ecuyer1988& RNG = *ptr; 
 		return RNG;
  }

  namespace traits {
    template <>
    struct input_parameter<boost::ecuyer1988&> {
      typedef typename Rcpp::ReferenceInputParameter<boost::ecuyer1988> type ;
    };
  }

}



