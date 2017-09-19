#include <boost/random/additive_combine.hpp>
#include <RcppCommon.h>

namespace Rcpp {
  template <> SEXP wrap<boost::ecuyer1988>(const boost::ecuyer1988 RNG);
  template <> boost::ecuyer1988 as(SEXP ptr_RNG);

  template <> class InputParameter<boost::ecuyer1988&>;

  namespace traits {
    template <> struct input_parameter<boost::ecuyer1988&>;
  }
}


#include <Rcpp.h>

namespace Rcpp {
  SEXP wrap(const boost::ecuyer1988 RNG){
    Rcpp::XPtr<boost::ecuyer1988, false> Xptr_RNG(RNG);
    return Rcpp::wrap(Xptr_RNG);
  }

  template<> boost::ecuyer1988 as(SEXP ptr_RNG) {
    Rcpp::XPtr<boost::ecuyer1988> ptr(ptr_RNG);
    boost::ecuyer1988& RNG = *ptr; 
 		return RNG;
  }

  template <>
  class InputParameter<boost::ecuyer1988&> {
    public:
      InputParameter(SEXP x_) : x(x_, false){}
      inline operator boost::ecuyer1988&() { 
        boost::ecuyer1988& RNG = *x;
        return RNG;
      }
    private:
      Rcpp::XPtr<boost::ecuyer1988, false> x ;
  };

  namespace traits {
    template <>
    struct input_parameter<boost::ecuyer1988&> {
      typedef typename Rcpp::InputParameter<boost::ecuyer1988&> type ;
    };
  }

}



