#ifndef RSTAN_EXPORTER_H
#define RSTAN_EXPORTER_H

#include <RcppCommon.h>
#include <boost/random/mixmax.hpp>
#include <iostream>

namespace Rcpp {
  SEXP wrap(boost::random::mixmax RNG);
  SEXP wrap(boost::random::mixmax& RNG);
  SEXP wrap(std::ostream stream);
  template <> boost::random::mixmax as(SEXP ptr_RNG);
  template <> boost::random::mixmax& as(SEXP ptr_RNG);
  template <> std::ostream* as(SEXP ptr_stream);
  namespace traits {
    template <> class Exporter<boost::random::mixmax&>;
    template <> struct input_parameter<boost::random::mixmax&>;
  }
}


#include <Rcpp.h>

namespace Rcpp {
  SEXP wrap(boost::random::mixmax RNG){
    boost::random::mixmax* ptr_RNG = &RNG;
    Rcpp::XPtr<boost::random::mixmax> Xptr_RNG(ptr_RNG);
    return Xptr_RNG;
  }

  SEXP wrap(boost::random::mixmax& RNG){
    boost::random::mixmax* ptr_RNG = &RNG;
    Rcpp::XPtr<boost::random::mixmax> Xptr_RNG(ptr_RNG);
    return Xptr_RNG;
  }

  SEXP wrap(std::ostream stream) {
    std::ostream* ptr_stream = &stream;
    Rcpp::XPtr<std::ostream> Xptr_stream(ptr_stream);
    return Xptr_stream;
  }

  template <> boost::random::mixmax as(SEXP ptr_RNG) {
    Rcpp::XPtr<boost::random::mixmax> ptr(ptr_RNG);
    boost::random::mixmax& RNG = *ptr; 
 		return RNG;
  }

  template <> boost::random::mixmax& as(SEXP ptr_RNG) {
    Rcpp::XPtr<boost::random::mixmax> ptr(ptr_RNG);
    boost::random::mixmax& RNG = *ptr; 
 		return RNG;
  }

  template <> std::ostream* as(SEXP ptr_stream) {
    Rcpp::XPtr<std::ostream> ptr(ptr_stream);
    return ptr;
  }


  namespace traits {
    template <> class Exporter<boost::random::mixmax&> {
    public:
      Exporter( SEXP x ) : t(Rcpp::as<boost::random::mixmax&>(x)) {}
      inline boost::random::mixmax& get(){ return t ; }
    private:
      boost::random::mixmax& t ;
    } ; 

    template <>
    struct input_parameter<boost::random::mixmax&> {
      typedef typename Rcpp::ConstReferenceInputParameter<boost::random::mixmax&> type ;
      //typedef typename boost::random::mixmax& type ;
    };
  }

}

#endif
