#ifndef RSTAN_EXPORTER_H
#define RSTAN_EXPORTER_H

#include <RcppCommon.h>
#include <stan/services/util/create_rng.hpp>
#include <iostream>

namespace Rcpp {
  SEXP wrap(stan::rng_t RNG);
  SEXP wrap(stan::rng_t& RNG);
  SEXP wrap(std::ostream stream);
  template <> stan::rng_t as(SEXP ptr_RNG);
  template <> stan::rng_t& as(SEXP ptr_RNG);
  template <> std::ostream* as(SEXP ptr_stream);
  namespace traits {
    template <> class Exporter<stan::rng_t&>;
    template <> struct input_parameter<stan::rng_t&>;
  }
}


#include <Rcpp.h>

namespace Rcpp {
  SEXP wrap(stan::rng_t RNG){
    stan::rng_t* ptr_RNG = &RNG;
    Rcpp::XPtr<stan::rng_t> Xptr_RNG(ptr_RNG);
    return Xptr_RNG;
  }

  SEXP wrap(stan::rng_t& RNG){
    stan::rng_t* ptr_RNG = &RNG;
    Rcpp::XPtr<stan::rng_t> Xptr_RNG(ptr_RNG);
    return Xptr_RNG;
  }

  SEXP wrap(std::ostream stream) {
    std::ostream* ptr_stream = &stream;
    Rcpp::XPtr<std::ostream> Xptr_stream(ptr_stream);
    return Xptr_stream;
  }

  template <> stan::rng_t as(SEXP ptr_RNG) {
    Rcpp::XPtr<stan::rng_t> ptr(ptr_RNG);
    stan::rng_t& RNG = *ptr;
 		return RNG;
  }

  template <> stan::rng_t& as(SEXP ptr_RNG) {
    Rcpp::XPtr<stan::rng_t> ptr(ptr_RNG);
    stan::rng_t& RNG = *ptr;
 		return RNG;
  }

  template <> std::ostream* as(SEXP ptr_stream) {
    Rcpp::XPtr<std::ostream> ptr(ptr_stream);
    return ptr;
  }


  namespace traits {
    template <> class Exporter<stan::rng_t&> {
    public:
      Exporter( SEXP x ) : t(Rcpp::as<stan::rng_t&>(x)) {}
      inline stan::rng_t& get(){ return t ; }
    private:
      stan::rng_t& t ;
    } ;

    template <>
    struct input_parameter<stan::rng_t&> {
      typedef typename Rcpp::ConstReferenceInputParameter<stan::rng_t&> type ;
      //typedef typename stan::rng_t& type ;
    };
  }

}

#endif
