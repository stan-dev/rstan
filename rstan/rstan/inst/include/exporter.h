#include <RcppCommon.h>
#include <boost/random/additive_combine.hpp>
#include <stan/math/prim/mat/fun/accumulator.hpp>
#include <iostream>

namespace Rcpp {
  SEXP wrap(boost::ecuyer1988 RNG);
  SEXP wrap(boost::ecuyer1988& RNG);
  SEXP wrap(stan::math::accumulator<double>& RNG);
  SEXP wrap(std::ostream stream);
  template <> boost::ecuyer1988 as(SEXP ptr_RNG);
  template <> boost::ecuyer1988& as(SEXP ptr_RNG);
  template <> stan::math::accumulator<double>& as(SEXP ptr_acc);
  template <> std::ostream* as(SEXP ptr_stream);
  namespace traits {
    template <> class Exporter<boost::ecuyer1988&>;
    template <> struct input_parameter<boost::ecuyer1988&>;
    template <> class Exporter<stan::math::accumulator<double>&>;
    template <> struct input_parameter<stan::math::accumulator<double>&>;
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

  SEXP wrap(stan::math::accumulator<double>& ACC){
    stan::math::accumulator<double>* ptr_ACC = &ACC;
    Rcpp::XPtr<stan::math::accumulator<double>> Xptr_ACC(ptr_ACC);
    return Xptr_ACC;
  }

  SEXP wrap(std::ostream stream) {
    std::ostream* ptr_stream = &stream;
    Rcpp::XPtr<std::ostream> Xptr_stream(ptr_stream);
    return Xptr_stream;
  }

  template <> boost::ecuyer1988 as(SEXP ptr_RNG) {
    Rcpp::XPtr<boost::ecuyer1988> ptr(ptr_RNG);
    boost::ecuyer1988& RNG = *ptr; 
    return RNG;
  }

  template <> boost::ecuyer1988& as(SEXP ptr_RNG) {
    Rcpp::XPtr<boost::ecuyer1988> ptr(ptr_RNG);
    boost::ecuyer1988& RNG = *ptr; 
    return RNG;
  }

  template <> stan::math::accumulator<double>& as(SEXP ptr_ACC) {
    Rcpp::XPtr<stan::math::accumulator<double>> ptr(ptr_ACC);
    stan::math::accumulator<double>& ACC = *ptr; 
    return ACC;
  }

  template <> std::ostream* as(SEXP ptr_stream) {
    Rcpp::XPtr<std::ostream> ptr(ptr_stream);
    return ptr;
  }


  namespace traits {
    template <> class Exporter<boost::ecuyer1988&> {
    public:
      Exporter( SEXP x ) : t(Rcpp::as<boost::ecuyer1988&>(x)) {}
      inline boost::ecuyer1988& get(){ return t ; }
    private:
      boost::ecuyer1988& t ;
    } ; 

    template <> class Exporter<stan::math::accumulator<double>&> {
    public:
      Exporter( SEXP x ) : t(Rcpp::as<stan::math::accumulator<double>&>(x)) {}
      inline stan::math::accumulator<double>& get(){ return t ; }
    private:
      stan::math::accumulator<double>& t ;
    } ; 

    template <>
    struct input_parameter<boost::ecuyer1988&> {
      typedef typename Rcpp::ConstReferenceInputParameter<boost::ecuyer1988&> type ;
    };

    template <>
    struct input_parameter<stan::math::accumulator<double>&> {
      typedef typename Rcpp::ConstReferenceInputParameter<stan::math::accumulator<double>&> type ;
    };
  }

}

