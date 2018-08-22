#include <Rcpp.h>
#include <boost/random/additive_combine.hpp>
#include <stan/math/prim/mat/fun/accumulator.hpp>

RcppExport SEXP get_stream_() {
  std::ostream* pstream(&Rcpp::Rcout);
  Rcpp::XPtr<std::ostream> ptr(pstream, false);
  return ptr;
}

RcppExport SEXP get_rng_(SEXP seed) {
  int seed_ = Rcpp::as<int>(seed);
  boost::ecuyer1988* rng = new boost::ecuyer1988(seed_);
  Rcpp::XPtr<boost::ecuyer1988> ptr(rng, true);
  return ptr;
}

RcppExport SEXP get_accumulator_() {
  stan::math::accumulator<double>* acc = new stan::math::accumulator<double>();
  Rcpp::XPtr<stan::math::accumulator<double>> ptr(acc, true);
  return ptr;
}

RcppExport SEXP check_accumulator_(SEXP ptr_ACC) {
  Rcpp::XPtr<stan::math::accumulator<double>> ptr(ptr_ACC);
  stan::math::accumulator<double>& acc = *ptr;
  double total = acc.sum();
  return Rcpp::wrap(total);
}


