#include <Rcpp.h>
#include <boost/random/additive_combine.hpp>
#include <stan/math/prim/math/fun/accumulator.hpp>

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

RcppExport SEXP get_accumulator_(SEXP start) {
  int start_ = Rcpp::as<int>(start);
  stan::math::accumulator<double>* acc = new stan::math::accumulator<double>(start_);
  Rcpp::XPtr<stan::math::accumulator<double>> ptr(acc, true);
  return ptr;
}


