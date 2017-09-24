#include <Rcpp.h>
#include <boost/random/additive_combine.hpp>

RcppExport SEXP get_stream() {
  std::ostream* pstream(&Rcpp::Rcout);
  Rcpp::XPtr<std::ostream> ptr(pstream, false);
  return ptr;
}

RcppExport SEXP get_rng(SEXP seed) {
  int seed_ = Rcpp::as<int>(seed);
  boost::ecuyer1988* rng = new boost::ecuyer1988(seed_);
  Rcpp::XPtr<boost::ecuyer1988> ptr(rng, true);
  return ptr;
}


