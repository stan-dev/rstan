#include <Rcpp.h>
#include <boost/random/mixmax.hpp>

RcppExport SEXP get_stream_() {
  std::ostream* pstream(&Rcpp::Rcout);
  Rcpp::XPtr<std::ostream> ptr(pstream, false);
  return ptr;
}

RcppExport SEXP get_rng_(SEXP seed) {
  int seed_ = Rcpp::as<int>(seed);
  boost::random::mixmax* rng = new boost::random::mixmax(seed_);
  Rcpp::XPtr<boost::random::mixmax> ptr(rng, true);
  return ptr;
}


