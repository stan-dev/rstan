#include <Rcpp.h>
#include <boost/random/additive_combine.hpp>

RcppExport SEXP get_rng(SEXP seed) {
  int seed_ = Rcpp::as<int>(seed);
  boost::ecuyer1988* base_rng__= new boost::ecuyer1988(seed_);
  Rcpp::XPtr<boost::ecuyer1988> ptr(base_rng__,true);
  return ptr;
}





