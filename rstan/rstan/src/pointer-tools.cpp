#include <Rcpp.h>
#include <stan/services/util/create_rng.hpp>

RcppExport SEXP get_stream_() {
  std::ostream* pstream(&Rcpp::Rcout);
  Rcpp::XPtr<std::ostream> ptr(pstream, false);
  return ptr;
}

RcppExport SEXP get_rng_(SEXP seed) {
  int seed_ = Rcpp::as<int>(seed);
  stan::rng_t* rng = new stan::rng_t(seed_);
  Rcpp::XPtr<stan::rng_t> ptr(rng, true);
  return ptr;
}
