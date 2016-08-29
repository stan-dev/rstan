#ifndef RSTAN__STAN_FIT_HPP
#define RSTAN__STAN_FIT_HPP

#include <Rcpp.h>
#include <rstan/io/rlist_ref_var_context.hpp>
#include <rstan/io/r_ostream.hpp>

namespace rstan {

  template <class Model>
  class stan_fit {
  private: 
    io::rlist_ref_var_context data_;
    Model model_;
    Rcpp::Function cxxfunction; // keep a reference to the cxxfun, no functional purpose.

    stan_fit(SEXP data, SEXP cxxf) :
      data_(data),
      model_(data_, &rstan::io::rcout),
      cxxfunction(cxxf) {
    }
  };

}
