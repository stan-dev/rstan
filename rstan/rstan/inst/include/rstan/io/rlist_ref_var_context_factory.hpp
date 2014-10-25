#ifndef __RSTAN__IO__RLIST_REF_VAR_CONTEXT_FACTORY_HPP__
#define __RSTAN__IO__RLIST_REF_VAR_CONTEXT_FACTORY_HPP__

#include <stan/common/context_factory.hpp>
#include <rstan/io/rlist_ref_var_context.hpp>

namespace rstan {
  namespace io {
    class rlist_ref_var_context_factory 
      : public stan::common::var_context_factory<rlist_ref_var_context> {
    public:
      rlist_ref_var_context_factory(SEXP in) :
	rlist_(in) {
      }
      
      rlist_ref_var_context operator()(const std::string source) {
	return rstan::io::rlist_ref_var_context(rlist_);
      }

      Rcpp::List rlist_;
    };
  }
}

#endif
