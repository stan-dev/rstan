#ifndef RSTAN__IO__RLIST_REF_VAR_CONTEXT_FACTORY_HPP
#define RSTAN__IO__RLIST_REF_VAR_CONTEXT_FACTORY_HPP

#include <stan/interface/var_context_factory/dump_factory.hpp>
#include <stan/interface/var_context_factory/var_context_factory.hpp>
#include <rstan/io/rlist_ref_var_context.hpp>


namespace rstan {
  namespace io {
    class rlist_ref_var_context_factory 
      : public stan::interface::var_context_factory::var_context_factory<rlist_ref_var_context> {
    public:
      rlist_ref_var_context_factory(SEXP in) : rlist_(in) { }
      rlist_ref_var_context operator()(const std::string source) {
        return rstan::io::rlist_ref_var_context(rlist_);
      }
      Rcpp::List rlist_;
    };
  }
}

#endif
