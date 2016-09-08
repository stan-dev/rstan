#ifndef RSTAN__STAN_MODEL_HPP
#define RSTAN__STAN_MODEL_HPP

#include <stan/services/diagnose/diagnose.hpp>
#include <stan/io/empty_var_context.hpp>
#include <stan/callbacks/noop_writer.hpp>
#include <stan/callbacks/stream_writer.hpp>
#include <stan/callbacks/noop_interrupt.hpp>
#include <rstan/io/rlist_ref_var_context.hpp>

namespace rstan {
  
  template <class Model>
  class stan_model_with_data {
  private:
    io::rlist_ref_var_context data_;    
    Model model_;
    
  public:
    stan_model_with_data(SEXP data)
      : data_(data),
        model_(data_, &rstan::io::rcout) {
    }
    
    int diagnose(const Rcpp::List& init_list, int random_seed, int chain, double init_radius, double epsilon, double error) {
      io::rlist_ref_var_context init(init_list);
      stan::callbacks::noop_interrupt interrupt;
      stan::callbacks::noop_writer message_writer;
      stan::callbacks::noop_writer init_writer;
      stan::callbacks::stream_writer parameter_writer(rstan::io::rcout);
      rstan::io::rcout << "init contains y: " << init.contains_r("y") << std::endl;
      return stan::services::diagnose::diagnose(model_,
                                                init,
                                                random_seed,
                                                chain,
                                                init_radius,
                                                epsilon,
                                                error,
                                                interrupt,
                                                message_writer,
                                                init_writer,
                                                parameter_writer);
    }
    
  };
}
#endif
