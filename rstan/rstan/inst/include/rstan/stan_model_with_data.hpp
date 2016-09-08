#ifndef RSTAN__STAN_MODEL_HPP
#define RSTAN__STAN_MODEL_HPP

#include <stan/services/diagnose/diagnose.hpp>
#include <stan/services/optimize/bfgs.hpp>
#include <stan/services/optimize/lbfgs.hpp>
#include <stan/services/optimize/newton.hpp>
#include <stan/services/sample/fixed_param.hpp>
#include <stan/services/sample/hmc_nuts_dense_e.hpp>
#include <stan/services/sample/hmc_nuts_dense_e_adapt.hpp>
#include <stan/services/sample/hmc_nuts_diag_e.hpp>
#include <stan/services/sample/hmc_nuts_diag_e_adapt.hpp>
#include <stan/services/sample/hmc_nuts_unit_e.hpp>
#include <stan/services/sample/hmc_nuts_unit_e_adapt.hpp>
#include <stan/services/sample/hmc_static_dense_e.hpp>
#include <stan/services/sample/hmc_static_dense_e_adapt.hpp>
#include <stan/services/sample/hmc_static_diag_e.hpp>
#include <stan/services/sample/hmc_static_diag_e_adapt.hpp>
#include <stan/services/sample/hmc_static_unit_e.hpp>
#include <stan/services/sample/hmc_static_unit_e_adapt.hpp>
#include <stan/services/experimental/advi/fullrank.hpp>
#include <stan/services/experimental/advi/meanfield.hpp>

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
    stan::callbacks::noop_interrupt interrupt_;
    stan::callbacks::noop_writer message_writer_;
    stan::callbacks::noop_writer init_writer_;
    stan::callbacks::stream_writer parameter_writer_;
    
  public:
    stan_model_with_data(SEXP data)
      : data_(data),
        model_(data_, &rstan::io::rcout),
        interrupt_(),
        message_writer_(),
        init_writer_(),
        parameter_writer_(rstan::io::rcout) {
    }
    
    int diagnose(const Rcpp::List& init_list, int random_seed, int id, double init_radius, double epsilon, double error) {
      io::rlist_ref_var_context init_context(init_list);
      return stan::services::diagnose::diagnose(model_,
                                                init_context,
                                                random_seed,
                                                id,
                                                init_radius,
                                                epsilon,
                                                error,
                                                interrupt_,
                                                message_writer_,
                                                init_writer_,
                                                parameter_writer_);
    }
    
    int optimize_bfgs(const Rcpp::List& init_list, int random_seed, int id, double init_radius,
                      double init_alpha, double tol_obj, double tol_rel_obj,
                      double tol_grad, double tol_rel_grad, double tol_param,
                      int num_iterations, bool save_iterations, int refresh) {
      io::rlist_ref_var_context init_context(init_list);
      return stan::services::optimize::bfgs(model_,
                                            init_context,
                                            random_seed,
                                            id,
                                            init_radius,
                                            init_alpha,
                                            tol_obj,
                                            tol_rel_obj,
                                            tol_grad,
                                            tol_rel_grad,
                                            tol_param,
                                            num_iterations,
                                            save_iterations,
                                            refresh,
                                            interrupt_,
                                            message_writer_,
                                            init_writer_,
                                            parameter_writer_);
      
    }

    int optimize_lbfgs(const Rcpp::List& init_list, int random_seed, int id, double init_radius,
                       int history_size, double init_alpha, double tol_obj, double tol_rel_obj,
                       double tol_grad, double tol_rel_grad, double tol_param,
                       int num_iterations, bool save_iterations, int refresh) {
      io::rlist_ref_var_context init_context(init_list);
      return stan::services::optimize::lbfgs(model_,
                                             init_context,
                                             random_seed,
                                             id,
                                             init_radius,
                                             history_size,
                                             init_alpha,
                                             tol_obj,
                                             tol_rel_obj,
                                             tol_grad,
                                             tol_rel_grad,
                                             tol_param,
                                             num_iterations,
                                             save_iterations,
                                             refresh,
                                             interrupt_,
                                             message_writer_,
                                             init_writer_,
                                             parameter_writer_);
    }

    int optimize_newton(const Rcpp::List& init_list, int random_seed, int id, double init_radius,
                        int num_iterations, bool save_iterations) {
      io::rlist_ref_var_context init_context(init_list);
      return stan::services::optimize::newton(model_,
                                              init_context,
                                              random_seed,
                                              id,
                                              init_radius,
                                              num_iterations,
                                              save_iterations,
                                              interrupt_,
                                              message_writer_,
                                              init_writer_,
                                              parameter_writer_);
    }
    
  };
}
#endif
