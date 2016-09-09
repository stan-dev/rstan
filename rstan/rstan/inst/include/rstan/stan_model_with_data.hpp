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
    stan::callbacks::noop_writer error_writer_;
    stan::callbacks::noop_writer init_writer_;
    stan::callbacks::stream_writer parameter_writer_;
    stan::callbacks::noop_writer diagnostic_writer_;

  public:
    stan_model_with_data(SEXP data)
      : data_(data),
        model_(data_, &rstan::io::rcout),
        interrupt_(),
        message_writer_(),
        error_writer_(),
        init_writer_(),
        parameter_writer_(rstan::io::rcout),
        diagnostic_writer_() {
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


    int sample_fixed_param(const Rcpp::List& init_list, int random_seed, int id, double init_radius,
                           int num_samples,
                           int num_thin,
                           int refresh) {
      io::rlist_ref_var_context init_context(init_list);
      return stan::services::sample::fixed_param(model_,
                                                 init_context,
                                                 random_seed,
                                                 id,
                                                 init_radius,
                                                 num_samples,
                                                 num_thin,
                                                 refresh,
                                                 interrupt_,
                                                 message_writer_,
                                                 error_writer_,
                                                 init_writer_,
                                                 parameter_writer_,
                                                 diagnostic_writer_);
    }

    int sample_hmc_nuts_dense_e(const Rcpp::List& init_list, int random_seed, int id, double init_radius,
                                int num_warmup,
                                int num_samples,
                                int num_thin,
                                bool save_warmup,
                                int refresh,
                                double stepsize,
                                double stepsize_jitter,
                                int max_depth) {
      io::rlist_ref_var_context init_context(init_list);
      return stan::services::sample::hmc_nuts_dense_e(model_,
                                                      init_context,
                                                      random_seed,
                                                      id,
                                                      init_radius,
                                                      num_warmup,
                                                      num_samples,
                                                      num_thin,
                                                      save_warmup,
                                                      refresh,
                                                      stepsize,
                                                      stepsize_jitter,
                                                      max_depth,
                                                      interrupt_,
                                                      message_writer_,
                                                      error_writer_,
                                                      init_writer_,
                                                      parameter_writer_,
                                                      diagnostic_writer_);
    }

    int sample_hmc_nuts_dense_e_adapt(const Rcpp::List& init_list, int random_seed, int id, double init_radius,
                                      int num_warmup,
                                      int num_samples,
                                      int num_thin,
                                      bool save_warmup,
                                      int refresh,
                                      double stepsize,
                                      double stepsize_jitter,
                                      int max_depth,
                                      double delta,
                                      double gamma,
                                      double kappa,
                                      double t0,
                                      unsigned int init_buffer,
                                      unsigned int term_buffer,
                                      unsigned int window) {
      io::rlist_ref_var_context init_context(init_list);
      return stan::services::sample::hmc_nuts_dense_e_adapt(model_,
                                                            init_context,
                                                            random_seed,
                                                            id,
                                                            init_radius,
                                                            num_warmup,
                                                            num_samples,
                                                            num_thin,
                                                            save_warmup,
                                                            refresh,
                                                            stepsize,
                                                            stepsize_jitter,
                                                            max_depth,
                                                            delta,
                                                            gamma,
                                                            kappa,
                                                            t0,
                                                            init_buffer,
                                                            term_buffer,
                                                            window,
                                                            interrupt_,
                                                            message_writer_,
                                                            error_writer_,
                                                            init_writer_,
                                                            parameter_writer_,
                                                            diagnostic_writer_);
    }

    int sample_hmc_nuts_diag_e(const Rcpp::List& init_list, int random_seed, int id, double init_radius,
                               int num_warmup,
                               int num_samples,
                               int num_thin,
                               bool save_warmup,
                               int refresh,
                               double stepsize,
                               double stepsize_jitter,
                               int max_depth) {
      io::rlist_ref_var_context init_context(init_list);
      return stan::services::sample::hmc_nuts_diag_e(model_,
                                                     init_context,
                                                     random_seed,
                                                     id,
                                                     init_radius,
                                                     num_warmup,
                                                     num_samples,
                                                     num_thin,
                                                     save_warmup,
                                                     refresh,
                                                     stepsize,
                                                     stepsize_jitter,
                                                     max_depth,
                                                     interrupt_,
                                                     message_writer_,
                                                     error_writer_,
                                                     init_writer_,
                                                     parameter_writer_,
                                                     diagnostic_writer_);
    }

    int sample_hmc_nuts_diag_e_adapt(const Rcpp::List& init_list, int random_seed, int id, double init_radius,
                                     int num_warmup,
                                     int num_samples,
                                     int num_thin,
                                     bool save_warmup,
                                     int refresh,
                                     double stepsize,
                                     double stepsize_jitter,
                                     int max_depth,
                                     double delta,
                                     double gamma,
                                     double kappa,
                                     double t0,
                                     unsigned int init_buffer,
                                     unsigned int term_buffer,
                                     unsigned int window) {
      io::rlist_ref_var_context init_context(init_list);
      return stan::services::sample::hmc_nuts_diag_e_adapt(model_,
                                                           init_context,
                                                           random_seed,
                                                           id,
                                                           init_radius,
                                                           num_warmup,
                                                           num_samples,
                                                           num_thin,
                                                           save_warmup,
                                                           refresh,
                                                           stepsize,
                                                           stepsize_jitter,
                                                           max_depth,
                                                           delta,
                                                           gamma,
                                                           kappa,
                                                           t0,
                                                           init_buffer,
                                                           term_buffer,
                                                           window,
                                                           interrupt_,
                                                           message_writer_,
                                                           error_writer_,
                                                           init_writer_,
                                                           parameter_writer_,
                                                           diagnostic_writer_);
    }

    int sample_hmc_nuts_unit_e(const Rcpp::List& init_list, int random_seed, int id, double init_radius,
                               int num_warmup,
                               int num_samples,
                               int num_thin,
                               bool save_warmup,
                               int refresh,
                               double stepsize,
                               double stepsize_jitter,
                               int max_depth) {
      io::rlist_ref_var_context init_context(init_list);
      return stan::services::sample::hmc_nuts_unit_e(model_,
                                                     init_context,
                                                     random_seed,
                                                     id,
                                                     init_radius,
                                                     num_warmup,
                                                     num_samples,
                                                     num_thin,
                                                     save_warmup,
                                                     refresh,
                                                     stepsize,
                                                     stepsize_jitter,
                                                     max_depth,
                                                     interrupt_,
                                                     message_writer_,
                                                     error_writer_,
                                                     init_writer_,
                                                     parameter_writer_,
                                                     diagnostic_writer_);
    }

    int sample_hmc_nuts_unit_e_adapt(const Rcpp::List& init_list, int random_seed, int id, double init_radius,
                                     int num_warmup,
                                     int num_samples,
                                     int num_thin,
                                     bool save_warmup,
                                     int refresh,
                                     double stepsize,
                                     double stepsize_jitter,
                                     int max_depth,
                                     double delta,
                                     double gamma,
                                     double kappa,
                                     double t0) {
      io::rlist_ref_var_context init_context(init_list);
      return stan::services::sample::hmc_nuts_unit_e_adapt(model_,
                                                           init_context,
                                                           random_seed,
                                                           id,
                                                           init_radius,
                                                           num_warmup,
                                                           num_samples,
                                                           num_thin,
                                                           save_warmup,
                                                           refresh,
                                                           stepsize,
                                                           stepsize_jitter,
                                                           max_depth,
                                                           delta,
                                                           gamma,
                                                           kappa,
                                                           t0,
                                                           interrupt_,
                                                           message_writer_,
                                                           error_writer_,
                                                           init_writer_,
                                                           parameter_writer_,
                                                           diagnostic_writer_);
    }

    int sample_hmc_static_dense_e(const Rcpp::List& init_list, int random_seed, int id, double init_radius,
                                  int num_warmup,
                                  int num_samples,
                                  int num_thin,
                                  bool save_warmup,
                                  int refresh,
                                  double stepsize,
                                  double stepsize_jitter,
                                  double int_time) {
      io::rlist_ref_var_context init_context(init_list);
      return stan::services::sample::hmc_static_dense_e(model_,
                                                        init_context,
                                                        random_seed,
                                                        id,
                                                        init_radius,
                                                        num_warmup,
                                                        num_samples,
                                                        num_thin,
                                                        save_warmup,
                                                        refresh,
                                                        stepsize,
                                                        stepsize_jitter,
                                                        int_time,
                                                        interrupt_,
                                                        message_writer_,
                                                        error_writer_,
                                                        init_writer_,
                                                        parameter_writer_,
                                                        diagnostic_writer_);
    }

    int sample_hmc_static_dense_e_adapt(const Rcpp::List& init_list, int random_seed, int id, double init_radius,
                                        int num_warmup,
                                        int num_samples,
                                        int num_thin,
                                        bool save_warmup,
                                        int refresh,
                                        double stepsize,
                                        double stepsize_jitter,
                                        double int_time,
                                        double delta,
                                        double gamma,
                                        double kappa,
                                        double t0,
                                        unsigned int init_buffer,
                                        unsigned int term_buffer,
                                        unsigned int window) {
      io::rlist_ref_var_context init_context(init_list);
      return stan::services::sample::hmc_static_dense_e_adapt(model_,
                                                              init_context,
                                                              random_seed,
                                                              id,
                                                              init_radius,
                                                              num_warmup,
                                                              num_samples,
                                                              num_thin,
                                                              save_warmup,
                                                              refresh,
                                                              stepsize,
                                                              stepsize_jitter,
                                                              int_time,
                                                              delta,
                                                              gamma,
                                                              kappa,
                                                              t0,
                                                              init_buffer,
                                                              term_buffer,
                                                              window,
                                                              interrupt_,
                                                              message_writer_,
                                                              error_writer_,
                                                              init_writer_,
                                                              parameter_writer_,
                                                              diagnostic_writer_);
    }

    int sample_hmc_static_diag_e(const Rcpp::List& init_list, int random_seed, int id, double init_radius,
                                 int num_warmup,
                                 int num_samples,
                                 int num_thin,
                                 bool save_warmup,
                                 int refresh,
                                 double stepsize,
                                 double stepsize_jitter,
                                 double int_time) {
      io::rlist_ref_var_context init_context(init_list);
      return stan::services::sample::hmc_static_diag_e(model_,
                                                       init_context,
                                                       random_seed,
                                                       id,
                                                       init_radius,
                                                       num_warmup,
                                                       num_samples,
                                                       num_thin,
                                                       save_warmup,
                                                       refresh,
                                                       stepsize,
                                                       stepsize_jitter,
                                                       int_time,
                                                       interrupt_,
                                                       message_writer_,
                                                       error_writer_,
                                                       init_writer_,
                                                       parameter_writer_,
                                                       diagnostic_writer_);
    }

    int sample_hmc_static_diag_e_adapt(const Rcpp::List& init_list, int random_seed, int id, double init_radius,
                                       int num_warmup,
                                       int num_samples,
                                       int num_thin,
                                       bool save_warmup,
                                       int refresh,
                                       double stepsize,
                                       double stepsize_jitter,
                                       double int_time,
                                       double delta,
                                       double gamma,
                                       double kappa,
                                       double t0,
                                       unsigned int init_buffer,
                                       unsigned int term_buffer,
                                       unsigned int window) {
      io::rlist_ref_var_context init_context(init_list);
      return stan::services::sample::hmc_static_diag_e_adapt(model_,
                                                             init_context,
                                                             random_seed,
                                                             id,
                                                             init_radius,
                                                             num_warmup,
                                                             num_samples,
                                                             num_thin,
                                                             save_warmup,
                                                             refresh,
                                                             stepsize,
                                                             stepsize_jitter,
                                                             int_time,
                                                             delta,
                                                             gamma,
                                                             kappa,
                                                             t0,
                                                             init_buffer,
                                                             term_buffer,
                                                             window,
                                                             interrupt_,
                                                             message_writer_,
                                                             error_writer_,
                                                             init_writer_,
                                                             parameter_writer_,
                                                             diagnostic_writer_);
    }

    int sample_hmc_static_unit_e(const Rcpp::List& init_list, int random_seed, int id, double init_radius,
                                 int num_warmup,
                                 int num_samples,
                                 int num_thin,
                                 bool save_warmup,
                                 int refresh,
                                 double stepsize,
                                 double stepsize_jitter,
                                 double int_time) {
      io::rlist_ref_var_context init_context(init_list);
      return stan::services::sample::hmc_static_unit_e(model_,
                                                       init_context,
                                                       random_seed,
                                                       id,
                                                       init_radius,
                                                       num_warmup,
                                                       num_samples,
                                                       num_thin,
                                                       save_warmup,
                                                       refresh,
                                                       stepsize,
                                                       stepsize_jitter,
                                                       int_time,
                                                       interrupt_,
                                                       message_writer_,
                                                       error_writer_,
                                                       init_writer_,
                                                       parameter_writer_,
                                                       diagnostic_writer_);
    }

    int sample_hmc_static_unit_e_adapt(const Rcpp::List& init_list, int random_seed, int id, double init_radius,
                                       int num_warmup,
                                       int num_samples,
                                       int num_thin,
                                       bool save_warmup,
                                       int refresh,
                                       double stepsize,
                                       double stepsize_jitter,
                                       double int_time,
                                       double delta,
                                       double gamma,
                                       double kappa,
                                       double t0) {
      io::rlist_ref_var_context init_context(init_list);
      return stan::services::sample::hmc_static_unit_e_adapt(model_,
                                                             init_context,
                                                             random_seed,
                                                             id,
                                                             init_radius,
                                                             num_warmup,
                                                             num_samples,
                                                             num_thin,
                                                             save_warmup,
                                                             refresh,
                                                             stepsize,
                                                             stepsize_jitter,
                                                             int_time,
                                                             delta,
                                                             gamma,
                                                             kappa,
                                                             t0,
                                                             interrupt_,
                                                             message_writer_,
                                                             error_writer_,
                                                             init_writer_,
                                                             parameter_writer_,
                                                             diagnostic_writer_);
    }

    int experimental_advi_fullrank(const Rcpp::List& init_list, int random_seed, int id, double init_radius,
                                   int grad_samples,
                                   int elbo_samples,
                                   int max_iterations,
                                   double tol_rel_obj,
                                   double eta,
                                   bool adapt_engaged,
                                   int adapt_iterations,
                                   int eval_elbo,
                                   int output_samples) {
      io::rlist_ref_var_context init_context(init_list);
      return stan::services::experimental::advi::fullrank(model_,
                                                          init_context,
                                                          random_seed,
                                                          id,
                                                          init_radius,
                                                          grad_samples,
                                                          elbo_samples,
                                                          max_iterations,
                                                          tol_rel_obj,
                                                          eta,
                                                          adapt_engaged,
                                                          adapt_iterations,
                                                          eval_elbo,
                                                          output_samples,
                                                          interrupt_,
                                                          message_writer_,
                                                          init_writer_,
                                                          parameter_writer_,
                                                          diagnostic_writer_);
    }

    int experimental_advi_meanfield(const Rcpp::List& init_list, int random_seed, int id, double init_radius,
                                    int grad_samples,
                                    int elbo_samples,
                                    int max_iterations,
                                    double tol_rel_obj,
                                    double eta,
                                    bool adapt_engaged,
                                    int adapt_iterations,
                                    int eval_elbo,
                                    int output_samples) {
      io::rlist_ref_var_context init_context(init_list);
      return stan::services::experimental::advi::meanfield(model_,
                                                           init_context,
                                                           random_seed,
                                                           id,
                                                           init_radius,
                                                           grad_samples,
                                                           elbo_samples,
                                                           max_iterations,
                                                           tol_rel_obj,
                                                           eta,
                                                           adapt_engaged,
                                                           adapt_iterations,
                                                           eval_elbo,
                                                           output_samples,
                                                           interrupt_,
                                                           message_writer_,
                                                           init_writer_,
                                                           parameter_writer_,
                                                           diagnostic_writer_);
    }
    
  };
}
#endif
