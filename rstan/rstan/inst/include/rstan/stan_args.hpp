#ifndef RSTAN__STAN_ARGS_HPP
#define RSTAN__STAN_ARGS_HPP


#include <Rcpp.h>
// #include <R.h>
// #include <Rinternals.h>

#include <algorithm>
#include <rstan/io/r_ostream.hpp>
#include <stan/version.hpp>
#include <boost/lexical_cast.hpp>

namespace rstan {

  namespace {
    /*
     * Get an element of Rcpp::List by name. If not found, set it to a
     * default value.
     * @param lst The list to look for elements
     * @param n The name of an element of interest
     * @param t Where to save the element
     * @param v0 The default value if not found in the list
     */
    template <class T>
    bool get_rlist_element(const Rcpp::List& lst, const char* n, T& t, const T& v0) {
      bool b = lst.containsElementNamed(n);
      if (b)  t = Rcpp::as<T>(const_cast<Rcpp::List&>(lst)[n]);
      else  t = T(v0);
      return b;
    }

    template <class T>
    bool get_rlist_element(const Rcpp::List& lst, const char* n, T& t) {
      bool b = lst.containsElementNamed(n);
      if (b) t = Rcpp::as<T>(const_cast<Rcpp::List&>(lst)[n]);
      return b;
    }

    template <>
    bool get_rlist_element(const Rcpp::List& lst, const char* n, SEXP& t) {
      bool b = lst.containsElementNamed(n);
      if (b) t = const_cast<Rcpp::List&>(lst)[n];
      return b;
    }

    inline unsigned int sexp2seed(SEXP seed) {
      if (TYPEOF(seed) == STRSXP)
        return boost::lexical_cast<unsigned int>(Rcpp::as<std::string>(seed));
      return Rcpp::as<unsigned int>(seed);
    }

    void write_comment(std::ostream& o) {
      o << "#" << std::endl;
    }

    template <typename M>
    void write_comment(std::ostream& o, const M& msg) {
      o << "# " << msg << std::endl;
    }

    template <typename K, typename V>
    void write_comment_property(std::ostream& o, const K& key, const V& val) {
      o << "# " << key << "=" << val << std::endl;
    }

    /**
     * Find the index of an element in a vector.
     * @param v the vector in which an element are searched.
     * @param e the element that we are looking for.
     * @return If e is in v, return the index (0 to size - 1);
     *  otherwise, return the size.
     */

    template <class T, class T2>
    size_t find_index(const std::vector<T>& v, const T2& e) {
      return std::distance(v.begin(), std::find(v.begin(), v.end(), T(e)));
    }
  }

  enum sampling_algo_t { NUTS = 1, HMC = 2, Metropolis = 3, Fixed_param = 4};
  enum optim_algo_t { Newton = 1, BFGS = 3, LBFGS = 4};
  enum variational_algo_t { MEANFIELD = 1, FULLRANK = 2};
  enum sampling_metric_t { UNIT_E = 1, DIAG_E = 2, DENSE_E = 3};
  enum stan_args_method_t { SAMPLING = 1, OPTIM = 2, TEST_GRADIENT = 3, VARIATIONAL = 4};

  /**
   *
   */
  class stan_args {
  private:
    unsigned int random_seed;
    unsigned int chain_id;
    std::string init;
    SEXP init_list;
    double init_radius;
    // FIXME(syclik): remove `enable_random_init`
    bool enable_random_init; // enable randomly partially specifying inits 
    std::string sample_file; // the file for outputting the samples
    bool append_samples;
    bool sample_file_flag; // true: write out to a file; false, do not
    stan_args_method_t method;
    std::string diagnostic_file;
    bool diagnostic_file_flag;
    union {
      struct {
        int iter;   // number of iterations
        int refresh;  //
        sampling_algo_t algorithm;
        int warmup; // number of warmup
        int thin;
        bool save_warmup; // weather to save warmup samples (true by default)
        int iter_save; // number of iterations saved
        int iter_save_wo_warmup; // number of iterations saved wo warmup
        bool adapt_engaged;
        double adapt_gamma;
        double adapt_delta;
        double adapt_kappa;
        unsigned int adapt_init_buffer;
        unsigned int adapt_term_buffer;
        unsigned int adapt_window;
        double adapt_t0;
        sampling_metric_t metric; // UNIT_E, DIAG_E, DENSE_E;
        double stepsize; // defaut to 1;
        double stepsize_jitter;
        int max_treedepth; // for NUTS, default to 10.
        double int_time; // for HMC, default to 2 * pi
      } sampling;
      struct {
        int iter; // default to 2000
        int refresh; // default to 100
        optim_algo_t algorithm; // Newton, (L)BFGS
        bool save_iterations; // default to false
        double init_alpha; // default to 0.001, for (L)BFGS
        double tol_obj; // default to 1e-12, for (L)BFGS
        double tol_grad; // default to 1e-8, for (L)BFGS
        double tol_param; // default to 1e-8, for (L)BFGS
        double tol_rel_obj; // default to 1e4, for (L)BFGS
        double tol_rel_grad; // default to 1e7, for (L)BFGS
        int history_size; // default to 5, for LBFGS only
      } optim;
      struct {
        int iter; // default to 10000
        variational_algo_t algorithm;  // MEANFIELD or FULLRANK
        int grad_samples; // default to 1
        int elbo_samples; // default to 100
        int eval_elbo;    // default to 100
        int output_samples; // default to 1000
        double eta; // defaults to 1.0
        bool adapt_engaged; // defaults to 1
        int adapt_iter; // defaults to 50
        double tol_rel_obj; // default to 0.01
        int refresh; // default to 1
      } variational;
      struct {
        double epsilon; // default to 1e-6, for test_grad
        double error;  // default to 1e-6, for test_grad
      } test_grad;
    } ctrl;

  private:
    void validate_args() {
      if (init_radius < 0) {
        std::stringstream msg;
        msg << "Invalid value for parameter init_r (found "
            << init_radius << "; require >= 0).";
        throw std::invalid_argument(msg.str());
      }
      switch (method) {
        case SAMPLING:
          if (ctrl.sampling.adapt_gamma < 0) {
            std::stringstream msg;
            msg << "Invalid adaptation parameter (found gamma="
                << ctrl.sampling.adapt_gamma << "; require >0).";
            throw std::invalid_argument(msg.str());
          }
          if (ctrl.sampling.adapt_delta <= 0 || ctrl.sampling.adapt_delta >= 1) {
            std::stringstream msg;
            msg << "Invalid adaptation parameter (found delta="
                << ctrl.sampling.adapt_delta << "; require 0<delta<1).";
            throw std::invalid_argument(msg.str());
          }
          if (ctrl.sampling.adapt_kappa < 0) {
            std::stringstream msg;
            msg << "Invalid adaptation parameter (found kappa="
                << ctrl.sampling.adapt_kappa << "; require >0).";
            throw std::invalid_argument(msg.str());
          }
          if (ctrl.sampling.adapt_t0 < 0) {
            std::stringstream msg;
            msg << "Invalid adaptation parameter (found t0="
                << ctrl.sampling.adapt_t0 << "; require >0).";
            throw std::invalid_argument(msg.str());
          }
          if (ctrl.sampling.stepsize < 0) {
            std::stringstream msg;
            msg << "Invalid adaptation parameter (found stepsize="
                << ctrl.sampling.stepsize << "; require stepsize > 0).";
            throw std::invalid_argument(msg.str());
          }
          if (ctrl.sampling.stepsize_jitter < 0 || ctrl.sampling.stepsize_jitter > 1) {
            std::stringstream msg;
            msg << "Invalid adaptation parameter (found stepsize_jitter="
                << ctrl.sampling.stepsize_jitter << "; require 0<=stepsize_jitter<=1).";
            throw std::invalid_argument(msg.str());
          }
          if (ctrl.sampling.algorithm == NUTS && ctrl.sampling.max_treedepth < 0) {
            std::stringstream msg;
            msg << "Invalid adaptation parameter (found max_treedepth="
                << ctrl.sampling.max_treedepth << "; require max_treedepth>0).";
            throw std::invalid_argument(msg.str());
          }
          if (ctrl.sampling.algorithm == HMC && ctrl.sampling.int_time < 0) {
            std::stringstream msg;
            msg << "Invalid adaptation parameter (found int_time="
                << ctrl.sampling.int_time << "; require int_time>0).";
            throw std::invalid_argument(msg.str());
          }
          break;
        case OPTIM:
          if (ctrl.optim.init_alpha < 0) {
            std::stringstream msg;
            msg << "Invalid adaptation parameter (found init_alpha="
                << ctrl.optim.init_alpha << "; require init_alpha > 0).";
            throw std::invalid_argument(msg.str());
          }
          break;
        case TEST_GRADIENT: break;
        case VARIATIONAL:
          if (ctrl.variational.grad_samples <= 0) {
            std::stringstream msg;
            msg << "Invalid parameter grad_samples (found grad_samples="
                << ctrl.variational.grad_samples << "; require 0 < grad_samples).";
            throw std::invalid_argument(msg.str());
          }
          if (ctrl.variational.elbo_samples <= 0) {
            std::stringstream msg;
            msg << "Invalid parameter elbo_samples (found elbo_samples="
                << ctrl.variational.elbo_samples << "; require 0 < elbo_samples).";
            throw std::invalid_argument(msg.str());
          }
          if (ctrl.variational.iter <= 0) {
            std::stringstream msg;
            msg << "Invalid parameter iter (found iter="
                << ctrl.variational.iter << "; require 0 < iter).";
            throw std::invalid_argument(msg.str());
          }
          if (ctrl.variational.tol_rel_obj <= 0) {
            std::stringstream msg;
            msg << "Invalid parameter tol_rel_obj (found tol_rel_obj="
                << ctrl.variational.tol_rel_obj << "; require 0 < tol_rel_obj).";
            throw std::invalid_argument(msg.str());
          }
          if (ctrl.variational.eta <= 0) {
            std::stringstream msg;
            msg << "Invalid parameter eta (found eta="
                << ctrl.variational.eta << "; require 0 < eta).";
            throw std::invalid_argument(msg.str());
          }
          if (ctrl.variational.eval_elbo <= 0) {
            std::stringstream msg;
            msg << "Invalid parameter eval_elbo (found eval_elbo="
                << ctrl.variational.eval_elbo << "; require 0 < eval_elbo).";
            throw std::invalid_argument(msg.str());
          }
          if (ctrl.variational.output_samples <= 0) {
            std::stringstream msg;
            msg << "Invalid parameter output_samples (found output_samples="
                << ctrl.variational.output_samples << "; require 0 < output_samples).";
            throw std::invalid_argument(msg.str());
          }
          if (ctrl.variational.adapt_iter <= 0) {
            std::stringstream msg;
            msg << "Invalid parameter adapt_iter (found adapt_iter="
                << ctrl.variational.adapt_iter << "; require 0 < adapt_iter).";
            throw std::invalid_argument(msg.str());
          }
          break;
      }
    }

  public:
    stan_args(const Rcpp::List& in) : init_list(R_NilValue) {

      std::string t_str;
      SEXP t_sexp;
      bool b;
      get_rlist_element(in, "chain_id", chain_id, static_cast<unsigned int>(1));
      get_rlist_element(in, "append_samples", append_samples, false);
      b = get_rlist_element(in, "method", t_str);
      if (!b) method = SAMPLING;
      else {
        if ("sampling" == t_str)  method = SAMPLING;
        else if ("optim" == t_str)  method = OPTIM;
        else if ("test_grad" == t_str)  method = TEST_GRADIENT;
        else if ("variational" == t_str) method = VARIATIONAL;
        else method = SAMPLING;
      }

      sample_file_flag = get_rlist_element(in, "sample_file", sample_file);
      diagnostic_file_flag = get_rlist_element(in, "diagnostic_file", diagnostic_file);
      b = get_rlist_element(in, "seed", t_sexp);
      if (b) random_seed = sexp2seed(t_sexp);
      else random_seed = std::time(0);

      int calculated_thin;
      get_rlist_element(in, "control", t_sexp, R_NilValue);
      Rcpp::List ctrl_lst(t_sexp);

      switch (method) {
        case VARIATIONAL:
          get_rlist_element(in, "iter", ctrl.variational.iter, 10000);
          get_rlist_element(in, "grad_samples", ctrl.variational.grad_samples, 1);
          get_rlist_element(in, "elbo_samples", ctrl.variational.elbo_samples, 100);
          get_rlist_element(in, "eval_elbo", ctrl.variational.eval_elbo, 100);
          get_rlist_element(in, "output_samples", ctrl.variational.output_samples, 1000);
          get_rlist_element(in, "adapt_iter", ctrl.variational.adapt_iter, 50);
          get_rlist_element(in, "eta", ctrl.variational.eta, 1.0);
          get_rlist_element(in, "adapt_engaged", ctrl.variational.adapt_engaged, true);
          get_rlist_element(in, "tol_rel_obj", ctrl.variational.tol_rel_obj, 0.01);
          get_rlist_element(in, "refresh", ctrl.variational.refresh, 1);
          ctrl.variational.algorithm = MEANFIELD;
          if (get_rlist_element(in, "algorithm", t_str)) {
            if (t_str == "fullrank") ctrl.variational.algorithm = FULLRANK;
          }
          break;
        case SAMPLING:
          get_rlist_element(in, "iter", ctrl.sampling.iter, 2000);
          get_rlist_element(in, "warmup", ctrl.sampling.warmup, ctrl.sampling.iter / 2);
          get_rlist_element(in, "save_warmup", ctrl.sampling.save_warmup, true);

          calculated_thin = (ctrl.sampling.iter - ctrl.sampling.warmup) / 1000;
          if (calculated_thin < 1) calculated_thin = 1;
          get_rlist_element(in, "thin", ctrl.sampling.thin, calculated_thin);

          ctrl.sampling.iter_save_wo_warmup
            = 1 + (ctrl.sampling.iter - ctrl.sampling.warmup - 1) / ctrl.sampling.thin;
          ctrl.sampling.iter_save
            = ctrl.sampling.iter_save_wo_warmup;
          if (ctrl.sampling.save_warmup)
            ctrl.sampling.iter_save += 
              1 + (ctrl.sampling.warmup - 1) / ctrl.sampling.thin;

          ctrl.sampling.refresh = (ctrl.sampling.iter >= 20) ?
                                  ctrl.sampling.iter / 10 : 1;
          get_rlist_element(in, "refresh", ctrl.sampling.refresh);

          get_rlist_element(ctrl_lst, "adapt_engaged", ctrl.sampling.adapt_engaged, true);
          get_rlist_element(ctrl_lst, "adapt_gamma", ctrl.sampling.adapt_gamma, 0.05);
          get_rlist_element(ctrl_lst, "adapt_delta", ctrl.sampling.adapt_delta, 0.8);
          get_rlist_element(ctrl_lst, "adapt_kappa", ctrl.sampling.adapt_kappa, 0.75);
          get_rlist_element(ctrl_lst, "adapt_t0", ctrl.sampling.adapt_t0, 10.0);
          get_rlist_element(ctrl_lst, "adapt_init_buffer", ctrl.sampling.adapt_init_buffer, 75U);
          get_rlist_element(ctrl_lst, "adapt_term_buffer", ctrl.sampling.adapt_term_buffer, 50U);
          get_rlist_element(ctrl_lst, "adapt_window", ctrl.sampling.adapt_window, 25U);
          get_rlist_element(ctrl_lst, "stepsize", ctrl.sampling.stepsize, 1.0);
          get_rlist_element(ctrl_lst, "stepsize_jitter", ctrl.sampling.stepsize_jitter, 0.0);

          if (get_rlist_element(in, "algorithm", t_str)) {
            if (t_str == "HMC") ctrl.sampling.algorithm = HMC;
            else if (t_str == "Metropolis") ctrl.sampling.algorithm = Metropolis;
            else if (t_str == "NUTS") ctrl.sampling.algorithm = NUTS;
            else if (t_str == "Fixed_param") {
              ctrl.sampling.algorithm = Fixed_param;
              ctrl.sampling.adapt_engaged = false;
              ctrl.sampling.warmup = 0;
              ctrl.sampling.iter_save_wo_warmup = 1 + (ctrl.sampling.iter - 1) / ctrl.sampling.thin;
              ctrl.sampling.iter_save = ctrl.sampling.iter_save_wo_warmup;
              ctrl.sampling.save_warmup = false;
            } else {
              std::stringstream msg;
              msg << "Invalid value for parameter algorithm (found "
                  << t_str << "; require HMC, Metropolis, Fixed_param, or NUTS).";
              throw std::invalid_argument(msg.str());
            }
          } else {
            ctrl.sampling.algorithm = NUTS;
          }

          if (get_rlist_element(ctrl_lst, "metric", t_str)) {
            if ("unit_e" == t_str) ctrl.sampling.metric = UNIT_E;
            else if ("diag_e" == t_str) ctrl.sampling.metric = DIAG_E;
            else if ("dense_e" == t_str) ctrl.sampling.metric = DENSE_E;
          } else ctrl.sampling.metric = DIAG_E;

          switch (ctrl.sampling.algorithm) {
            case NUTS:
              get_rlist_element(ctrl_lst, "max_treedepth", ctrl.sampling.max_treedepth, 10);
              break;
             case HMC:
              get_rlist_element(ctrl_lst, "int_time", ctrl.sampling.int_time,
                                6.283185307179586476925286766559005768e+00);
               break;
             case Metropolis: break;
             case Fixed_param: break;
          }
          break;

        case OPTIM:
          get_rlist_element(in, "iter", ctrl.optim.iter, 2000);
          if (get_rlist_element(in, "algorithm", t_str)) {
            if ("BFGS" == t_str)  ctrl.optim.algorithm = BFGS;
            else if ("Newton" == t_str)  ctrl.optim.algorithm = Newton;
            else if ("LBFGS" == t_str)  ctrl.optim.algorithm = LBFGS;
            else {
              std::stringstream msg;
              msg << "Invalid value for parameter algorithm (found "
                  << t_str << "; require (L)BFGS or Newton).";
              throw std::invalid_argument(msg.str());
            }
          } else {
            ctrl.optim.algorithm = LBFGS;
          }

          if (!get_rlist_element(in, "refresh", ctrl.optim.refresh)) {
            ctrl.optim.refresh = ctrl.optim.iter / 100;
            if (ctrl.optim.refresh < 1) ctrl.optim.refresh = 1;
          }

          get_rlist_element(in, "init_alpha", ctrl.optim.init_alpha, 0.001);
          get_rlist_element(in, "tol_obj", ctrl.optim.tol_obj, 1e-12);
          get_rlist_element(in, "tol_grad", ctrl.optim.tol_grad, 1e-8);
          get_rlist_element(in, "tol_param", ctrl.optim.tol_param, 1e-8);
          get_rlist_element(in, "tol_rel_obj", ctrl.optim.tol_rel_obj, 1e4);
          get_rlist_element(in, "tol_rel_grad", ctrl.optim.tol_rel_grad, 1e7);
          get_rlist_element(in, "save_iterations", ctrl.optim.save_iterations, true);
          get_rlist_element(in, "history_size", ctrl.optim.history_size, static_cast<int>(5));
          break;

        case TEST_GRADIENT:
          get_rlist_element(ctrl_lst, "epsilon", ctrl.test_grad.epsilon,1e-6);
          get_rlist_element(ctrl_lst, "error", ctrl.test_grad.error,1e-6);
          break;
      }

      if (get_rlist_element(in, "init", t_sexp)) {
        switch (TYPEOF(t_sexp)) {
          case STRSXP: init = Rcpp::as<std::string>(t_sexp); break;
          case VECSXP: init = "user"; init_list = t_sexp; break;
          default: init = "random";
        }
      } else {
        init = "random";
      }
      get_rlist_element(in, "init_r", init_radius, 2.0);
      if (0 >= init_radius)  init = "0";
      if (init == "0") init_radius = 0;
      get_rlist_element(in, "enable_random_init", enable_random_init, true);
      validate_args();
    }

    /**
     * return all the arguments used as an R list
     * @return An R list containing all the arguments for a chain.
     */
    SEXP stan_args_to_rlist() const {
      std::map<std::string, SEXP> args;
      std::map<std::string, SEXP> ctrl_args;
      std::stringstream ss;
      ss << random_seed;
      args["random_seed"] = Rcpp::wrap(ss.str());
      args["chain_id"] = Rcpp::wrap(chain_id);
      args["init"] = Rcpp::wrap(init);
      args["init_list"] = init_list;
      args["init_radius"] = Rcpp::wrap(init_radius);
      args["enable_random_init"] = Rcpp::wrap(enable_random_init);
      args["append_samples"] = Rcpp::wrap(append_samples);
      if (sample_file_flag)
        args["sample_file"] = Rcpp::wrap(sample_file);
      if (diagnostic_file_flag)
        args["diagnostic_file_flag"] = Rcpp::wrap(diagnostic_file);

      std::string sampler_t;
      switch (method) {
        case VARIATIONAL:
          args["method"] = Rcpp::wrap("variational");
          args["iter"] = Rcpp::wrap(ctrl.variational.iter);
          args["grad_samples"] = Rcpp::wrap(ctrl.variational.grad_samples);
          args["elbo_samples"] = Rcpp::wrap(ctrl.variational.elbo_samples);
          args["eval_elbo"] = Rcpp::wrap(ctrl.variational.eval_elbo);
          args["output_samples"] = Rcpp::wrap(ctrl.variational.output_samples);
          args["eta"] = Rcpp::wrap(ctrl.variational.eta);
          args["adapt_engaged"] = Rcpp::wrap(ctrl.variational.adapt_engaged);
          args["tol_rel_obj"] = Rcpp::wrap(ctrl.variational.tol_rel_obj);
          args["adapt_iter"] = Rcpp::wrap(ctrl.variational.adapt_iter);
          switch (ctrl.variational.algorithm) {
            case MEANFIELD: args["algorithm"] = Rcpp::wrap("meanfield"); break;
            case FULLRANK: args["algorithm"] = Rcpp::wrap("fullrank"); break;
          }
          break;
        case SAMPLING:
          args["method"] = Rcpp::wrap("sampling");
          args["iter"] = Rcpp::wrap(ctrl.sampling.iter);
          args["warmup"] = Rcpp::wrap(ctrl.sampling.warmup);
          args["thin"] = Rcpp::wrap(ctrl.sampling.thin);
          args["refresh"] = Rcpp::wrap(ctrl.sampling.refresh);
          args["test_grad"] = Rcpp::wrap(false);
          args["save_warmup"] = Rcpp::wrap(ctrl.sampling.save_warmup);
          ctrl_args["adapt_engaged"] = Rcpp::wrap(ctrl.sampling.adapt_engaged);
          ctrl_args["adapt_gamma"] = Rcpp::wrap(ctrl.sampling.adapt_gamma);
          ctrl_args["adapt_delta"] = Rcpp::wrap(ctrl.sampling.adapt_delta);
          ctrl_args["adapt_kappa"] = Rcpp::wrap(ctrl.sampling.adapt_kappa);
          ctrl_args["adapt_t0"] = Rcpp::wrap(ctrl.sampling.adapt_t0);
          ctrl_args["adapt_init_buffer"] = Rcpp::wrap(ctrl.sampling.adapt_init_buffer);
          ctrl_args["adapt_term_buffer"] = Rcpp::wrap(ctrl.sampling.adapt_term_buffer);
          ctrl_args["adapt_window"] = Rcpp::wrap(ctrl.sampling.adapt_window);
          ctrl_args["stepsize"] = Rcpp::wrap(ctrl.sampling.stepsize);
          ctrl_args["stepsize_jitter"] = Rcpp::wrap(ctrl.sampling.stepsize_jitter);
          switch (ctrl.sampling.algorithm) {
            case NUTS:
              ctrl_args["max_treedepth"] = Rcpp::wrap(ctrl.sampling.max_treedepth);
              sampler_t.append("NUTS");
              break;
            case HMC:
              ctrl_args["int_time"] = Rcpp::wrap(ctrl.sampling.int_time);
              sampler_t.append("HMC");
              break;
            case Metropolis:
              sampler_t.append("Metropolis");
              break;
            default: break;
          }
          if (ctrl.sampling.algorithm != Metropolis) {
            switch (ctrl.sampling.metric) {
              case UNIT_E:
                ctrl_args["metric"] = Rcpp::wrap("unit_e");
                sampler_t.append("(unit_e)");
                break;
              case DIAG_E:
                ctrl_args["metric"] = Rcpp::wrap("diag_e");
                sampler_t.append("(diag_e)");
                break;
              case DENSE_E:
                ctrl_args["metric"] = Rcpp::wrap("dense_e");
                sampler_t.append("(dense_e)");
                break;
            }
          }
          args["sampler_t"] = Rcpp::wrap(sampler_t);
          args["control"] = Rcpp::wrap(ctrl_args);
          break;
        case OPTIM:
          args["method"] = Rcpp::wrap("optim");
          args["iter"] = Rcpp::wrap(ctrl.optim.iter);
          args["refresh"] = Rcpp::wrap(ctrl.optim.refresh);
          args["save_iterations"] = Rcpp::wrap(ctrl.optim.save_iterations);
          switch (ctrl.optim.algorithm) {
            case Newton: args["algorithm"] = Rcpp::wrap("Newton"); break;
            case LBFGS: args["algorithm"] = Rcpp::wrap("LBFGS");
                        args["init_alpha"] = Rcpp::wrap(ctrl.optim.init_alpha);
                        args["tol_param"] = Rcpp::wrap(ctrl.optim.tol_param);
                        args["tol_obj"] = Rcpp::wrap(ctrl.optim.tol_obj);
                        args["tol_grad"] = Rcpp::wrap(ctrl.optim.tol_grad);
                        args["tol_rel_obj"] = Rcpp::wrap(ctrl.optim.tol_rel_obj);
                        args["tol_rel_grad"] = Rcpp::wrap(ctrl.optim.tol_rel_grad);
                        args["history_size"] = Rcpp::wrap(ctrl.optim.history_size);
                        break;
            case BFGS: args["algorithm"] = Rcpp::wrap("BFGS");
                       args["init_alpha"] = Rcpp::wrap(ctrl.optim.init_alpha);
                       args["tol_param"] = Rcpp::wrap(ctrl.optim.tol_param);
                       args["tol_obj"] = Rcpp::wrap(ctrl.optim.tol_obj);
                       args["tol_grad"] = Rcpp::wrap(ctrl.optim.tol_grad);
                       args["tol_rel_obj"] = Rcpp::wrap(ctrl.optim.tol_rel_obj);
                       args["tol_rel_grad"] = Rcpp::wrap(ctrl.optim.tol_rel_grad);
                       break;
          }
          break;
        case TEST_GRADIENT:
          args["method"] = Rcpp::wrap("test_grad");
          args["test_grad"] = Rcpp::wrap(true);
          ctrl_args["epsilon"] = Rcpp::wrap(ctrl.test_grad.epsilon);
          ctrl_args["error"] = Rcpp::wrap(ctrl.test_grad.error);
          args["control"] = Rcpp::wrap(ctrl_args);
      }
      return Rcpp::wrap(args);
    }

    inline const std::string& get_sample_file() const {
      return sample_file;
    }
    inline bool get_sample_file_flag() const {
      return sample_file_flag;
    }
    inline bool get_diagnostic_file_flag() const {
      return diagnostic_file_flag;
    }
    inline const std::string& get_diagnostic_file() const {
      return diagnostic_file;
    }

    void set_random_seed(unsigned int seed) {
      random_seed = seed;
    }

    inline unsigned int get_random_seed() const {
      return random_seed;
    }

    inline int get_ctrl_variational_grad_samples() const {
      return ctrl.variational.grad_samples;
    }
    inline int get_ctrl_variational_elbo_samples() const {
      return ctrl.variational.elbo_samples;
    }
    inline int get_ctrl_variational_output_samples() const {
      return ctrl.variational.output_samples;
    }
    inline int get_ctrl_variational_eval_elbo() const {
      return ctrl.variational.eval_elbo;
    }
    inline double get_ctrl_variational_eta() const {
      return ctrl.variational.eta;
    }
    inline bool get_ctrl_variational_adapt_engaged() const {
      return ctrl.variational.adapt_engaged;
    }
    inline double get_ctrl_variational_tol_rel_obj() const {
      return ctrl.variational.tol_rel_obj;
    }
    inline variational_algo_t get_ctrl_variational_algorithm() const {
      return ctrl.variational.algorithm;
    }
    inline int get_ctrl_variational_adapt_iter() const {
      return ctrl.variational.adapt_iter;
    }
    
    inline int get_ctrl_sampling_refresh() const {
      return ctrl.sampling.refresh;
    }
    inline sampling_metric_t get_ctrl_sampling_metric() const {
      return ctrl.sampling.metric;
    }
    inline sampling_algo_t get_ctrl_sampling_algorithm() const {
      return ctrl.sampling.algorithm;
    }
    inline int get_ctrl_sampling_warmup() const {
      return ctrl.sampling.warmup;
    }
    void set_ctrl_sampling_warmup(int n) {
      ctrl.sampling.warmup = n;
    }
    inline int get_ctrl_sampling_thin() const {
      return ctrl.sampling.thin;
    }
    inline double get_ctrl_sampling_int_time() const {
      return ctrl.sampling.int_time;
    }
    inline bool get_append_samples() const {
      return append_samples;
    }
    inline stan_args_method_t get_method() const {
      return method;
    }
    inline int get_refresh() const {
      switch (method) {
        case SAMPLING: return ctrl.sampling.refresh;
        case OPTIM: return ctrl.optim.refresh;
        case VARIATIONAL: return ctrl.variational.refresh;
        case TEST_GRADIENT: return 0;
      }
      return 0;
    }
    inline int get_iter() const {
      switch (method) {
        case SAMPLING: return ctrl.sampling.iter;
        case OPTIM: return ctrl.optim.iter;
        case VARIATIONAL: return ctrl.variational.iter;
        case TEST_GRADIENT: return 0;
      }
      return 0;
    }
    inline bool get_ctrl_sampling_adapt_engaged() const {
      return ctrl.sampling.adapt_engaged;
    }
    inline double get_ctrl_sampling_adapt_gamma() const {
      return ctrl.sampling.adapt_gamma;
    }
    inline double get_ctrl_sampling_adapt_delta() const {
      return ctrl.sampling.adapt_delta;
    }
    inline double get_ctrl_sampling_adapt_kappa() const {
      return ctrl.sampling.adapt_kappa;
    }
    inline double get_ctrl_sampling_adapt_t0() const {
      return ctrl.sampling.adapt_t0;
    }
    inline unsigned int get_ctrl_sampling_adapt_init_buffer() const {
      return ctrl.sampling.adapt_init_buffer;
    }
    inline unsigned int get_ctrl_sampling_adapt_term_buffer() const {
      return ctrl.sampling.adapt_term_buffer;
    }
    inline unsigned int get_ctrl_sampling_adapt_window() const {
      return ctrl.sampling.adapt_window;
    }
    inline double get_ctrl_sampling_stepsize() const {
       return ctrl.sampling.stepsize;
    }
    inline double get_ctrl_sampling_stepsize_jitter() const {
       return ctrl.sampling.stepsize_jitter;
    }
    inline int get_ctrl_sampling_max_treedepth() const {
       return ctrl.sampling.max_treedepth;
    }
    inline int get_ctrl_sampling_iter_save_wo_warmup() const {
       return ctrl.sampling.iter_save_wo_warmup;
    }
    inline int get_ctrl_sampling_iter_save() const {
       return ctrl.sampling.iter_save;
    }
    inline bool get_ctrl_sampling_save_warmup() const {
       return ctrl.sampling.save_warmup; // was true
    }
    inline optim_algo_t get_ctrl_optim_algorithm() const {
      return ctrl.optim.algorithm;
    }
    inline int get_ctrl_optim_refresh() const {
      return ctrl.optim.refresh;
    }
    inline bool get_ctrl_optim_save_iterations() const {
      return ctrl.optim.save_iterations;
    }
    inline double get_ctrl_optim_init_alpha() const {
      return ctrl.optim.init_alpha;
    }
    inline double get_ctrl_optim_tol_obj() const {
      return ctrl.optim.tol_obj;
    }
    inline double get_ctrl_optim_tol_grad() const {
      return ctrl.optim.tol_grad;
    }
    inline double get_ctrl_optim_tol_param() const {
      return ctrl.optim.tol_param;
    }
    inline double get_ctrl_optim_tol_rel_obj() const {
      return ctrl.optim.tol_rel_obj;
    }
    inline double get_ctrl_optim_tol_rel_grad() const {
      return ctrl.optim.tol_rel_grad;
    }
    inline int get_ctrl_optim_history_size() const {
      return ctrl.optim.history_size;
    }
    inline double get_ctrl_test_grad_epsilon() const {
      return ctrl.test_grad.epsilon;
    }
    inline double get_ctrl_test_grad_error() const {
      return ctrl.test_grad.error;
    }
    inline unsigned int get_chain_id() const {
      return chain_id;
    }
    inline double get_init_radius() const {
      return init_radius;
    }
    inline bool get_enable_random_init() const {
      return enable_random_init;
    }
    const std::string& get_init() const {
      return init;
    }
    SEXP get_init_list() const {
      return init_list;
    }

    void write_args_as_comment(std::ostream& ostream) const {
      write_comment_property(ostream,"init",init);
      write_comment_property(ostream,"enable_random_init",enable_random_init);
      write_comment_property(ostream,"seed",random_seed);
      write_comment_property(ostream,"chain_id",chain_id);
      write_comment_property(ostream,"iter",get_iter());
      switch (method) {
        case VARIATIONAL:
          write_comment_property(ostream,"grad_samples", ctrl.variational.grad_samples);
          write_comment_property(ostream,"elbo_samples", ctrl.variational.elbo_samples);
          write_comment_property(ostream,"output_samples", ctrl.variational.output_samples);
          write_comment_property(ostream,"eval_elbo", ctrl.variational.eval_elbo);
          write_comment_property(ostream,"eta", ctrl.variational.eta);
          write_comment_property(ostream,"tol_rel_obj", ctrl.variational.tol_rel_obj);
          switch (ctrl.variational.algorithm) {
            case MEANFIELD: write_comment_property(ostream,"algorithm", "meanfield"); break;
            case FULLRANK: write_comment_property(ostream,"algorithm", "fullrank"); break;
          }
          break;
        case SAMPLING:
          write_comment_property(ostream,"warmup",ctrl.sampling.warmup);
          write_comment_property(ostream,"save_warmup",ctrl.sampling.save_warmup);
          write_comment_property(ostream,"thin",ctrl.sampling.thin);
          write_comment_property(ostream,"refresh",ctrl.sampling.refresh);
          write_comment_property(ostream,"stepsize",ctrl.sampling.stepsize);
          write_comment_property(ostream,"stepsize_jitter",ctrl.sampling.stepsize_jitter);
          write_comment_property(ostream,"adapt_engaged",ctrl.sampling.adapt_engaged);
          write_comment_property(ostream,"adapt_gamma",ctrl.sampling.adapt_gamma);
          write_comment_property(ostream,"adapt_delta",ctrl.sampling.adapt_delta);
          write_comment_property(ostream,"adapt_kappa",ctrl.sampling.adapt_kappa);
          write_comment_property(ostream,"adapt_t0",ctrl.sampling.adapt_t0);
          switch (ctrl.sampling.algorithm) {
            case NUTS:
              write_comment_property(ostream,"max_treedepth",ctrl.sampling.max_treedepth);
              switch (ctrl.sampling.metric) {
                case UNIT_E: write_comment_property(ostream,"sampler_t","NUTS(unit_e)"); break;
                case DIAG_E: write_comment_property(ostream,"sampler_t","NUTS(diag_e)"); break;
                case DENSE_E: write_comment_property(ostream,"sampler_t","NUTS(dense_e)"); break;
              }
              break;
            case HMC: write_comment_property(ostream,"sampler_t", "HMC");
                      write_comment_property(ostream,"int_time", ctrl.sampling.int_time);
                      break;
            case Metropolis: write_comment_property(ostream,"sampler_t", "Metropolis"); break;
            case Fixed_param: write_comment_property(ostream, "sampler_t", "Fixed_param"); break;
            default: break;
          }
          break;

        case OPTIM:
          write_comment_property(ostream,"refresh",ctrl.optim.refresh);
          write_comment_property(ostream,"save_iterations",ctrl.optim.save_iterations);
          switch (ctrl.optim.algorithm) {
            case Newton: write_comment_property(ostream,"algorithm", "Newton"); break;
            case BFGS: write_comment_property(ostream,"algorithm", "BFGS");
                       write_comment_property(ostream,"init_alpha", ctrl.optim.init_alpha);
                       write_comment_property(ostream,"tol_obj", ctrl.optim.tol_obj);
                       write_comment_property(ostream,"tol_grad", ctrl.optim.tol_grad);
                       write_comment_property(ostream,"tol_param", ctrl.optim.tol_param);
                       write_comment_property(ostream,"tol_rel_obj", ctrl.optim.tol_rel_obj);
                       write_comment_property(ostream,"tol_rel_grad", ctrl.optim.tol_rel_grad);
                       break;
            case LBFGS: write_comment_property(ostream,"algorithm", "LBFGS");
                       write_comment_property(ostream,"init_alpha", ctrl.optim.init_alpha);
                       write_comment_property(ostream,"tol_obj", ctrl.optim.tol_obj);
                       write_comment_property(ostream,"tol_grad", ctrl.optim.tol_grad);
                       write_comment_property(ostream,"tol_param", ctrl.optim.tol_param);
                       write_comment_property(ostream,"tol_rel_obj", ctrl.optim.tol_rel_obj);
                       write_comment_property(ostream,"tol_rel_grad", ctrl.optim.tol_rel_grad);
                       write_comment_property(ostream,"history_size", ctrl.optim.history_size);
                       break;
          }
        case TEST_GRADIENT: break;
      }
      if (sample_file_flag)
        write_comment_property(ostream,"sample_file",sample_file);
      if (diagnostic_file_flag)
        write_comment_property(ostream,"diagnostic_file",diagnostic_file);
      write_comment_property(ostream,"append_samples",append_samples);
      write_comment(ostream);
    }
  };
}

#endif

