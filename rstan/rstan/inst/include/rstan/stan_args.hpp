
#ifndef __RSTAN__STAN_ARGS_HPP__
#define __RSTAN__STAN_ARGS_HPP__


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

  enum sampling_algo_t { NUTS = 1, HMC = 2, Metropolis = 3};
  enum optim_algo_t { Newton = 1, Nesterov = 2, BFGS = 3};
  enum sampling_metric_t { UNIT_E = 1, DIAG_E = 2, DENSE_E = 3};
  enum stan_args_method_t { SAMPLING = 1, OPTIM = 2, TEST_GRADIENT = 3};

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
        bool save_warmup; // weather to save warmup samples (always true now)
        int iter_save; // number of iterations saved 
        int iter_save_wo_warmup; // number of iterations saved wo warmup
        bool adapt_engaged; 
        double adapt_gamma;
        double adapt_delta;
        double adapt_kappa;
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
        optim_algo_t algorithm; // Newton, Nesterov, BFGS
        bool save_iterations; // default to false
        double stepsize; // default to 1, for Nesterov
        double init_alpha; // default to 0.001, for BFGS
        double tol_obj; // default to 1e-8, for BFGS
        double tol_grad; // default to 1e-8, for BFGS
        double tol_param; // default to 1e-8, for BFGS
      } optim; 
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
            msg << "Invalid adaptation parameter (found detal="
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
          if (ctrl.sampling.max_treedepth < 0) {
            std::stringstream msg;
            msg << "Invalid adaptation parameter (found max_treedepth="
                << ctrl.sampling.max_treedepth << "; require max_treedepth>0).";
            throw std::invalid_argument(msg.str());
          } 
          if (ctrl.sampling.int_time < 0) {
            std::stringstream msg;
            msg << "Invalid adaptation parameter (found int_time="
                << ctrl.sampling.int_time << "; require int_time>0).";
            throw std::invalid_argument(msg.str());
          } 
          break;
        case OPTIM:
          if (ctrl.optim.stepsize < 0) {
            std::stringstream msg; 
            msg << "Invalid adaptation parameter (found stepsize="
                << ctrl.optim.stepsize << "; require stepsize > 0).";
            throw std::invalid_argument(msg.str());
          } 
          if (ctrl.optim.init_alpha < 0) {  
            std::stringstream msg; 
            msg << "Invalid adaptation parameter (found init_alpha="
                << ctrl.optim.init_alpha << "; require init_alpha > 0).";
            throw std::invalid_argument(msg.str());
          } 
          break;
        case TEST_GRADIENT: break;
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
        else method = SAMPLING;
      } 

      sample_file_flag = get_rlist_element(in, "sample_file", sample_file);
      diagnostic_file_flag = get_rlist_element(in, "diagnostic_file", diagnostic_file);

      int calculated_thin;
      switch (method) { 
        case SAMPLING: 
          get_rlist_element(in, "iter", ctrl.sampling.iter, 2000);
          get_rlist_element(in, "warmup", ctrl.sampling.warmup, ctrl.sampling.iter / 2);
   
          calculated_thin = (ctrl.sampling.iter - ctrl.sampling.warmup) / 1000;
          if (calculated_thin < 1) calculated_thin = 1;
          get_rlist_element(in, "thin", ctrl.sampling.thin, calculated_thin);
  
          ctrl.sampling.iter_save_wo_warmup 
            = 1 + (ctrl.sampling.iter - ctrl.sampling.warmup - 1) / ctrl.sampling.thin; 
          ctrl.sampling.iter_save 
            = ctrl.sampling.iter_save_wo_warmup
              + 1 + (ctrl.sampling.warmup - 1) / ctrl.sampling.thin;
  
          ctrl.sampling.refresh = (ctrl.sampling.iter >= 20) ? 
                                  ctrl.sampling.iter / 10 : 1; 
          get_rlist_element(in, "refresh", ctrl.sampling.refresh);
         
          b = get_rlist_element(in, "seed", t_sexp);
          if (b) random_seed = sexp2seed(t_sexp);
          else random_seed = std::time(0);
  
          if (get_rlist_element(in, "algorithm", t_str)) {
            if (t_str == "HMC") ctrl.sampling.algorithm = HMC;
            else if (t_str == "Metropolis") ctrl.sampling.algorithm = Metropolis;
            else if (t_str == "NUTS") ctrl.sampling.algorithm = NUTS;
            else {
              std::stringstream msg;
              msg << "Invalid value for parameter algorithm (found "
                  << t_str << "; require HMC, Metropolis, or NUTS).";
              throw std::invalid_argument(msg.str());
            } 
          } else {
            ctrl.sampling.algorithm = NUTS;
          }
  
          if (get_rlist_element(in, "control", t_sexp)) {
            Rcpp::List ctrl_lst(t_sexp);
            get_rlist_element(ctrl_lst, "adapt_engaged", ctrl.sampling.adapt_engaged, true);
            get_rlist_element(ctrl_lst, "adapt_gamma", ctrl.sampling.adapt_gamma, 0.05);
            get_rlist_element(ctrl_lst, "adapt_delta", ctrl.sampling.adapt_delta, 0.65);
            get_rlist_element(ctrl_lst, "adapt_kappa", ctrl.sampling.adapt_kappa, 0.75);
            get_rlist_element(ctrl_lst, "adapt_t0", ctrl.sampling.adapt_t0, 10.0);
            get_rlist_element(ctrl_lst, "stepsize", ctrl.sampling.stepsize, 1.0);
            get_rlist_element(ctrl_lst, "stepsize_jitter", ctrl.sampling.stepsize_jitter, 0.0);
    
            switch (ctrl.sampling.algorithm) { 
              case NUTS: 
                get_rlist_element(ctrl_lst, "max_treedepth", ctrl.sampling.max_treedepth, 10);
                if (get_rlist_element(ctrl_lst, "metric", t_str)) { 
                  if ("unit_e" == t_str) ctrl.sampling.metric = UNIT_E;
                  else if ("diag_e" == t_str) ctrl.sampling.metric = DIAG_E;
                  else if ("dense_e" == t_str) ctrl.sampling.metric = DENSE_E;
                } else ctrl.sampling.metric = DIAG_E;
                break;
              case HMC: 
                get_rlist_element(ctrl_lst, "int_time", ctrl.sampling.int_time, 
                                  6.283185307179586476925286766559005768e+00);
                break;
              case Metropolis: break;
            }
          } else { 
            ctrl.sampling.adapt_engaged = true;
            ctrl.sampling.adapt_gamma = 0.05;
            ctrl.sampling.adapt_delta = 0.65;
            ctrl.sampling.adapt_kappa = 0.75;
            ctrl.sampling.adapt_t0  = 10;
            ctrl.sampling.max_treedepth = 10;
            ctrl.sampling.metric = DIAG_E;
            ctrl.sampling.stepsize = 1;
            ctrl.sampling.stepsize_jitter = 0;
            ctrl.sampling.int_time = 6.283185307179586476925286766559005768e+00;
          }
          break;

        case OPTIM: 
          get_rlist_element(in, "iter", ctrl.optim.iter, 2000);
          if (get_rlist_element(in, "algorithm", t_str)) {
            if ("BFGS" == t_str)  ctrl.optim.algorithm = BFGS;
            else if ("Newton" == t_str)  ctrl.optim.algorithm = Newton;
            else if ("Nesterov" == t_str)  ctrl.optim.algorithm = Nesterov;
            else {
              std::stringstream msg;
              msg << "Invalid value for parameter algorithm (found "
                  << t_str << "; require BFGS, Newton, or Nesterov).";
              throw std::invalid_argument(msg.str());
            }
          } else {
            ctrl.optim.algorithm = BFGS;
          } 
          
          if (!get_rlist_element(in, "refresh", ctrl.optim.refresh)) {
            ctrl.optim.refresh = ctrl.optim.iter / 100; 
            if (ctrl.optim.refresh < 1) ctrl.optim.refresh = 1;
          } 
  
          get_rlist_element(in, "stepsize", ctrl.optim.stepsize, 1.0);
          get_rlist_element(in, "init_alpha", ctrl.optim.init_alpha, 0.001);
          get_rlist_element(in, "tol_obj", ctrl.optim.tol_obj, 1e-8);
          get_rlist_element(in, "tol_grad", ctrl.optim.tol_grad, 1e-8);
          get_rlist_element(in, "tol_param", ctrl.optim.tol_param, 1e-8);
          get_rlist_element(in, "save_iterations", ctrl.optim.save_iterations, true);
          break;

        case TEST_GRADIENT: break;
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
      validate_args();
    } 

    /**
     * return all the arguments used as an R list
     * @return An R list containing all the arguments for a chain. 
     */ 
    SEXP stan_args_to_rlist() const {
      Rcpp::List lst;
      Rcpp::List ctrl_list;

      std::stringstream ss; 
      ss << random_seed; 
      lst["random_seed"] = ss.str();
      lst["chain_id"] = chain_id;
      lst["init"] = init;
      lst["init_list"] = init_list;
      lst["init_radius"] = init_radius;
      lst["append_samples"] = append_samples;
      if (sample_file_flag) 
        lst["sample_file"] = sample_file;
      if (diagnostic_file_flag) 
        lst["diagnostic_file_flag"] = diagnostic_file;

      switch (method) { 
        case SAMPLING: 
          lst["method"] = "sampling";
          lst["iter"] = ctrl.sampling.iter;
          lst["warmup"] = ctrl.sampling.warmup;
          lst["thin"] = ctrl.sampling.thin;
          lst["refresh"] = ctrl.sampling.refresh;
          lst["test_grad"] = false;
          ctrl_list["adapt_engaged"] = ctrl.sampling.adapt_engaged;
          ctrl_list["adapt_gamma"] = ctrl.sampling.adapt_gamma;
          ctrl_list["adapt_delta"] = ctrl.sampling.adapt_delta;
          ctrl_list["adapt_kappa"] = ctrl.sampling.adapt_kappa;
          ctrl_list["adapt_t0"] = ctrl.sampling.adapt_t0;
          switch (ctrl.sampling.algorithm) {  
            case NUTS:
              switch (ctrl.sampling.metric) { 
                case UNIT_E: 
                  ctrl_list["metric"] = Rcpp::wrap("unit_e"); 
                  lst["sampler_t"] = "NUTS(unit_e)";
                  break;
                case DIAG_E: 
                  ctrl_list["metric"] = Rcpp::wrap("diag_e");
                  lst["sampler_t"] = "NUTS(diag_e)";
                  break;
                case DENSE_E: 
                  ctrl_list["metric"] = Rcpp::wrap("dense_e");
                  lst["sampler_t"] = "NUTS(dense_e)";
                  break;
              }
              ctrl_list["stepsize"] = ctrl.sampling.stepsize;
              ctrl_list["stepsize_jitter"] = ctrl.sampling.stepsize_jitter;
              ctrl_list["max_treedepth"] = ctrl.sampling.max_treedepth;
              break;
            case HMC: 
              ctrl_list["int_time"] = ctrl.sampling.int_time;
              lst["sampler_t"] = "HMC";
              break;
            case Metropolis: 
              lst["sampler_t"] = "Metropolis";
              break;
          } 
          lst["control"] = ctrl_list;
          break;
        case OPTIM: 
          lst["method"] = "optim";
          lst["iter"] = ctrl.optim.iter;
          lst["refresh"] = ctrl.optim.refresh;
          lst["save_iterations"] = ctrl.optim.save_iterations;
          switch (ctrl.optim.algorithm) {
            case Newton: lst["algorithm"] = "Newton"; break;
            case Nesterov: lst["algorithm"] = "Nesterov"; 
                           lst["stepsize"] = ctrl.optim.stepsize;
                           break;
            case BFGS: lst["algorithm"] = "BFGS"; 
                       lst["init_alpha"] = ctrl.optim.init_alpha;
                       lst["tol_obj"] = ctrl.optim.tol_obj;
                       lst["tol_grad"] = ctrl.optim.tol_grad;
                       lst["tol_param"] = ctrl.optim.tol_param;
                       break;
          } 
          break;
        case TEST_GRADIENT:
          lst["method"] = "test_grad";
          lst["test_grad"] = true;
      } 
      return lst;
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

    inline int get_ctrl_sampling_refresh() const { 
      return ctrl.sampling.refresh; 
    } 
    const inline sampling_metric_t get_ctrl_sampling_metric() const { 
      return ctrl.sampling.metric;
    } 
    const inline sampling_algo_t get_ctrl_sampling_algorithm() const {
      return ctrl.sampling.algorithm;
    }
    inline int get_ctrl_sampling_warmup() const { 
      return ctrl.sampling.warmup;
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
    inline int get_iter() const {
      switch (method) {
        case SAMPLING: return ctrl.sampling.iter;
        case OPTIM: return ctrl.optim.iter;
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
       return true;
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
    inline double get_ctrl_optim_stepsize() const { 
      return ctrl.optim.stepsize;
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
    inline unsigned int get_chain_id() const {
      return chain_id;
    } 
    inline double get_init_radius() const {
      return init_radius;
    } 
    const std::string& get_init() const {
      return init;
    } 
    SEXP get_init_list() const {
      return init_list; 
    } 
    
    void write_args_as_comment(std::ostream& ostream) const { 
      write_comment_property(ostream,"init",init);
      write_comment_property(ostream,"seed",random_seed);
      write_comment_property(ostream,"chain_id",chain_id);
      write_comment_property(ostream,"iter",get_iter()); 
      switch (method) {
        case SAMPLING: 
          write_comment_property(ostream,"warmup",ctrl.sampling.warmup);
          write_comment_property(ostream,"save_warmup",1);
          write_comment_property(ostream,"thin",ctrl.sampling.thin);
          write_comment_property(ostream,"refresh",ctrl.sampling.refresh);
          write_comment_property(ostream,"max_treedepth",ctrl.sampling.max_treedepth);
          write_comment_property(ostream,"stepsize",ctrl.sampling.stepsize);
          write_comment_property(ostream,"stepsize_jitter",ctrl.sampling.stepsize_jitter);
          write_comment_property(ostream,"adapt_engaged",ctrl.sampling.adapt_engaged);
          write_comment_property(ostream,"adapt_gamma",ctrl.sampling.adapt_gamma);
          write_comment_property(ostream,"adapt_delta",ctrl.sampling.adapt_delta);
          write_comment_property(ostream,"adapt_kappa",ctrl.sampling.adapt_kappa);
          write_comment_property(ostream,"adapt_t0",ctrl.sampling.adapt_t0);
          switch (ctrl.sampling.algorithm) {
            case NUTS: 
              switch (ctrl.sampling.metric) {
                case UNIT_E: write_comment_property(ostream,"sampler_t","NUTS(unit_e)"); break;
                case DIAG_E: write_comment_property(ostream,"sampler_t","NUTS(diag_e)"); break;
                case DENSE_E: write_comment_property(ostream,"sampler_t","NUTS(dense_e)"); break;
              } 
              break;
            case HMC: write_comment_property(ostream,"sampler_t", "HMC"); break;
            case Metropolis: write_comment_property(ostream,"sampler_t", "Metropolis"); break;
          } 
          break;

        case OPTIM: 
          write_comment_property(ostream,"refresh",ctrl.optim.refresh);
          write_comment_property(ostream,"save_iterations",ctrl.optim.save_iterations);
          switch (ctrl.optim.algorithm) {
            case Newton: write_comment_property(ostream,"algorithm", "Newton"); break;
            case Nesterov: write_comment_property(ostream,"algorithm", "Nesterov");
                           write_comment_property(ostream,"stepsize", ctrl.optim.stepsize);
                           break;
            case BFGS: write_comment_property(ostream,"algorithm", "BFGS");
                       write_comment_property(ostream,"init_alpha", ctrl.optim.init_alpha);
                       write_comment_property(ostream,"tol_obj", ctrl.optim.tol_obj);
                       write_comment_property(ostream,"tol_grad", ctrl.optim.tol_grad);
                       write_comment_property(ostream,"tol_param", ctrl.optim.tol_param);
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

