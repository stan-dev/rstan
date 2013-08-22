
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
    inline unsigned int sexp2seed(SEXP seed) { 
      if (TYPEOF(seed) == STRSXP)  
        return boost::lexical_cast<unsigned int>(Rcpp::as<std::string>(seed));
      return Rcpp::as<unsigned int>(seed); 
    }

    void write_comment(std::ostream& o) {
      o << "#" << std::endl;
    }
  
    template <typename M>
    void write_comment(std::ostream& o,
                       const M& msg) {
      o << "# " << msg << std::endl;
    }
  
    template <typename K, typename V>
    void write_comment_property(std::ostream& o,
                                const K& key,
                                const V& val) {
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
    std::string random_seed_src; // "user" or "default" 
    unsigned int chain_id; 
    std::string chain_id_src; // "user" or "default" 
    std::string init; 
    SEXP init_list;  
    std::string sample_file; // the file for outputting the samples
    bool append_samples; 
    bool sample_file_flag; // true: write out to a file; false, do not 
    stan_args_method_t method; 
    std::string diagnostic_file; 
    bool diagnostic_file_flag; // 
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
        int iter;
        int refresh;
        optim_algo_t algorithm;
        bool save_iterations;
        double stepsize;
      } optim; 
      struct {
        int refresh2;
      } test_grad;
    } ctrl; 

  public:
    stan_args(const Rcpp::List& in) : init_list(R_NilValue) {
      std::vector<std::string> args_names 
        = Rcpp::as<std::vector<std::string> >(in.names()); 

      size_t idx = find_index(args_names, std::string("method")); 
      if (idx == args_names.size()) method = SAMPLING;
      else {
        std::string t = Rcpp::as<std::string>(in[idx]);
        if ("sampling" == t)  method = SAMPLING;
        else if ("optim" == t)  method = OPTIM;
        else if ("test_grad" == t)  method = TEST_GRADIENT;
        else method = SAMPLING;
      } 

      if (method == SAMPLING) {
        size_t idx = find_index(args_names, std::string("sample_file")); 
        if (idx == args_names.size()) sample_file_flag = false; 
        else {
          sample_file = Rcpp::as<std::string>(in[idx]); 
          sample_file_flag = true; 
        }

        idx = find_index(args_names, std::string("diagnostic_file"));
        if (idx == args_names.size()) diagnostic_file_flag = false;
        else {
          diagnostic_file = Rcpp::as<std::string>(in[idx]);
          diagnostic_file_flag = true;
        } 

        idx = find_index(args_names, std::string("iter")); 
        if (idx == args_names.size()) ctrl.sampling.iter = 2000;  
        else ctrl.sampling.iter = Rcpp::as<int>(in[idx]); 
  
        idx = find_index(args_names, std::string("warmup")); 
        if (idx == args_names.size()) ctrl.sampling.warmup = ctrl.sampling.iter / 2; 
        else ctrl.sampling.warmup = Rcpp::as<int>(in[idx]); 

        idx = find_index(args_names, std::string("thin")); 
        int calculated_thin = (ctrl.sampling.iter - ctrl.sampling.warmup) / 1000;
        if (idx == args_names.size()) 
          ctrl.sampling.thin = (calculated_thin > 1) ? calculated_thin : 1;
        else ctrl.sampling.thin = Rcpp::as<int>(in[idx]); 

        ctrl.sampling.iter_save_wo_warmup 
          = 1 + (ctrl.sampling.iter - ctrl.sampling.warmup - 1) / ctrl.sampling.thin; 
        ctrl.sampling.iter_save 
          = ctrl.sampling.iter_save_wo_warmup
            + 1 + (ctrl.sampling.warmup - 1) / ctrl.sampling.thin;

        ctrl.sampling.refresh = 1;
        idx = find_index(args_names, std::string("refresh"));
        if (idx == args_names.size()) {
          if (ctrl.sampling.iter >= 20) ctrl.sampling.refresh = ctrl.sampling.iter / 10; 
        } else ctrl.sampling.refresh = Rcpp::as<int>(in[idx]);

        idx = find_index(args_names, std::string("seed")); 
        if (idx == args_names.size()) {
          random_seed = std::time(0); 
          random_seed_src = "random"; 
        } else {
          random_seed = sexp2seed(in[idx]);
          random_seed_src = "user or from R"; 
        }

        idx = find_index(args_names, std::string("algorithm"));
        if (idx == args_names.size()) { 
          ctrl.sampling.algorithm = NUTS;
        } else {
          std::string algo_str = Rcpp::as<std::string>(in[idx]);
          if (algo_str == "HMC") ctrl.sampling.algorithm = HMC;
          else if (algo_str == "NUTS") ctrl.sampling.algorithm = NUTS;
          else if (algo_str == "Metropolis") ctrl.sampling.algorithm = Metropolis;
        } 

        idx = find_index(args_names, std::string("control"));
        if (idx == args_names.size()) { 
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
        } else { 
          Rcpp::List ctrl_lst(static_cast<SEXP>(in[idx]));
          std::vector<std::string> ctrl_names 
            = Rcpp::as<std::vector<std::string> >(ctrl_lst.names());
          size_t idx2 = find_index(ctrl_names, std::string("adapt_engaged"));
          if (idx2 == ctrl_names.size()) {
            ctrl.sampling.adapt_engaged = true;
          } else {
            ctrl.sampling.adapt_engaged = Rcpp::as<bool>(ctrl_lst[idx2]);
          } 
          idx2 = find_index(ctrl_names, std::string("adapt_gamma"));
          if (idx2 == ctrl_names.size()) ctrl.sampling.adapt_gamma = 0.05;
          else ctrl.sampling.adapt_gamma = Rcpp::as<double>(ctrl_lst[idx2]);
          idx2 = find_index(ctrl_names, std::string("adapt_delta"));
          if (idx2 == ctrl_names.size()) ctrl.sampling.adapt_delta = 0.65;
          else ctrl.sampling.adapt_delta = Rcpp::as<double>(ctrl_lst[idx2]);
          idx2 = find_index(ctrl_names, std::string("adapt_kappa"));
          if (idx2 == ctrl_names.size()) ctrl.sampling.adapt_kappa = 0.75;
          else ctrl.sampling.adapt_kappa = Rcpp::as<double>(ctrl_lst[idx2]);
          idx2 = find_index(ctrl_names, std::string("adapt_t0"));
          if (idx2 == ctrl_names.size()) ctrl.sampling.adapt_t0 = 10;
          else ctrl.sampling.adapt_t0 = Rcpp::as<double>(ctrl_lst[idx2]);
  
          switch (ctrl.sampling.algorithm) { 
            case NUTS: 
              idx2 = find_index(ctrl_names, "max_treedepth");
              if (idx2 == ctrl_names.size()) ctrl.sampling.max_treedepth = 10;
              else ctrl.sampling.max_treedepth = Rcpp::as<int>(ctrl_lst[idx2]);
              idx2 = find_index(ctrl_names, "metric");
              if (idx2 != ctrl_names.size()) {
                std::string s1 = Rcpp::as<std::string>(ctrl_lst[idx2]);
                if ("unit_e" == s1) ctrl.sampling.metric = UNIT_E;
                else if ("diag_e" == s1) ctrl.sampling.metric = DIAG_E;
                else if ("dense_e" == s1) ctrl.sampling.metric = DENSE_E;
              } else ctrl.sampling.metric = DIAG_E;
              idx2 = find_index(ctrl_names, "stepsize");
              if (idx2 == ctrl_names.size()) ctrl.sampling.stepsize = 1; 
              else ctrl.sampling.stepsize = Rcpp::as<double>(ctrl_lst[idx2]);
              idx2 = find_index(ctrl_names, "stepsize_jitter");
              if (idx2 == ctrl_names.size()) ctrl.sampling.stepsize_jitter = 0; 
              else ctrl.sampling.stepsize_jitter = Rcpp::as<double>(ctrl_lst[idx2]);
              break;
            case HMC: 
              idx2 = find_index(ctrl_names, "int_time");
              if (idx2 == ctrl_names.size()) 
                ctrl.sampling.int_time = 6.283185307179586476925286766559005768e+00;
              else ctrl.sampling.int_time = Rcpp::as<double>(ctrl_lst[idx2]);
              break;
            case Metropolis: break;
          }
        } 
      } else if (method == OPTIM) {
        idx = find_index(args_names, "algorithm");
        if (idx == args_names.size()) {
          ctrl.optim.algorithm = BFGS;
        } else {
          std::string t = Rcpp::as<std::string>(in[idx]);
          if ("BFGS" == t)  ctrl.optim.algorithm = BFGS;
          else if ("Newton" == t)  ctrl.optim.algorithm = Newton;
          else if ("Nesterov" == t)  ctrl.optim.algorithm = Nesterov;
        } 
        idx = find_index(args_names, "refresh");
        if (idx == args_names.size())
          ctrl.optim.refresh = ctrl.optim.iter / 100; 
        else 
          ctrl.optim.refresh = Rcpp::as<int>(in[idx]);
        idx = find_index(args_names, "stepsize");
        if (idx == args_names.size()) ctrl.optim.stepsize = 1;
        else ctrl.optim.stepsize = Rcpp::as<double>(in[idx]);
        idx = find_index(args_names, "save_iterations");
        if (idx == args_names.size()) ctrl.optim.save_iterations = true;
        else ctrl.optim.save_iterations = Rcpp::as<bool>(in[idx]);
      } else if (method == TEST_GRADIENT) {
      } 

      idx = find_index(args_names, std::string("chain_id")); 
      if (idx == args_names.size()) { 
        chain_id = 1; 
        chain_id_src = "default"; 
      } else {
        chain_id = Rcpp::as<unsigned int>(in[idx]); 
        chain_id_src = "user"; 
      }
      
      idx = find_index(args_names, std::string("init")); 
      if (idx == args_names.size()) {
        init = "random"; 
      } else {
        switch (TYPEOF(in[idx])) {
          case STRSXP: init = Rcpp::as<std::string>(in[idx]); break; 
          case VECSXP: init = "user"; init_list = in[idx]; break; 
          default: init = "random"; 
        } 
      }
      idx = find_index(args_names, std::string("append_samples")); 
      if (idx == args_names.size()) append_samples = false; 
      else append_samples = Rcpp::as<bool>(in[idx]); 
    } 

    /**
     * return all the arguments used as an R list
     * @return An R list containing all the arguments for a chain. 
     */ 
    SEXP stan_args_to_rlist() const {
      Rcpp::List lst;

      std::stringstream ss; 
      ss << random_seed; 
      lst["random_seed"] = ss.str();
      lst["chain_id"] = chain_id;
      lst["init"] = init;
      lst["init_list"] = init_list;

      if (method == SAMPLING) { // sampling 
        lst["method"] = "SAMPLING";
        Rcpp::List ctrl_list;
        ctrl_list["adapt_engaged"] = ctrl.sampling.adapt_engaged;
        ctrl_list["adapt_gamma"] = ctrl.sampling.adapt_gamma;
        ctrl_list["adapt_delta"] = ctrl.sampling.adapt_delta;
        ctrl_list["adapt_kappa"] = ctrl.sampling.adapt_kappa;
        ctrl_list["adapt_t0"] = ctrl.sampling.adapt_t0;
        if (ctrl.sampling.algorithm == 1) {
          switch (ctrl.sampling.metric) { 
            case UNIT_E : ctrl_list["metric"] = Rcpp::wrap("unit_e"); break;
            case DIAG_E : ctrl_list["metric"] = Rcpp::wrap("diag_e"); break;
            case DENSE_E : ctrl_list["metric"] = Rcpp::wrap("dense_e"); break;
          }
          ctrl_list["stepsize"] = ctrl.sampling.stepsize;
          ctrl_list["stepsize_jitter"] = ctrl.sampling.stepsize_jitter;
          ctrl_list["max_treedepth"] = ctrl.sampling.max_treedepth;
        } else if (ctrl.sampling.algorithm == 2) {
          ctrl_list["int_time"] = ctrl.sampling.int_time;
        } 

        lst["control"] = ctrl_list;
      } else if (method == OPTIM)  { //  
        lst["method"] = "OPTIM";
      } else if (method == TEST_GRADIENT)  { //  
        lst["method"] = "TEST_GRADIENT";
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

    // 1: optimization (point_estimate), 2: test_grad, 3: sampling 
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
    inline bool get_ctrl_optim_stepsize() const { 
      return ctrl.optim.stepsize;
    }
    inline unsigned int get_chain_id() const {
      return chain_id;
    } 
    const std::string& get_init() const {
      return init;
    } 
    SEXP get_init_list() const {
      return init_list; 
    } 
    void write_args_as_comment(std::ostream& ostream) const { 
     // TODO
    } 
    
/*
    const std::string& get_random_seed_src() const {
      return random_seed_src; 
    } 
    const std::string& get_chain_id_src() const {
      return chain_id_src; 
    } 

    int get_leapfrog_steps() const {
      return leapfrog_steps; 
    } 
    double get_epsilon() const {
      return epsilon; 
    } 
    double get_epsilon_pm() const {
      return epsilon_pm; 
    } 
    double get_delta() const {  
      return delta;
    } 
    double get_gamma() const { 
      return gamma;
    } 
    bool get_test_grad() const {
      return test_grad; 
    } 
    inline int get_point_estimate() const {
      return point_estimate;
    }
    unsigned int get_chain_id() const {
      return chain_id; 
    } 
    bool get_equal_step_sizes() const {
      return equal_step_sizes; 
    } 
    bool get_nondiag_mass() const {
      return nondiag_mass;
    } 
    void write_args_as_comment(std::ostream& ostream) const { 
      write_comment_property(ostream,"init",init);
      write_comment_property(ostream,"append_samples",append_samples);
      write_comment_property(ostream,"seed",random_seed);
      write_comment_property(ostream,"chain_id",chain_id);
      write_comment_property(ostream,"chain_id_src",chain_id_src);
      write_comment_property(ostream,"iter",iter); 
      write_comment_property(ostream,"warmup",warmup);
      write_comment_property(ostream,"save_warmup",1);
      write_comment_property(ostream,"thin",thin);
      write_comment_property(ostream,"leapfrog_steps",leapfrog_steps);
      write_comment_property(ostream,"max_treedepth",max_treedepth);
      write_comment_property(ostream,"epsilon",epsilon);
      write_comment_property(ostream,"equal_step_sizes",equal_step_sizes); 
      write_comment_property(ostream,"epsilon_pm",epsilon_pm);
      write_comment_property(ostream,"delta",delta);
      write_comment_property(ostream,"gamma",gamma);
      write_comment_property(ostream,"nondiag_mass",nondiag_mass);
      if (sample_file_flag) 
        write_comment_property(ostream,"sample_file",sample_file);
      if (diagnostic_file_flag)
        write_comment_property(ostream,"diagnostic_file",diagnostic_file);
      write_comment(ostream);
    }
*/
  }; 
} 

#endif 

