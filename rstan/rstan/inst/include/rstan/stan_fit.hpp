
#ifndef __RSTAN__STAN_FIT_HPP__
#define __RSTAN__STAN_FIT_HPP__

#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <stan/version.hpp>
// #include <stan/io/cmd_line.hpp>
// #include <stan/io/dump.hpp>

#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/random/additive_combine.hpp> // L'Ecuyer RNG
#include <boost/random/uniform_real_distribution.hpp>

#include <stan/model/util.hpp>

#include <stan/optimization/newton.hpp>
#include <stan/optimization/nesterov_gradient.hpp>
#include <stan/optimization/bfgs.hpp>

#include <stan/mcmc/hmc/static/adapt_unit_e_static_hmc.hpp>
#include <stan/mcmc/hmc/static/adapt_diag_e_static_hmc.hpp>
#include <stan/mcmc/hmc/static/adapt_dense_e_static_hmc.hpp>
#include <stan/mcmc/hmc/nuts/adapt_unit_e_nuts.hpp>
#include <stan/mcmc/hmc/nuts/adapt_diag_e_nuts.hpp>
#include <stan/mcmc/hmc/nuts/adapt_dense_e_nuts.hpp>

#include <stan/common/do_print.hpp>
#include <stan/common/print_progress.hpp>


#include <rstan/io/rlist_ref_var_context.hpp> 
#include <rstan/io/r_ostream.hpp> 
#include <rstan/stan_args.hpp> 
#include <rstan/mcmc_output.hpp>
// #include <stan/mcmc/chains.hpp>
#include <Rcpp.h>
// #include <Rinternals.h>

//http://cran.r-project.org/doc/manuals/R-exts.html#Allowing-interrupts
#include <R_ext/Utils.h>
// void R_CheckUserInterrupt(void);


// REF: stan/gm/command.hpp 

namespace rstan {

  namespace {
    /**
     *@tparam T The type by which we use for dimensions. T could be say size_t
     * or unsigned int. This whole business (not using size_t) is due to that
     * Rcpp::wrap/as does not support size_t on some platforms and R could not
     * deal with 64bits integers. 
     *
     */ 
    template <class T> 
    size_t calc_num_params(const std::vector<T>& dim) {
      T num_params = 1;
      for (size_t i = 0;  i < dim.size(); ++i)
        num_params *= dim[i];
      return num_params;
    }

    template <class T> 
    void calc_starts(const std::vector<std::vector<T> >& dims,
                     std::vector<T>& starts) { 
      starts.resize(0); 
      starts.push_back(0); 
      for (size_t i = 1; i < dims.size(); ++i)
        starts.push_back(starts[i - 1] + calc_num_params(dims[i - 1]));
    }

    template <class T> 
    T calc_total_num_params(const std::vector<std::vector<T> >& dims) {
      T num_params = 0;
      for (size_t i = 0; i < dims.size(); ++i)
        num_params += calc_num_params(dims[i]);
      return num_params;
    }

    /**
     *  Get the parameter indexes for a vector(array) parameter.
     *  For example, we have parameter beta, which has 
     *  dimension [2,3]. Then this function gets 
     *  the indexes as (if col_major = false)
     *  [0,0], [0,1], [0,2] 
     *  [1,0], [1,1], [1,2] 
     *  or (if col_major = true) 
     *  [0,0], [1,0] 
     *  [0,1], [1,1] 
     *  [0,2], [121] 
     *
     *  @param dim[in] the dimension of parameter
     *  @param idx[out] for keeping all the indexes
     *
     *  <p> when idx is empty (size = 0), idx 
     *  would contains an empty vector. 
     * 
     *
     */
    
    template <class T>
    void expand_indices(std::vector<T> dim,
                        std::vector<std::vector<T> >& idx,
                        bool col_major = false) {
      size_t len = dim.size();
      idx.resize(0);
      size_t total = calc_num_params(dim);
      if (0 >= total) return;
      std::vector<size_t> loopj;
      for (size_t i = 1; i <= len; ++i)
        loopj.push_back(len - i);
    
      if (col_major)
        for (size_t i = 0; i < len; ++i)
          loopj[i] = len - 1 - loopj[i];
    
      idx.push_back(std::vector<T>(len, 0));
      for (size_t i = 1; i < total; i++) {
        std::vector<T>  v(idx.back());
        for (size_t j = 0; j < len; ++j) {
          size_t k = loopj[j];
          if (v[k] < dim[k] - 1) {
            v[k] += 1;
            break;
          }
          v[k] = 0;
        }
        idx.push_back(v);
      }
    }

    /**
     * Get the names for an array of given dimensions 
     * in the way of column majored. 
     * For example, if we know an array named `a`, with
     * dimensions of [2, 3, 4], the names then are (starting
     * from 0):
     * a[0,0,0]
     * a[1,0,0]
     * a[0,1,0]
     * a[1,1,0]
     * a[0,2,0]
     * a[1,2,0]
     * a[0,0,1]
     * a[1,0,1]
     * a[0,1,1]
     * a[1,1,1]
     * a[0,2,1]
     * a[1,2,1]
     * a[0,0,2]
     * a[1,0,2]
     * a[0,1,2]
     * a[1,1,2]
     * a[0,2,2]
     * a[1,2,2]
     * a[0,0,3]
     * a[1,0,3]
     * a[0,1,3]
     * a[1,1,3]
     * a[0,2,3]
     * a[1,2,3]
     *
     * @param name The name of the array variable 
     * @param dim The dimensions of the array 
     * @param fnames[out] Where the names would be pushed. 
     * @param first_is_one[true] Where to start for the first index: 0 or 1. 
     *
     */
    template <class T> void
    get_flatnames(const std::string& name,
                  const std::vector<T>& dim,
                  std::vector<std::string>& fnames,
                  bool col_major = true,
                  bool first_is_one = true) {

      fnames.clear(); 
      if (0 == dim.size()) {
        fnames.push_back(name);
        return;
      }

      std::vector<std::vector<T> > idx;
      expand_indices(dim, idx, col_major); 
      size_t first = first_is_one ? 1 : 0;
      for (typename std::vector<std::vector<T> >::const_iterator it = idx.begin();
           it != idx.end();
           ++it) {
        std::stringstream stri;
        stri << name << "[";

        size_t lenm1 = it -> size() - 1;
        for (size_t i = 0; i < lenm1; i++)
          stri << ((*it)[i] + first) << ",";
        stri << ((*it)[lenm1] + first) << "]";
        fnames.push_back(stri.str());
      }
    }

    // vectorize get_flatnames 
    template <class T> 
    void get_all_flatnames(const std::vector<std::string>& names, 
                           const std::vector<T>& dims, 
                           std::vector<std::string>& fnames, 
                           bool col_major = true) {
      fnames.clear(); 
      for (size_t i = 0; i < names.size(); ++i) {
        std::vector<std::string> i_names; 
        get_flatnames(names[i], dims[i], i_names, col_major, true); // col_major = true
        fnames.insert(fnames.end(), i_names.begin(), i_names.end());
      } 
    } 

    /* To facilitate transform an array variable ordered by col-major index
     * to row-major index order by providing the transforming indices.
     * For example, we have "x[2,3]", then if ordered by col-major, we have
     * 
     * x[1,1], x[2,1], x[1,2], x[2,2], x[1,3], x[3,1]
     * 
     * Then the indices for transforming to row-major order are 
     * [0, 2, 4, 1, 3, 5] + start. 
     *
     * @param dim[in] the dimension of the array variable, empty means a scalar
     * @param midx[out] store the indices for mapping col-major to row-major
     * @param start shifts the indices with a starting point
     *
     */ 
    template <typename T, typename T2>
    void get_indices_col2row(const std::vector<T>& dim, std::vector<T2>& midx,
                             T start = 0) {
      size_t len = dim.size();
      if (len < 1) { 
        midx.push_back(start); 
        return; 
      }
    
      std::vector<T> z(len, 1);
      for (size_t i = 1; i < len; i++) {
        z[i] *= z[i - 1] * dim[i - 1];
      } 
    
      T total = calc_num_params(dim);
      midx.resize(total);
      std::fill_n(midx.begin(), total, start);
      std::vector<T> v(len, 0);
      for (T i = 1; i < total; i++) {
        for (size_t j = 0; j < len; ++j) {
          size_t k = len - j - 1;
          if (v[k] < dim[k] - 1) {
            v[k] += 1;
            break; 
          }
          v[k] = 0; 
        } 
        // v is the index of the ith element by row-major, for example v=[0,1,2]. 
        // obtain the position for v if it is col-major indexed. 
        T pos = 0;
        for (size_t j = 0; j < len; j++) 
          pos += z[j] * v[j];
        midx[i] += pos;
      } 
    } 
   
    template <class T>
    void get_all_indices_col2row(const std::vector<std::vector<T> >& dims,
                                 std::vector<size_t>& midx) {
      midx.clear();
      std::vector<T> starts; 
      calc_starts(dims, starts);
      for (size_t i = 0; i < dims.size(); ++i) {
        std::vector<size_t> midxi;
        get_indices_col2row(dims[i], midxi, starts[i]);
        midx.insert(midx.end(), midxi.begin(), midxi.end());
      } 
    } 
    
    template <class Model>
    std::vector<std::string> get_param_names(Model& m) { 
      std::vector<std::string> names;
      m.get_param_names(names);
      names.push_back("lp__");
      return names; 
    }

    template <class T>
    void print_vector(const std::vector<T>& v, std::ostream& o, 
                      const std::vector<size_t>& midx, 
                      const std::string& sep = ",") {
      if (v.size() > 0)
        o << v[0];
      for (size_t i = 1; i < v.size(); i++)
        o << sep << v[midx.at(i)];
      o << std::endl;
    }

    template <class T>
    void print_vector(const std::vector<T>& v, std::ostream& o, 
                      const std::string& sep = ",") {
      if (v.size() > 0)
        o << v[0];
      for (size_t i = 1; i < v.size(); i++)
        o << sep << v[i];
      o << std::endl;
    }

    template <class Model, class RNG_t>
    void run_markov_chain(stan::mcmc::base_mcmc* sampler_ptr,
                          const stan_args& args, 
                          bool is_warmup,
                          mcmc_output<Model>& outputter, 
                          stan::mcmc::sample& init_s,
                          Model& model,
                          std::vector<Rcpp::NumericVector>& chains, 
                          int& iter_save_i,
                          const std::vector<size_t>& qoi_idx,
                          std::vector<double>& sum_pars,
                          double& sum_lp,
                          std::vector<Rcpp::NumericVector>& sampler_params, 
                          std::vector<Rcpp::NumericVector>& iter_params, 
                          std::string& adaptation_info, 
                          RNG_t& base_rng) {
      int start = 0;
      int end = args.get_ctrl_sampling_warmup();
      if (!is_warmup) { 
        start = end;
        end = args.get_iter();
      } 
      for (int m = start; m < end; ++m) {
        stan::common::print_progress(m, 0, args.get_iter(), 
                                     args.get_ctrl_sampling_refresh(), is_warmup,
                                     "\r", "", rstan::io::rcout);
        R_CheckUserInterrupt();
        init_s = sampler_ptr -> transition(init_s);
        if (args.get_ctrl_sampling_save_warmup() && (((m - start) % args.get_ctrl_sampling_thin()) == 0)) {
          outputter.output_sample_params(base_rng, init_s, sampler_ptr, model,
                                        chains, is_warmup,
                                        sampler_params, iter_params,
                                        sum_pars, sum_lp, qoi_idx,
                                        iter_save_i, &rstan::io::rcout);
          iter_save_i++;
          outputter.output_diagnostic_params(init_s, sampler_ptr);
        }
      }
    }

    template <class Model, class RNG_t>
    void warmup_phase(stan::mcmc::base_mcmc* sampler_ptr,
                      stan_args& args, 
                      mcmc_output<Model>& outputter,
                      stan::mcmc::sample& init_s,
                      Model& model,
                      std::vector<Rcpp::NumericVector>& chains, 
                      int& iter_save_i,
                      const std::vector<size_t>& qoi_idx,
                      std::vector<double>& sum_pars,
                      double& sum_lp,
                      std::vector<Rcpp::NumericVector>& sampler_params, 
                      std::vector<Rcpp::NumericVector>& iter_params, 
                      std::string& adaptation_info, 
                      RNG_t& base_rng) {
      run_markov_chain<Model, RNG_t>(sampler_ptr, args, true, outputter,
                                     init_s, model, chains, iter_save_i, qoi_idx,
                                     sum_pars, sum_lp, sampler_params, iter_params,
                                     adaptation_info, base_rng);
    }

    template <class Model, class RNG_t>
    void sampling_phase(stan::mcmc::base_mcmc* sampler_ptr,
                        const stan_args& args,
                        mcmc_output<Model>& outputter,
                        stan::mcmc::sample& init_s,
                        Model& model,
                        std::vector<Rcpp::NumericVector>& chains, 
                        int& iter_save_i,
                        const std::vector<size_t>& qoi_idx,
                        std::vector<double>& sum_pars,
                        double& sum_lp,
                        std::vector<Rcpp::NumericVector>& sampler_params, 
                        std::vector<Rcpp::NumericVector>& iter_params, 
                        std::string& adaptation_info, 
                        RNG_t& base_rng) {
      run_markov_chain<Model, RNG_t>(sampler_ptr, args, false, outputter,
                                     init_s, model, chains, iter_save_i, qoi_idx,
                                     sum_pars, sum_lp, sampler_params, iter_params,
                                     adaptation_info,
                                     base_rng);
    }
    

    /**
     * Cast a size_t vector to an unsigned int vector. 
     * The reason is that first Rcpp::wrap/as does not 
     * support size_t on some platforms; second R 
     * could not deal with 64bits integers.  
     */ 

    std::vector<unsigned int> 
    sizet_to_uint(std::vector<size_t> v1) {
      std::vector<unsigned int> v2(v1.size());
      for (size_t i = 0; i < v1.size(); ++i) 
        v2[i] = static_cast<unsigned int>(v1[i]);
      return v2;
    } 

    template <class Model>
    std::vector<std::vector<unsigned int> > get_param_dims(Model& m) {
      std::vector<std::vector<size_t> > dims; 
      m.get_dims(dims); 

      std::vector<std::vector<unsigned int> > uintdims; 
      for (std::vector<std::vector<size_t> >::const_iterator it = dims.begin();
           it != dims.end(); 
           ++it) 
        uintdims.push_back(sizet_to_uint(*it)); 

      std::vector<unsigned int> scalar_dim; // for lp__
      uintdims.push_back(scalar_dim); 
      return uintdims; 
    } 

    template<class Sampler>
    void init_static_hmc(stan::mcmc::base_mcmc* sampler_ptr, const stan_args& args) {
      double epsilon = args.get_ctrl_sampling_stepsize();
      double epsilon_jitter = args.get_ctrl_sampling_stepsize_jitter();
      double int_time = args.get_ctrl_sampling_int_time();

      Sampler* sampler_ptr2 = dynamic_cast<Sampler*>(sampler_ptr); 
      sampler_ptr2->set_nominal_stepsize_and_T(epsilon, int_time);
      sampler_ptr2->set_stepsize_jitter(epsilon_jitter);
    } 

    template<class Sampler>
    void init_nuts(stan::mcmc::base_mcmc* sampler_ptr, const stan_args& args) {
      double epsilon = args.get_ctrl_sampling_stepsize();
      double epsilon_jitter = args.get_ctrl_sampling_stepsize_jitter();
      int max_depth = args.get_ctrl_sampling_max_treedepth();

      Sampler* sampler_ptr2 = dynamic_cast<Sampler*>(sampler_ptr); 
      sampler_ptr2->set_nominal_stepsize(epsilon);
      sampler_ptr2->set_stepsize_jitter(epsilon_jitter);
      sampler_ptr2->set_max_depth(max_depth);
    }
    
    template<class Sampler>
    void init_adapt(stan::mcmc::base_mcmc* sampler_ptr, const stan_args& args, 
                    const Eigen::VectorXd& cont_params) {

      if (!args.get_ctrl_sampling_adapt_engaged()) return;

      double delta = args.get_ctrl_sampling_adapt_delta();
      double gamma = args.get_ctrl_sampling_adapt_gamma();
      double kappa = args.get_ctrl_sampling_adapt_kappa();
      double t0 = args.get_ctrl_sampling_adapt_t0();
      double epsilon = args.get_ctrl_sampling_stepsize();

      Sampler* sampler_ptr2 = dynamic_cast<Sampler*>(sampler_ptr); 
      sampler_ptr2->get_stepsize_adaptation().set_mu(log(10 * epsilon));
      sampler_ptr2->get_stepsize_adaptation().set_delta(delta);
      sampler_ptr2->get_stepsize_adaptation().set_gamma(gamma);
      sampler_ptr2->get_stepsize_adaptation().set_kappa(kappa);
      sampler_ptr2->get_stepsize_adaptation().set_t0(t0);
      sampler_ptr2->engage_adaptation();
      sampler_ptr2->z().q = cont_params;
      sampler_ptr2->init_stepsize();
    }

    template<class Sampler>
    bool init_windowed_adapt(stan::mcmc::base_mcmc* sampler_ptr, const stan_args& args, 
                             const Eigen::VectorXd& cont_params) {
 
      init_adapt<Sampler>(sampler_ptr, args, cont_params);
      Sampler* sampler_ptr2 = dynamic_cast<Sampler*>(sampler_ptr); 
      sampler_ptr2->set_window_params(args.get_ctrl_sampling_warmup(),
                                      args.get_ctrl_sampling_adapt_init_buffer(),
                                      args.get_ctrl_sampling_adapt_term_buffer(),
                                      args.get_ctrl_sampling_adapt_window());
      return true;
    }

    template <class T>
    std::ostream& print(std::ostream& o, const std::string name, 
                        const std::vector<T> x) {
      o << name << " (" << x.size() << ")";
      if (x.size() > 0) {
        o << ":" << std::endl 
          << "  {" << x[0];
        for (size_t n = 1; n < x.size(); n++) {
          o << "," << x[n];
        }
        o << "}";
      }
      o << std::endl;
      return o;
    }

    template <class RNG_t>
    void print_execute_sampling_input(stan_args& args, 
                                      stan::mcmc::sample& s,
                                      const std::vector<size_t>& qoi_idx, 
                                      std::vector<double>& initv,
                                      const std::vector<std::string>& fnames_oi,
                                      RNG_t& base_rng) {
      Rcpp::Rcout << "---------- execute_sampling input ----------" << std::endl;
      args.write_args_as_comment(Rcpp::Rcout);
      Rcpp::Rcout << "s:" << std::endl
                  << "  cont_parms(): " << s.cont_params() << std::endl
                  << "  log_prob:     " << s.log_prob() << std::endl
                  << "  accept_stat:  " << s.accept_stat() << std::endl;
      print(Rcpp::Rcout, "qoi_idx", qoi_idx);
      print(Rcpp::Rcout, "initv", initv);
      print(Rcpp::Rcout, "fnames_oi", fnames_oi);
      Rcpp::Rcout << "base_rng: " << base_rng << std::endl;
      Rcpp::Rcout << "--------------------------------------------" << std::endl;
      return;
    }

    void print_execute_sampling_output(Rcpp::List& holder) {
      Rcpp::Rcout << "---------- execute_sampling output ----------" << std::endl
                  << "holder:" << std::endl;
      ::Rf_PrintValue(holder);
      Rcpp::Rcout << "---------------------------------------------" << std::endl;
    }

    template <class Model, class RNG_t> 
    void execute_sampling(stan_args& args, Model& model, Rcpp::List& holder,
                          stan::mcmc::base_mcmc* sampler_ptr, 
                          stan::mcmc::sample& s,
                          const std::vector<size_t>& qoi_idx, 
                          std::vector<double>& initv,
                          std::fstream& sample_stream,
                          std::fstream& diagnostic_stream,
                          const std::vector<std::string>& fnames_oi, RNG_t& base_rng) {  
      print_execute_sampling_input(args, s, qoi_idx, initv, fnames_oi, base_rng);

      int iter_save_i = 0;
      double mean_lp(0);
      std::string adaptation_info;
      
      std::vector<Rcpp::NumericVector> chains; 
      for (unsigned int i = 0; i < qoi_idx.size(); i++) 
        chains.push_back(Rcpp::NumericVector(args.get_ctrl_sampling_iter_save())); 
      std::vector<double> mean_pars;
      mean_pars.resize(initv.size(), 0);
      std::vector<Rcpp::NumericVector> sampler_params;
      std::vector<Rcpp::NumericVector> iter_params;
      std::vector<std::string> sampler_param_names;
      std::vector<std::string> iter_param_names;

      mcmc_output<Model> outputter(&sample_stream, &diagnostic_stream);
      outputter.set_output_names(s, sampler_ptr, model, iter_param_names, sampler_param_names);
      outputter.init_sampler_params(sampler_params, args.get_ctrl_sampling_iter_save());
      outputter.init_iter_params(iter_params, args.get_ctrl_sampling_iter_save());

      if (!args.get_append_samples()) {
        outputter.print_sample_names();
        outputter.output_diagnostic_names(s, sampler_ptr, model);
      } 
      // Warm-Up
      clock_t start = clock();
      warmup_phase<Model, RNG_t>(sampler_ptr, args, outputter,
                                 s, model, chains, iter_save_i,
                                 qoi_idx, mean_pars, mean_lp,
                                 sampler_params, iter_params, adaptation_info,
                                 base_rng); 
      clock_t end = clock();
      double warmDeltaT = (double)(end - start) / CLOCKS_PER_SEC;
      if (args.get_ctrl_sampling_adapt_engaged()) { 
        dynamic_cast<stan::mcmc::base_adapter*>(sampler_ptr)->disengage_adaptation();
        outputter.output_adapt_finish(sampler_ptr, adaptation_info);
      }
      // Sampling
      start = clock();
      sampling_phase<Model, RNG_t>(sampler_ptr, args, outputter,
                                   s, model, chains, iter_save_i,
                                   qoi_idx, mean_pars, mean_lp, 
                                   sampler_params, iter_params, adaptation_info,  
                                   base_rng); 
      end = clock();
      double sampleDeltaT = (double)(end - start) / CLOCKS_PER_SEC;

      if (args.get_ctrl_sampling_iter_save_wo_warmup() > 0) {
        mean_lp /= args.get_ctrl_sampling_iter_save_wo_warmup();
        for (std::vector<double>::iterator it = mean_pars.begin();
             it != mean_pars.end(); 
             ++it) 
          (*it) /= args.get_ctrl_sampling_iter_save_wo_warmup();
      } 
      if (args.get_ctrl_sampling_refresh() > 0) { 
        rstan::io::rcout << std::endl;
        outputter.print_timing(warmDeltaT, sampleDeltaT, &rstan::io::rcout);
      }
      
      outputter.output_timing(warmDeltaT, sampleDeltaT);
      if (args.get_sample_file_flag()) {
        rstan::io::rcout << "Sample of chain " 
                         << args.get_chain_id() 
                         << " is written to file " << args.get_sample_file() << "."
                         << std::endl;
        sample_stream.close();
      }
      if (args.get_diagnostic_file_flag()) 
        diagnostic_stream.close();
     
      holder = Rcpp::List(chains.begin(), chains.end());
      holder.attr("test_grad") = Rcpp::wrap(false); 
      holder.attr("args") = args.stan_args_to_rlist(); 
      holder.attr("inits") = initv; 
      holder.attr("mean_pars") = mean_pars; 
      holder.attr("mean_lp__") = mean_lp; 
      holder.attr("adaptation_info") = adaptation_info;
      // put sampler parameters such as treedepth together with iter_params 
      iter_params.insert(iter_params.end(), sampler_params.begin(), sampler_params.end());
      // put iter_param_names and sampler_param_name tegether for returning to R
      iter_param_names.insert(iter_param_names.end(),
                              sampler_param_names.begin(),
                              sampler_param_names.end());
      Rcpp::List slst(iter_params.begin(), iter_params.end());
      slst.names() = iter_param_names;
      slst.erase(outputter.get_index_for_lp());
      holder.attr("sampler_params") = slst;
      holder.names() = fnames_oi;

      print_execute_sampling_output(holder);
    } 


    /**
     * @tparam Model 
     * @tparam RNG 
     *
     * @param args: the instance that wraps the arguments passed for sampling.
     * @param model: the model instance.
     * @param holder[out]: the object to hold all the information returned to R. 
     * @param qoi_idx: the indexes for all parameters of interest.  
     * @param fnames_oi: the parameter names of interest.  
     * @param base_rng: the boost RNG instance. 
     */
    template <class Model, class RNG_t> 
    int sampler_command(stan_args& args, Model& model, Rcpp::List& holder,
                        const std::vector<size_t>& qoi_idx, 
                        const std::vector<std::string>& fnames_oi, RNG_t& base_rng) {

      base_rng.seed(args.get_random_seed());
      // (2**50 = 1T samples, 1000 chains)
      static boost::uintmax_t DISCARD_STRIDE = 
        static_cast<boost::uintmax_t>(1) << 50;
      // rstan::io::rcout << "DISCARD_STRIDE=" << DISCARD_STRIDE << std::endl;
      base_rng.discard(DISCARD_STRIDE * (args.get_chain_id() - 1));
      
      std::vector<double> cont_vector = std::vector<double>(model.num_params_r(), 0);
      std::vector<int> disc_vector = std::vector<int>(model.num_params_i(),0);
      std::vector<double> params_inr_etc; // cont, disc, and others
      std::vector<double> init_grad;
      std::string init_val = args.get_init();
      double init_log_prob;
      int num_init_tries = 0;
      // parameter initialization
      if (init_val == "0") {
        try {
          init_log_prob
            = stan::model::log_prob_grad<true,true>(model,
                                                    cont_vector, 
                                                    disc_vector, 
                                                    init_grad, 
                                                    &rstan::io::rcout);
        } catch (const std::domain_error& e) {
          std::string msg("Domain error during initialization with 0:\n"); 
          msg += e.what();
          throw std::runtime_error(msg);
        }
        if (!boost::math::isfinite(init_log_prob))  
          throw std::runtime_error("Error during initialization with 0: vanishing density.");
        for (size_t i = 0; i < init_grad.size(); i++) {
          if (!boost::math::isfinite(init_grad[i])) 
            throw std::runtime_error("Error during initialization with 0: divergent gradient.");
        }
        Rcpp::Rcout << "Inside sampler_command, init_val == \"0\", base_rng looks like: " 
                    << base_rng << std::endl;
      } else if (init_val == "user") {
        try { 
          rstan::io::rlist_ref_var_context init_var_context(args.get_init_list()); 
          model.transform_inits(init_var_context,disc_vector,cont_vector);
          init_log_prob
            = stan::model::log_prob_grad<true,true>(model,
                                                    cont_vector, 
                                                    disc_vector, 
                                                    init_grad, 
                                                    &rstan::io::rcout);
        } catch (const std::exception& e) {
          std::string msg("Error during user-specified initialization:\n"); 
          msg += e.what(); 
          throw std::runtime_error(msg);
        } 
        if (!boost::math::isfinite(init_log_prob))  
          throw std::runtime_error("Rejecting user-specified initialization because of vanishing density.");
        for (size_t i = 0; i < init_grad.size(); i++) {
          if (!boost::math::isfinite(init_grad[i])) 
            throw std::runtime_error("Rejecting user-specified initialization because of divergent gradient.");
        }
      } else {
        init_val = "random"; 
        double r = args.get_init_radius();
        boost::random::uniform_real_distribution<double> 
          init_range_distribution(-r, r);
        boost::variate_generator<RNG_t&, boost::random::uniform_real_distribution<double> >
          init_rng(base_rng,init_range_distribution);
        
        static int MAX_INIT_TRIES = 100;
        for (; num_init_tries < MAX_INIT_TRIES; ++num_init_tries) {
          for (size_t i = 0; i < cont_vector.size(); ++i)
            cont_vector[i] = init_rng();
          try {
            init_log_prob 
              = stan::model::log_prob_grad<true,true>(model,cont_vector,disc_vector,init_grad,&rstan::io::rcout);
          } catch (const std::domain_error& e) {
            // write_error_msg(&rstan::io::rcout, e);
            rstan::io::rcout << e.what(); 
            rstan::io::rcout << "Rejecting proposed initial value with zero density." << std::endl;
            continue;
          } 
          if (!boost::math::isfinite(init_log_prob))
            continue;
          for (size_t i = 0; i < init_grad.size(); ++i)
            if (!boost::math::isfinite(init_grad[i]))
              continue;
          break;
        }
        if (num_init_tries == MAX_INIT_TRIES) {
          rstan::io::rcout << "Initialization failed after " << MAX_INIT_TRIES 
                           << " attempts. "
                           << " Try specifying initial values,"
                           << " reducing ranges of constrained values,"
                           << " or reparameterizing the model."
                           << std::endl;
          return -1;
        }
      }
      // keep a record of the initial values 
      std::vector<double> initv;
      model.write_array(base_rng,cont_vector,disc_vector,initv); 

      Rcpp::Rcout << "--- Inside sampler_command, after model.write_array(): " 
                  << base_rng << std::endl;

      if (TEST_GRADIENT == args.get_method()) {
        rstan::io::rcout << std::endl << "TEST GRADIENT MODE" << std::endl;
        std::stringstream ss; 
        double epsilon = args.get_ctrl_test_grad_epsilon();
        double error = args.get_ctrl_test_grad_error();
        int num_failed = stan::model::test_gradients<true,true>(model,cont_vector,disc_vector,epsilon,error,ss);
        rstan::io::rcout << ss.str() << std::endl; 
        holder = Rcpp::List::create(Rcpp::_["num_failed"] = num_failed);
        holder.attr("test_grad") = Rcpp::wrap(true);
        holder.attr("inits") = initv; 
        return 0;
      } 

      std::fstream sample_stream; 
      std::fstream diagnostic_stream;
      bool append_samples(args.get_append_samples());
      if (args.get_sample_file_flag()) {
        std::ios_base::openmode samples_append_mode
          = append_samples ? (std::fstream::out | std::fstream::app)
                           : std::fstream::out;
        sample_stream.open(args.get_sample_file().c_str(), samples_append_mode);
      }

      if (OPTIM == args.get_method()) { // point estimation
        if (BFGS == args.get_ctrl_optim_algorithm()) {
          rstan::io::rcout << "STAN OPTIMIZATION COMMAND (BFGS)" << std::endl;
          rstan::io::rcout << "init = " << init_val << std::endl;
          if (num_init_tries > 0)
            rstan::io::rcout << "init tries = " << num_init_tries << std::endl;
          if (args.get_sample_file_flag()) 
            rstan::io::rcout << "output = " << args.get_sample_file() << std::endl;
          rstan::io::rcout << "save_iterations = " << args.get_ctrl_optim_save_iterations() << std::endl;
          rstan::io::rcout << "init_alpha = " << args.get_ctrl_optim_init_alpha() << std::endl;
          rstan::io::rcout << "tol_obj = " << args.get_ctrl_optim_tol_obj() << std::endl;
          rstan::io::rcout << "tol_grad = " << args.get_ctrl_optim_tol_grad() << std::endl;
          rstan::io::rcout << "tol_param = " << args.get_ctrl_optim_tol_param() << std::endl;
          rstan::io::rcout << "seed = " << args.get_random_seed() << std::endl;
          
          if (args.get_sample_file_flag()) { 
            write_comment(sample_stream,"Point Estimate Generated by Stan (BFGS)");
            write_comment(sample_stream);
            write_comment_property(sample_stream,"stan_version_major",stan::MAJOR_VERSION);
            write_comment_property(sample_stream,"stan_version_minor",stan::MINOR_VERSION);
            write_comment_property(sample_stream,"stan_version_patch",stan::PATCH_VERSION);
            write_comment_property(sample_stream,"init",init_val);
            write_comment_property(sample_stream,"save_iterations",args.get_ctrl_optim_save_iterations());
            write_comment_property(sample_stream,"init_alpha",args.get_ctrl_optim_init_alpha());
            write_comment_property(sample_stream,"tol_obj",args.get_ctrl_optim_tol_obj());
            write_comment_property(sample_stream,"tol_grad",args.get_ctrl_optim_tol_grad());
            write_comment_property(sample_stream,"tol_param",args.get_ctrl_optim_tol_param());
            write_comment_property(sample_stream,"seed",args.get_random_seed());
            write_comment(sample_stream);
            
            sample_stream << "lp__,"; // log probability first
            model.write_csv_header(sample_stream);
          } 
          
          stan::optimization::BFGSLineSearch<Model> bfgs(model, cont_vector, disc_vector,
                                                         &rstan::io::rcout);
          bfgs._opts.alpha0 = args.get_ctrl_optim_init_alpha();
          bfgs._opts.tolF = args.get_ctrl_optim_tol_obj();
          bfgs._opts.tolGrad = args.get_ctrl_optim_tol_grad();
          bfgs._opts.tolX = args.get_ctrl_optim_tol_param();
          bfgs._opts.maxIts = args.get_iter();
          
          double lp = bfgs.logp();
          
          rstan::io::rcout << "initial log joint probability = " << lp << std::endl;
          int ret = 0;
          while (0 == ret) {
            int i = bfgs.iter_num();
            
            if (stan::common::do_print(i, i==0, 50 * args.get_ctrl_optim_refresh())) {
              rstan::io::rcout << "    Iter ";
              rstan::io::rcout << "     log prob ";
              rstan::io::rcout << "       ||dx|| ";
              rstan::io::rcout << "     ||grad|| ";
              rstan::io::rcout << "      alpha ";
              rstan::io::rcout << " # evals ";
              rstan::io::rcout << " Notes " << std::endl;
            }
            ret = bfgs.step();
            lp = bfgs.logp();
            bfgs.params_r(cont_vector);

            if (stan::common::do_print(i, i==0, args.get_ctrl_optim_refresh()) 
                || ret != 0 || !bfgs.note().empty()) {
              rstan::io::rcout << " " << std::setw(7) << i << " ";
              rstan::io::rcout << " " << std::setw(12) << std::setprecision(6) << lp << " ";
              rstan::io::rcout << " " << std::setw(12) << std::setprecision(6) << bfgs.prev_step_size() << " ";
              rstan::io::rcout << " " << std::setw(12) << std::setprecision(6) << bfgs.curr_g().norm() << " ";
              rstan::io::rcout << " " << std::setw(10) << std::setprecision(4) << bfgs.alpha() << " ";
              rstan::io::rcout << " " << std::setw(7) << bfgs.grad_evals() << " ";
              rstan::io::rcout << " " << bfgs.note() << " ";
              rstan::io::rcout << std::endl;
              rstan::io::rcout.flush();
            }

            if (args.get_ctrl_optim_save_iterations()) {
              sample_stream << lp << ',';
              model.write_csv(base_rng,cont_vector,disc_vector,sample_stream);
              sample_stream.flush();
            }
          }
          if (args.get_ctrl_optim_refresh() > 0) {
            if (ret >= 0) {
              rstan::io::rcout << "Optimization terminated normally: ";
            } else {
              rstan::io::rcout << "Optimization terminated with error: ";
            }
            rstan::io::rcout << bfgs.get_code_string(ret) << std::endl;
          } 
          
          model.write_array(base_rng,cont_vector,disc_vector, params_inr_etc);
          holder = Rcpp::List::create(Rcpp::_["par"] = params_inr_etc, 
                                      Rcpp::_["value"] = lp);
          if (args.get_sample_file_flag()) { 
            sample_stream << lp << ',';
            print_vector(params_inr_etc, sample_stream);
            sample_stream.close();
          }
        } else if (Newton == args.get_ctrl_optim_algorithm()) {
          rstan::io::rcout << "STAN OPTIMIZATION COMMAND (Newton)" << std::endl;
          if (args.get_sample_file_flag()) {
            write_comment(sample_stream,"Point Estimate Generated by Stan (Newton)");
            write_comment(sample_stream);
            write_comment_property(sample_stream,"stan_version_major",stan::MAJOR_VERSION);
            write_comment_property(sample_stream,"stan_version_minor",stan::MINOR_VERSION);
            write_comment_property(sample_stream,"stan_version_patch",stan::PATCH_VERSION);
            write_comment_property(sample_stream,"init",init_val);
            write_comment_property(sample_stream,"seed",args.get_random_seed());
            write_comment(sample_stream);
  
            sample_stream << "lp__,"; // log probability first
            model.write_csv_header(sample_stream);
          }
          std::vector<double> gradient;
          double lp = stan::model::log_prob_grad<true,true>(model, cont_vector, disc_vector, gradient);
          
          double lastlp = lp - 1;
          rstan::io::rcout << "initial log joint probability = " << lp << std::endl;
          int m = 0;
          while ((lp - lastlp) / fabs(lp) > 1e-8) {
            R_CheckUserInterrupt(); 
            lastlp = lp;
            lp = stan::optimization::newton_step(model, cont_vector, disc_vector);
            if (args.get_ctrl_optim_refresh() > 0) {
                rstan::io::rcout << "Iteration ";
                rstan::io::rcout << std::setw(2) << (m + 1) << ". ";
                rstan::io::rcout << "Log joint probability = " << std::setw(10) << lp;
                rstan::io::rcout << ". Improved by " << (lp - lastlp) << ".";
                rstan::io::rcout << std::endl;
                rstan::io::rcout.flush();
            }
            m++;
            if (args.get_sample_file_flag()) { 
              sample_stream << lp << ',';
              model.write_csv(base_rng,cont_vector,disc_vector,sample_stream);
            }
          }
          model.write_array(base_rng, cont_vector, disc_vector, params_inr_etc);
          holder = Rcpp::List::create(Rcpp::_["par"] = params_inr_etc, 
                                      Rcpp::_["value"] = lp);
  
          if (args.get_sample_file_flag()) { 
            sample_stream << lp << ',';
            print_vector(params_inr_etc, sample_stream);
            sample_stream.close();
          }
        } else if (Nesterov == args.get_ctrl_optim_algorithm()) {
          rstan::io::rcout << "STAN OPTIMIZATION COMMAND (Nesterov)" << std::endl;
          rstan::io::rcout << "stepsize = " << args.get_ctrl_optim_stepsize() << std::endl;
          if (args.get_sample_file_flag()) {
            write_comment(sample_stream,"Point Estimate Generated by Stan (Nesterov)");
            write_comment(sample_stream);
            write_comment_property(sample_stream,"stan_version_major",stan::MAJOR_VERSION);
            write_comment_property(sample_stream,"stan_version_minor",stan::MINOR_VERSION);
            write_comment_property(sample_stream,"stan_version_patch",stan::PATCH_VERSION);
            write_comment_property(sample_stream,"init",init_val);
            write_comment_property(sample_stream,"seed",args.get_random_seed());
            write_comment_property(sample_stream,"stepsize",args.get_ctrl_optim_stepsize());
            write_comment(sample_stream);
  
            sample_stream << "lp__,"; // log probability first
            model.write_csv_header(sample_stream);
          }
  
          stan::optimization::NesterovGradient<Model> ng(model, cont_vector, disc_vector,
                                                         args.get_ctrl_optim_stepsize(),
                                                         &rstan::io::rcout);
          double lp = ng.logp();
          double lastlp = lp - 1;
          rstan::io::rcout << "initial log joint probability = " << lp << std::endl;
          int m = 0;
          for (int i = 0; i < args.get_iter(); i++) {
            R_CheckUserInterrupt(); 
            lastlp = lp;
            lp = ng.step();
            ng.params_r(cont_vector);
            if (stan::common::do_print(i, i==0, args.get_ctrl_optim_refresh())) {
              rstan::io::rcout << "Iteration ";
              rstan::io::rcout << std::setw(2) << (m + 1) << ". ";
              rstan::io::rcout << "Log joint probability = " << std::setw(10) << lp;
              rstan::io::rcout << ". Improved by " << (lp - lastlp) << ".";
              rstan::io::rcout << std::endl;
              rstan::io::rcout.flush();
            }
            m++;
            if (args.get_sample_file_flag()) {
              sample_stream << lp << ',';
              model.write_csv(base_rng,cont_vector,disc_vector,sample_stream);
            }
          }
  
          sample_stream << lp << ',';
          model.write_array(base_rng,cont_vector,disc_vector,params_inr_etc);
          holder = Rcpp::List::create(Rcpp::_["par"] = params_inr_etc, 
                                      Rcpp::_["value"] = lp);
          if (args.get_sample_file_flag()) { 
            sample_stream << lp << ',';
            print_vector(params_inr_etc, sample_stream);
            sample_stream.close();
          }
        } 
        return 0;
      }  
      Eigen::VectorXd cont_params = Eigen::VectorXd::Zero(model.num_params_r());
      for (size_t i = 0; i < cont_vector.size(); i++) cont_params(i) = cont_vector[i];
  
      // method = 3 //sampling 
      if (args.get_diagnostic_file_flag())
        diagnostic_stream.open(args.get_diagnostic_file().c_str(), std::fstream::out);

      if (args.get_sample_file_flag()) {
        write_comment(sample_stream,"Samples Generated by Stan");
        write_comment_property(sample_stream,"stan_version_major",stan::MAJOR_VERSION);
        write_comment_property(sample_stream,"stan_version_minor",stan::MINOR_VERSION);
        write_comment_property(sample_stream,"stan_version_patch",stan::PATCH_VERSION);
        args.write_args_as_comment(sample_stream); 
      } 
      if (args.get_diagnostic_file_flag()) {
        write_comment(diagnostic_stream,"Samples Generated by Stan");
        write_comment_property(diagnostic_stream,"stan_version_major",stan::MAJOR_VERSION);
        write_comment_property(diagnostic_stream,"stan_version_minor",stan::MINOR_VERSION);
        write_comment_property(diagnostic_stream,"stan_version_patch",stan::PATCH_VERSION);
        args.write_args_as_comment(diagnostic_stream);
      } 
     
      int engine_index = 0;
       
      sampling_metric_t metric = args.get_ctrl_sampling_metric(); // unit_e, diag_e, dense_e;
      int metric_index = 0;
      if (UNIT_E == metric) metric_index = 0;
      else if (DIAG_E == metric) metric_index = 1;
      else if(DENSE_E == metric) metric_index = 2;
      sampling_algo_t algorithm = args.get_ctrl_sampling_algorithm();
      switch (algorithm) {
         case Metropolis: engine_index = 3; break;
         case HMC: engine_index = 0; break;
         case NUTS: engine_index = 1; break;
      } 

      stan::mcmc::sample s(cont_params, 0, 0);

      int sampler_select = engine_index + 10 * metric_index;
      if (args.get_ctrl_sampling_adapt_engaged())  sampler_select += 100;
      switch (sampler_select) {
        case 0: {
          typedef stan::mcmc::unit_e_static_hmc<Model, RNG_t> sampler_t;
          sampler_t sampler(model, base_rng);
          init_static_hmc<sampler_t>(&sampler, args);
          execute_sampling(args, model, holder, &sampler, s, qoi_idx, initv,
                           sample_stream, diagnostic_stream, fnames_oi,
                           base_rng);
          break;
        }
        case 1: {
          typedef stan::mcmc::unit_e_nuts<Model, RNG_t> sampler_t;
          sampler_t sampler(model, base_rng, &rstan::io::rcout, &rstan::io::rcerr);
          init_nuts<sampler_t>(&sampler, args);
          execute_sampling(args, model, holder, &sampler, s, qoi_idx, initv,
                           sample_stream, diagnostic_stream, fnames_oi,
                           base_rng);
          break;
        }
        case 10: {
          typedef stan::mcmc::diag_e_static_hmc<Model, RNG_t> sampler_t;
          sampler_t sampler(model, base_rng, &rstan::io::rcout, &rstan::io::rcerr);
          init_static_hmc<sampler_t>(&sampler, args);
          execute_sampling(args, model, holder, &sampler, s, qoi_idx, initv,
                           sample_stream, diagnostic_stream, fnames_oi,
                           base_rng);
          break;
        }
        case 11: {
          typedef stan::mcmc::diag_e_nuts<Model, RNG_t> sampler_t;
          sampler_t sampler(model, base_rng, &rstan::io::rcout, &rstan::io::rcerr);
          init_nuts<sampler_t>(&sampler, args);
          execute_sampling(args, model, holder, &sampler, s, qoi_idx, initv,
                           sample_stream, diagnostic_stream, fnames_oi,
                           base_rng);
          break;
        }
        case 20: {
          typedef stan::mcmc::dense_e_static_hmc<Model, RNG_t> sampler_t;
          sampler_t sampler(model, base_rng, &rstan::io::rcout, &rstan::io::rcerr);
          init_static_hmc<sampler_t>(&sampler, args);
          execute_sampling(args, model, holder, &sampler, s, qoi_idx, initv,
                           sample_stream, diagnostic_stream, fnames_oi,
                           base_rng);
          break;
        }
        case 21: {
          typedef stan::mcmc::dense_e_nuts<Model, RNG_t> sampler_t;
          sampler_t sampler(model, base_rng, &rstan::io::rcout, &rstan::io::rcerr);
          init_nuts<sampler_t>(&sampler, args);
          execute_sampling(args, model, holder, &sampler, s, qoi_idx, initv,
                           sample_stream, diagnostic_stream, fnames_oi,
                           base_rng);
          break;
        }
        case 100: {
          typedef stan::mcmc::adapt_unit_e_static_hmc<Model, RNG_t> sampler_t;
          sampler_t sampler(model, base_rng, &rstan::io::rcout, &rstan::io::rcerr);
          init_static_hmc<sampler_t>(&sampler, args);
          init_adapt<sampler_t>(&sampler, args, cont_params);
          execute_sampling(args, model, holder, &sampler, s, qoi_idx, initv,
                           sample_stream, diagnostic_stream, fnames_oi,
                           base_rng);
          break;
        }
        case 101: {
          typedef stan::mcmc::adapt_unit_e_nuts<Model, RNG_t> sampler_t;
          sampler_t sampler(model, base_rng, &rstan::io::rcout, &rstan::io::rcerr);
          init_nuts<sampler_t>(&sampler, args);
          init_adapt<sampler_t>(&sampler, args, cont_params);
          execute_sampling(args, model, holder, &sampler, s, qoi_idx, initv,
                           sample_stream, diagnostic_stream, fnames_oi,
                           base_rng);
          break;
        }
        case 110: {
          typedef stan::mcmc::adapt_diag_e_static_hmc<Model, RNG_t> sampler_t;
          sampler_t sampler(model, base_rng, &rstan::io::rcout, &rstan::io::rcerr);
          init_static_hmc<sampler_t>(&sampler, args);
          init_windowed_adapt<sampler_t>(&sampler, args, cont_params);
          execute_sampling(args, model, holder, &sampler, s, qoi_idx, initv,
                           sample_stream, diagnostic_stream, fnames_oi,
                           base_rng);
          break;
        }
        case 111: {
          Rcpp::Rcout << "*** Inside sampler, adapt_diag_e_nuts: " 
                      << base_rng << std::endl;

          typedef stan::mcmc::adapt_diag_e_nuts<Model, RNG_t> sampler_t;
          sampler_t sampler(model, base_rng, &rstan::io::rcout, &rstan::io::rcerr);
          Rcpp::Rcout << "*** created sampler: "
                      << base_rng << std::endl;

          init_nuts<sampler_t>(&sampler, args);

          Rcpp::Rcout << "*** after init_nuts: "
                      << base_rng << std::endl;


          init_windowed_adapt<sampler_t>(&sampler, args, cont_params);

          Rcpp::Rcout << "*** after init_windowed_adapt: "
                      << base_rng << std::endl;
          
          execute_sampling(args, model, holder, &sampler, s, qoi_idx, initv,
                           sample_stream, diagnostic_stream, fnames_oi,
                           base_rng);
          break;
        }
        case 120: {
          typedef stan::mcmc::adapt_dense_e_static_hmc<Model, RNG_t> sampler_t;
          sampler_t sampler(model, base_rng, &rstan::io::rcout, &rstan::io::rcerr);
          init_static_hmc<sampler_t>(&sampler, args);
          init_windowed_adapt<sampler_t>(&sampler, args, cont_params);
          execute_sampling(args, model, holder, &sampler, s, qoi_idx, initv,
                           sample_stream, diagnostic_stream, fnames_oi,
                           base_rng);
          break;
        }
        case 121: {
          typedef stan::mcmc::adapt_dense_e_nuts<Model, RNG_t> sampler_t;
          sampler_t sampler(model, base_rng, &rstan::io::rcout, &rstan::io::rcerr);
          init_nuts<sampler_t>(&sampler, args);
          init_windowed_adapt<sampler_t>(&sampler, args, cont_params);
          execute_sampling(args, model, holder, &sampler, s, qoi_idx, initv,
                           sample_stream, diagnostic_stream, fnames_oi,
                           base_rng);
          break;
        }
        default: 
          throw std::invalid_argument("No sampler matching HMC specification!");
      } 
      return 0;
    }
  }



  template <class Model, class RNG_t> 
  class stan_fit {

  private:
    io::rlist_ref_var_context data_;
    Model model_;
    RNG_t base_rng; 
    const std::vector<std::string> names_;
    const std::vector<std::vector<unsigned int> > dims_; 
    const unsigned int num_params_; 

    std::vector<std::string> names_oi_; // parameters of interest 
    std::vector<std::vector<unsigned int> > dims_oi_; 
    std::vector<size_t> names_oi_tidx_;  // the total indexes of names2.
    // std::vector<size_t> midx_for_col2row; // indices for mapping col-major to row-major
    std::vector<unsigned int> starts_oi_;  
    unsigned int num_params2_;  // total number of POI's.   
    std::vector<std::string> fnames_oi_; 
    Rcpp::Function cxxfunction; // keep a reference to the cxxfun, no functional purpose.

  private: 
    /**
     * Tell if a parameter name is an element of an array parameter. 
     * Note that it only supports full specified name; slicing 
     * is not supported. The test only tries to see if there 
     * are brackets. 
     */
    bool is_flatname(const std::string& name) {
      return name.find('[') != name.npos && name.find(']') != name.npos; 
    } 
  
    /*
     * Update the parameters we are interested for the model. 
     * As well, the dimensions vector for the parameters are 
     * updated. 
     */
    void update_param_oi0(const std::vector<std::string>& pnames) {
      names_oi_.clear(); 
      dims_oi_.clear(); 
      names_oi_tidx_.clear(); 

      std::vector<unsigned int> starts; 
      calc_starts(dims_, starts);
      for (std::vector<std::string>::const_iterator it = pnames.begin(); 
           it != pnames.end(); 
           ++it) { 
        size_t p = find_index(names_, *it); 
        if (p != names_.size()) {
          names_oi_.push_back(*it); 
          dims_oi_.push_back(dims_[p]); 
          if (*it == "lp__") {
            names_oi_tidx_.push_back(-1); // -1 for lp__ as it is not really a parameter  
            continue;
          } 
          size_t i_num = calc_num_params(dims_[p]); 
          size_t i_start = starts[p]; 
          for (size_t j = i_start; j < i_start + i_num; j++)
            names_oi_tidx_.push_back(j);
        } 
      }
      calc_starts(dims_oi_, starts_oi_);
      num_params2_ = names_oi_tidx_.size(); 
    } 

  public:
    SEXP update_param_oi(SEXP pars) {
      std::vector<std::string> pnames = 
        Rcpp::as<std::vector<std::string> >(pars);  
      if (std::find(pnames.begin(), pnames.end(), "lp__") == pnames.end()) 
        pnames.push_back("lp__"); 
      update_param_oi0(pnames); 
      get_all_flatnames(names_oi_, dims_oi_, fnames_oi_, true); 
      return Rcpp::wrap(true); 
    } 

    stan_fit(SEXP data, SEXP cxxf) : 
      data_(data),
      model_(data_, &rstan::io::rcout),  
      base_rng(static_cast<boost::uint32_t>(std::time(0))),
      names_(get_param_names(model_)), 
      dims_(get_param_dims(model_)), 
      num_params_(calc_total_num_params(dims_)), 
      names_oi_(names_), 
      dims_oi_(dims_),
      num_params2_(num_params_),
      cxxfunction(cxxf)  
    {
      for (size_t j = 0; j < num_params2_ - 1; j++) 
        names_oi_tidx_.push_back(j);
      names_oi_tidx_.push_back(-1); // lp__
      calc_starts(dims_oi_, starts_oi_);
      get_all_flatnames(names_oi_, dims_oi_, fnames_oi_, true); 
      // get_all_indices_col2row(dims_, midx_for_col2row);
    }             

    /**
     * Transform the parameters from its defined support
     * to unconstrained space 
     * 
     * @param par An R list as for specifying the initial values
     *  for a chain 
     */
    SEXP unconstrain_pars(SEXP par) {
      BEGIN_RCPP
      rstan::io::rlist_ref_var_context par_context(par);
      std::vector<int> params_i;
      std::vector<double> params_r;
      model_.transform_inits(par_context, params_i, params_r);
      SEXP __sexp_result;
      PROTECT(__sexp_result = Rcpp::wrap(params_r));
      UNPROTECT(1);
      return __sexp_result;
      END_RCPP
    } 

    SEXP unconstrained_param_names(SEXP include_tparams, SEXP include_gqs) {
      BEGIN_RCPP
      std::vector<std::string> n;
      model_.unconstrained_param_names(n, Rcpp::as<bool>(include_tparams), 
                                       Rcpp::as<bool>(include_gqs));
      SEXP __sexp_result;
      PROTECT(__sexp_result = Rcpp::wrap(n));
      UNPROTECT(1);
      return __sexp_result;
      END_RCPP
    } 

    SEXP constrained_param_names(SEXP include_tparams, SEXP include_gqs) {
      BEGIN_RCPP
      std::vector<std::string> n;
      model_.unconstrained_param_names(n, Rcpp::as<bool>(include_tparams), 
                                       Rcpp::as<bool>(include_gqs));
      SEXP __sexp_result;
      PROTECT(__sexp_result = Rcpp::wrap(n));
      UNPROTECT(1);
      return __sexp_result;
      END_RCPP
    } 

    /**
     * Contrary to unconstrain_pars, transform parameters
     * from unconstrained support to the constrained. 
     * 
     * @param upar The parameter values on the unconstrained 
     *  space 
     */ 
    SEXP constrain_pars(SEXP upar) {
      BEGIN_RCPP
      std::vector<double> par;
      std::vector<double> params_r = Rcpp::as<std::vector<double> >(upar);
      if (params_r.size() != model_.num_params_r()) {
        std::stringstream msg; 
        msg << "Number of unconstrained parameters does not match " 
               "that of the model (" 
            << params_r.size() << " vs " 
            << model_.num_params_r() 
            << ").";
        throw std::domain_error(msg.str()); 
      } 
      std::vector<int> params_i(model_.num_params_i());
      model_.write_array(base_rng, params_r, params_i, par);
      SEXP __sexp_result;
      PROTECT(__sexp_result = Rcpp::wrap(par));
      UNPROTECT(1);
      return __sexp_result;
      END_RCPP
    } 
  
    /**
     * Expose the log_prob of the model to stan_fit so R user
     * can call this function. 
     * 
     * @param upar The real parameters on the unconstrained 
     *  space. 
     */
    SEXP log_prob(SEXP upar, SEXP jacobian_adjust_transform, SEXP gradient) {
      BEGIN_RCPP
      using std::vector;
      vector<double> par_r = Rcpp::as<vector<double> >(upar);
      if (par_r.size() != model_.num_params_r()) {
        std::stringstream msg; 
        msg << "Number of unconstrained parameters does not match " 
               "that of the model (" 
            << par_r.size() << " vs " 
            << model_.num_params_r() 
            << ").";
        throw std::domain_error(msg.str()); 
      } 
      vector<int> par_i(model_.num_params_i(), 0);
      if (!Rcpp::as<bool>(gradient)) { 
        if (Rcpp::as<bool>(jacobian_adjust_transform)) {
          return Rcpp::wrap(stan::model::log_prob_propto<true>(model_, par_r, par_i, &rstan::io::rcout));
        } else {
          return Rcpp::wrap(stan::model::log_prob_propto<false>(model_, par_r, par_i, &rstan::io::rcout));
        } 
      } 

      std::vector<double> gradient; 
      double lp;
      if (Rcpp::as<bool>(jacobian_adjust_transform)) 
        lp = stan::model::log_prob_grad<true,true>(model_, par_r, par_i, gradient, &rstan::io::rcout);
      else 
        lp = stan::model::log_prob_grad<true,false>(model_, par_r, par_i, gradient, &rstan::io::rcout);
      Rcpp::NumericVector lp2 = Rcpp::wrap(lp);
      lp2.attr("gradient") = gradient;
      SEXP __sexp_result;
      PROTECT(__sexp_result = Rcpp::wrap(lp2));
      UNPROTECT(1);
      return __sexp_result;
      END_RCPP
    } 

    /**
     * Expose the grad_log_prob of the model to stan_fit so R user
     * can call this function. 
     * 
     * @param upar The real parameters on the unconstrained 
     *  space. 
     * @param jacobian_adjust_transform TRUE/FALSE, whether
     *  we add the term due to the transform from constrained
     *  space to unconstrained space implicitly done in Stan.
     */
    SEXP grad_log_prob(SEXP upar, SEXP jacobian_adjust_transform) {
      BEGIN_RCPP
      std::vector<double> par_r = Rcpp::as<std::vector<double> >(upar);
      if (par_r.size() != model_.num_params_r()) {
        std::stringstream msg; 
        msg << "Number of unconstrained parameters does not match " 
               "that of the model (" 
            << par_r.size() << " vs " 
            << model_.num_params_r() 
            << ").";
        throw std::domain_error(msg.str()); 
      } 
      std::vector<int> par_i(model_.num_params_i(), 0);
      std::vector<double> gradient; 
      double lp;
      if (Rcpp::as<bool>(jacobian_adjust_transform)) 
        lp = stan::model::log_prob_grad<true,true>(model_, par_r, par_i, gradient, &rstan::io::rcout);
      else 
        lp = stan::model::log_prob_grad<true,false>(model_, par_r, par_i, gradient, &rstan::io::rcout);
      Rcpp::NumericVector grad = Rcpp::wrap(gradient); 
      grad.attr("log_prob") = lp;
      SEXP __sexp_result;
      PROTECT(__sexp_result = Rcpp::wrap(grad));
      UNPROTECT(1);
      return __sexp_result;
      END_RCPP
    } 

    /**
     * Return the number of unconstrained parameters 
     */ 
    SEXP num_pars_unconstrained() {
      BEGIN_RCPP
      int n = model_.num_params_r();
      SEXP __sexp_result;
      PROTECT(__sexp_result = Rcpp::wrap(n));
      UNPROTECT(1);
      return __sexp_result;
      END_RCPP
    } 
    
    SEXP call_sampler(SEXP args_) { 
      BEGIN_RCPP 
      Rcpp::List lst_args(args_); 
      stan_args args(lst_args); 
      Rcpp::List holder;

      int ret;
      ret = sampler_command(args, model_, holder, names_oi_tidx_, 
                            fnames_oi_, base_rng);
      if (ret != 0) {
        return R_NilValue;  // indicating error happened 
      } 
      SEXP __sexp_result;
      PROTECT(__sexp_result = Rcpp::wrap(holder));
      UNPROTECT(1);
      return __sexp_result;
      END_RCPP 
    } 

    SEXP param_names() const {
      BEGIN_RCPP 
      SEXP __sexp_result;
      PROTECT(__sexp_result = Rcpp::wrap(names_));
      UNPROTECT(1);
      return __sexp_result;
      END_RCPP 
    } 

    SEXP param_names_oi() const {
      BEGIN_RCPP 
      SEXP __sexp_result;
      PROTECT(__sexp_result = Rcpp::wrap(names_oi_));
      UNPROTECT(1);
      return __sexp_result;
      END_RCPP 
    } 

    /**
     * tidx (total indexes) 
     * the index is among those parameters of interest, not 
     * all the parameters. 
     */ 
    SEXP param_oi_tidx(SEXP pars) {
      BEGIN_RCPP 
      std::vector<std::string> names = 
        Rcpp::as<std::vector<std::string> >(pars); 
      std::vector<std::string> names2; 
      std::vector<std::vector<unsigned int> > indexes; 
      for (std::vector<std::string>::const_iterator it = names.begin();
           it != names.end(); 
           ++it) {
        if (is_flatname(*it)) { // an element of an array  
          size_t ts = std::distance(fnames_oi_.begin(),
                                    std::find(fnames_oi_.begin(), 
                                              fnames_oi_.end(), *it));       
          if (ts == fnames_oi_.size()) // not found 
            continue; 
          names2.push_back(*it); 
          indexes.push_back(std::vector<unsigned int>(1, ts)); 
          continue;
        }
        size_t j = std::distance(names_oi_.begin(), 
                                 std::find(names_oi_.begin(),    
                                           names_oi_.end(), *it)); 
        if (j == names_oi_.size()) // not found 
          continue; 
        unsigned int j_size = calc_num_params(dims_oi_[j]); 
        unsigned int j_start = starts_oi_[j]; 
        std::vector<unsigned int> j_idx; 
        for (unsigned int k = 0; k < j_size; k++) {
          j_idx.push_back(j_start + k); 
        } 
        names2.push_back(*it); 
        indexes.push_back(j_idx); 
      }
      Rcpp::List lst = Rcpp::wrap(indexes); 
      lst.names() = names2; 
      SEXP __sexp_result;
      PROTECT(__sexp_result = Rcpp::wrap(lst));
      UNPROTECT(1);
      return __sexp_result;
      END_RCPP
    } 


    SEXP param_dims() const {
      BEGIN_RCPP 
      Rcpp::List lst = Rcpp::wrap(dims_); 
      lst.names() = names_; 
      SEXP __sexp_result;
      PROTECT(__sexp_result = Rcpp::wrap(lst));
      UNPROTECT(1);
      return __sexp_result;
      END_RCPP
    } 

    SEXP param_dims_oi() const {
      BEGIN_RCPP 
      Rcpp::List lst = Rcpp::wrap(dims_oi_); 
      lst.names() = names_oi_; 
      SEXP __sexp_result;
      PROTECT(__sexp_result = Rcpp::wrap(lst));
      UNPROTECT(1);
      return __sexp_result;
      END_RCPP
    } 
    
    SEXP param_fnames_oi() const {
      BEGIN_RCPP 
      std::vector<std::string> fnames; 
      get_all_flatnames(names_oi_, dims_oi_, fnames, true); 
      SEXP __sexp_result;
      PROTECT(__sexp_result = Rcpp::wrap(fnames));
      UNPROTECT(1);
      return __sexp_result;
      END_RCPP
    } 
  };
} 

#endif 

/*
 * compile to check syntax error
 */
/*
STAN= ../../../../../ 
RCPPINC=`Rscript -e "cat(system.file('include', package='Rcpp'))"`
RINC=`Rscript -e "cat(R.home('include'))"` 
g++ -Wall -I${RINC} -I"${STAN}/lib/boost_1.51.0" -I"${STAN}/lib/eigen_3.1.1"  -I"${STAN}/src" -I"${RCPPINC}" -I"../" stan_fit.hpp 
*/

