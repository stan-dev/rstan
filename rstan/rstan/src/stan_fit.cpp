#include <rstan_next/stan_fit.hpp>

#include <cstring>
#include <iomanip>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>

//#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/version.hpp>

#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/random/additive_combine.hpp> // L'Ecuyer RNG
#include <boost/random/uniform_real_distribution.hpp>

//#include <Rcpp.h>
//#include <RcppEigen.h>

//http://cran.r-project.org/doc/manuals/R-exts.html#Allowing-interrupts
#include <R_ext/Utils.h>
// void R_CheckUserInterrupt(void);


// REF: cmdstan: src/cmdstan/command.hpp
#include <stan/callbacks/interrupt.hpp>
#include <stan/callbacks/stream_logger.hpp>
#include <stan/callbacks/stream_writer.hpp>
#include <stan/callbacks/writer.hpp>
#include <stan/io/empty_var_context.hpp>
#include <stan/model/model_base.hpp>
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
#include <stan/services/sample/standalone_gqs.hpp>
#include <stan/services/experimental/advi/fullrank.hpp>
#include <stan/services/experimental/advi/meanfield.hpp>


#include <rstan/io/rlist_ref_var_context.hpp>
#include <rstan/io/r_ostream.hpp>
#include <rstan/stan_args.hpp>
#include <rstan/filtered_values.hpp>
#include <rstan/sum_values.hpp>
#include <rstan/value.hpp>
#include <rstan/values.hpp>
#include <rstan/rstan_writer.hpp>
#include <rstan/logger.hpp>


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
    for (auto&& it : idx) {
      std::stringstream stri;
      stri << name << "[";

      const size_t lenm1 = it.size() - 1;
      for (size_t i = 0; i < lenm1; i++) {
        stri << (it[i] + first) << ",";
      }
      stri << (it[lenm1] + first) << "]";
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

std::vector<std::string> get_param_names(stan::model::model_base* m) {
  std::vector<std::string> names;
  m->get_param_names(names);
  names.push_back("lp__");
  return names;
}

template <class T>
void print_vector(const std::vector<T>& v, std::ostream& o,
                  const std::vector<size_t>& midx,
                  const std::string& sep = ",") {
  if (v.size() > 0)
    o << v[midx.at(0)];
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

void write_stan_version_as_comment(std::ostream& output) {
  write_comment_property(output,"stan_version_major",stan::MAJOR_VERSION);
  write_comment_property(output,"stan_version_minor",stan::MINOR_VERSION);
  write_comment_property(output,"stan_version_patch",stan::PATCH_VERSION);
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

std::vector<std::vector<unsigned int> > get_param_dims(stan::model::model_base* m) {
  std::vector<std::vector<size_t> > dims;
  m->get_dims(dims);

  std::vector<std::vector<unsigned int> > uintdims;
  for (std::vector<std::vector<size_t> >::const_iterator it = dims.begin();
       it != dims.end();
       ++it)
    uintdims.push_back(sizet_to_uint(*it));

  std::vector<unsigned int> scalar_dim; // for lp__
  uintdims.push_back(scalar_dim);
  return uintdims;
}

struct R_CheckUserInterrupt_Functor : public stan::callbacks::interrupt {
  void operator()() {
    R_CheckUserInterrupt();
  }
};

std::vector<double> unconstrained_to_constrained(stan::model::model_base* model,
                                                 unsigned int random_seed,
                                                 unsigned int id,
                                                 std::vector<double>& params) {
  std::vector<int> params_i;
  std::vector<double> constrained_params;
  boost::ecuyer1988 rng = stan::services::util::create_rng(random_seed, id);
  model->write_array(rng, params, params_i, constrained_params);
  return constrained_params;
}
/**
*
* @param args: the instance that wraps the arguments passed for sampling.
* @param model: pointer to the model instance.
* @param holder[out]: the object to hold all the information returned to R.
* @param qoi_idx: the indexes for all parameters of interest.
* @param fnames_oi: the parameter names of interest.
* @param base_rng: the PRNG
*/
int command(stan_args& args,
            stan::model::model_base* model,
            Rcpp::List& holder,
            const std::vector<size_t>& qoi_idx,
            const std::vector<std::string>& fnames_oi,
            boost::ecuyer1988& base_rng) {

  if (args.get_method() == SAMPLING
        && model->num_params_r() == 0
        && args.get_ctrl_sampling_algorithm() != Fixed_param)
        throw std::runtime_error("Must use algorithm=\"Fixed_param\" for "
                                   "model that has no parameters.");
  int refresh = args.get_refresh();
  unsigned int id = args.get_chain_id();

  std::ostream nullout(nullptr);
  std::ostream& c_out = refresh ? Rcpp::Rcout : nullout;
  std::ostream& c_err = refresh ? rstan::io::rcerr : nullout;

  stan::callbacks::stream_logger_with_chain_id
    logger(c_out, c_out, c_out, c_err, c_err, id);

  R_CheckUserInterrupt_Functor interrupt;

  std::fstream sample_stream;
  std::fstream diagnostic_stream;
  std::stringstream comment_stream;
  bool append_samples(args.get_append_samples());
  if (args.get_sample_file_flag()) {
    std::ios_base::openmode samples_append_mode
    = append_samples ? (std::fstream::out | std::fstream::app)
      : std::fstream::out;
    sample_stream.open(args.get_sample_file().c_str(), samples_append_mode);

    if (args.get_method() == TEST_GRADIENT)
      write_comment(sample_stream, "Output generated by Stan (test_grad)");
    else if (args.get_method() == OPTIM)
      write_comment(sample_stream, "Point Estimate Generated by Stan");
    else if (args.get_method() == SAMPLING)
      write_comment(sample_stream, "Sample generated by Stan");
    else if (args.get_method() == VARIATIONAL)
      write_comment(sample_stream, "Sample generated by Stan (Variational Bayes)");
    write_stan_version_as_comment(sample_stream);
    args.write_args_as_comment(sample_stream);
  }
  if (args.get_diagnostic_file_flag()) {
    diagnostic_stream.open(args.get_diagnostic_file().c_str(), std::fstream::out);

    if (args.get_method() == TEST_GRADIENT)
      write_comment(diagnostic_stream, "Output generated by Stan (test_grad)");
    else if (args.get_method() == OPTIM)
      write_comment(diagnostic_stream, "Point Estimate Generated by Stan");
    else if (args.get_method() == SAMPLING)
      write_comment(diagnostic_stream, "Sample generated by Stan");
    else if (args.get_method() == VARIATIONAL)
      write_comment(diagnostic_stream, "Sample generated by Stan (Variational Bayes)");
    write_stan_version_as_comment(diagnostic_stream);
    args.write_args_as_comment(diagnostic_stream);
  }

  stan::callbacks::stream_writer diagnostic_writer(diagnostic_stream, "# ");
  std::unique_ptr<stan::io::var_context> init_context_ptr;
  if (args.get_init() == "user")
    init_context_ptr.reset(new io::rlist_ref_var_context(args.get_init_list()));
  else
    init_context_ptr.reset(new stan::io::empty_var_context());

  std::vector<std::string> constrained_param_names;
  model->constrained_param_names(constrained_param_names);
  rstan::value init_writer;

  int return_code = stan::services::error_codes::CONFIG;


  unsigned int random_seed = args.get_random_seed();
  double init_radius = args.get_init_radius();

  if (args.get_method() == TEST_GRADIENT) {
    double epsilon = args.get_ctrl_test_grad_epsilon();
    double error = args.get_ctrl_test_grad_error();
    stan::callbacks::writer sample_writer;
    return_code = stan::services::diagnose::diagnose(*model,
                                                     *init_context_ptr,
                                                     random_seed, id,
                                                     init_radius,
                                                     epsilon, error,
                                                     interrupt,
                                                     logger,
                                                     init_writer,
                                                     sample_writer);
    holder = Rcpp::List::create(Rcpp::_["num_failed"] = return_code);
    holder.attr("test_grad") = Rcpp::wrap(true);
    holder.attr("inits") = unconstrained_to_constrained(model, random_seed, id,
                                                        init_writer.x());
  }
  if (args.get_method() == OPTIM) {
    rstan::value sample_writer;
    bool save_iterations = args.get_ctrl_optim_save_iterations();
    int num_iterations = args.get_iter();

    if (args.get_ctrl_optim_algorithm() == Newton) {
      return_code
      = stan::services::optimize::newton(*model, *init_context_ptr,
                                         random_seed, id, init_radius,
                                         num_iterations,
                                         save_iterations,
                                         interrupt, logger,
                                         init_writer, sample_writer);
    }
    if (args.get_ctrl_optim_algorithm() == BFGS) {
      double init_alpha = args.get_ctrl_optim_init_alpha();
      double tol_obj= args.get_ctrl_optim_tol_obj();
      double tol_rel_obj = args.get_ctrl_optim_tol_rel_obj();
      double tol_grad = args.get_ctrl_optim_tol_grad();
      double tol_rel_grad = args.get_ctrl_optim_tol_rel_grad();
      double tol_param = args.get_ctrl_optim_tol_param();
      return_code
        = stan::services::optimize::bfgs(*model, *init_context_ptr,
                                         random_seed, id, init_radius,
                                         init_alpha,
                                         tol_obj,
                                         tol_rel_obj,
                                         tol_grad,
                                         tol_rel_grad,
                                         tol_param,
                                         num_iterations,
                                         save_iterations,
                                         refresh,
                                         interrupt, logger,
                                         init_writer, sample_writer);
    }
    if (args.get_ctrl_optim_algorithm() == LBFGS) {
      int history_size = args.get_ctrl_optim_history_size();
      double init_alpha = args.get_ctrl_optim_init_alpha();
      double tol_obj= args.get_ctrl_optim_tol_obj();
      double tol_rel_obj = args.get_ctrl_optim_tol_rel_obj();
      double tol_grad = args.get_ctrl_optim_tol_grad();
      double tol_rel_grad = args.get_ctrl_optim_tol_rel_grad();
      double tol_param = args.get_ctrl_optim_tol_param();
      return_code
        = stan::services::optimize::lbfgs(*model, *init_context_ptr,
                                          random_seed, id, init_radius,
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
                                          interrupt, logger,
                                          init_writer, sample_writer);
    }
    std::vector<double> params = sample_writer.x();
    double lp = params.front();
    params.erase(params.begin());
    holder = Rcpp::List::create(Rcpp::_["par"] = params,
                                Rcpp::_["value"] = lp);

  }
  if (args.get_method() == SAMPLING) {
    std::vector<std::string> sample_names;
    stan::mcmc::sample::get_sample_param_names(sample_names);
    std::vector<std::string> sampler_names;

    std::unique_ptr<rstan_sample_writer> sample_writer_ptr;
    size_t sample_writer_offset;

    int num_warmup = args.get_ctrl_sampling_warmup();
    int num_samples = args.get_iter() - num_warmup;
    int num_thin = args.get_ctrl_sampling_thin();
    bool save_warmup = args.get_ctrl_sampling_save_warmup();
    int num_iter_save = args.get_ctrl_sampling_iter_save();
    int num_warmup_save = num_iter_save - args.get_ctrl_sampling_iter_save_wo_warmup();

    if (args.get_ctrl_sampling_algorithm() == Fixed_param) {
      sampler_names.resize(0);
      sample_writer_ptr.reset(sample_writer_factory(&sample_stream,
                                                    comment_stream, "# ",
                                                    sample_names.size(),
                                                    sampler_names.size(),
                                                    constrained_param_names.size(),
                                                    num_iter_save,
                                                    num_warmup_save,
                                                    qoi_idx).release());
      return_code
        = stan::services::sample::fixed_param(*model, *init_context_ptr,
                                              random_seed, id, init_radius,
                                              num_samples,
                                              num_thin,
                                              refresh,
                                              interrupt,
                                              logger, init_writer,
                                              *sample_writer_ptr, diagnostic_writer);
    } else if (args.get_ctrl_sampling_algorithm() == NUTS) {
      sampler_names.resize(5);
      sampler_names[0] = "stepsize__";
      sampler_names[1] = "treedepth__";
      sampler_names[2] = "n_leapfrog__";
      sampler_names[3] = "divergent__";
      sampler_names[4] = "energy__";
      sample_writer_offset = sample_names.size() + sampler_names.size();

      sample_writer_ptr.reset(sample_writer_factory(&sample_stream,
                                                    comment_stream, "# ",
                                                    sample_names.size(),
                                                    sampler_names.size(),
                                                    constrained_param_names.size(),
                                                    num_iter_save,
                                                    num_warmup_save,
                                                    qoi_idx).release());

      double stepsize = args.get_ctrl_sampling_stepsize();
      double stepsize_jitter = args.get_ctrl_sampling_stepsize_jitter();
      int max_depth = args.get_ctrl_sampling_max_treedepth();

      if (args.get_ctrl_sampling_metric() == DENSE_E) {
        if (!args.get_ctrl_sampling_adapt_engaged()) {
          return_code = stan::services::sample
          ::hmc_nuts_dense_e(*model, *init_context_ptr,
                             random_seed, id, init_radius,
                             num_warmup, num_samples,
                             num_thin, save_warmup, refresh,
                             stepsize, stepsize_jitter, max_depth,
                             interrupt, logger, init_writer,
                             *sample_writer_ptr, diagnostic_writer);
        } else {
          double delta = args.get_ctrl_sampling_adapt_delta();
          double gamma = args.get_ctrl_sampling_adapt_gamma();
          double kappa = args.get_ctrl_sampling_adapt_kappa();
          double t0 = args.get_ctrl_sampling_adapt_t0();
          unsigned int init_buffer = args.get_ctrl_sampling_adapt_init_buffer();
          unsigned int term_buffer = args.get_ctrl_sampling_adapt_term_buffer();
          unsigned int window = args.get_ctrl_sampling_adapt_window();

          return_code = stan::services::sample
            ::hmc_nuts_dense_e_adapt(*model, *init_context_ptr,
                                     random_seed, id, init_radius,
                                     num_warmup, num_samples,
                                     num_thin, save_warmup, refresh,
                                     stepsize, stepsize_jitter, max_depth,
                                     delta, gamma, kappa,
                                     t0, init_buffer, term_buffer, window,
                                     interrupt, logger, init_writer,
                                     *sample_writer_ptr, diagnostic_writer);
        }
      } else if (args.get_ctrl_sampling_metric() == DIAG_E) {
        if (!args.get_ctrl_sampling_adapt_engaged()) {
          return_code = stan::services::sample
          ::hmc_nuts_diag_e(*model, *init_context_ptr,
                            random_seed, id, init_radius,
                            num_warmup, num_samples,
                            num_thin, save_warmup, refresh,
                            stepsize, stepsize_jitter, max_depth,
                            interrupt, logger, init_writer,
                            *sample_writer_ptr, diagnostic_writer);
        } else {
          double delta = args.get_ctrl_sampling_adapt_delta();
          double gamma = args.get_ctrl_sampling_adapt_gamma();
          double kappa = args.get_ctrl_sampling_adapt_kappa();
          double t0 = args.get_ctrl_sampling_adapt_t0();
          unsigned int init_buffer = args.get_ctrl_sampling_adapt_init_buffer();
          unsigned int term_buffer = args.get_ctrl_sampling_adapt_term_buffer();
          unsigned int window = args.get_ctrl_sampling_adapt_window();

          return_code = stan::services::sample
            ::hmc_nuts_diag_e_adapt(*model, *init_context_ptr,
                                    random_seed, id, init_radius,
                                    num_warmup, num_samples,
                                    num_thin, save_warmup, refresh,
                                    stepsize, stepsize_jitter, max_depth,
                                    delta, gamma, kappa,
                                    t0, init_buffer, term_buffer, window,
                                    interrupt, logger, init_writer,
                                    *sample_writer_ptr, diagnostic_writer);
        }
      } else if (args.get_ctrl_sampling_metric() == UNIT_E) {
        if (!args.get_ctrl_sampling_adapt_engaged()) {
          return_code = stan::services::sample
          ::hmc_nuts_unit_e(*model, *init_context_ptr,
                            random_seed, id, init_radius,
                            num_warmup, num_samples,
                            num_thin, save_warmup, refresh,
                            stepsize, stepsize_jitter, max_depth,
                            interrupt, logger, init_writer,
                            *sample_writer_ptr, diagnostic_writer);
        } else {
          double delta = args.get_ctrl_sampling_adapt_delta();
          double gamma = args.get_ctrl_sampling_adapt_gamma();
          double kappa = args.get_ctrl_sampling_adapt_kappa();
          double t0 = args.get_ctrl_sampling_adapt_t0();

          return_code = stan::services::sample
            ::hmc_nuts_unit_e_adapt(*model, *init_context_ptr,
                                    random_seed, id, init_radius,
                                    num_warmup, num_samples,
                                    num_thin, save_warmup, refresh,
                                    stepsize, stepsize_jitter, max_depth,
                                    delta, gamma, kappa, t0,
                                    interrupt, logger, init_writer,
                                    *sample_writer_ptr, diagnostic_writer);
        }
      }
    } else if (args.get_ctrl_sampling_algorithm() == HMC) {
      sampler_names.resize(3);
      sampler_names[0] = "stepsize__";
      sampler_names[1] = "int_time__";
      sampler_names[2] = "energy__";
      sample_writer_offset = sample_names.size() + sampler_names.size();

      sample_writer_ptr.reset(sample_writer_factory(&sample_stream,
                                                    comment_stream, "# ",
                                                    sample_names.size(),
                                                    sampler_names.size(),
                                                    constrained_param_names.size(),
                                                    num_iter_save,
                                                    num_warmup_save,
                                                    qoi_idx).release());

      double stepsize = args.get_ctrl_sampling_stepsize();
      double stepsize_jitter = args.get_ctrl_sampling_stepsize_jitter();
      double int_time = args.get_ctrl_sampling_int_time();

      if (args.get_ctrl_sampling_metric() == DENSE_E) {
        if (!args.get_ctrl_sampling_adapt_engaged()) {
          return_code = stan::services::sample
          ::hmc_static_dense_e(*model, *init_context_ptr,
                               random_seed, id, init_radius,
                               num_warmup, num_samples,
                               num_thin, save_warmup, refresh,
                               stepsize, stepsize_jitter, int_time,
                               interrupt, logger, init_writer,
                               *sample_writer_ptr, diagnostic_writer);
        } else {
          double delta = args.get_ctrl_sampling_adapt_delta();
          double gamma = args.get_ctrl_sampling_adapt_gamma();
          double kappa = args.get_ctrl_sampling_adapt_kappa();
          double t0 = args.get_ctrl_sampling_adapt_t0();
          unsigned int init_buffer = args.get_ctrl_sampling_adapt_init_buffer();
          unsigned int term_buffer = args.get_ctrl_sampling_adapt_term_buffer();
          unsigned int window = args.get_ctrl_sampling_adapt_window();

          return_code = stan::services::sample
            ::hmc_static_dense_e_adapt(*model, *init_context_ptr,
                                       random_seed, id, init_radius,
                                       num_warmup, num_samples,
                                       num_thin, save_warmup, refresh,
                                       stepsize, stepsize_jitter, int_time,
                                       delta, gamma, kappa, t0,
                                       init_buffer, term_buffer, window,
                                       interrupt, logger, init_writer,
                                       *sample_writer_ptr, diagnostic_writer);

        }
      } else if (args.get_ctrl_sampling_metric() == DIAG_E) {
        if (!args.get_ctrl_sampling_adapt_engaged()) {
          return_code = stan::services::sample
          ::hmc_static_diag_e(*model, *init_context_ptr,
                              random_seed, id, init_radius,
                              num_warmup, num_samples,
                              num_thin, save_warmup, refresh,
                              stepsize, stepsize_jitter, int_time,
                              interrupt, logger, init_writer,
                              *sample_writer_ptr, diagnostic_writer);

        } else {
          double delta = args.get_ctrl_sampling_adapt_delta();
          double gamma = args.get_ctrl_sampling_adapt_gamma();
          double kappa = args.get_ctrl_sampling_adapt_kappa();
          double t0 = args.get_ctrl_sampling_adapt_t0();
          unsigned int init_buffer = args.get_ctrl_sampling_adapt_init_buffer();
          unsigned int term_buffer = args.get_ctrl_sampling_adapt_term_buffer();
          unsigned int window = args.get_ctrl_sampling_adapt_window();

          return_code = stan::services::sample
            ::hmc_static_diag_e_adapt(*model, *init_context_ptr,
                                      random_seed, id, init_radius,
                                      num_warmup, num_samples,
                                      num_thin, save_warmup, refresh,
                                      stepsize, stepsize_jitter, int_time,
                                      delta, gamma, kappa, t0,
                                      init_buffer, term_buffer, window,
                                      interrupt, logger, init_writer,
                                      *sample_writer_ptr, diagnostic_writer);
        }
      } else if (args.get_ctrl_sampling_metric() == UNIT_E) {
        if (args.get_ctrl_sampling_adapt_engaged()) {
          return_code = stan::services::sample
          ::hmc_static_unit_e(*model, *init_context_ptr,
                              random_seed, id, init_radius,
                              num_warmup, num_samples,
                              num_thin, save_warmup, refresh,
                              stepsize, stepsize_jitter, int_time,
                              interrupt, logger, init_writer,
                              *sample_writer_ptr, diagnostic_writer);

        } else  {
          double delta = args.get_ctrl_sampling_adapt_delta();
          double gamma = args.get_ctrl_sampling_adapt_gamma();
          double kappa = args.get_ctrl_sampling_adapt_kappa();
          double t0 = args.get_ctrl_sampling_adapt_t0();

          return_code = stan::services::sample
            ::hmc_static_unit_e_adapt(*model, *init_context_ptr,
                                      random_seed, id, init_radius,
                                      num_warmup, num_samples,
                                      num_thin, save_warmup, refresh,
                                      stepsize, stepsize_jitter, int_time,
                                      delta, gamma, kappa, t0,
                                      interrupt, logger, init_writer,
                                      *sample_writer_ptr, diagnostic_writer);
        }
      }
    }
    double mean_lp(0);
    std::vector<double> mean_pars;
    mean_pars.resize(constrained_param_names.size(), 0);

    sample_writer_offset = sample_names.size() + sampler_names.size();
    if (args.get_ctrl_sampling_iter_save_wo_warmup() > 0) {
      double inverse_saved = 1.0 / args.get_ctrl_sampling_iter_save_wo_warmup();
      mean_lp = sample_writer_ptr->sum_.sum()[0] * inverse_saved;
      for (size_t n = 0; n < mean_pars.size(); n++) {
        mean_pars[n] = sample_writer_ptr->sum_.sum()[sample_writer_offset + n] * inverse_saved;
      }
    }

    holder = Rcpp::List(sample_writer_ptr->values_.x().begin(),
                        sample_writer_ptr->values_.x().end());
    holder.attr("test_grad") = Rcpp::wrap(false);
    holder.attr("args") = args.stan_args_to_rlist();
    holder.attr("inits") = unconstrained_to_constrained(model, random_seed, id,
                                                        init_writer.x());
    holder.attr("mean_pars") = mean_pars;
    holder.attr("mean_lp__") = mean_lp;

    std::string comments = comment_stream.str();
    size_t start = 0;
    size_t end = 0;
    std::string adaptation_info;
    if ((start = comments.find("# Adaptation")) != std::string::npos) {
      end = comments.find("# \n", start);
      adaptation_info = comments.substr(start,
                                        end - start);
    }
    double warmDeltaT = 0;
    double sampleDeltaT = 0;
    if ((start = comments.find("Elapsed Time: ")) != std::string::npos) {
      start += strlen("Elapsed Time: ");
      end = comments.find("seconds", start + 1);
      std::stringstream ss(comments.substr(start, end));
      ss >> warmDeltaT;

      start = comments.find("# ", end) + strlen("# ");
      end = comments.find("seconds (Sampling)", start + 1);
      ss.str(comments.substr(start, end));
      ss >> sampleDeltaT;
    }
    holder.attr("adaptation_info") = adaptation_info;
    holder.attr("elapsed_time") =
      Rcpp::NumericVector::create(Rcpp::_["warmup"] = warmDeltaT,
                                  Rcpp::_["sample"] = sampleDeltaT);

    Rcpp::List slst(sample_writer_ptr->sampler_values_.x().begin()+1,
                    sample_writer_ptr->sampler_values_.x().end());

    std::vector<std::string> slst_names(sample_names.begin()+1, sample_names.end());
    slst_names.insert(slst_names.end(), sampler_names.begin(), sampler_names.end());
    slst.names() = slst_names;
    holder.attr("sampler_params") = slst;
    holder.names() = fnames_oi;
    sample_writer_ptr.reset();
  }
  if (args.get_method() == VARIATIONAL) {
    int grad_samples = args.get_ctrl_variational_grad_samples();
    int elbo_samples = args.get_ctrl_variational_elbo_samples();
    int max_iterations = args.get_iter();
    double tol_rel_obj = args.get_ctrl_variational_tol_rel_obj();
    double eta = args.get_ctrl_variational_eta();
    bool adapt_engaged = args.get_ctrl_variational_adapt_engaged();
    int adapt_iterations = args.get_ctrl_variational_adapt_iter();
    int eval_elbo = args.get_ctrl_variational_eval_elbo();
    int output_samples = args.get_ctrl_variational_output_samples();

    stan::callbacks::stream_writer sample_writer(sample_stream, "# ");

    if (args.get_ctrl_variational_algorithm() == FULLRANK) {
      return_code = stan::services::experimental::advi
      ::fullrank(*model, *init_context_ptr,
                 random_seed, id, init_radius,
                 grad_samples, elbo_samples,
                 max_iterations, tol_rel_obj, eta,
                 adapt_engaged, adapt_iterations,
                 eval_elbo, output_samples,
                 interrupt, logger, init_writer,
                 sample_writer, diagnostic_writer);
    } else {
      return_code = stan::services::experimental::advi
      ::meanfield(*model, *init_context_ptr,
                  random_seed, id, init_radius,
                  grad_samples, elbo_samples,
                  max_iterations, tol_rel_obj, eta,
                  adapt_engaged, adapt_iterations,
                  eval_elbo, output_samples,
                  interrupt, logger, init_writer,
                  sample_writer, diagnostic_writer);
    }
    holder = Rcpp::List::create(Rcpp::_["samples"] = R_NilValue);
    holder.attr("args") = args.stan_args_to_rlist();
    holder.attr("inits") = unconstrained_to_constrained(model, random_seed, id,
                                                        init_writer.x());
  }

  init_context_ptr.reset();
  if (sample_stream.is_open())
    sample_stream.close();
  if (diagnostic_stream.is_open())
    diagnostic_stream.close();

  return return_code;
}
}

bool stan_fit::is_flatname(const std::string& name) {
    return name.find('[') != name.npos && name.find(']') != name.npos;
  }

  /*
   * Update the parameters we are interested for the model->
   * As well, the dimensions vector for the parameters are
   * updated.
   */
void stan_fit::update_param_oi0(const std::vector<std::string>& pnames) {
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

bool stan_fit::update_param_oi(std::vector<std::string> pnames) {
    if (std::find(pnames.begin(), pnames.end(), "lp__") == pnames.end())
      pnames.push_back("lp__");
    update_param_oi0(pnames);
    get_all_flatnames(names_oi_, dims_oi_, fnames_oi_, true);
    return true;
  }

stan_fit::stan_fit(SEXP model_sexp, int seed) :
    model_sexp_(model_sexp),
    model_xptr_(model_sexp),
    model_(model_xptr_.get()),
    base_rng(static_cast<boost::uint32_t>(seed)),
    names_(get_param_names(model_)),
    dims_(get_param_dims(model_)),
    num_params_(calc_total_num_params(dims_)),
    names_oi_(names_),
    dims_oi_(dims_),
    num_params2_(num_params_)
{
  for (size_t j = 0; j < num_params2_ - 1; j++)
    names_oi_tidx_.push_back(j);
  names_oi_tidx_.push_back(-1); // lp__
  calc_starts(dims_oi_, starts_oi_);
  get_all_flatnames(names_oi_, dims_oi_, fnames_oi_, true);
  Rcpp::Rcpp_PreserveObject(model_sexp_);
}

stan_fit::~stan_fit() {
   Rcpp::Rcpp_ReleaseObject(model_sexp_);
}

  /**
   * Transform the parameters from its defined support
   * to unconstrained space
   *
   * @param par An R list as for specifying the initial values
   *  for a chain
   * @return A standard vector of doubles
   */
std::vector<double> stan_fit::unconstrain_pars(Rcpp::List par) {
    rstan::io::rlist_ref_var_context par_context(par);
    std::vector<int> params_i;
    std::vector<double> params_r;
    model_->transform_inits(par_context, params_i, params_r, &rstan::io::rcout);
    return params_r;
  }


  /**
   * Contrary to unconstrain_pars, transform parameters
   * from unconstrained support to the constrained.
   *
   * @param upar A standard vector of doubles on the unconstrained space
   * @return A standard vector of doubles on the constrained space
   */
  std::vector<double> stan_fit::constrain_pars(std::vector<double> upar) {
    if (upar.size() != model_->num_params_r()) {
      std::stringstream msg;
      msg << "Number of unconstrained parameters does not match "
      "that of the model ("
      << upar.size() << " vs "
      << model_->num_params_r()
      << ").";
      throw std::domain_error(msg.str());
    }
    std::vector<int> params_i(model_->num_params_i());
    std::vector<double> par;
    model_->write_array(base_rng, upar, params_i, par);
    return par;
  }

  /**
   * Get the unconstrained or constrained parameter names
   *
   * @param include_tparams Flag to include transformed parameter names
   * @param include_gqs Flag to include generated quantitiy names
   * @return A standard vector of standard strings
   */

  std::vector<std::string> stan_fit::unconstrained_param_names(bool include_tparams,
                                                     bool include_gqs) {
    std::vector<std::string> n;
    model_->unconstrained_param_names(n, include_tparams, include_gqs);
    return n;
  }

  std::vector<std::string> stan_fit::constrained_param_names(bool include_tparams,
                                                   bool include_gqs) {
    std::vector<std::string> n;
    model_->constrained_param_names(n, include_tparams, include_gqs);
    return n;
  }

  /**
   * Expose the log_prob of the model to stan_fit so R users
   * can call this function.
   *
   * @param upar The real parameters on the unconstrained
   *  space.
   * @param jacobian_adjust_transform A flag to indicate whether
   *   the Jacobian adjustment is included
   * @param gradient A flag to indicate whether to return the
   *   gradient as an attribute
   * @param A numeric vector of size 1, possibly with a grad attribute
   */
  Rcpp::NumericVector stan_fit::log_prob(std::vector<double> upar,
                               bool jacobian_adjust_transform,
                               bool gradient) {
    if (upar.size() != model_->num_params_r()) {
      std::stringstream msg;
      msg << "Number of unconstrained parameters does not match "
      "that of the model ("
      << upar.size() << " vs "
      << model_->num_params_r()
      << ").";
      throw std::domain_error(msg.str());
    }
    std::vector<int> par_i(model_->num_params_i(), 0);
    if (!gradient) {
      double lp = 0;
      if (jacobian_adjust_transform) {
        lp = model_->log_prob_jacobian(upar, par_i, &rstan::io::rcout);
      } else {
        lp = model_->log_prob_propto(upar, par_i, &rstan::io::rcout);
      }
      Rcpp::NumericVector lp2 = Rcpp::wrap(lp);
      return lp2;
    }

    std::vector<double> grad;
    double lp = jacobian_adjust_transform ?
      stan::model::log_prob_grad<true,true >(*model_, upar, par_i, grad, &rstan::io::rcout) :
      stan::model::log_prob_grad<true,false>(*model_, upar, par_i, grad, &rstan::io::rcout);
    Rcpp::NumericVector lp2 = Rcpp::wrap(lp);
    lp2.attr("gradient") = grad;
    return lp2;
  }

  /**
   * Expose the grad_log_prob of the model to stan_fit so R user
   * can call this function.
   *
   * @param upar The real parameters on the unconstrained
   *  space.
   * @param jacobian_adjust_transform A flag to indicate whether
   *   the Jacobian adjustment is included
   * @return A numeric vector whose size is equal to the size of upar
   */
  Rcpp::NumericVector stan_fit::grad_log_prob(std::vector<double> upar,
                                    bool jacobian_adjust_transform) {
    if (upar.size() != model_->num_params_r()) {
      std::stringstream msg;
      msg << "Number of unconstrained parameters does not match "
      "that of the model ("
      << upar.size() << " vs "
      << model_->num_params_r()
      << ").";
      throw std::domain_error(msg.str());
    }
    std::vector<int> par_i(model_->num_params_i(), 0);
    std::vector<double> gradient;
    double lp = jacobian_adjust_transform ?
      stan::model::log_prob_grad<true,true >(*model_, upar, par_i, gradient, &rstan::io::rcout) :
      stan::model::log_prob_grad<true,false>(*model_, upar, par_i, gradient, &rstan::io::rcout);
    Rcpp::NumericVector grad = Rcpp::wrap(gradient);
    grad.attr("log_prob") = lp;
    return grad;
  }

  /**
   * Return the number of unconstrained parameters
   */
  double stan_fit::num_pars_unconstrained() {
    return static_cast<double>(model_->num_params_r());
  }

  /**
   * Drive the sampler / optimizer / approximator
   *
   * @param args_ A R(cpp) list of arguments
   * @return A R(cpp) list of fit stuff
   */

  Rcpp::List stan_fit::call_sampler(Rcpp::List args_) {
    stan_args args(args_);
    Rcpp::List holder;

    int ret = command(args, model_, holder, names_oi_tidx_,
                      fnames_oi_, base_rng);
    holder.attr("return_code") = ret;
    return holder;
  }

  /**
   * Drive the generated quantities
   *
   * @param draws A matrix of posterior draws
   * @param seed An unsigned integer to seed the PRNG
   * @return A R(cpp) list of realizations from generated quantities
   */

  Rcpp::List stan_fit::standalone_gqs(const Eigen::Map<Eigen::MatrixXd> draws,
                            unsigned int seed) {
    Rcpp::List holder;

    R_CheckUserInterrupt_Functor interrupt;
    stan::callbacks::stream_logger logger(Rcpp::Rcout, Rcpp::Rcout, Rcpp::Rcout,
                                          rstan::io::rcerr, rstan::io::rcerr);

    std::unique_ptr<rstan_sample_writer> sample_writer_ptr;
    std::fstream sample_stream;
    std::stringstream comment_stream;
    std::vector<std::string> all_names;
    model_->constrained_param_names(all_names, true, true);
    std::vector<std::string> some_names;
    model_->constrained_param_names(some_names, true, false);
    int gq_size = all_names.size() - some_names.size();
    std::vector<size_t> gq_idx(gq_size);
    for (int i = 0; i < gq_size; i++) gq_idx[i] = i;
    sample_writer_ptr.reset(sample_writer_factory(&sample_stream,
                                                  comment_stream, "# ",
                                                  0, 0,
                                                  gq_size,
                                                  draws.rows(), 0,
                                                  gq_idx).release());

    int ret = stan::services::error_codes::CONFIG;
    ret = stan::services::standalone_generate(*model_, draws,
                                              seed, interrupt,
                                              logger, *sample_writer_ptr);

    holder = Rcpp::List(sample_writer_ptr->values_.x().begin(),
                        sample_writer_ptr->values_.x().end());

    return holder;
  }

  /**
   * Return names (of interest)
   *
   * @return A standard vector of standard strings
   */
  std::vector<std::string> stan_fit::param_names() const {
    return names_;
  }

  std::vector<std::string> stan_fit::param_names_oi() const {
    return names_oi_;
  }

  /**
   * Return the indices among those parameters of interest,
   * rather than all the parameters
   *
   * @param names A standard vector of standard strings naming POIs
   * @return A R(cpp) list of indices thereof
   */
  Rcpp::List stan_fit::param_oi_tidx(const std::vector<std::string>& names) {
    std::vector<std::string> names2;
    std::vector<std::vector<unsigned int> > indexes;
    for (auto& it : names) {
      if (is_flatname(it)) { // an element of an array
        size_t ts = std::distance(fnames_oi_.begin(),
                                  std::find(fnames_oi_.begin(),
                                            fnames_oi_.end(), it));
        if (ts == fnames_oi_.size()) // not found
          continue;
        names2.push_back(it);
        indexes.push_back(std::vector<unsigned int>(1, ts));
        continue;
      }
      size_t j = std::distance(names_oi_.begin(),
                               std::find(names_oi_.begin(),
                                         names_oi_.end(), it));
      if (j == names_oi_.size()) // not found
        continue;
      unsigned int j_size = calc_num_params(dims_oi_[j]);
      unsigned int j_start = starts_oi_[j];
      std::vector<unsigned int> j_idx;
      for (unsigned int k = 0; k < j_size; k++) {
        j_idx.push_back(j_start + k);
      }
      names2.push_back(it);
      indexes.push_back(j_idx);
    }
    Rcpp::List lst = Rcpp::wrap(indexes);
    lst.names() = names2;
    return lst;
  }

  /**
   * Get dimensions
   *
   * @return A R(cpp) list of dimensions (of interest)
   */

  Rcpp::List stan_fit::param_dims() const {
    Rcpp::List lst = Rcpp::wrap(dims_);
    lst.names() = names_;
    return lst;
  }

  Rcpp::List stan_fit::param_dims_oi() const {
    Rcpp::List lst = Rcpp::wrap(dims_oi_);
    lst.names() = names_oi_;
    return lst;
  }

  /**
   * Get flatnames of interest
   *
   * @return A standard vector of standard strings of FOIs
   */

  std::vector<std::string> stan_fit::param_fnames_oi() const {
    std::vector<std::string> fnames;
    get_all_flatnames(names_oi_, dims_oi_, fnames, true);
    return fnames;
  }

}
