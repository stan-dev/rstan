#include <rstan/io/rlist_ref_var_context.hpp>
#include <rstan/boost_random_R.hpp>
#include <stan/math/prim/mat.hpp>
#include <stan/io/array_var_context.hpp>
#include <stan/io/dump.hpp>
#include <stan/io/reader.hpp>
#include <stan/io/writer.hpp>
#include <stan/lang/rethrow_located.hpp>
#include <stan/model/indexing.hpp>
#include <boost/exception/all.hpp>
#include <Rcpp.h>
#include <vector>
#include <string>

class standalone_generate {
protected:
  const bool propto__ = true;
  const bool jacobian__ = false;
  const boost_random_R base_rng__ = boost_random_R();
  size_t num_params_r__ = 0;
  std::vector<std::pair<int, int> > param_ranges_i__ =
    std::vector<std::pair<int, int> >(0);
  std::ostream* pstream__ = &Rcpp::Rcout;
  
public:
  virtual void transform_inits(const stan::io::var_context& context__,
                               std::vector<int>& params_i__,
                               std::vector<double>& params_r__,
                               std::ostream* pstream__) = 0;
  virtual void write_array( // RNG& base_rng__,
      std::vector<double>& params_r__,
      std::vector<int>& params_i__,
      std::vector<double>& vars__,
      bool include_tparams__ = true,
      bool include_gqs__ = true,
      std::ostream* pstream__ = 0) = 0;
  virtual void constrained_param_names(std::vector<std::string>& param_names__,
                               bool include_tparams__ = true,
                               bool include_gqs__ = true) = 0;
    
  virtual void get_dims(std::vector<std::vector<size_t> >& dimss__) = 0;
  virtual ~standalone_generate() { }
    
  Rcpp::NumericMatrix gqs(Rcpp::NumericMatrix draws) {
    
    std::vector<std::string> param_names;
    constrained_param_names(param_names, false, true);
    size_t upper_bound = param_names.size();
    param_names.clear();
    constrained_param_names(param_names, false, false);
    size_t lower_bound = param_names.size();
    size_t num_gqs = upper_bound - lower_bound;
    std::vector<std::vector<size_t> > param_dimss;
    get_dims(param_dimss); // does this include lp__?
    param_dimss.erase(param_dimss.begin() + lower_bound, param_dimss.end());
    
    if (lower_bound != draws.cols()) {
      throw std::runtime_error("draws has the wrong number of columns");
    }
    
    std::vector<int> dummy_params_i;
    std::vector<double> unconstrained_params_r;
    std::vector<double> gqs;
    std::vector<double> draws_i(draws.cols());
    std::stringstream msg;
    Rcpp::NumericMatrix output(draws.rows(), num_gqs);
    for (size_t i = 0; i < draws.rows(); ++i) {
      dummy_params_i.clear();
      unconstrained_params_r.clear();
      for (size_t j = 0; j < draws_i.size(); j++) 
        draws_i[j] = draws(i, j);
      try {
        stan::io::array_var_context context(param_names, draws_i,
                                            param_dimss);
        transform_inits(context, dummy_params_i, unconstrained_params_r,
                        &msg);
      } catch (const std::exception& e) {
        throw std::runtime_error(e.what());
      }
      // call out to interrupt and fail
      write_array(unconstrained_params_r, dummy_params_i, gqs, 
                  false, true, pstream__);
      for (size_t j = 0; j < num_gqs; j++)
        output(i, j) = gqs[j];
    }
    return output;
  }
};
