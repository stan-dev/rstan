#ifndef RSTAN_FIT_HPP
#define RSTAN_FIT_HPP


#include "stan_fit_base.hpp"

#include <stan/model/model_base.hpp>

namespace rstan {

class stan_fit : public stan_fit_base {
  
private:
  stan::model::model_base* model_;
  boost::ecuyer1988 base_rng;
  const std::vector<std::string> names_;
  const std::vector<std::vector<unsigned int> > dims_;
  const unsigned int num_params_;
  
  std::vector<std::string> names_oi_;                // parameters of interest
  std::vector<std::vector<unsigned int> > dims_oi_;  // and their dimensions
  std::vector<size_t> names_oi_tidx_;                // total indexes of names2
  std::vector<unsigned int> starts_oi_;              // do not know what this is
  unsigned int num_params2_;                         // total number of POI's.
  std::vector<std::string> fnames_oi_;               // flatnames of
                                                     // interest

  stan::math::ChainableStack ad_stack_;

private:
  /**
   * Tell if a parameter name is an element of an array parameter.
   * Note that it only supports full specified name; slicing
   * is not supported. The test only tries to see if there
   * are brackets.
   */
  bool is_flatname(const std::string& name);
  /*
   * Update the parameters we are interested for the model->
   * As well, the dimensions vector for the parameters are
   * updated.
   */
  void update_param_oi0(const std::vector<std::string>& pnames);
  
public:
  bool update_param_oi(std::vector<std::string> pnames);
  
  stan_fit(Rcpp::XPtr<stan::model::model_base> model, int seed);
  
  /**
   * Transform the parameters from its defined support
   * to unconstrained space
   *
   * @param par An R list as for specifying the initial values
   *  for a chain
   * @return A standard vector of doubles
   */
  std::vector<double> unconstrain_pars(Rcpp::List par);
  
  /**
   * Contrary to unconstrain_pars, transform parameters
   * from unconstrained support to the constrained.
   *
   * @param upar A standard vector of doubles on the unconstrained space
   * @return A standard vector of doubles on the constrained space
   */
  std::vector<double> constrain_pars(std::vector<double> upar);
  
  /**
   * Get the unconstrained or constrained parameter names
   * 
   * @param include_tparams Flag to include transformed parameter names
   * @param include_gqs Flag to include generated quantitiy names
   * @return A standard vector of standard strings
   */
  
  std::vector<std::string> unconstrained_param_names(bool include_tparams, 
                                                     bool include_gqs);
  
  std::vector<std::string> constrained_param_names(bool include_tparams, 
                                                   bool include_gqs);
  
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
  Rcpp::NumericVector log_prob(std::vector<double> upar, 
                               bool jacobian_adjust_transform, 
                               bool gradient);
  
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
  Rcpp::NumericVector grad_log_prob(std::vector<double> upar, 
                                    bool jacobian_adjust_transform);
  
  /**
   * Return the number of unconstrained parameters
   */
  int num_pars_unconstrained();

  /**
   * Drive the sampler / optimizer / approximator
   * 
   * @param args_ A R(cpp) list of arguments
   * @return A R(cpp) list of fit stuff
   */
  
  Rcpp::List call_sampler(Rcpp::List args_);
  
  /**
   * Drive the generated quantities
   * 
   * @param draws A matrix of posterior draws
   * @param seed An unsigned integer to seed the PRNG 
   * @return A R(cpp) list of realizations from generated quantities
   */
  
  Rcpp::List standalone_gqs(const Eigen::Map<Eigen::MatrixXd> draws, 
                            unsigned int seed);
  
  /**
   * Return names (of interest)
   * 
   * @return A standard vector of standard strings
   */
  std::vector<std::string> param_names() const;
  
  std::vector<std::string> param_names_oi() const;
  
  /**
   * Return the indices among those parameters of interest, 
   * rather than all the parameters
   * 
   * @param names A standard vector of standard strings naming POIs
   * @return A R(cpp) list of indices thereof
   */
  Rcpp::List param_oi_tidx(std::vector<std::string> names);
  
  /**
   * Get dimensions
   * 
   * @return A R(cpp) list of dimensions (of interest)
   */
  
  Rcpp::List param_dims() const;
  
  Rcpp::List param_dims_oi() const;

  /**
   * Get flatnames of interest
   * 
   * @return A standard vector of standard strings of FOIs
   */
  
  std::vector<std::string> param_fnames_oi() const;
};


}

#endif
