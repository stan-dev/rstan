#ifndef STAN_FIT_BASE_HPP
#define STAN_FIT_BASE_HPP

#include <vector>
#include <stan/math/prim/mat/fun/Eigen.hpp>

#include <Rcpp.h>
#include <RcppEigen.h>

namespace rstan {

class stan_fit_base {
 public:
  virtual ~stan_fit_base();
  virtual bool update_param_oi(std::vector<std::string> pnames) = 0;
  virtual std::vector<double> unconstrain_pars(Rcpp::List par) = 0;
  virtual std::vector<double> constrain_pars(std::vector<double> upar) = 0;
  virtual std::vector<std::string> unconstrained_param_names(bool include_tparams,
                                                             bool include_gqs) = 0;
  virtual std::vector<std::string> constrained_param_names(bool include_tparams,
                                                           bool include_gqs) = 0;
  virtual Rcpp::NumericVector log_prob(std::vector<double> upar,
                                       bool jacobian_adjust_transform,
                                       bool gradient) = 0;
  virtual Rcpp::NumericVector grad_log_prob(std::vector<double> upar,
                                            bool jacobian_adjust_transform) = 0;
  virtual double num_pars_unconstrained() = 0;
  virtual Rcpp::List call_sampler(Rcpp::List args_) = 0;
  virtual Rcpp::List standalone_gqs(const Eigen::Map<Eigen::MatrixXd> draws,
                                    unsigned int seed) = 0;
  virtual std::vector<std::string> param_names() const = 0;
  virtual std::vector<std::string> param_names_oi() const = 0;
  virtual Rcpp::List param_oi_tidx(const std::vector<std::string>& names) = 0;
  virtual Rcpp::List param_dims() const = 0;
  virtual Rcpp::List param_dims_oi() const = 0;
  virtual std::vector<std::string> param_fnames_oi() const = 0;
};

}

#endif
