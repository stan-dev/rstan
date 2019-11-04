#include <rstan_next/stan_fit_base.hpp>

namespace rstan {

class stan_fit_proxy : public stan_fit_base {

  stan_fit_base& fit_;

 public:
  stan_fit_proxy(Rcpp::XPtr<rstan::stan_fit_base> fit) : fit_(*fit.get()) {}
  ~stan_fit_proxy() {}


  bool update_param_oi(std::vector<std::string> pnames) {
    return fit_.update_param_oi(pnames);
  }
  
  std::vector<double> unconstrain_pars(Rcpp::List par) {
    return fit_.unconstrain_pars(par);
  }
  
  std::vector<double> constrain_pars(std::vector<double> upar) {
    return fit_.constrain_pars(upar);
  }
  
  std::vector<std::string> unconstrained_param_names(bool include_tparams, 
                                                     bool include_gqs) {
    return fit_.unconstrained_param_names(include_tparams, include_gqs);
  }
  
  std::vector<std::string> constrained_param_names(bool include_tparams, 
                                                   bool include_gqs) {
    return fit_.constrained_param_names(include_tparams, include_gqs);
  }
  
  Rcpp::NumericVector log_prob(std::vector<double> upar, 
                               bool jacobian_adjust_transform, 
                               bool gradient) {
    return fit_.log_prob(upar, jacobian_adjust_transform, gradient);
  }
  
  Rcpp::NumericVector grad_log_prob(std::vector<double> upar, 
                                    bool jacobian_adjust_transform) {
    return fit_.grad_log_prob(upar, jacobian_adjust_transform);
  }
  
  int num_pars_unconstrained() {
    return fit_.num_pars_unconstrained();
  }

  Rcpp::List call_sampler(Rcpp::List args_) {
    return fit_.call_sampler(args_);
  }
  
  Rcpp::List standalone_gqs(const Eigen::Map<Eigen::MatrixXd> draws, 
                            unsigned int seed) {
    return fit_.standalone_gqs(draws, seed);
  }
  
  std::vector<std::string> param_names() const {
    return fit_.param_names();
  }
  
  std::vector<std::string> param_names_oi() const {
    return fit_.param_names_oi();
  }
  
  Rcpp::List param_oi_tidx(std::vector<std::string> names) {
    return fit_.param_oi_tidx(names);
  }
  
  Rcpp::List param_dims() const {
    return fit_.param_dims();
  }
  
  Rcpp::List param_dims_oi() const {
    return fit_.param_dims_oi();
  }

  std::vector<std::string> param_fnames_oi() const {
    return fit_.param_fnames_oi();
  }
  
};

}
RCPP_MODULE(class_stan_fit){
  Rcpp::class_<rstan::stan_fit_proxy>("stan_fit")
      .constructor<Rcpp::XPtr<rstan::stan_fit_base>>()
      .method("update_param_oi",&rstan::stan_fit_proxy::update_param_oi)
      .method("unconstrain_pars", &rstan::stan_fit_proxy::unconstrain_pars)
      .method("constrain_pars", &rstan::stan_fit_proxy::constrain_pars)
      .method("unconstrained_param_names", &rstan::stan_fit_proxy::unconstrained_param_names)
      .method("constrained_param_names", &rstan::stan_fit_proxy::constrained_param_names)
      .method("log_prob", &rstan::stan_fit_proxy::log_prob)
      .method("grad_log_prob", &rstan::stan_fit_proxy::grad_log_prob)
      .method("num_pars_unconstrained", &rstan::stan_fit_proxy::num_pars_unconstrained)
      .method("call_sampler", &rstan::stan_fit_proxy::call_sampler)
      .method("standalone_gqs", &rstan::stan_fit_proxy::standalone_gqs)
      .method("param_names", &rstan::stan_fit_proxy::param_names)
      .method("param_names_oi", &rstan::stan_fit_proxy::param_names_oi)
      .method("param_oi_tidx", &rstan::stan_fit_proxy::param_oi_tidx)
      .method("param_dims", &rstan::stan_fit_proxy::param_dims)
      .method("param_dims_oi", &rstan::stan_fit_proxy::param_dims_oi)
      .method("param_fnames_oi", &rstan::stan_fit_proxy::param_fnames_oi)
  ;
}
