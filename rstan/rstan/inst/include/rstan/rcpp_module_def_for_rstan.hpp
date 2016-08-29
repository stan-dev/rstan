#include <Rcpp.h>
#include <rstan/stan_model_with_data.hpp>

/**
 * Define Rcpp Module to expose stan_fit's functions to R.
 */
RCPP_MODULE(stan_fit4%model_name%_mod){
  Rcpp::class_<rstan::stan_fit<%model_name%_namespace::%model_name%> >("stan_fit4%model_name%")
               .constructor<SEXP, SEXP>()
    .method("diagnose",
            &rstan::stan_fit<%model_name%_namespace::%model_name%>::diagnose);
  // .constructor<Rcpp::List>()
  // // .constructor<SEXP, SEXP>()
  // .method("call_sampler",
  //         &rstan::stan_fit<%model_name%_namespace::%model_name%, boost::random::ecuyer1988>::call_sampler)
  // .method("param_names",
  //         &rstan::stan_fit<%model_name%_namespace::%model_name%, boost::random::ecuyer1988>::param_names)
  // .method("param_names_oi",
  //         &rstan::stan_fit<%model_name%_namespace::%model_name%, boost::random::ecuyer1988>::param_names_oi)
  // .method("param_fnames_oi",
  //         &rstan::stan_fit<%model_name%_namespace::%model_name%, boost::random::ecuyer1988>::param_fnames_oi)
  // .method("param_dims",
  //         &rstan::stan_fit<%model_name%_namespace::%model_name%, boost::random::ecuyer1988>::param_dims)
  // .method("param_dims_oi",
  //         &rstan::stan_fit<%model_name%_namespace::%model_name%, boost::random::ecuyer1988>::param_dims_oi)
  // .method("update_param_oi",
  //         &rstan::stan_fit<%model_name%_namespace::%model_name%, boost::random::ecuyer1988>::update_param_oi)
  // .method("param_oi_tidx",
  //         &rstan::stan_fit<%model_name%_namespace::%model_name%, boost::random::ecuyer1988>::param_oi_tidx)
  // .method("grad_log_prob",
  //         &rstan::stan_fit<%model_name%_namespace::%model_name%, boost::random::ecuyer1988>::grad_log_prob)
  // .method("log_prob",
  //         &rstan::stan_fit<%model_name%_namespace::%model_name%, boost::random::ecuyer1988>::log_prob)
  // .method("unconstrain_pars",
  //         &rstan::stan_fit<%model_name%_namespace::%model_name%, boost::random::ecuyer1988>::unconstrain_pars)
  // .method("constrain_pars",
  //         &rstan::stan_fit<%model_name%_namespace::%model_name%, boost::random::ecuyer1988>::constrain_pars)
  // .method("num_pars_unconstrained",
  //         &rstan::stan_fit<%model_name%_namespace::%model_name%, boost::random::ecuyer1988>::num_pars_unconstrained)
  // .method("unconstrained_param_names",
  //         &rstan::stan_fit<%model_name%_namespace::%model_name%, boost::random::ecuyer1988>::unconstrained_param_names)
  // .method("constrained_param_names",
  //         &rstan::stan_fit<%model_name%_namespace::%model_name%, boost::random::ecuyer1988>::constrained_param_names)
}

RCPP_MODULE(stan_model4%model_name%_mod){
  Rcpp::class_<rstan::stan_model_with_data<%model_name%_namespace::%model_name%> >("stan_model4%model_name%")
  .constructor<SEXP>()
  .method("diagnose",
    &rstan::stan_model_with_data<%model_name%_namespace::%model_name%>::diagnose);
}
