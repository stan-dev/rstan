/**
 * Define Rcpp Module to expose stan_fit's functions to R.
 */

/*
RCPP_MODULE(stan_fit4%model_name%_mod){
  Rcpp::class_<rstan::stan_fit<stan_model,
               boost::random::mixmax> >("stan_fit4%model_name%")
    // .constructor<Rcpp::List>()
    .constructor<SEXP, SEXP, SEXP>()
    // .constructor<SEXP, SEXP>()
    .method("call_sampler", &rstan::stan_fit<stan_model, boost::random::mixmax>::call_sampler)
    .method("param_names", &rstan::stan_fit<stan_model, boost::random::mixmax>::param_names)
    .method("param_names_oi", &rstan::stan_fit<stan_model, boost::random::mixmax>::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<stan_model, boost::random::mixmax>::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<stan_model, boost::random::mixmax>::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<stan_model, boost::random::mixmax>::param_dims_oi)
    .method("update_param_oi",&rstan::stan_fit<stan_model, boost::random::mixmax>::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<stan_model, boost::random::mixmax>::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<stan_model, boost::random::mixmax>::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<stan_model, boost::random::mixmax>::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<stan_model, boost::random::mixmax>::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<stan_model, boost::random::mixmax>::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<stan_model, boost::random::mixmax>::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<stan_model, boost::random::mixmax>::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<stan_model, boost::random::mixmax>::constrained_param_names)
    .method("standalone_gqs", &rstan::stan_fit<stan_model, boost::random::mixmax>::standalone_gqs)
  ;
}
*/
