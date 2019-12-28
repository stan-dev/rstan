#include <Rcpp.h>
#include <stan/model/model_base.hpp>
#include <stan/services/util/create_rng.hpp>
#include <rstan/io/rlist_ref_var_context.hpp>
#include <rstan/io/r_ostream.hpp>
#include <rstan/stan_args.hpp>

/*
RCPP_MODULE(class_model_base) {
    using namespace Rcpp ;

    class_<stan::model::model_base>("model_base")

    .method("model_name", &stan::model::model_base::model_name)
    ;
}
*/

inline std::vector<std::string> 
get_param_names(stan::model::model_base* user_model) {
  std::vector<std::string> names;
  user_model->get_param_names(names);
  return names;
}

inline Rcpp::List
get_dims(stan::model::model_base* user_model) {
  std::vector<std::vector<size_t> > dimss;
  user_model->get_dims(dimss);
  dimss.push_back({});
  Rcpp::List lst(dimss.begin(), dimss.end());
  std::vector<std::string> names;
  user_model->get_param_names(names);
  names.push_back("lp__");
  lst.names() = names;
  return lst;
}

inline std::vector<std::string>
constrained_param_names(stan::model::model_base* user_model,
                        bool include_tparams = true,
                        bool include_gqs = true) {
  std::vector<std::string> param_names;
  user_model->constrained_param_names(param_names, 
                                      include_tparams, include_gqs);
  return param_names;
}

inline std::vector<std::string>
unconstrained_param_names(stan::model::model_base* user_model,
                          bool include_tparams = true,
                          bool include_gqs = true) {
  std::vector<std::string> param_names;
  user_model->unconstrained_param_names(param_names, 
                                        include_tparams, include_gqs);
  return param_names;
}

inline double
log_prob(stan::model::model_base* user_model, 
         std::vector<double>& params_r) {
  std::vector<int> params_i;
  return user_model->log_prob(params_r, params_i, &Rcpp::Rcout);
}

inline double
log_prob_jacobian(stan::model::model_base* user_model, 
                  std::vector<double>& params_r) {
  std::vector<int> params_i;
  return user_model->log_prob_jacobian(params_r, params_i, &Rcpp::Rcout);
}

inline double // this is always 0 so should we not expose it
log_prob_propto(stan::model::model_base* user_model, 
                std::vector<double>& params_r) {
  std::vector<int> params_i;
  return user_model->log_prob_propto(params_r, params_i, &Rcpp::Rcout);
}

inline double
log_prob_propto_jacobian(stan::model::model_base* user_model, 
                         std::vector<double>& params_r) {
  std::vector<int> params_i;
  return user_model->log_prob_propto_jacobian(params_r, params_i, 
                                              &Rcpp::Rcout);
}

inline std::vector<double>
transform_inits(stan::model::model_base* user_model,
                const rstan::io::rlist_ref_var_context context) {
  std::vector<int> params_i;
  std::vector<double> params_r;
  user_model->transform_inits(context, params_i, params_r, 
                              &Rcpp::Rcout);
  return params_r;
}

inline std::vector<double>
write_array(stan::model::model_base* user_model,
            std::vector<double>& params_r,
            bool include_tparams = true, bool include_gqs = true,
            unsigned int random_seed = 0, unsigned int id = 0) {
  std::vector<int> params_i;
  std::vector<double> constrained_params;
  boost::ecuyer1988 rng = stan::services::util::create_rng(random_seed, id);
  user_model->write_array(rng, params_r, params_i, constrained_params, 
                          include_tparams, include_gqs, &Rcpp::Rcout);
  return constrained_params;
}

inline Rcpp::List
write_list(stan::model::model_base* user_model,
           std::vector<double>& params_r,
           bool include_tparams = true, bool include_gqs = true,
             unsigned int random_seed = 0, unsigned int id = 0) {
  std::vector<double> params = 
    write_array(user_model, params_r, include_tparams,
                include_gqs, random_seed, id);
  std::vector<std::vector<size_t> > dims;
  user_model->get_dims(dims);
  unsigned int K = dims.size();
  Rcpp::List out(K);
  unsigned int pos = 0;
  for (unsigned int k = 0; k < K; k++) {
    if (dims[k].empty()) {
      out[k] = params[pos++];
    } else {
      std::vector<size_t> d = dims[k];
      size_t p = std::accumulate(d.begin(), d.end(),
                                 1, std::multiplies<size_t>());
      Rcpp::NumericVector x(p);
      for (size_t j = 0; j < p; j++) {
        x[j] = params[pos++];
      }
      x.attr("dim") = d;
      out[k] = x;
    }
  }
  std::vector<std::string> names;
  user_model->get_param_names(names);
  out.names() = names;
  return out;
}


inline stan::model::model_base*
new_model(Rcpp::XPtr<stan::model::model_base> user_model) {
  return user_model;
}

RCPP_MODULE(class_model_base) {
  using namespace Rcpp ;
  
  class_<stan::model::model_base>("model_base")
  .factory<Rcpp::XPtr<stan::model::model_base> >(new_model)
  .method("model_name", &stan::model::model_base::model_name,
          "takes no arguments and returns the MD5 hashed model")
  .method("get_param_names", &get_param_names,
          "takes no arguments and returns a character vector of "
          "parameter names")
  .method("get_dims", &get_dims,
          "takes no arguments and returns a list of dimensions")
  .method("constrained_param_names", &constrained_param_names,
          "takes flags for include_tparams and include_gqs and "
          "returns a character vector of names of unknown quantities")
  .method("unconstrained_param_names", &unconstrained_param_names,
          "takes flags for include_tparams and include_gqs and "
          "returns a character vector of names of unknown quantities "
          "in the unconstrained space")
  .method("log_prob", &log_prob,
          "takes a numeric vector of parameters and returns the "
          "log of the unnormalized density with constants but "
          "without a Jacobian correction")
  .method("log_prob_jacobian", &log_prob_jacobian,
          "takes a numeric vector of parameters and returns the "
          "log of the unnormalized density with constants and "
          "a Jacobian correction")
  .method("log_prob_propto", &log_prob_propto,
          "takes a numeric vector of parameters and returns the "
          "log of the unnormalized density without constants and "
          "without a Jacobian correction")
  .method("log_prob_propto_jacobian", &log_prob_propto_jacobian,
          "takes a numeric vector of parameters and returns the "
          "log of the unnormalized density without constants but "
          "with a Jacobian correction")
  .method("transform_inits", &transform_inits,
          "takes a list of constrained parameters and returns a "
          "numeric vector of unconstrained parameters")
  .method("write_array", &write_array,
          "takes a vector of unconstrained parameters, flags for "
          "include_tparams and include_gqs, as well as integers "
          "for id and seed and returns a vector of constrained "
          "parameters")
  .method("write_list", &write_list,
          "takes a vector of unconstrained parameters, flags for "
          "include_tparams and include_gqs, as well as integers "
          "for id and seed and returns a list of constrained "
          "parameters with the appropriate dimensions")
  ;
}

