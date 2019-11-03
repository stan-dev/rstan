# This file is part of RStan
# Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019 Trustees of Columbia University
#
# RStan is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# RStan is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# get_Rcpp_module_def_code <- function(model_name) {
#   def_Rcpp_module_hpp_file <-
#     system.file('include', '/rstan/rcpp_module_def_for_rstan.hpp', package = 'rstan')
#   if (def_Rcpp_module_hpp_file == '')
#     stop("Rcpp module definition file for rstan is not found.\n")
#   src <- paste(readLines(def_Rcpp_module_hpp_file), collapse = '\n')
#   gsub("%model_name%", model_name, src)
# }

get_Rcpp_module_def_code <- function(model_name) {
  RCPP_MODULE <-
'

#include <rstan_next/stan_fit.hpp>

struct stan_model_holder {
    stan_model_holder(rstan::io::rlist_ref_var_context rcontext,
                      unsigned int random_seed)
    : rcontext_(rcontext), random_seed_(random_seed),
      model_(new stan_model(rcontext_, random_seed_), true),
      fit_(new rstan::stan_fit(model_, random_seed_), true) {
    }

   //stan::math::ChainableStack ad_stack;
   rstan::io::rlist_ref_var_context rcontext_;
   unsigned int random_seed_;
   Rcpp::XPtr<stan::model::model_base> model_;
   Rcpp::XPtr<rstan::stan_fit_base> fit_;
};

Rcpp::XPtr<stan::model::model_base> model_ptr(stan_model_holder* smh) {
  //Rcpp::XPtr<stan::model::model_base> ptr(&static_cast<stan::model::model_base&>(*sm), false);
  //Rcpp::XPtr<stan_model> ptr(sm, false);
  return smh->model_;
}

Rcpp::XPtr<rstan::stan_fit_base> fit_ptr(stan_model_holder* smh) {
  //Rcpp::XPtr<stan::model::model_base> ptr(&static_cast<stan::model::model_base&>(*sm), false);
  //Rcpp::XPtr<stan_model> ptr(sm, false);
  return smh->fit_;
}

std::string model_name(stan_model_holder* smh) {
  return smh->model_.get()->model_name();
}


/*

Rcpp::XPtr<stan_model> ptr(stan_model* sm) {
  Rcpp::XPtr<stan_model> ptr(sm, true);
  //Rcpp::XPtr<stan_model> ptr(sm, false);
  return ptr;
}

Rcpp::XPtr<stan_model> ptr(SEXP rsm) {
  Rcpp::class_<stan_model> mod(r
  Rcpp::XPtr<stan_model> ptr(sm, true);
  //Rcpp::XPtr<stan_model> ptr(sm, false);
  return ptr;
}


Rcpp::XPtr<stan::model::model_base>
make_model(const rstan::io::rlist_ref_var_context context,
           unsigned int random_seed) {
   Rcpp::XPtr<stan::model::model_base> ptr(new stan_model(context, random_seed));
   return ptr;
}
*/

RCPP_MODULE(stan_fit4%model_name%_mod){
  Rcpp::class_<stan_model_holder>("stan_fit4%model_name%")
  .constructor<rstan::io::rlist_ref_var_context, unsigned int>()
  //.constructor<rstan::io::rlist_ref_var_context&, unsigned int>()
  .method("model_ptr", &model_ptr)
  .method("fit_ptr", &fit_ptr)
  .method("model_name", &model_name)
//  .method("log_prob", &stan_model::log_prob)
//  .method("write_array", &stan_model::write_array)
  ;
}
'
gsub("%model_name%", model_name, RCPP_MODULE)
}

