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

get_Rcpp_module_def_code <- function(model_name) {
  RCPP_MODULE <-
'
#include <rstan_next/stan_fit.hpp>

struct stan_model_holder {
    stan_model_holder(rstan::io::rlist_ref_var_context rcontext,
                      unsigned int random_seed)
    : rcontext_(rcontext), random_seed_(random_seed),
      model_raw_(new stan_model(rcontext_, random_seed_)),
      model_(model_raw_, false)
     {
       R_PreserveObject(model_);
     }

   //stan::math::ChainableStack ad_stack;
   rstan::io::rlist_ref_var_context rcontext_;
   unsigned int random_seed_;
   stan::model::model_base* model_raw_;
   Rcpp::XPtr<stan::model::model_base> model_;
};

Rcpp::XPtr<stan::model::model_base> model_ptr(stan_model_holder* smh) {
  return smh->model_;
}

Rcpp::XPtr<rstan::stan_fit_base> fit_ptr(stan_model_holder* smh) {
  return Rcpp::XPtr<rstan::stan_fit_base>(new rstan::stan_fit(smh->model_, smh->random_seed_), true);
}

std::string model_name(stan_model_holder* smh) {
  return smh->model_.get()->model_name();
}

void finalize_stan_model_holder(stan_model_holder* smh) {
//   R_ReleaseObject(smh->model_raw_);
   delete smh->model_raw_;
}

RCPP_MODULE(stan_fit4%model_name%_mod){
  Rcpp::class_<stan_model_holder>("stan_fit4%model_name%")
  .constructor<rstan::io::rlist_ref_var_context, unsigned int>()
  .method("model_ptr", &model_ptr)
  .method("fit_ptr", &fit_ptr)
  .method("model_name", &model_name)
  .finalizer(&finalize_stan_model_holder)
  ;
}
'
gsub("%model_name%", model_name, RCPP_MODULE)
}
