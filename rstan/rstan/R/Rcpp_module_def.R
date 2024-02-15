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
RCPP_MODULE(stan_fit4%model_name%_mod) {
  class_<rstan::stan_fit<stan_model, boost::random::mixmax> >(
      "stan_fit4%model_name%")

      .constructor<SEXP, SEXP, SEXP>()

      .method(
          "call_sampler",
          &rstan::stan_fit<stan_model, boost::random::mixmax>::call_sampler)
      .method(
          "param_names",
          &rstan::stan_fit<stan_model, boost::random::mixmax>::param_names)
      .method("param_names_oi",
              &rstan::stan_fit<stan_model,
                               boost::random::mixmax>::param_names_oi)
      .method("param_fnames_oi",
              &rstan::stan_fit<stan_model,
                               boost::random::mixmax>::param_fnames_oi)
      .method(
          "param_dims",
          &rstan::stan_fit<stan_model, boost::random::mixmax>::param_dims)
      .method("param_dims_oi",
              &rstan::stan_fit<stan_model,
                               boost::random::mixmax>::param_dims_oi)
      .method("update_param_oi",
              &rstan::stan_fit<stan_model,
                               boost::random::mixmax>::update_param_oi)
      .method("param_oi_tidx",
              &rstan::stan_fit<stan_model,
                               boost::random::mixmax>::param_oi_tidx)
      .method("grad_log_prob",
              &rstan::stan_fit<stan_model,
                               boost::random::mixmax>::grad_log_prob)
      .method("log_prob",
              &rstan::stan_fit<stan_model, boost::random::mixmax>::log_prob)
      .method("unconstrain_pars",
              &rstan::stan_fit<stan_model,
                               boost::random::mixmax>::unconstrain_pars)
      .method("constrain_pars",
              &rstan::stan_fit<stan_model,
                               boost::random::mixmax>::constrain_pars)
      .method(
          "num_pars_unconstrained",
          &rstan::stan_fit<stan_model,
                           boost::random::mixmax>::num_pars_unconstrained)
      .method(
          "unconstrained_param_names",
          &rstan::stan_fit<
              stan_model, boost::random::mixmax>::unconstrained_param_names)
      .method(
          "constrained_param_names",
          &rstan::stan_fit<stan_model,
                           boost::random::mixmax>::constrained_param_names)
      .method("standalone_gqs",
              &rstan::stan_fit<stan_model,
                               boost::random::mixmax>::standalone_gqs);
}
'
gsub("%model_name%", model_name, RCPP_MODULE)
}
