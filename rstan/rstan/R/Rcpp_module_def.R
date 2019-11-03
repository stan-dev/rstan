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

auto self_ptr(stan_model* sm) {
  Rcpp::XPtr<stan_model> ptr(sm, true);
  return ptr;  
}

double log_prob(stan_model* sm, std::vector<double>& params_r__) {
  std::vector<int> params_i__;
  return sm->log_prob<false, false, double>(params_r__,
                                            params_i__, &Rcpp::Rcout);
}
  
RCPP_MODULE(stan_fit4%model_name%_mod){
  Rcpp::class_<stan_model>("stan_fit4%model_name%")
  .constructor<const rstan::io::rlist_ref_var_context, unsigned int>()
  .method("this", &self_ptr)
  .method("log_prob", &log_prob)
// .method("log_prob", &stan_model::log_prob<false, false, double>)
//  .method("write_array", &stan_model::write_array)
  ;
}
'
gsub("%model_name%", model_name, RCPP_MODULE)
}

get_Rcpp_module_def_code <- function(model_name) {
  
  code <-
    '
// [[Rcpp::export]]
inline
Rcpp::XPtr<stan_model>
new_model(const rstan::io::rlist_ref_var_context context__,
          const unsigned int seed) {
  Rcpp::XPtr<stan_model> ptr(new stan_model(context__, seed, &Rcpp::Rcout), true);
  return ptr;
}
'
return(code)
}
