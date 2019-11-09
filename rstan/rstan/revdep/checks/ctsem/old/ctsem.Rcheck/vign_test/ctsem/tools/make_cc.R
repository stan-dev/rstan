# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
options(warn = 3L)
options("useFancyQuotes" = FALSE)

make_cc <- function(file) {
  file <- sub("\\.cc$", ".stan", file)
  cppcode <- rstan::stanc(file, allow_undefined = TRUE,
                          obfuscate_model_name = FALSE)$cppcode
  #disabled include meta header
  # cppcode <- sub("(class[[:space:]]+[A-Za-z_][A-Za-z0-9_]*[[:space:]]*: public prob_grad \\{)",
                 # paste("#include <meta_header.hpp>\n", "\\1"), cppcode)
# 
  # w32 <- .Machine$sizeof.pointer == 4
  stan_files<-paste0('stan_files',ifelse(.Machine$sizeof.pointer == 4,'32',''))
  
  cat(readLines(dir(stan_files, pattern = "license.stan", recursive = TRUE, full.names = TRUE)),
      "#ifndef MODELS_HPP", "#define MODELS_HPP", "#define STAN__SERVICES__COMMAND_HPP",
      "#include <rstan/rstaninc.hpp>",
      cppcode, "#endif", file = sub("\\.stan$", ".hpp", file),
      sep = "\n", append = FALSE)
  
  f <- sub("\\.stan$", "", basename(file))
  Rcpp::exposeClass(class = paste0("model_", f),
                    constructors = list(c("SEXP", "SEXP", "SEXP")), fields = character(),
                    methods = c("call_sampler", 
                                "param_names", "param_names_oi", "param_fnames_oi", 
                                "param_dims",  "param_dims_oi", "update_param_oi", "param_oi_tidx", 
                                "grad_log_prob", "log_prob", 
                                "unconstrain_pars", "constrain_pars", "num_pars_unconstrained", 
                                "unconstrained_param_names", "constrained_param_names"), 
                    file = file.path(stan_files, paste0(f, ".cc")), 
                    header = paste0('#include "', f, '.hpp"'),
                    module = paste0("stan_fit4", f, "_mod"), 
                    CppClass = "rstan::stan_fit<stan_model, boost::random::ecuyer1988> ",
                    Rfile = FALSE)
  return(invisible(NULL))
}
