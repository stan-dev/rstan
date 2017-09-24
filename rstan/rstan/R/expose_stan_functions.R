# This file is part of RStan
# Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017 Trustees of Columbia University
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

expose_stan_functions_hacks <- function(code) {
  code <- sub("_header.hpp>\n", "_header.hpp>\n#include <stan/model/model_header.hpp>", 
              code, fixed = TRUE)
  code <- sub("// [[Rcpp::depends(rstan)]]", 
              "// [[Rcpp::depends(rstan)]]\n#include <exporter.h>\n#include <RcppEigen.h>", code, fixed = TRUE)
  return(code)
}

expose_stan_functions <- function(file, ...) {
  model_code <- get_model_strcode(file, NULL)
  model_cppname <- legitimate_model_name(basename(file), obfuscate_name = TRUE)
  r <- .Call("stanfuncs", model_code, model_cppname, allow_undefined = FALSE)
  code <- expose_stan_functions_hacks(r$cppcode)
  compiled <- Rcpp::sourceCpp(code = paste(code, collapse = "\n"), ...)
  return(invisible(compiled$functions))
}
