# This file is part of RStan
# Copyright (C) 2012, 2013, 2014, 2015 Trustees of Columbia University
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
  def_Rcpp_module_hpp_file <- 
    system.file('include', '/rstan/rcpp_module_def_for_rstan.hpp', package = 'rstan') 
  if (def_Rcpp_module_hpp_file == '') 
    stop("Rcpp module definition file for rstan is not found.\n") 
  src <- paste(readLines(def_Rcpp_module_hpp_file), collapse = '\n')
  gsub("%model_name%", model_name, src)
} 

