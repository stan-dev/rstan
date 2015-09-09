# This file is part of RStan
# Copyright (C) 2012, 2013, 2014, 2015 Jiqiang Guo and Benjamin Goodrich
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

## Use an environment to keep some options, especially, 
## for plotting. 

.rstan_opt_env <- new.env() 

init_rstan_opt_env <- function(e) {
  tmat <- matrix(c(254, 237, 222, 
                   253, 208, 162, 
                   253, 174, 107, 
                   253, 141, 60, 
                   230, 85, 13, 
                   166, 54, 3), 
                 byrow = TRUE, ncol = 3)

  rhat_cols <- rgb(tmat, alpha = 150, names = paste(1:nrow(tmat)),
                   maxColorValue = 255)

  assign("plot_rhat_breaks", c(1.1, 1.2, 1.5, 2), e)
  # in this default setting, 
  # if rhat < rhat.breaks[i], the color is rhat_cols[i]
  assign("plot_rhat_cols", rhat_cols[1:4], e)

  # when R-hat is NA, NaN, or Inf
  assign("plot_rhat_nan_col", rhat_cols[6] , e)
  # when R-hat is large than max(rhat.breaks) 
  assign("plot_rhat_large_col", rhat_cols[6], e)

  # color for indicating or important info. 
  # for example, the color of star and text saying
  # the variable is truncated in stan_plot_inferences
  assign("rstan_alert_col", rgb(230, 85, 13, maxColorValue = 255), e)

  # color for plot chains in trace plot and stan_plot_inferences 
  assign("rstan_chain_cols", rstancolc, e)
   
  # set the default number of parameters we are considered 
  # for plot of stanfit: when the number of parameters in a 
  # vector/array parameter is less than 
  # what is set here, we would have empty space. But when the 
  # number of parameters is larger than max, they are truncated. 
  assign('plot_standard_npar', 30, e)
  assign('plot_max_npar', 40, e)

  # color for shading the area of warmup trace plot
  assign("rstan_warmup_bg_col", rstancolgrey[3], e)

#   stan_lib_path  <- system.file('include', 'stanlib', package = 'rstan')
#   boost_dir <- dir(stan_lib_path, pattern = 'boost.*')
  # boost lib path 
#   boost_lib_path <- file.path(stan_lib_path, boost_dir)
  boost_lib_path <- system.file('include', package = 'BH')
  eigen_lib_path <- system.file('include', package = 'RcppEigen')
#   eigen_dir <- dir(stan_lib_path, pattern = 'eigen.*')
#   eigen_lib_path <- file.path(stan_lib_path, eigen_dir)
  assign("eigen_lib", eigen_lib_path, e) 
  assign("boost_lib", boost_lib_path, e) 

  ya_boost  <- system.file('include', 'boost_not_in_BH', package = 'rstan')
  assign('boost_lib2', ya_boost, e)

  assign('auto_write', FALSE, e)
  assign("iter", 2000, e)
  assign("chains", 4, e)
  
  # cat("init_rstan_opt_env called.\n")
  invisible(e)
}

# init_rstan_opt_env(.rstan_opt_env)

rstan_options <- function(...) { 
  # Set/get options in RStan
  # Args: any options can be defined, using 'name = value' 
  #
  # e <- rstan:::.rstan_opt_env 
  e <- .rstan_opt_env 
  if (length(as.list(e)) == 0) 
    init_rstan_opt_env(e) 

  a <-  list(...)
  len <- length(a) 
  if (len < 1) return(NULL) 
  a_names <- names(a) 
  if (is.null(a_names)) { # case like rstan_options("a", "b")
    empty <- rep(TRUE, len)
    empty_len <- len 
  } else { # case like rstan_options(a = 3, b = 4, "c")
    empty <- (a_names == '') 
    empty_len <- sum(empty) 
  } 
  for (i in which(empty)) {
    if (!is.character(a[[i]])) stop("rstan_options only accepts arguments as 'name=value' or 'name'") 
  } 
  
  r <- if (empty_len < len) mget(a_names[!empty], envir = e, ifnotfound = NA) 
  if (empty_len > 0) 
    r <- c(r, mget(unlist(a[empty]), envir = e, 
                   ifnotfound = list(function(x) { warning("rstan option '", x, "' not found"); NA })))

  # set options 
  for (n in a_names[!empty]) {
    if (n == 'plot_rhat_breaks') { assign(n, sort(a[[n]]), e); next }
    assign(n, a[[n]], e)
  } 

  if (len == 1) return(invisible(r[[1]])) 
  invisible(r)
} 
