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

## 
stan_model <- function(file, 
                       model_name = "anon_model", 
                       model_code = '', 
                       stanc_ret = NULL, 
                       boost_lib = NULL, 
                       eigen_lib = NULL, 
                       save_dso = TRUE,
                       verbose = FALSE, 
                       auto_write = rstan_options("auto_write"), ...) { 

  # Construct a stan model from stan code 
  # 
  # Args: 
  #   file: the file that has the model in Stan model language.
  #   model_name: a character for naming the model. 
  #   stanc_ret: An alternative way to specify the model
  #     by using returned results from stanc. 
  #   model_code: if file is not specified, we can used 
  #     a character to specify the model.   

  if (is.null(stanc_ret)) {
    model_name2 <- deparse(substitute(model_code))
    if (is.null(attr(model_code, "model_name2")))
      attr(model_code, "model_name2") <- model_name2
    if (missing(model_name)) model_name <- NULL 
    
    if(missing(file)) {
      tf <- tempfile()
      writeLines(model_code, con = tf)
      file <- file.path(dirname(tf), paste0(tools::md5sum(tf), ".stan"))
      if(!file.exists(file)) file.rename(from = tf, to = file)
    }
    mtime <- file.info(file)$mtime
    file.rda <- gsub("stan$", "rda", file)
    if (!file.exists(file.rda)) {
      file.rda <- file.path(tempdir(), paste0(tools::md5sum(file), ".rda"))
    }
    if( mtime < as.POSIXct(packageDescription("rstan")$Date) ||
       !file.exists(file.rda) ||
       file.info(file.rda)$mtime <  mtime ||
       !is(obj <- readRDS(file.rda), "stanmodel") ||
       !is_sm_valid(obj)) {
         
          stanc_ret <- stanc(file = file, model_code = model_code, 
                             model_name = model_name, verbose, ...)
    }
    else return(invisible(obj))
  }
  if (!is.list(stanc_ret)) {
    stop("stanc_ret needs to be the returned object from stanc.")
  } 
  m <- match(c("cppcode", "model_name", "status"), names(stanc_ret)) 
  if (any(is.na(m))) {
    stop("stanc_ret does not have element `cppcode', `model_name', and `status'") 
  } else {
    if (!stanc_ret$status)
      stop("stanc_ret is not a successfully returned list from stanc")
  }
  
  # check for compilers
  check <- system2(file.path(R.home(component = "bin"), "R"), 
                   args = "CMD config CXX", 
                   stdout = TRUE, stderr = FALSE)
  if(identical(check, "")) {
    WIKI <- "https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started"
    warning(paste("C++ compiler not found on system. If absent, see", WIKI))
  }
  
  model_cppname <- stanc_ret$model_cppname 
  model_name <- stanc_ret$model_name 
  model_code <- stanc_ret$model_code 
  inc <- paste("#define STAN__SERVICES__COMMAND_HPP",
               stanc_ret$cppcode,
               "#include <rstan/rstaninc.hpp>\n", 
               get_Rcpp_module_def_code(model_cppname), 
               sep = '')  

  cat("COMPILING THE C++ CODE FOR MODEL '", model_name, "' NOW.\n", sep = '') 
  if (verbose) cat(system_info(), "\n")
  if (!is.null(boost_lib)) { 
    old.boost_lib <- rstan_options(boost_lib = boost_lib) 
    on.exit(rstan_options(boost_lib = old.boost_lib)) 
  } 

  if (!is.null(eigen_lib)) { 
    old.eigen_lib <- rstan_options(eigen_lib = eigen_lib) 
    on.exit(rstan_options(eigen_lib = old.eigen_lib), add = TRUE) 
  }

  dso <- cxxfunctionplus(signature(), body = paste(" return Rcpp::wrap(\"", model_name, "\");", sep = ''), 
                         includes = inc, plugin = "rstan", save_dso = save_dso | auto_write,
                         module_name = paste('stan_fit4', model_cppname, '_mod', sep = ''), 
                         verbose = verbose) 
               
  obj <- new("stanmodel", model_name = model_name, 
             model_code = model_code, 
             dso = dso, # keep a reference to dso
             model_cpp = list(model_cppname = model_cppname, 
                              model_cppcode = stanc_ret$cppcode))
  
  if(missing(file) || (file.access(dirname(file), mode = 2) != 0) || !isTRUE(auto_write)) {
    tf <- tempfile()
    writeLines(model_code, con = tf)
    file <- file.path(tempdir(), paste0(tools::md5sum(tf), ".stan"))
    if(!file.exists(file)) file.rename(from = tf, to = file)
    saveRDS(obj, file = gsub("stan$", "rda", file))
  }
  else if(isTRUE(auto_write)) saveRDS(obj, file = gsub("stan$", "rda", file))
  
  invisible(obj) 
  ## We keep a reference to *dso* above to avoid dso to be 
  ## deleted by R's garbage collection. Note that if dso 
  ## is freed, we can lose the compiled shared object, which
  ## can cause segfault later. 
}

is_sm_valid <- function(sm) {
  # Test if a stan model (compiled object) is still valid. 
  # It could become invalid when the user do not specify
  # save_dso when calling stan_model. So when the user
  # use the model created in another R session, the dso
  # is lost. 
  # 
  # Args:
  #   sm: the stanmodel object 
  # 
  if (is_dso_loaded(sm@dso)) return(TRUE)
  sm@dso@dso_saved && identical(sm@dso@system, R.version$system)
} 

##
##
## 

stan <- function(file, model_name = "anon_model", 
                 model_code = '', 
                 fit = NA, 
                 data = list(), 
                 pars = NA, 
                 chains = 4, iter = 2000, 
                 warmup = floor(iter / 2), 
                 thin = 1, 
                 init = "random", 
                 seed = sample.int(.Machine$integer.max, 1), 
                 algorithm = c("NUTS", "HMC", "Fixed_param"), #, "Metropolis"),
                 control = NULL,
                 sample_file = NULL, # the file to which the samples are written
                 diagnostic_file = NULL, # the file to which diagnostics are written 
                 save_dso = TRUE,
                 verbose = FALSE, 
                 cores = getOption("mc.cores", 1L),
                 open_progress = interactive() && !isatty(stdout()) &&
                                 !identical(Sys.getenv("RSTUDIO"), "1"), 
                 ...,
                 boost_lib = NULL, 
                 eigen_lib = NULL) {
  # Return a fitted model (stanfit object)  from a stan model, data, etc.  
  # A wrap of method stan_model and sampling of class stanmodel. 
  # 
  # Args:
  # 
  # Returns: 
  #   A S4 class stanfit object  

  dot_arg_names <- names(list(...))
  is_arg_recognizable(dot_arg_names, 
                      c("chain_id", "init_r", "test_grad", 
                        "append_samples", "refresh", "control",
                        "enable_random_init",
                        "obfuscate_model_name"),
                      pre_msg = "passing unknown arguments: ", 
                      call. = FALSE)

  if (is(fit, "stanfit")) sm <- get_stanmodel(fit)
  else { 
    attr(model_code, "model_name2") <- deparse(substitute(model_code))  
    if (missing(model_name)) model_name <- NULL
    if (cores == 1) {
      sr <- stanc(file, model_name = model_name, model_code = model_code,
                  verbose = verbose)
    }
    else sr <- NULL
    sm <- stan_model(file, model_name = model_name, 
                     model_code = model_code, stanc_ret = sr,
                     boost_lib = boost_lib, eigen_lib = eigen_lib, 
                     save_dso = save_dso, verbose = verbose, ...)
  }

  if (is.null(sample_file))  sample_file <- NA 
  if (is.null(diagnostic_file))  diagnostic_file <- NA 

  sampling(sm, data, pars, chains, iter, warmup, thin, seed, init, 
           check_data = TRUE, sample_file = sample_file, 
           diagnostic_file = diagnostic_file,
           verbose = verbose, algorithm = match.arg(algorithm), 
           control = control, check_unknown_args = FALSE, 
           cores = cores, open_progress = open_progress, ...)
} 
