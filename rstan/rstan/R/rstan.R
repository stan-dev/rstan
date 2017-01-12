# This file is part of RStan
# Copyright (C) 2012, 2013, 2014, 2015, 2016 Trustees of Columbia University
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
                       auto_write = rstan_options("auto_write"), 
                       obfuscate_model_name = TRUE,
                       allow_undefined = FALSE, 
                       includes = NULL) {

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
      else file.remove(tf)
    }
    else file <- normalizePath(file)
    
    stanc_ret <- stanc(file = file, model_code = model_code, 
                       model_name = model_name, verbose = verbose,
                       obfuscate_model_name = obfuscate_model_name, 
                       allow_undefined = allow_undefined)
    
    # find possibly identical stanmodels
    model_re <- "(^[[:alnum:]]{2,}.*$)|(^[A-E,G-S,U-Z,a-z].*$)|(^[F,T].+)"
    if(!is.null(model_name))
      if(!grepl(model_re, model_name))
        stop("model name must match ", model_re)
    S4_objects <- apropos(model_re, mode="S4", ignore.case=FALSE)
    if (length(S4_objects) > 0) {
      e <- environment()
      stanfits <- sapply(mget(S4_objects, envir = e, inherits = TRUE), 
                         FUN = is, class2 = "stanfit")
      stanmodels <- sapply(mget(S4_objects, envir = e, inherits = TRUE), 
                           FUN = is, class2 = "stanmodel")
      if (any(stanfits)) for (i in names(which(stanfits))) {
        obj <- get_stanmodel(get(i, envir = e, inherits = TRUE))
        if (identical(obj@model_code[1], stanc_ret$model_code[1])) return(obj)
      }
      if (any(stanmodels)) for (i in names(which(stanmodels))) {
        obj <- get(i, envir = e, inherits = TRUE)
        if (identical(obj@model_code[1], stanc_ret$model_code[1])) return(obj)
      }
    }
    
    mtime <- file.info(file)$mtime
    file.rds <- gsub("stan$", "rds", file)
    md5 <- tools::md5sum(file)
    if (!file.exists(file.rds)) {
      file.rds <- file.path(tempdir(), paste0(md5, ".rds"))
    }
    if(!file.exists(file.rds) ||
       (mtime.rds <- file.info(file.rds)$mtime) < 
       as.POSIXct(packageDescription("rstan")$Date) ||
       !is(obj <- readRDS(file.rds), "stanmodel") ||
       !is_sm_valid(obj) ||
       (!identical(stanc_ret$model_code, obj@model_code) && is.null(
        message("hash mismatch so recompiling; make sure Stan code ends with a blank line"))) ||
       dirname(file.rds) == gsub("\\", "/", tempdir(), fixed = TRUE) &&
       avoid_crash(obj@dso@.CXXDSOMISC$module) && is.null(
        message("recompiling to avoid crashing R session"))) {

        # do nothing
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
  if (.Platform$OS.type == "windows") find_rtools()
  else {
    CXX <- get_CXX()
    if (nchar(CXX) == 0) {
      WIKI <- "https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started"
      warning(paste("C++ compiler not found on system. If absent, see\n", WIKI))
    }
    else if (grepl("69", CXX, fixed = TRUE))
      warning("You may need to launch Xcode once to accept its license")
  }
  
  model_cppname <- stanc_ret$model_cppname 
  model_name <- stanc_ret$model_name 
  model_code <- stanc_ret$model_code 
  inc <- paste("#define STAN__SERVICES__COMMAND_HPP",
               # include, stanc_ret$cppcode,
               if(is.null(includes)) stanc_ret$cppcode else
                 sub("(class.*: public prob_grad \\{)", 
                     paste(includes, "\\1"), stanc_ret$cppcode),
               "#include <rstan/rstaninc.hpp>\n", 
               get_Rcpp_module_def_code(model_cppname), 
               sep = '')  

  if (verbose && interactive())
    cat("COMPILING THE C++ CODE FOR MODEL '", model_name, "' NOW.\n", sep = '')
  if (verbose) cat(system_info(), "\n")
  if (!is.null(boost_lib)) { 
    old.boost_lib <- rstan_options(boost_lib = boost_lib) 
    on.exit(rstan_options(boost_lib = old.boost_lib)) 
  } 
  if (!file.exists(rstan_options("boost_lib")))
    stop("Boost not found; call install.packages('BH')")
  
  if (!is.null(eigen_lib)) {
    old.eigen_lib <- rstan_options(eigen_lib = eigen_lib) 
    on.exit(rstan_options(eigen_lib = old.eigen_lib), add = TRUE) 
  }
  if (!file.exists(rstan_options("eigen_lib")))
    stop("Eigen not found; call install.packages('RcppEigen')")
  
  if (packageVersion("StanHeaders") > packageVersion("rstan"))
    stop("StanHeaders version is ahead of rstan version; ",
         "see https://github.com/stan-dev/rstan/wiki/RStan-Transition-Periods")
    

  
  dso <- cxxfunctionplus(signature(), body = paste(" return Rcpp::wrap(\"", model_name, "\");", sep = ''), 
                         includes = inc, plugin = "rstan", save_dso = save_dso | auto_write,
                         module_name = paste('stan_fit4', model_cppname, '_mod', sep = ''), 
                         verbose = verbose) 
               
  obj <- new("stanmodel", model_name = model_name, 
             model_code = model_code, 
             dso = dso, # keep a reference to dso
             mk_cppmodule = mk_cppmodule,  # mk_cppmodule function is defined in file stanmodel-class.R
             model_cpp = list(model_cppname = model_cppname, 
                              model_cppcode = stanc_ret$cppcode))
  
  if(missing(file) || (file.access(dirname(file), mode = 2) != 0) || !isTRUE(auto_write)) {
    tf <- tempfile()
    writeLines(model_code, con = tf)
    file <- file.path(tempdir(), paste0(tools::md5sum(tf), ".stan"))
    if(!file.exists(file)) file.rename(from = tf, to = file)
    else file.remove(tf)
    saveRDS(obj, file = gsub("stan$", "rds", file))
  }
  else if(isTRUE(auto_write)) {
    file <- gsub("stan$", "rds", file)
    if (file.exists(file)) {
      rds <- readRDS(file)
      if (!is(rds, "stanmodel"))
        warning(rds, " exists but is not a 'stanmodel' so not overwriting")
      saveRDS(obj, file = file)
    }
    else saveRDS(obj, file = file)
  }
  
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
  # This function works only for a stanmodel that is
  # created from function stan_model in package rstan.
  # 
  # Args:
  #   sm: the stanmodel object 
  # 
  if (isS4(sm@dso) && is_dso_loaded(sm@dso)) return(TRUE)
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
                 verbose = FALSE, include = TRUE,
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
                        "enable_random_init", "save_warmup",
                        "obfuscate_model_name"),
                      pre_msg = "passing unknown arguments: ", 
                      call. = FALSE)

  if (!is.list(data) && !is.environment(data) && !is.character(data))
    stop("'data' must be a list, environment, or character vector")
  
  if (is(fit, "stanfit")) sm <- get_stanmodel(fit)
  else { 
    attr(model_code, "model_name2") <- deparse(substitute(model_code))  
    if (missing(model_name)) model_name <- NULL
    sm <- stan_model(file, model_name = model_name, 
                     model_code = model_code, stanc_ret = NULL,
                     boost_lib = boost_lib, eigen_lib = eigen_lib, 
                     save_dso = save_dso, verbose = verbose)
  }

  if (is.null(sample_file))  sample_file <- NA 
  if (is.null(diagnostic_file))  diagnostic_file <- NA 

  sampling(sm, data, pars, chains, iter, warmup, thin, seed, init, 
           check_data = TRUE, sample_file = sample_file, 
           diagnostic_file = diagnostic_file,
           verbose = verbose, algorithm = match.arg(algorithm), 
           control = control, check_unknown_args = FALSE, 
           cores = cores, open_progress = open_progress, include = include, ...)
} 
