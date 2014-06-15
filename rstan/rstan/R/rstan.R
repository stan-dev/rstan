
.stanmodel_list <- create_sm_list()

rstan_options(cache_stanmodel = TRUE)

rstan_options(cache_path = file.path(path.expand("~"), ".stan"),
              cache_folder = "stanmodel") 

cache_full_path <- function() {
  file.path(rstan_options("cache_path"), rstan_options("cache_folder"))
}

load_stanmodel <- function(md5hash) {
  sm_saved <- file.path(cache_full_path(), paste0(md5hash, '.RData'))
  if (file.exists(sm_saved)) {
    obj <- local({load(sm_saved); obj}) 
    .stanmodel_list$set(md5hash, obj)
    return(obj)
  }
  NULL
}

get_stanmodel_by_md5 <- function(md5hash) {
  sm <- .stanmodel_list$get(md5hash) 
  if (is.null(sm)) {
    sm <- load_stanmodel(md5hash)
  }
  sm 
}

rm_stanmodel_cache <- function() {
  x <- cache_full_path()
  if (0 == unlink(x, recursive = TRUE)) { 
    .stanmodel_list <- create_sm_list()
    message("Succeeded in removing folder ", x, " that caches stanmodel.\n")
  } else {
    message("Failed to remove folder ", x, " that caches stanmodel files.\n")
  }
  invisible(NULL)
}

## 
stan_model <- function(file, 
                       model_name = "anon_model", 
                       model_code = '', 
                       stanc_ret = NULL, 
                       boost_lib = NULL, 
                       eigen_lib = NULL, 
                       save_dso = TRUE,
                       cache_stanmodel = rstan_options('cache_stanmodel'),
                       verbose = FALSE, ...) { 

  # Construct a stan model from stan code 
  # 
  # Args: 
  #   file: the file that has the model in Stan model language.
  #   model_name: a character for naming the model. 
  #   stanc_ret: An alternative way to specify the model
  #     by using returned results from stanc. 
  #   model_code: if file is not specified, we can use
  #     a character to specify the model.   

  if (is.null(stanc_ret)) {
    model_name2 <- deparse(substitute(model_code))
    if (is.null(attr(model_code, "model_name2")))
      attr(model_code, "model_name2") <- model_name2
    if (missing(model_name)) model_name <- NULL 
    stanc_ret <- stanc(file = file, model_code = model_code, 
                       model_name = model_name, verbose, ...)
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

  model_cppname <- stanc_ret$model_cppname 
  model_name <- stanc_ret$model_name 
  model_code <- stanc_ret$model_code 
  if (cache_stanmodel) {
    md5hash <- md5sum_stan_code(model_name, model_code) 
    sm <- get_stanmodel_by_md5(md5hash)
    if (!is.null(sm)) return(sm)
  }
  
  inc <- paste(stanc_ret$cppcode, 
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
                         includes = inc, plugin = "rstan", save_dso = save_dso,
                         module_name = paste('stan_fit4', model_cppname, '_mod', sep = ''), 
                         verbose = verbose) 
               
  obj <- new("stanmodel", model_name = model_name, 
             model_code = model_code, 
             dso = dso, # keep a reference to dso
             model_cpp = list(model_cppname = model_cppname, 
                              model_cppcode = stanc_ret$cppcode),
             md5hash = md5hash) 
  if (cache_stanmodel) {
    sm_cache_path <- cache_full_path()
    sm_cache_file <- file.path(sm_cache_path, paste0(md5hash, '.RData'))
    .stanmodel_list$set(md5hash, obj)
    if (!file.exists(sm_cache_path)) 
      dir.create(sm_cache_path, recursive = TRUE, mode = '0777')
    save(obj, file = sm_cache_file)
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
                 algorithm = c("NUTS", "HMC"), #, "Metropolis"),
                 control = NULL,
                 sample_file, # the file to which the samples are written
                 diagnostic_file, # the file to which diagnostics are written 
                 save_dso = TRUE,
                 cache_stanmodel = rstan_options('cache_stanmodel'),
                 verbose = FALSE, ...,
                 boost_lib = NULL, 
                 eigen_lib = NULL) {
  # Return a fitted model (stanfit object)  from a stan model, data, etc.  
  # A wrap of method stan_model and sampling of class stanmodel. 
  # 
  # Args:
  # 
  # Returns: 
  #   A S4 class stanfit object  

  if (is(fit, "stanfit")) sm <- get_stanmodel(fit)
  else { 
    attr(model_code, "model_name2") <- deparse(substitute(model_code))  
    if (missing(model_name)) model_name <- NULL 
    sm <- stan_model(file, model_name = model_name, model_code = model_code,
                     boost_lib = boost_lib, eigen_lib = eigen_lib, 
                     save_dso = save_dso, cache_stanmodel = cache_stanmodel, 
                     verbose = verbose, 
                     ...)
  }

  if (missing(sample_file))  sample_file <- NA 
  if (missing(diagnostic_file))  diagnostic_file <- NA 

  sampling(sm, data, pars, chains, iter, warmup, thin, seed, init, 
           check_data = TRUE, sample_file = sample_file, 
           diagnostic_file = diagnostic_file,
           verbose = verbose, algorithm = match.arg(algorithm), 
           control = control, ...) 
} 
