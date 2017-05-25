# This file is part of RStan
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
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

help_from_instance <- function(topic) {
  "Returns documentation string for the $method given by 'topic'"
  CLASS <- .self$getRefClass()
  if (missing(topic)) print(CLASS$help())
  else {
    if (is.name(substitute(topic))) topic <- substitute(topic)
    print(CLASS$help(parse(text = topic)))
  }
  return(invisible(NULL))
}

#' A ReferenceClass representing a Stan module
#' @importFrom methods as callGeneric new setRefClass
#' @importFrom Rcpp .__C__C++Object
#' @importFrom stats fivenum
#' @importFrom utils head str tail
StanModule <- setRefClass("StanModule", fields = list(module = "C++Object"), 
                          methods = list(
help = help_from_instance, 
hmc_nuts_diag_e_adapt = function(chains = 4, iter = 2000, warmup = 1000,
                                 thin = 1, init = "random", 
                                 sample_file = NULL, diagnostic_file = NULL, 
                                 adapt_gamma = 0.05, adapt_delta = 0.8,
                                 adapt_kappa = 0.75, adapt_t0 = 10, 
                                 adapt_init_buffer = 75,
                                 adapt_term_buffer = 50,
                                 adapt_window = 25, stepsize = 1,
                                 stepsize_jitter = 0, max_treedepth = 10, 
                                 ...) {
  "No U-Turn Sampling (NUTS) with adaptive integration time and a Euclidean metric
  'chains' = 4: number of Markov chains to execute
  'iter' = 2000: total number of iterations to execute per chain, including 'warmup'
  'warmup' = 1000: number of warmup iterations
  'thin' = 1: thinning interval
  'sample_file' = NULL: optional path specifying where to store draws on disk
  'diagnostic_file' = NULL: optional path specifying where to store diagnostics on disk
  'adapt_gamma' = 0.05: gamma tuning parameter
  'adapt_delta' = 0.8: delta tuning parameter
  'adapt_kappa' = 0.75: kappa tuning parameter
  'adapt_t0' = 10: t0 tuning parameter
  'adapt_init_buffer' = 75: pre-adaptation buffer
  'adapt_term_buffer' = 50: final adaptation buffer
  'adapt_window' = 25: interior phase buffer
  'stepsize' = 1: initial stepsize, final stepsize is adapted
  'stepsize_jitter' = 0: bad idea, do not change this
  'max_treedepth' = 10: maximum number of leapfrog steps per iteration'
  '...' = various other arguments such as:
    chain_id (integer >= 0)
    init_r (numeric > 0)
    refresh (integer >= 0)
    save_warmup (logical)"
  
  control <- list(adapt_engaged = TRUE, adapt_gamma = adapt_gamma,
                  adapt_kappa = adapt_kappa, adapt_t0 = adapt_t0,
                  adapt_init_buffer = as.integer(adapt_init_buffer),
                  adapt_term_buffer = as.integer(adapt_term_buffer),
                  adapt_window = as.integer(adapt_window), 
                  stepsize = stepsize, stepsize_jitter = stepsize_jitter,
                  max_treedepth = as.integer(max_treedepth))
  args_list <- rstan:::config_argss(chains = as.integer(chains), 
                                    iter = as.integer(iter),
                                    warmup = as.integer(warmup), 
                                    thin = as.integer(thin), init,
                                    seed = 12345L, # FIXME
                                    sample_file, diagnostic_file, 
                                    algorithm = "NUTS", 
                                    control = control, ...)
  samples <- lapply(args_list, FUN = function(x) .self$module$call_sampler(x))
  return(StanFitMCMC(sample_draws = samples))
},
sample = function() {
  "calls the 'hmc_nuts_diag_e_adapt' method with its default arguments"
  return(.self$hmc_nuts_diag_e_adapt())
} 
)
)

# A ReferenceClass representing the fixed aspects of a Stan program
# Question: Is seed a field?
setRefClass("StanConfig",
            contains = "VIRTUAL",
            fields = list(program.name = "character", 
                          shared.object = "character"),
            # in child classes there is one field for the C++ class and 
            # fields for each element of the data block in a Stan program
            methods = list(
finalize = function() {
  try(dyn.unload(.self$shared.object), silent = TRUE)
  return(invisible(NULL))
},
data = function(...) {
  if (nargs() == 0) {
    dotlist <- sapply(tail(names(.self$getRefClass()$fields()), -3L),
                      FUN = rstan:::dynGet, ifnotfound = NULL,
                      inherits = TRUE, simplify = FALSE)
  }
  else dotlist <- list(...)
  nms <- names(dotlist)
  if (is.null(nms)) nms <- as.character(match.call())[-1]
  for (i in seq_along(dotlist)) if (!is.null(dotlist[[i]])) {
    .self$field(nms[i], as(dotlist[[i]], class(.self$field(nms[i]))))
  }
},
# Question: Should build be a method or a standalone function()
build = function(seed = sample.int(.Machine$integer.max)) { 
  nms <- names(.self$getRefClass()$fields())[-c(1:2)]
  dummy <- function(...) stop("This function should not get called")
  datlist <- sapply(nms[-1], simplify = FALSE, 
                    FUN = function(nm) .self$field(nm)@values)
  mod <- new(.self$field(nms[1]), datlist, dummy)
  # FIXME: why the hell does this require rstan3::: ?
  return(rstan3:::StanModule(module = mod))
},
show = function() {
  nms <- tail(names(.self$getRefClass()$fields()), -3L)
  cat("Stan program", sQuote(.self$program.name), "with data", "\n", sep = " ")
  for (nm in nms) {
    cat(nm, ": ", sep = "")
    print(.self$field(nm))
    cat("\n")
  }
},
help = help_from_instance
))

#' Constructor
#' 
#' @param file A character vector of length one specifying a path to a file 
#'   written in the Stan language, which is ignored if \code{model_code} is
#'   specified
#' @param model_name A character vector of length one giving a name to the
#'   Stan program, which is autogenerated by default
#' @param model_code An optional character vector written in the Stan language
#' @param includes An optional character vector specifying includes and stuff
#' @export
stan_config <- function(file = file.choose(), 
                        model_name = NA_character_,
                        model_code, 
                        includes = character()) {
  allow_undefined <- !is.null(includes)
  if (missing(model_code)) {
    md5 <- tools::md5sum(file)
    tf <- file.path(tempdir(), paste0(md5, ".stan"))
    file.copy(from = file, to = tf)
    if (is.na(model_name))
      model_name <- sub("\\..*$", "", basename(file))
  }
  else {
    tf <- tempfile(fileext = ".stan")
    writeLines(model_code, tf)
    md5 <- tools::md5sum(tf)
    new_name <- file.path(tempdir(), paste0(md5, ".stan"))
    file.rename(from = tf, to = new_name)
    tf <- new_name
    if (is.na(model_name))
      model_name <- "anon_model"
  }
  on.exit(unlink(tf))
  ret <- rstan::stanc_builder(tf, verbose = FALSE,
                              allow_undefined = allow_undefined)
  class_name <- paste("model", ret$model_name, sep = "_")
  inc <- paste("#define STAN__SERVICES__COMMAND_HPP\n",
               "// [[Rcpp::depends(rstan)]]\n",
               if(length(includes)) ret$cppcode else
                 sub("(class[[:space:]]+[A-Za-z_][A-Za-z0-9_]*[[:space:]]*: public prob_grad \\{)",
                     paste(includes, "\\1"), ret$cppcode),
               "#include <rstan/rstaninc.hpp>\n",
               rstan:::get_Rcpp_module_def_code(class_name),
               sep = '')
  stuff <- Rcpp::sourceCpp(code = inc, env = environment())

  cppcode <- scan(what = character(), sep = "\n", quiet = TRUE, 
                  text = ret$cppcode)
  private <- grep("^private:$", cppcode) + 1L
  public <- grep("^public:$", cppcode) - 1L
  objects <- gsub("^.* ([0-9A-Za-z_]+).*;.*$", "\\1", cppcode[private:public])
  in_data <- grep("context__.vals_", cppcode, fixed = TRUE, 
                  value = TRUE)
  in_data <- gsub("^.*\\(\"(.*)\"\\).*;$", "\\1", in_data)
  objects <- intersect(objects, in_data)
  
  fields <- make_fields(ret$model_code, objects)
  fields <- c("C++Class", fields)
  obj_name <- paste0("stan_fit4", class_name)
  names(fields)[1] <- obj_name
  SPWD <- setRefClass(md5, contains = "StanConfig", fields = fields, 
                      where = .GlobalEnv)
  SPWD$lock(c("shared.object", obj_name))
  out <- SPWD(program.name = model_name,
              shared.object = dir(stuff$buildDirectory, pattern = "\\.so$", 
                                  full.names = TRUE))
  out$field(obj_name, get(obj_name, inherits = FALSE))
  return(out)
}

# Dynamically creates field names corresponding to objects in the data block
make_fields <- function(stan_code, objects) {
  if (length(objects) == 0) return(list())
  tokens <- sourcetools::tokenize_string(stan_code)
  matches <- sapply(objects, FUN = function(x) which(tokens$value == x)[1])
  rows <- tokens$row[matches]
  tokens <- tokens[tokens$row %in% rows,]
  semis <- which(tokens$type == "semi")
  arrays <- sapply(semis, FUN = function(s) {
    while(TRUE) {
      s <- s - 1L
      if (tokens[s, "value"] == "]") return(TRUE)
      if (tokens[s, "type"] == "symbol") return(FALSE)
    }
  })
  if (length(arrays) == 0) arrays <- integer()
  StanTypes <- names(getClass("StanType")@subclasses)
  types <- sapply(tokens$value, FUN = match, table = StanTypes)
  types <- types[!is.na(types)]
  nms <- sapply(1:length(semis), FUN = function(i) {
    s <- semis[i]
    if (!is.na(arrays[i]) && arrays[i]) while(TRUE) {
      s <- s - 1L
      if (tokens[s, "value"] == "[") return(tokens[s - 1L, "value"])
    }
    else while(TRUE) {
      s <- s - 1L
      if (tokens[s, "type"] == "symbol") return(tokens[s, "value"])
    }
  })
  names(types) <- nms
  return(sapply(types, simplify = FALSE, FUN = function(j) StanTypes[j]))
}

setClassUnion("scalarORarray", c("numeric", "array"))

#' A S4 class tree mirroring the types in the Stan language
#' 
#' These classes are not intended to be generated by users, but there are
#' associated S4 methods that can be called by users.
#' 
#' The \code{StanType-class} is virtual and the inheritance logic is as 
#' follows:
#' \itemize{
#'   \item real: inherits from \code{StanType-class} and represents 
#'     real scalars
#'   \itemize{
#'     \item int: inherits from \code{real-class} and represents
#'       integer scalars with optional slots for 'levels' and 'ordered'
#'       that can be used when coercing a \R \code{\link{factor}}
#'   }
#'   \item mat: inherits from  \code{StanType-class} and represents
#'     real matrices
#'   \itemize{
#'     \item cov_matrix: inherits from \code{mat-class} and
#'       represents covariance matrices
#'      \itemize{
#'        \item corr_matrix: inherits from \code{mat-class} and
#'          represents correlation matrices
#'      }
#'      \item cholesky_factor_cov: inherits from \code{mat-class} and
#'         represents a Cholesky factor of a covariance matrix
#'      \itemize{
#'         \item cholesky_factor_corr: inherits from 
#'           \code{cholesky_factor_cov-class} and represents a Cholesky
#'           factor of a correlation matrix
#'      }
#'      \item vec: inherits from \code{mat-class} and represents
#'         a vector, even though in the Stan language a vector and a
#'         matrix are both primitive types
#'      \itemize{
#'        \item row_vector: inherits from \code{vec-class} and
#'          represents a row vector, even though in the Stan language a
#'          vector and a row_vector are both primitive types
#'        \item simplex: inherits from \code{vec-class} and
#'          represents a simplex
#'        \item unit_vector: inherits from \code{vec-class} and represents
#'          a unit vector on a hypersphere
#'        \item ordered: inherits from \code{vec-class} and represents
#'          an ordered vector
#'        \itemize{
#'          \item positive_ordered: inherits from \code{ordered-class} and
#'            represents an ordered vector whose elements are positive
#'        }
#'      }
#'    }
#' }
#' 
#' @slot values A numeric (possibly integer) scalar or \code{link[base]{array}}
#' @slot lower A numeric scalar lower bound
#' @slot upper A numeric scalar upper bound
#' @slot type A length-one character vector among
#'   \itemize{
#'     \item data
#'     \item parameter
#'     \item transformed parameter
#'     \item generated quantity
#'     \item diagnostic
#'   }
#' @slot constrained A logical scalar indicating whether an unknown is in
#'   its constrained state   
#'
setClass("StanType", contains = "VIRTUAL",
         slots = list(values = "scalarORarray",
                      lower = "numeric", upper = "numeric",
                      type = "character", constrained = "logical"),
         prototype = list(values = array(NA_integer_, dim = 0L),
                          lower = -Inf, upper = Inf,
                          type = "data",
                          constrained = TRUE),
         validity = function(object) {
           if (!is.numeric(object@values))
             return("'values' must be numeric")
           if (length(object@type) != 1)
             return("'type' must be of length one")
           if (length(object@lower) != 1)
             return("'lower' must be of length one")
           if (length(object@upper) != 1)
             return("'upper' must be of length one")
           if (object@lower > object@upper)
             return("'lower' must be less than or equal to upper")
           if (any(object@values < object@lower))
             return("'lower' bound violated")
           if (any(object@values > object@upper))
             return("'upper' bound violated")
           types <- c("data", "parameter", "transformed parameter",
                      "generated quantity", "diagnostic")
           if (!(object@type %in% types))
             return(paste("'type' must be among:", 
                          paste(types, collapse = "\n"), sep = "\n"))
           return(TRUE)
         }
)

#' @rdname StanType-class
setClass("real", contains = "StanType")

#' @rdname StanType-class
setClass("int", contains = "real",
         slots = list(labels = "character", ordered = "logical"),
         prototype = list(lower = -.Machine$integer.max,
                          upper =  .Machine$integer.max),
         validity = function(object) {
           if (!is.integer(object@values)) return("'values' must contain integers")
           if (identical(object@type, "paramter") ||
               identical(object@type, "transformed parameters"))
                 return("'int' cannot be used in parameters or transformed parameters")
           return(TRUE)
         }
)

#' @rdname StanType-class
setClass("mat", contains = "StanType")

#' @rdname StanType-class
setClass("cov_matrix", contains = "mat")
#' @rdname StanType-class
setClass("corr_matrix", contains = "cov_matrix")
#' @rdname StanType-class
setClass("cholesky_factor_cov", contains = "matrix")
#' @rdname StanType-class
setClass("cholesky_factor_corr", contains = "cholesky_factor_cov")
#' @rdname StanType-class
setClass("vec", contains = "mat")

#' @rdname StanType-class
setClass("row_vector", contains = "vec")
#' @rdname StanType-class
setClass("simplex", contains = "vec", validity = function(object) {
          if (any(object@values < 0))
            return("'values' must be non-negative")
          if (!(sum(object@values) %in% 0:1))
            return("'values' must sum to 1")
          return(TRUE)
        }
)
#' @rdname StanType-class
setClass("unit_vector", contains = "vec", validity = function(object) {
          if (any(object@values < -1))
            return("'values' must >= -1")
          if (any(object@values >  1))
            return("'values' must <= 1")
          if (!(sum(object@values ^ 2) %in% 0:1))
            return("sum of squared 'values' must be 1")
          return(TRUE)
        }
)
#' @rdname StanType-class
setClass("ordered", contains = "vec", validity = function(object) {
          if (any(diff(object@values) < 0))
            return("'values' must be in increasing order")
          return(TRUE)
        }
)
#' @rdname StanType-class
setClass("positive_ordered", contains = "ordered", validity = function(object) {
          if (any(object@values < 0)) return("'values' must be non-negative")
          return(TRUE)
        }
)

not_implemented <- function() stop("this is not implemented yet")

#' A ReferenceClass representing MCMC output
#' 
#' @field warmup_draws A list of \code{\link{StanType}}s during warmup
#' @field sample_draws A list of \code{\link{StanType}}s after warmup
#' 
#' @seealso \code{\link{StanFitOptimization-methods}}
StanFitMCMC <- 
  setRefClass("StanFitMCMC", fields = list(warmup_draws = "list",
                                           sample_draws = "list"),
            # FIXME: add other extraction methods but without the get_ prefixes
            methods = list(show = not_implemented,
                           summary = not_implemented,
                           as.mcmc.list = not_implemented,
                           extract = not_implemented,
                           add_params = not_implemented,
                           help = help_from_instance)
)
StanFitMCMC$lock(setdiff(names(StanFitMCMC$fields()), "added_draws"))

