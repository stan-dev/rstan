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

setMethod("show", "stanmodel",
          function(object) {
            cat("S4 class stanmodel '", object@model_name, "' coded as follows:\n" ,sep = '')
            cat(object@model_code, "\n")
          })

setGeneric(name = 'optimizing',
           def = function(object, ...) {

             if (is.sparc()) {
               msg <- "optimizing() will likely crash on SPARC."
               if (interactive()) stop(msg, " Run in batch mode to test.")
               else message(msg)
             }
             standardGeneric("optimizing")
})

setGeneric(name = 'vb',
           def = function(object, ...) {

             if (is.sparc()) {
               msg <- "vb() will likely crash on SPARC."
               if (interactive()) stop(msg, " Run in batch mode to test.")
               else message(msg)
             }
             standardGeneric("vb")
})

setGeneric(name = "sampling",
           def = function(object, ...) { standardGeneric("sampling")})

setMethod('get_stancode', signature = "stanmodel",
          function(object, print = FALSE) {
            code <- object@model_code
            if (print) cat(code, "\n")
            return(code)
          })

setGeneric(name = "get_cppcode",
           def = function(object, ...) { standardGeneric("get_cppcode") })

setMethod("get_cppcode", "stanmodel",
          function(object) {
            object@model_cpp$model_cppcode
          })

setGeneric(name = "get_cxxflags",
           def = function(object, ...) { standardGeneric("get_cxxflags") })
setMethod("get_cxxflags", "stanmodel", function(object) { object@dso@cxxflags })

new_empty_stanfit <- function(stanmodel, miscenv = new.env(parent = emptyenv()),
                              model_pars = character(0), par_dims = list(),
                              mode = 2L, sim = list(),
                              inits = list(), stan_args = list()) {
  new("stanfit",
      model_name = stanmodel@model_name,
      model_pars = model_pars,
      par_dims = par_dims,
      mode = mode,
      sim = sim,
      inits = inits,
      stan_args = stan_args,
      stanmodel = stanmodel,
      date = date(),
      .MISC = miscenv)
}


prep_call_sampler <- function(object) {
  if (!is_sm_valid(object))
    stop(paste("the compiled object from C++ code for this model is invalid, possible reasons:\n",
               "  - compiled with save_dso=FALSE;\n",
               "  - compiled on a different platform;\n",
               "  - does not exist (created from reading csv files).", sep = ''))
  if (!is_dso_loaded(object@dso)) {
    # load the dso if available
    grab_cxxfun(object@dso)
  }
}

# This function gets used when a stanmodel instance is created
# in function stan_model.
mk_cppmodule <- function(object) {
  prep_call_sampler(object)
  model_cppname <- object@model_cpp$model_cppname
  mod <- get("module", envir = object@dso@.CXXDSOMISC, inherits = FALSE)
  eval(call("$", mod, paste('stan_fit4', model_cppname, sep = '')))
}

setMethod("vb", "stanmodel",
          function(object, data = list(), pars = NA, include = TRUE,
                   seed = sample.int(.Machine$integer.max, 1),
                   init = 'random',
                   check_data = TRUE, sample_file = tempfile(fileext = '.csv'),
                   algorithm = c("meanfield", "fullrank"),
                   importance_resampling = FALSE,
                   keep_every = 1,
                   ...) {
            if (isTRUE(rstan_options("threads_per_chain") > 1L)) {
              Sys.setenv("STAN_NUM_THREADS" = rstan_options("threads_per_chain"))
            }

            if (is.list(data) & !is.data.frame(data)) {
              parsed_data <- with(data, parse_data(get_cppcode(object)))
              if (!is.list(parsed_data)) {
                message("failed to get names of data from the model; sampling not done")
                return(invisible(new_empty_stanfit(object)))
              }
              # for (nm in names(data)) parsed_data[[nm]] <- data[[nm]]
              # parsed_data <- parsed_data[!sapply(parsed_data, is.null)]
              data <- parsed_data
            } else if (is.character(data)) { # names of objects
              data <- try(mklist(data))
              if (is(data, "try-error")) {
                message("failed to create the data; sampling not done")
                return(invisible(new_empty_stanfit(object)))
              }
            }

            if (check_data) {
              data <- try(force(data))
              if (is(data, "try-error")) {
                message("failed to evaluate the data; sampling not done")
                return(invisible(new_empty_stanfit(object)))
              }

              if (!missing(data) && length(data) > 0) {
                data <- try(data_preprocess(data))
                if (is(data, "try-error")) {
                  message("failed to preprocess the data; variational Bayes not done")
                  return(invisible(new_empty_stanfit(object)))
                }
              } else data <- list()
            }
            cxxfun <- grab_cxxfun(object@dso)
            stan_fit_cpp_module <- object@mk_cppmodule(object)
            if (test_221(object@model_cpp$model_cppcode) &&
                stan_fit_cpp_module@constructors[[1]]$nargs == 2L) {
              mod <- try(new(stan_fit_cpp_module, data, as.integer(seed)))
              if (is(mod, "try-error")) {
                message('failed to create the sampler; sampling not done')
                return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv)))
              }
              sampler <- try(new(stan_fit, mod$fit_ptr()))
              if (is(sampler, "try-error")) {
                message('failed to create the sampler; sampling not done')
                return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv)))
              }
            } else {
              sampler <- try(new(stan_fit_cpp_module, data, as.integer(seed), cxxfun))
              if (is(sampler, "try-error")) {
                message('failed to create the optimizer; optimization not done')
                return(invisible(list(stanmodel = object)))
              }
            }

            if (is.numeric(init)) init <- as.character(init)
            if (is.function(init)) init <- init()
            if (!is.list(init) && !is.character(init)) {
              message("wrong specification of initial values")
              return(invisible(new_empty_stanfit(object)))
            }
            seed <- check_seed(seed, warn = 1)
            if (is.null(seed))
              return(invisible(list(stanmodel = object)))
            args <- list(init = init, seed = seed, chain_id = 1L,
                         method = "variational",
                         algorithm = match.arg(algorithm))

            if (!is.null(sample_file) && !is.na(sample_file))
              args$sample_file <- writable_sample_file(sample_file)
            dotlist <- list(...)
            is_arg_recognizable(names(dotlist),
                                c("iter", "init_r",
                                  "append_samples",
                                  "elbo_samples",
                                  "eta",
                                  "adapt_engaged",
                                  "eval_elbo",
                                  "grad_samples",
                                  "output_samples",
                                  "adapt_iter",
                                  "tol_rel_obj",
                                  "refresh"),
                                 pre_msg = "passing unknown arguments: ",
                                 call. = FALSE)
            if (!is.null(dotlist$method))  dotlist$method <- NULL
            if (is.null(dotlist$eta)) dotlist$eta <- 1.0
            else {
              if (dotlist$eta <= 0) stop("'eta' must be > 0")
            }

            sfmiscenv <- new.env(parent = emptyenv())
            assign("stan_fit_instance", sampler, envir = sfmiscenv)

            m_pars <- sampler$param_names()
            p_dims <- sampler$param_dims()
            if(!include) {
              if (length(pars) == 1 && is.na(pars)) pars <- "lp__"
              else pars <- setdiff(m_pars, pars)
              if (length(pars) == 0) pars <- "lp__"
            }

            if (!missing(pars) && !any(is.na(pars)) && length(pars) > 0) {
              sampler$update_param_oi(pars)
              m <- which(match(pars, m_pars, nomatch = 0) == 0)
              if (length(m) > 0) {
                message("no parameter ", paste(pars[m], collapse = ', '), "; sampling not done")
                return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv, m_pars, p_dims, 2L)))
              }
            }
            else pars <- m_pars

            skeleton <- create_skeleton(m_pars, p_dims)
            cC <- sapply(names(skeleton), simplify = FALSE, FUN = function(x) {
              param <- skeleton[[x]]
              if (x == "lp__") TRUE
              else if (x %in% pars) rep(TRUE, length(param))
              else rep(FALSE, length(param))
            })

            vbres <- sampler$call_sampler(c(args, dotlist))
            samples <- read_one_stan_csv(attr(vbres, "args")$sample_file)
            diagnostic_columns <- which(grepl('__$',colnames(samples)))[-1]
            if (length(diagnostic_columns)>0) {
              diagnostics <- samples[-1,diagnostic_columns]
              samples <- samples[,-diagnostic_columns]
            } else {
              diagnostics <- NULL
            }
            pest <- rstan_relist(as.numeric(samples[1,-1]), skeleton[-length(skeleton)])
            means <- sapply(samples[-1,], mean)
            means <- as.matrix(c(means[-1], means[1]))
            colnames(means) <- "chain:1"
            assign("posterior_mean_4all", means, envir = sfmiscenv)
            inits_used <- rstan_relist(attr(vbres, "inits"), skeleton)
            if ("lp__" %in% names(inits_used))  inits_used$lp__ <- NULL
            samples <- cbind(samples[-1,-1,drop=FALSE],
                             "lp__" = samples[-1,1])[,unlist(cC)]
            cC <- cC[sapply(cC, any)] # any(logical()) is FALSE
            count <- 1L
            for (i in seq_along(cC)) {
              len <- length(cC[[i]])
              cC[[i]] <- count:(count + len - 1L)
              count <- count + len
            }
            samples <- samples[,unlist(cC[c(pars, "lp__")]), drop = FALSE]
            fnames_oi <- sampler$param_fnames_oi()
            n_flatnames <- length(fnames_oi)
            iter <- nrow(samples)
            if ("log_g__" %in% colnames(diagnostics)) {
              if (length(extralp <- which(grepl('lp__.1', colnames(samples)))) > 0)
                  samples <- samples[,-extralp]
              lr <- diagnostics$log_p-diagnostics$log_g
              lr[lr == -Inf] <- -800
              p <- suppressWarnings(loo::psis(lr, r_eff = 1))
              p$log_weights <- p$log_weights - log_sum_exp(p$log_weights)
              theta_pareto_k <- suppressWarnings(apply(samples, 2L, function(col) {
                if (all(is.finite(col))) loo::psis(log1p(col^2) / 2 + lr, r_eff = 1)$diagnostics$pareto_k
                else NaN
              }))
              ## todo: change fixed threshold to an option
              if (p$diagnostics$pareto_k > 1) {
                warning("Pareto k diagnostic value is ",
                        round(p$diagnostics$pareto_k,2),
                        ". Resampling is disabled.",
                        " Decreasing tol_rel_obj may help if variational algorithm has terminated prematurely.",
                        " Otherwise consider using sampling instead.", call. = FALSE, immediate. = TRUE)
                #importance_resampling <- FALSE
              } else if (p$diagnostics$pareto_k > 0.7) {
                warning("Pareto k diagnostic value is ",
                        round(p$diagnostics$pareto_k,2),
                        ". Resampling is unreliable.",
                        " Increasing the number of draws or decreasing tol_rel_obj may help.",
                        call. = FALSE, immediate. = TRUE)
              }
              psis <- loo::nlist(pareto_k = p$diagnostics$pareto_k,
                                 n_eff = p$diagnostics$n_eff / keep_every)
              ## importance_resampling
              if (importance_resampling) {
                iter <- ceiling(dim(samples)[1] / keep_every)
                ir_idx <- sample_indices(exp(p$log_weights),
                                         n_draws = iter)
                samples <- samples[ir_idx,]
                ## SIR mcse and n_eff
                w_sir <- as.numeric(table(ir_idx)) / length(ir_idx)
                mcse <- apply(samples[!duplicated(ir_idx),], 2L, function(col) {
                  if (all(is.finite(col))) sqrt(sum(w_sir^2*(col-mean(col))^2))
                  else NaN
                })
                n_eff <- round(apply(samples[!duplicated(ir_idx),], 2L, var) / (mcse^2), digits = 0)
              } else {
                ir_idx <- NULL
                mcse <- rep(NaN, length(theta_pareto_k))
                n_eff <- rep(NaN, length(theta_pareto_k))
              }
              diagnostics <- list(as.list(diagnostics), psis = psis,
                                  theta_pareto_k = theta_pareto_k,
                                  ir_idx = ir_idx, mcse = mcse, n_eff = n_eff)
            } else {
              diagnostics <- list(as.list(diagnostics))
            }
            sim <- list(samples = list(as.list(samples)),
                        diagnostics = diagnostics,
                        iter = iter, thin = keep_every,
                        warmup = 0L,
                        chains = 1L,
                        n_save = iter,
                        est = pest, # point estimate
                        warmup2 = 0L, # number of warmup iters in n_save
                        permutation = list(1:iter),
                        pars_oi = sampler$param_names_oi(),
                        dims_oi = sampler$param_dims_oi(),
                        fnames_oi = fnames_oi,
                        n_flatnames = n_flatnames)
            nfit <- new("stanfit",
                        model_name = object@model_name,
                        model_pars = m_pars,
                        par_dims = p_dims,
                        mode = 0L,
                        sim = sim,
                        # keep a record of the initial values
                        inits = inits_used,
                        stan_args = list(args),
                        stanmodel = object,
                        # keep a ref to avoid garbage collection
                        # (see comments in fun stan_model)
                        date = date(),
                        .MISC = sfmiscenv)
            return(nfit)
          })

setMethod("optimizing", "stanmodel",
          function(object, data = list(),
                   seed = sample.int(.Machine$integer.max, 1),
                   init = 'random', check_data = TRUE, sample_file = NULL,
                   algorithm = c("LBFGS", "BFGS", "Newton"),
                   verbose = FALSE, hessian = FALSE, as_vector = TRUE,
                   draws = 0, constrained = TRUE,
                   importance_resampling = FALSE, ...) {
            if (isTRUE(rstan_options("threads_per_chain") > 1L)) {
              Sys.setenv("STAN_NUM_THREADS" = rstan_options("threads_per_chain"))
            }

            if (is.list(data) & !is.data.frame(data)) {
              parsed_data <- with(data, parse_data(get_cppcode(object)))
              if (!is.list(parsed_data)) {
                message("failed to get names of data from the model; sampling not done")
                return(invisible(new_empty_stanfit(object)))
              }
              # for (nm in names(data)) parsed_data[[nm]] <- data[[nm]]
              # parsed_data <- parsed_data[!sapply(parsed_data, is.null)]
              data <- parsed_data
            } else if (is.character(data)) { # names of objects
              data <- try(mklist(data))
              if (is(data, "try-error")) {
                message("failed to create the data; sampling not done")
                return(invisible(new_empty_stanfit(object)))
              }
            }

            if (check_data) {
              data <- try(force(data))
              if (is(data, "try-error")) {
                message("failed to evaluate the data; sampling not done")
                return(invisible(NULL))
              }

              if (!missing(data) && length(data) > 0) {
                data <- try(data_preprocess(data))
                if (is(data, "try-error")) {
                  message("failed to preprocess the data; optimization not done")
                  return(invisible(list(stanmodel = object)))
                }
              } else data <- list()
            }
            cxxfun <- grab_cxxfun(object@dso)
            sfmiscenv <- new.env(parent = emptyenv())
            stan_fit_cpp_module <- object@mk_cppmodule(object)
            if (test_221(object@model_cpp$model_cppcode) &&
                stan_fit_cpp_module@constructors[[1]]$nargs == 2L) {
              mod <- try(new(stan_fit_cpp_module, data, as.integer(seed)))
              if (is(mod, "try-error")) {
                message('failed to create the sampler; sampling not done')
                return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv)))
              }
              sampler <- try(new(stan_fit, mod$fit_ptr()))
              if (is(sampler, "try-error")) {
                message('failed to create the sampler; sampling not done')
                return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv)))
              }
            } else {
              sampler <- try(new(stan_fit_cpp_module, data, as.integer(seed), cxxfun))
              if (is(sampler, "try-error")) {
                message('failed to create the optimizer; optimization not done')
                return(invisible(list(stanmodel = object)))
              }
            }

            m_pars <- sampler$param_names()
            idx_wo_lp <- which(m_pars != "lp__")
            m_pars <- m_pars[idx_wo_lp]
            p_dims <- sampler$param_dims()[idx_wo_lp]
            if (is.numeric(init)) init <- as.character(init)
            if (is.function(init)) init <- init()
            if (!is.list(init) && !is.character(init)) {
              message("wrong specification of initial values")
              return(invisible(list(stanmodel = object)))
            }
            seed <- check_seed(seed, warn = 1)
            if (is.null(seed))
              return(invisible(list(stanmodel = object)))
            args <- list(init = init,
                         seed = seed,
                         method = "optim",
                         algorithm = match.arg(algorithm))

            if (!is.null(sample_file) && !is.na(sample_file))
              args$sample_file <- writable_sample_file(sample_file)
            dotlist <- list(...)
            is_arg_recognizable(names(dotlist),
                                c("iter",
                                  "save_iterations",
                                  "refresh",
                                  "init_alpha",
                                  "tol_obj",
                                  "tol_grad",
                                  "tol_param",
                                  "tol_rel_obj",
                                  "tol_rel_grad",
                                  "history_size"),
                                 pre_msg = "passing unknown arguments: ",
                                 call. = FALSE)
            if (!is.null(dotlist$method))  dotlist$method <- NULL
            if (!verbose && is.null(dotlist$refresh)) dotlist$refresh <- 0L
            optim <- sampler$call_sampler(c(args, dotlist))
            optim$return_code <- attr(optim, "return_code")
            if (optim$return_code != 0) warning("non-zero return code in optimizing")
            attr(optim, "return_code") <- NULL
            fnames <- sampler$param_fnames_oi()
            names(optim$par) <- fnames[-length(fnames)]
            skeleton <- create_skeleton(m_pars, p_dims)
            theta <- rstan_relist(optim$par, skeleton)
            theta <- sampler$unconstrain_pars(theta)
            if (hessian || draws) {
              fn <- function(theta) {
                sampler$log_prob(theta, FALSE, FALSE)
              }
              gr <- function(theta) {
                sampler$grad_log_prob(theta, FALSE)
              }
              H <- optimHess(theta, fn, gr, control = list(fnscale = -1))
              colnames(H) <- rownames(H) <- sampler$unconstrained_param_names(FALSE, FALSE)
              if (hessian) optim$hessian <- H
            }
            if (draws > 0 && is.matrix(R <- try(chol(-H)))) {
                K <- ncol(R)
                R_inv <- backsolve(R, diag(K))
                Z <- matrix(rnorm(K * draws), K, draws)
                theta_tilde <- t(theta + R_inv %*% Z)
                if (importance_resampling) {
                  log_p <- apply(theta_tilde, 1, FUN = function(theta) {
                    sampler$log_prob(theta, adjust_transform = TRUE, gradient = FALSE)
                  })
                  log_g <- colSums(dnorm(Z, log = TRUE)) - sum(log(diag(R_inv)))
                  optim$log_p <- log_p
                  optim$log_g <- log_g
                } else {
                  optim$log_p <- rep(NaN, length(theta))
                  optim$log_g <- rep(NaN, length(theta))
                }
                colnames(theta_tilde) <- colnames(H)
                optim$log_prob <- sampler$log_prob
                optim$grad_log_prob <- sampler$grad_log_prob
            } else {
              theta_tilde <- t(theta)
            }
            if (constrained) {
              theta_tilde <- t(apply(theta_tilde, 1, FUN = function(theta) {
                sampler$constrain_pars(theta)
              }))
              if (NCOL(theta_tilde) != length(optim$par)) theta_tilde <- t(theta_tilde)
            }
            colnames(theta_tilde) <- names(optim$par)
            optim$theta_tilde <- theta_tilde
            if (!as_vector) optim$par <- rstan_relist(optim$par, skeleton)
            return(optim)
          })

setMethod("sampling", "stanmodel",
          function(object, data = list(), pars = NA, chains = 4, iter = 2000,
                   warmup = floor(iter / 2),
                   thin = 1, seed = sample.int(.Machine$integer.max, 1),
                   init = "random", check_data = TRUE,
                   sample_file = NULL, diagnostic_file = NULL, verbose = FALSE,
                   algorithm = c("NUTS", "HMC", "Fixed_param"), #, "Metropolis"),
                   control = NULL, include = TRUE,
                   cores = getOption("mc.cores", 1L),
                   open_progress = interactive() && !isatty(stdout()) &&
                     !identical(Sys.getenv("RSTUDIO"), "1"),
                   show_messages = TRUE, ...) {
            is_arg_deprecated(names(list(...)),
                              c("enable_random_init"),
                              pre_msg = "passing deprecated arguments: ")

            if (isTRUE(rstan_options("threads_per_chain") > 1L)) {
              Sys.setenv("STAN_NUM_THREADS" = rstan_options("threads_per_chain"))
            }

            objects <- ls()
            if (is.list(data) & !is.data.frame(data)) {
              parsed_data <- with(data, parse_data(get_cppcode(object)))
              if (!is.list(parsed_data)) {
                message("failed to get names of data from the model; sampling not done")
                return(invisible(new_empty_stanfit(object)))
              }
              # for (nm in names(data)) parsed_data[[nm]] <- data[[nm]]
              # parsed_data <- parsed_data[!sapply(parsed_data, is.null)]
              data <- parsed_data
            } else if (is.character(data)) { # names of objects
              data <- try(mklist(data))
              if (is(data, "try-error")) {
                message("failed to create the data; sampling not done")
                return(invisible(new_empty_stanfit(object)))
              }
            }
            # check data and preprocess
            if (verbose)
              cat('\n', "CHECKING DATA AND PREPROCESSING FOR MODEL '", object@model_name,
                    "' NOW.\n", sep = '')
            if (check_data) {
              data <- try(force(data))
              if (is(data, "try-error")) {
                message("failed to evaluate the data; sampling not done")
                return(invisible(new_empty_stanfit(object)))
              }
              if (!missing(data) && length(data) > 0) {
                data <- try(data_preprocess(data))
                if (is(data, "try-error")) {
                  message("failed to preprocess the data; sampling not done")
                  return(invisible(new_empty_stanfit(object)))
                }
              } else data <- list()
            }
            if (verbose)
              cat('\n', "COMPILING MODEL '", object@model_name,
                  "' NOW.\n", sep = '')
            dots <- list(...)
            data$CHAIN_ID <- dots$chain_id
            if (verbose)
              cat('\n', "STARTING SAMPLER FOR MODEL '", object@model_name,
                  "' NOW.\n", sep = '')
            sfmiscenv <- new.env(parent = emptyenv())
            stan_fit_cpp_module <- object@mk_cppmodule(object)
            cxxfun <- grab_cxxfun(object@dso)
            if (test_221(object@model_cpp$model_cppcode) &&
                stan_fit_cpp_module@constructors[[1]]$nargs == 2L) {
              mod <- try(new(stan_fit_cpp_module, data, as.integer(seed)))
              if (is(mod, "try-error")) {
                message('failed to create the sampler; sampling not done')
                return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv)))
              }
              sampler <- try(new(stan_fit, mod$fit_ptr()))
              if (is(sampler, "try-error")) {
                message('failed to create the sampler; sampling not done')
                return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv)))
              }
            } else {
              sampler <- try(new(stan_fit_cpp_module, data, as.integer(seed), cxxfun))
              if (is(sampler, "try-error")) {
                message('failed to create the sampler; sampling not done')
                return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv)))
              }
            }

            assign("stan_fit_instance", sampler, envir = sfmiscenv)
            m_pars = sampler$param_names()
            p_dims = sampler$param_dims()
            mode <- if (!is.null(dots$test_grad) && dots$test_grad)
              "TESTING GRADIENT" else "SAMPLING"

            if (is.numeric(init)) init <- as.character(init)
            if (is.function(init)) {
              if ("chain_id" %in% names(formals(init)))
                init <- lapply(1:chains, FUN = init)
              else
                init <- lapply(1:chains, function(id) init())
            }
            if (!is.list(init) && !is.character(init)) {
              message("wrong specification of initial values")
              return(invisible(new_empty_stanfit(object)))
            }
            if (is.list(list)) init <- lapply(init, function(x) x)

            if (cores > 1 && mode == "SAMPLING" && chains > 1) {
              .dotlist <- c(sapply(objects, simplify = FALSE, FUN = get,
                                  envir = environment()), list(...))
              .dotlist$chains <- 1L
              .dotlist$cores <- 0L
              .dotlist$open_progress <- FALSE
              callFun <- function(i) {
                .dotlist$chain_id <- i
                if(is.list(.dotlist$init)) .dotlist$init <- .dotlist$init[i]
                if(is.character(.dotlist$sample_file)) {
                  if (grepl("\\.csv$", .dotlist$sample_file))
                    .dotlist$sample_file <- sub("\\.csv$", paste0("_", i, ".csv"),
                                                .dotlist$sample_file)
                  else .dotlist$sample_file <- paste0(.dotlist$sample_file,
                                                      "_", i, ".csv")
                }
                if(is.character(.dotlist$diagnostic_file)) {
                  if (grepl("\\.csv$", .dotlist$diagnostic_file))
                    .dotlist$diagnostic_file <- sub("\\.csv$", paste0("_", i, ".csv"),
                                                    .dotlist$diagnostic_file)
                  else
                    .dotlist$diagnostic_file <- paste0(.dotlist$diagnostic_file,
                                                       "_", i, ".csv")
                }
                out <- do.call(rstan::sampling, args = .dotlist)
                return(out)
              }
              if ( .Platform$OS.type == "unix" &&
                   (!interactive() || isatty(stdout())) ) {
                nfits <- parallel::mclapply(1:chains, FUN = callFun,
                                            mc.preschedule = FALSE,
                                            mc.cores = min(chains, cores))
              }
              else {
                tfile <- tempfile()
                sinkfile <- paste0(tfile, "_StanProgress.txt")
                cat("Click the Refresh button to see progress of the chains\n", file = sinkfile)
                if (open_progress &&
                    !identical(browser <- getOption("browser"), "false")) {
                  if (identical(Sys.getenv("RSTUDIO"), "1"))
                    stop("you cannot specify 'open_progress = TRUE' when using RStudio")
                  sinkfile_html <- paste0(tfile, "_StanProgress.html")
                  create_progress_html_file(sinkfile_html, sinkfile)
                  utils::browseURL(paste0("file://", sinkfile_html))
                }
                else if (identical(Sys.getenv("RSTUDIO"), "1") &&
                         .Platform$OS.type == "windows" && interactive()) {
                  if (!requireNamespace("rstudioapi"))
                    stop("must install the rstudioapi package when using RStan in parallel via RStudio")
                  if (rstudioapi::isAvailable("0.98.423")) {
                    v <- rstudioapi::getVersion()
                    if (v > "0.99.100") rstudioapi::viewer(sinkfile, height = "maximize")
                    else if (v >= "0.98.423") rstudioapi::viewer(sinkfile)
                  }
                  else sinkfile <- ""
                }
                else sinkfile <- ""
                if (!is.null(dots$refresh) && dots$refresh == 0)
                  cl <- parallel::makeCluster(min(cores, chains), useXDR = FALSE,
                                              setup_strategy = "sequential")
                else
                  cl <- parallel::makeCluster(min(cores, chains), outfile = sinkfile, 
                                              useXDR = FALSE, setup_strategy = "sequential")
                on.exit(parallel::stopCluster(cl))
                dependencies <- c("rstan", "Rcpp", "ggplot2")
                .paths <- unique(c(.libPaths(), sapply(dependencies, FUN = function(d) {
                  dirname(system.file(package = d))
                })))
                .paths <- .paths[.paths != ""]
                parallel::clusterExport(cl, varlist = ".paths", envir = environment())
                parallel::clusterEvalQ(cl, expr = .libPaths(.paths))
                parallel::clusterEvalQ(cl, expr =
                                      suppressPackageStartupMessages(require(rstan, quietly = TRUE)))
                parallel::clusterExport(cl, varlist = ".dotlist", envir = environment())
                data_e <- as.environment(data)
                parallel::clusterExport(cl, varlist = names(data_e), envir = data_e)
                nfits <- parallel::parLapplyLB(cl, X = 1:chains, fun = callFun)
              }
              valid <- sapply(nfits, is, class2 = "stanfit") &
                       sapply(nfits, FUN = function(x) x@mode == 0)
              if(all(valid)) {
                nfits <- sflist2stanfit(nfits)
                nfits@.MISC <- sfmiscenv
                throw_sampler_warnings(nfits)
                return(nfits)
              }
              else {
                warning("some chains had errors; consider specifying chains = 1 to debug")
                message("here are whatever error messages were returned")
                print(nfits[!valid])
                if(any(valid)) {
                  nfits <- sflist2stanfit(nfits[valid])
                  nfits@.MISC <- sfmiscenv
                  throw_sampler_warnings(nfits)
                  return(nfits)
                }
                return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv,
                                                   m_pars, p_dims, 2L)))
              }
            }
            check_unknown_args <- dots$check_unknown_args
            if (is.null(check_unknown_args) || check_unknown_args) {
              is_arg_recognizable(names(dots),
                                  c("chain_id", "init_r", "test_grad",
                                    "obfuscate_model_name",
                                    "enable_random_init",
                                    "append_samples", "refresh", "control",
                                    "include", "cores", "open_progress",
                                    "save_warmup"),
                                  pre_msg = "passing unknown arguments: ",
                                  call. = FALSE)
            }
            if(!include) {
              if (length(pars) == 1 && is.na(pars)) pars <- "lp__"
              else pars <- setdiff(m_pars, pars)
              if (length(pars) == 0) pars <- "lp__"
            }

            if (!missing(pars) && !any(is.na(pars)) && length(pars) > 0) {
              sampler$update_param_oi(pars)
              m <- which(match(pars, m_pars, nomatch = 0) == 0)
              if (length(m) > 0) {
                message("no parameter ", paste(pars[m], collapse = ', '), "; sampling not done")
                return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv, m_pars, p_dims, 2L)))
              }
            }

            if (chains < 1) {
              message("the number of chains is less than 1; sampling not done")
              return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv, m_pars, p_dims, 2L)))
            }

            args_list <- try(config_argss(chains = chains, iter = iter,
                                          warmup = warmup, thin = thin,
                                          init = init, seed = seed, sample_file = sample_file,
                                          diagnostic_file = diagnostic_file,
                                          algorithm = match.arg(algorithm), control = control, ...))

            if (is(args_list, "try-error")) {
              message('error in specifying arguments; sampling not done')
              return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv, m_pars, p_dims, 2L)))
            }

            # number of samples saved after thinning
            warmup2 <- 1 + (warmup - 1) %/% thin
            if (!is.null(dots$save_warmup) && !dots$save_warmup)
              warmup2 <- 0L
            n_kept <- 1 + (iter - warmup - 1) %/% thin
            n_save <- n_kept + warmup2

            samples <- vector("list", chains)

            for (i in 1:chains) {
              cid <- args_list[[i]]$chain_id
              if (is.null(dots$refresh) || dots$refresh > 0) {
                cat('\n', mode, " FOR MODEL '", object@model_name,
                    "' NOW (CHAIN ", cid, ").\n", sep = '')
              }
              if (is.character(show_messages))
                messages <- normalizePath(show_messages, mustWork = FALSE)
              else messages <- tempfile()
              mfile <- file(messages, open = "wt")
              sink(mfile, type = "message")
              samples_i <- try(sampler$call_sampler(args_list[[i]]))
              sink(NULL, type = "message")
              if (!file.exists(messages)) {
                report <- "log file disappeared"
              }
              else {
                close(mfile)
                report <- scan(file = messages, what = character(),
                               sep = "\n", quiet = TRUE)
                unlink(messages)
              }
              if (is(samples_i, "try-error") || is.null(samples_i)) {
                print(report)
                msg <- "error occurred during calling the sampler; sampling not done"
                if (.Platform$OS.type == "windows") print(msg)
                else message(msg)
                return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv,
                                                   m_pars, p_dims, 2L)))
              }
              if (length(report) > 0 && isTRUE(show_messages)) {
                end <- unique(grep("if ", report, ignore.case = TRUE, value = TRUE))
                report <- grep("if ", report, ignore.case = TRUE, value = TRUE, invert = TRUE)
                report <- gsub(" because of the following issue:", "", report, fixed = TRUE)
                report <- grep("^Exception thrown at line", report, value = TRUE)
                report <- gsub("stan::math::", "", report, fixed = TRUE)
                report <- strtrim(report, width = 100)
                if (length(report) > 0) {
                  tab <- sort(table(report), decreasing = TRUE)
                  msg <- paste("The following numerical problems occurred",
                               "the indicated number of times on chain",
                               cid)
                  if (.Platform$OS.type == "windows") print(msg)
                  else message(msg)
                  mat <- as.matrix(tab)
                  colnames(mat) <- "count"
                  if (.Platform$OS.type == "windows") {
                    print(mat)
                    print("When a numerical problem occurs, the Hamiltonian proposal gets rejected.")
                    print("See https://mc-stan.org/misc/warnings.html#exception-hamiltonian-proposal-rejected")
                    print(paste("If the number in the 'count' column is small, ",
                                "there is no need to ask about this message on stan-users."))
                  }
                  else {
                    message(paste(capture.output(print(mat)), collapse = "\n"))
                    message("When a numerical problem occurs, the Hamiltonian proposal gets rejected.")
                    message("See https://mc-stan.org/misc/warnings.html#exception-hamiltonian-proposal-rejected")
                    message("If the number in the 'count' column is small, ",
                            "there is no need to ask about this message on stan-users.")
                  }
                }
              }
              samples[[i]] <- samples_i
            }

            idx_wo_lp <- which(m_pars != 'lp__')
            skeleton <- create_skeleton(m_pars[idx_wo_lp], p_dims[idx_wo_lp])
            inits_used = lapply(lapply(samples, function(x) attr(x, "inits")),
                                function(y) rstan_relist(y, skeleton))

            # test_gradient mode: no sample
            if (attr(samples[[1]], 'test_grad')) {
              sim = list(num_failed = sapply(samples, function(x) x$num_failed))
              return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv,
                                                 m_pars, p_dims, 1L, sim = sim,
                                                 inits = inits_used,
                                                 stan_args = args_list)))
            }

            # perm_lst <- lapply(1:chains, function(id) rstan_seq_perm(n_kept, chains, seed, chain_id = id))
            ## sample_int is a little bit faster than our own rstan_seq_perm (one
            ## reason is that the RNG used in R is faster),
            ## but without controlling the seed
            perm_lst <- lapply(1:chains, function(id) sample.int(n_kept))

            fnames_oi <- sampler$param_fnames_oi()
            n_flatnames <- length(fnames_oi)
            sim = list(samples = samples,
                       iter = iter, thin = thin,
                       warmup = warmup,
                       chains = chains,
                       n_save = rep(n_save, chains),
                       warmup2 = rep(warmup2, chains), # number of warmup iters in n_save
                       permutation = perm_lst,
                       pars_oi = sampler$param_names_oi(),
                       dims_oi = sampler$param_dims_oi(),
                       fnames_oi = fnames_oi,
                       n_flatnames = n_flatnames)
            nfit <- new("stanfit",
                        model_name = object@model_name,
                        model_pars = m_pars,
                        par_dims = p_dims,
                        mode = 0L,
                        sim = sim,
                        # keep a record of the initial values
                        inits = inits_used,
                        stan_args = args_list,
                        stanmodel = object,
                          # keep a ref to avoid garbage collection
                          # (see comments in fun stan_model)
                        date = date(),
                        .MISC = sfmiscenv)
            if (cores > 0) throw_sampler_warnings(nfit)
            return(nfit)
          })


setGeneric(name = "gqs",
           def = function(object, ...) { standardGeneric("gqs") })
setMethod("gqs", "stanmodel",
          function(object, data = list(), draws,
                seed = sample.int(.Machine$integer.max, size = 1L)) {
  draws <- as.matrix(draws)
  objects <- ls()
  if (is.list(data) & !is.data.frame(data)) {
    parsed_data <- with(data, parse_data(get_cppcode(object)))
    if (!is.list(parsed_data)) {
      message("failed to get names of data from the model; sampling not done")
      return(invisible(new_empty_stanfit(object)))
    }
    data <- parsed_data
  } else if (is.character(data)) { # names of objects
    data <- try(mklist(data))
    if (is(data, "try-error")) {
      message("failed to create the data; sampling not done")
      return(invisible(new_empty_stanfit(object)))
    }
  }
  if (TRUE) { # check_data
    data <- try(force(data))
    if (is(data, "try-error")) {
      message("failed to evaluate the data; sampling not done")
      return(invisible(new_empty_stanfit(object)))
    }
    if (!missing(data) && length(data) > 0) {
      data <- try(data_preprocess(data))
      if (is(data, "try-error")) {
        message("failed to preprocess the data; sampling not done")
        return(invisible(new_empty_stanfit(object)))
      }
    } else data <- list()
  }
  stan_fit_cpp_module <- object@mk_cppmodule(object)
  cxxfun <- grab_cxxfun(object@dso)
  sfmiscenv <- new.env(parent = emptyenv())
  stan_fit_cpp_module <- object@mk_cppmodule(object)
  cxxfun <- grab_cxxfun(object@dso)
  if (test_221(object@model_cpp$model_cppcode) &&
      stan_fit_cpp_module@constructors[[1]]$nargs == 2L) {
    mod <- try(new(stan_fit_cpp_module, data, as.integer(seed)))
    if (is(mod, "try-error")) {
      message('failed to create the sampler; sampling not done')
      return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv)))
    }
    sampler <- try(new(stan_fit, mod$fit_ptr()))
    if (is(sampler, "try-error")) {
      message('failed to create the sampler; sampling not done')
      return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv)))
    }
  } else {
    sampler <- try(new(stan_fit_cpp_module, data, as.integer(seed), cxxfun))
    if (is(sampler, "try-error")) {
      message('failed to create the sampler; sampling not done')
      return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv)))
    }
  }
  
  assign("stan_fit_instance", sampler, envir = sfmiscenv)
  m_pars = sampler$param_names()
  p_dims = sampler$param_dims()
  p_names <- unique(sub("\\..*$", "", sampler$constrained_param_names(TRUE, FALSE)))
  all_names <- sampler$constrained_param_names(TRUE, TRUE)
  some_names <- sampler$constrained_param_names(TRUE, FALSE)
  draws_colnames <- sub("\\.", "[", some_names)
  draws_colnames <- gsub("\\.", ",", draws_colnames)
  draws_colnames[grep("\\[", draws_colnames)] <- 
    paste0(draws_colnames[grep("\\[", draws_colnames)], "]")
  draws <- draws[, draws_colnames, drop = FALSE]
  gq_names <- unique(sub("\\..*$", "", setdiff(all_names, some_names)))
  sampler$update_param_oi(gq_names)
  samples <- try(sampler$standalone_gqs(draws, as.integer(seed)))
  if (is(samples, "try-error") || is.null(samples)) {
    msg <- "error occurred during calling the sampler; sampling not done"
    message(msg)
    return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv,
                                       m_pars, p_dims, 2L)))
  }

  skeleton <- create_skeleton(gq_names, p_dims[gq_names])

  perm_lst <- list(1:nrow(draws)) # not actually permuted

  fnames_oi <- setdiff(sampler$param_fnames_oi(), "lp__")
  n_flatnames <- length(fnames_oi)
  sim = list(samples = list(samples),
             iter = nrow(draws), thin = 1L,
             warmup = 0L,
             chains = 1L,
             n_save = nrow(draws),
             warmup2 = 0L, # number of warmup iters in n_save
             permutation = perm_lst,
             pars_oi = gq_names,
             dims_oi = sampler$param_dims_oi()[gq_names],
             fnames_oi = fnames_oi,
             n_flatnames = n_flatnames)
  nfit <- new("stanfit",
              model_name = object@model_name,
              model_pars = gq_names,
              par_dims = p_dims[gq_names],
              mode = 0L,
              sim = sim,
              # keep a record of the initial values
              inits = list(),
              stan_args = list(list(seed = seed)),
              stanmodel = object,
              # keep a ref to avoid garbage collection
              # (see comments in fun stan_model)
              date = date(),
              .MISC = sfmiscenv)
  return(nfit)
})
