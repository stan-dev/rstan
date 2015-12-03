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

setMethod("show", "stanmodel",
          function(object) {
            cat("S4 class stanmodel '", object@model_name, "' coded as follows:\n" ,sep = '') 
            cat(object@model_code, "\n")
          }) 

setGeneric(name = 'optimizing',
           def = function(object, ...) { standardGeneric("optimizing")})

setGeneric(name = 'vb',
           def = function(object, ...) { standardGeneric("vb")})

setGeneric(name = "sampling",
           def = function(object, ...) { standardGeneric("sampling")})

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
               "  - not existed (created from reading csv files).", sep = '')) 
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
                   check_data = TRUE, sample_file = tempfile(fileext = '.csv'),
                   algorithm = c("meanfield", "fullrank"), ...) {
            stan_fit_cpp_module <- object@mk_cppmodule(object)
            if (is.list(data) & !is.data.frame(data)) {
              parsed_data <- parse_data(get_cppcode(object))
              if (!is.list(parsed_data)) {
                message("failed to get names of data from the model; sampling not done")
                return(invisible(new_empty_stanfit(object)))
              }
              for (nm in names(data)) parsed_data[[nm]] <- data[[nm]]
              parsed_data <- parsed_data[!sapply(parsed_data, is.null)]
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
            sampler <- try(new(stan_fit_cpp_module, data, object@dso@.CXXDSOMISC$cxxfun))
            if (is(sampler, "try-error")) {
              message('failed to create the model; variational Bayes not done')
              return(invisible(new_empty_stanfit(object)))
            }
            seed <- check_seed(seed, warn = 1)
            if (is.null(seed))
              return(invisible(list(stanmodel = object)))
            args <- list(seed = seed, chain_id = 1L,
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
                                  "tol_rel_obj"),
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
            if(!include) pars <- setdiff(m_pars, pars)
            
            if (!missing(pars) && !is.na(pars) && length(pars) > 0) {
              sampler$update_param_oi(pars)
              m <- which(match(pars, m_pars, nomatch = 0) == 0)
              if (length(m) > 0) {
                message("no parameter ", paste(pars[m], collapse = ', '), "; sampling not done") 
                return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv, m_pars, p_dims, 2L))) 
              }
            }
            else pars <- m_pars

            skeleton <- create_skeleton(m_pars, p_dims)
            cC <- unlist(lapply(names(skeleton), FUN = function(x) {
              param <- skeleton[[x]]
              if (x == "lp__") TRUE
              else if (x %in% pars) rep(TRUE, length(param))
              else rep(FALSE, length(param))
            }))

            vbres <- sampler$call_sampler(c(args, dotlist))
            samples <- read_one_stan_csv(attr(vbres, "args")$sample_file)
            means <- sapply(samples, mean)
            means <- as.matrix(c(means[-1], "lp__" = means[1]))
            colnames(means) <- "chain:1"
            assign("posterior_mean_4all", means, envir = sfmiscenv)
            idx_wo_lp <- which(m_pars != "lp__")
            skeleton <- create_skeleton(m_pars[idx_wo_lp], p_dims[idx_wo_lp])
            inits_used <- rstan_relist(as.numeric(samples[1,]), skeleton)
            samples <- cbind(samples[,-1,drop=FALSE], "lp__" = samples[,1])[,cC]
            
            fnames_oi <- sampler$param_fnames_oi()
            n_flatnames <- length(fnames_oi)
            iter <- nrow(samples)
            sim = list(samples = list(as.list(samples)),
                       iter = iter, thin = 1L,
                       warmup = 0L,
                       chains = 1L,
                       n_save = iter,
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
                   draws = 0, constrained = TRUE, ...) {
            stan_fit_cpp_module <- object@mk_cppmodule(object)

            if (is.list(data) & !is.data.frame(data)) {
              parsed_data <- parse_data(get_cppcode(object))
              if (!is.list(parsed_data)) {
                message("failed to get names of data from the model; sampling not done")
                return(invisible(new_empty_stanfit(object)))
              }
              for (nm in names(data)) parsed_data[[nm]] <- data[[nm]]
              parsed_data <- parsed_data[!sapply(parsed_data, is.null)]
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
            sampler <- try(new(stan_fit_cpp_module, data, cxxfun))
            if (is(sampler, "try-error")) {
              message('failed to create the optimizer; optimization not done') 
              return(invisible(list(stanmodel = object)))
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
            names(optim$par) <- flatnames(m_pars, p_dims, col_major = TRUE)
            skeleton <- create_skeleton(m_pars, p_dims)
            if (hessian || draws) {
              fn <- function(theta) {
                sampler$log_prob(theta, FALSE, FALSE)
              }
              gr <- function(theta) {
                sampler$grad_log_prob(theta, FALSE)
              }
              theta <- rstan_relist(optim$par, skeleton)
              theta <- sampler$unconstrain_pars(theta)
              H <- optimHess(theta, fn, gr, control = list(fnscale = -1))
              colnames(H) <- rownames(H) <- sampler$unconstrained_param_names(FALSE, FALSE)
              if (hessian) optim$hessian <- H
              if (draws > 0 && is.matrix(R <- try(chol(-H)))) {
                K <- ncol(R)
                R_inv <- backsolve(R, diag(K))
                Z <- matrix(rnorm(K * draws), K, draws)
                theta_tilde <- t(theta + R_inv %*% Z)
                if (constrained) {
                  theta_tilde <- t(apply(theta_tilde, 1, FUN = function(theta) {
                    sampler$constrain_pars(theta)  
                  }))
                  colnames(theta_tilde) <- names(optim$par)
                }
                else {
                  log_p <- apply(theta_tilde, 1, FUN = function(theta) {
                    sampler$log_prob(theta, adjust_transform = TRUE, gradient = FALSE)
                  })
                  log_g <- colSums(dnorm(Z, log = TRUE)) - sum(log(diag(R_inv)))
                  optim$log_p <- log_p
                  optim$log_g <- log_g
                  colnames(theta_tilde) <- colnames(H)
                  optim$log_prob <- sampler$log_prob
                  optim$grad_log_prob <- sampler$grad_log_prob
                }
                optim$theta_tilde <- theta_tilde
              }
            }
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
            objects <- ls()
            if (is.list(data) & !is.data.frame(data)) {
              parsed_data <- try(parse_data(get_cppcode(object)))
              if (!is.list(parsed_data)) {
                message("failed to get names of data from the model; sampling not done")
                return(invisible(new_empty_stanfit(object)))
              }
              for (nm in names(data)) parsed_data[[nm]] <- data[[nm]]
              parsed_data <- parsed_data[!sapply(parsed_data, is.null)]
              data <- parsed_data
            } else if (is.character(data)) { # names of objects
              data <- try(mklist(data))
              if (is(data, "try-error")) {
                message("failed to create the data; sampling not done")
                return(invisible(new_empty_stanfit(object)))
              }
            }
            # check data and preprocess
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
            stan_fit_cpp_module <- object@mk_cppmodule(object)
            cxxfun <- grab_cxxfun(object@dso)
            sampler <- try(new(stan_fit_cpp_module, data, cxxfun))
            sfmiscenv <- new.env(parent = emptyenv())
            if (is(sampler, "try-error")) {
              message('failed to create the sampler; sampling not done') 
              return(invisible(new_empty_stanfit(object, miscenv = sfmiscenv)))
            } 
            assign("stan_fit_instance", sampler, envir = sfmiscenv)
            m_pars = sampler$param_names()
            p_dims = sampler$param_dims()
            dots <- list(...)
            mode <- if (!is.null(dots$test_grad) && dots$test_grad) 
              "TESTING GRADIENT" else "SAMPLING"
            
            if (cores > 1 && mode == "SAMPLING" && chains > 1) {
              .dotlist <- c(sapply(objects, simplify = FALSE, FUN = get,
                                  envir = environment()), list(...))
              .dotlist$chains <- 1L
              .dotlist$cores <- 1L
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
              cl <- parallel::makeCluster(min(cores, chains), 
                                          outfile = sinkfile, useXDR = FALSE)
              on.exit(parallel::stopCluster(cl))
              dependencies <- read.dcf(file = system.file("DESCRIPTION", package = "rstan"), 
                                       fields = "Imports")[1,]
              dependencies <- scan(what = character(), sep = ",", strip.white = TRUE, 
                                   quiet = TRUE, text = dependencies)
              dependencies <- c("Rcpp", "rstan", "rstanarm", dependencies)
              .paths <- unique(sapply(dependencies, FUN = function(d) {
                dirname(system.file(package = d))
              }))
              .paths <- .paths[.paths != ""]
              parallel::clusterExport(cl, varlist = ".paths", envir = environment())
              parallel::clusterEvalQ(cl, expr = .libPaths(.paths))
              parallel::clusterEvalQ(cl, expr = 
                                    suppressPackageStartupMessages(require(rstan, quietly = TRUE)))
              callFun <- function(i) {
                .dotlist$chain_id <- i
                if(is.list(.dotlist$init)) .dotlist$init <- .dotlist$init[i]
                if(is.character(.dotlist$sample_file)) {
                  .dotlist$sample_file <- paste0(.dotlist$sample_file, i)
                }
                if(is.character(.dotlist$diagnostic_file)) {
                  .dotlist$diagnostic_file <- paste0(.dotlist$diagnostic_file, i)
                }
                Sys.sleep(0.5 * i)
                out <- do.call(rstan::sampling, args = .dotlist)
                return(out)
              }
              parallel::clusterExport(cl, varlist = ".dotlist", envir = environment())
              data_e <- as.environment(data)
              parallel::clusterExport(cl, varlist = names(data_e), envir = data_e)
              nfits <- parallel::parLapply(cl, X = 1:chains, fun = callFun)
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
                                    "include", "cores", "open_progress"), 
                                  pre_msg = "passing unknown arguments: ",
                                  call. = FALSE)
            }

            if(!include) pars <- setdiff(m_pars, pars)
            
            if (!missing(pars) && !is.na(pars) && length(pars) > 0) {
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
              close(mfile)
              report <- scan(file = messages, what = character(),
                             sep = "\n", quiet = TRUE)
              if (is(samples_i, "try-error") || is.null(samples_i)) {
                print(tail(report, n = 10))
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
                report <- gsub("stan::math::", "", report, fixed = TRUE)
                report <- grep("^Informational", report, value = TRUE, invert = TRUE)
                report <- grep("^[[:digit:]]+", report, value = TRUE, invert = TRUE)
                report <- strtrim(report, width = 100)
                if (length(report) > 0) {
                  tab <- sort(table(report), decreasing = TRUE)
                  msg <- paste("The following numerical problems occured",
                               "the indicated number of times after warmup on chain", 
                               cid)
                  if (.Platform$OS.type == "windows") print(msg)
                  else message(msg)
                  mat <- as.matrix(tab)
                  colnames(mat) <- "count"
                  if (.Platform$OS.type == "windows") {
                    print(mat)
                    print("When a numerical problem occurs, the Metropolis proposal gets rejected.")
                    print(paste("However, by design Metropolis proposals sometimes get rejected ", 
                                "even when there are no numerical problems."))
                    print(paste("Thus, if the number in the 'count' column is small, ",
                                "do not ask about this message on stan-users."))
                  }
                  else {
                    message(paste(capture.output(print(mat)), collapse = "\n"))
                    message("When a numerical problem occurs, the Metropolis proposal gets rejected.")
                    message("However, by design Metropolis proposals sometimes get rejected ", 
                            "even when there are no numerical problems.")
                    message("Thus, if the number in the 'count' column is small, ",
                            "do not ask about this message on stan-users.")
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
            if (interactive() && cores <= 1) throw_sampler_warnings(nfit)
            return(nfit)
          }) 

