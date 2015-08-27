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

# require(methods)
# source('AllClass.R')

setMethod("show", "stanfit", 
          function(object) {
            print.stanfit(x = object, pars = object@sim$pars_oi)
          })  

print.stanfit <- function(x, pars = x@sim$pars_oi, 
                          probs = c(0.025, 0.25, 0.5, 0.75, 0.975), 
                          digits_summary = 2, include = TRUE, ...) { 
  if (x@mode == 1L) { 
    cat("Stan model '", x@model_name, "' is of mode 'test_grad';\n",
        "sampling is not conducted.\n", sep = '')
    return(invisible(NULL)) 
  } else if (x@mode == 2L) {
    cat("Stan model '", x@model_name, "' does not contain samples.\n", sep = '') 
    return(invisible(NULL)) 
  } 

  if(!include) pars <- setdiff(x@sim$pars_oi, pars)
  s <- summary(x, pars, probs, ...)  
  if (is.null(s)) return(invisible(NULL))
  n_kept <- x@sim$n_save - x@sim$warmup2
  cat("Inference for Stan model: ", x@model_name, '.\n', sep = '')
  cat(x@sim$chains, " chains, each with iter=", x@sim$iter, 
      "; warmup=", x@sim$warmup, "; thin=", x@sim$thin, "; \n", 
      "post-warmup draws per chain=", n_kept[1], ", ", 
      "total post-warmup draws=", sum(n_kept), ".\n\n", sep = '')

  if (!is.null(x@stan_args[[1]]$method) && 
               x@stan_args[[1]]$method == "variational") {
    print(round(s$summary, digits_summary), ...) 
    cat("\nApproximate samples were drawn using VB(", x@stan_args[[1]]$algorithm, ") at ", x@date, 
        ".\n", sep = '')
    message("We recommend genuine 'sampling' from the posterior distribution for final inferences!")
    return(invisible(NULL))
  }
  
  # round n_eff to integers
  s$summary[, 'n_eff'] <- round(s$summary[, 'n_eff'], 0)

  print(round(s$summary, digits_summary), ...) 

  sampler <- attr(x@sim$samples[[1]], "args")$sampler_t

  cat("\nSamples were drawn using ", sampler, " at ", x@date, ".\n",
      "For each parameter, n_eff is a crude measure of effective sample size,\n", 
      "and Rhat is the potential scale reduction factor on split chains (at \n",
      "convergence, Rhat=1).\n", sep = '')
  return(invisible(NULL)) 
}  

setMethod("plot", signature(x = "stanfit", y = "missing"), 
          function(x, ..., pars, include = TRUE, unconstrain = FALSE) {
            if (x@mode == 1L) {
              cat("Stan model '", x@model_name, "' is of mode 'test_grad';\n",
                  "sampling is not conducted.\n", sep = '')
              return(invisible(NULL)) 
            } else if (x@mode == 2L) {
              cat("Stan model '", x@model_name, "' does not contain samples.\n", sep = '') 
              return(invisible(NULL)) 
            } 

            if (isTRUE(all.equal(x@sim$n_save, x@sim$warmup2,
                check.attributes = FALSE, check.names = FALSE))) {
              cat("Stan model '", x@model_name, "' does not contain samples after warmup.\n", sep = '')
              return(invisible(NULL))
            }
            args <- list(object = x, include = include, unconstrain = unconstrain, ...)
            if (!missing(pars)) { 
              if ("log-posterior" %in% pars)
                pars[which(pars == "log-posterior")] <- "lp__"
              args$pars <- pars
            }
            do.call("stan_plot", args)
          }) 

setGeneric(name = "get_stancode",
           def = function(object, ...) { standardGeneric("get_stancode")}) 

setGeneric(name = "get_cppo_mode", 
           def = function(object, ...) { standardGeneric("get_cppo_mode") }) 

setMethod('get_cppo_mode', signature = "stanfit", 
           function(object) { 
             cxxf <- get_cxxflags(object@stanmodel)
             if (identical(cxxf, character(0))) return(NA)
             l <- get_cxxo_level(cxxf)
             if ("" == l) l <- "0" 
             p <- match(l, c("3", "2", "1", "0")) 
             c("fast", "presentation2", "presentation1", "debug")[p]
           }) 

setMethod('get_stancode', signature = "stanfit", 
          function(object, print = FALSE) {
            code <- object@stanmodel@model_code
            if (print) cat(code, "\n") 
            return(code)
          }) 

setGeneric(name = 'get_stanmodel', 
           def = function(object, ...) { standardGeneric("get_stanmodel")})

setMethod("get_stanmodel", signature = "stanfit", 
          function(object) { 
            return(object@stanmodel)
          }) 

setGeneric(name = 'get_inits', 
           def = function(object, ...) { standardGeneric("get_inits")})

setMethod("get_inits", signature = "stanfit", 
          function(object) { return(object@inits) })

setGeneric(name = 'get_seed', 
           def = function(object, ...) { standardGeneric("get_seed")})

setMethod("get_seed", signature = "stanfit", 
          function(object) { 
            if (length(object@stan_args) < 1L) return(NULL) 
            object@stan_args[[1]]$seed })

setGeneric(name = 'get_seeds', 
           def = function(object, ...) { standardGeneric("get_seeds")})

setMethod("get_seeds", signature = "stanfit", 
          function(object) { 
            if (length(object@stan_args) < 1L) return(NULL) 
            sapply(object@stan_args, function(x) x$seed) })

### HELPER FUNCTIONS
### 


get_kept_samples <- function(n, sim) {
  #
  # Args:
  #  sim: the sim slot in object stanfit 
  #  n: the nth parameter (starting from 1) 
  # Note: 
  #  samples from different chains are merged. 

  # get chain kept samples (gcks) 
  gcks <- function(s, nw, permutation) {
    s <- s[[n]][-(1:nw)] 
    s[permutation] 
  } 
  ss <- mapply(gcks, sim$samples, sim$warmup2, sim$permutation,
               SIMPLIFY = FALSE, USE.NAMES = FALSE) 
  do.call(c, ss) 
} 

get_kept_samples2 <- function(n, sim) {
  # a different implementation of get_kept_samples 
  # It seems this one is faster than get_kept_samples 
  # TODO: to understand why it is faster? 
  lst <- vector("list", sim$chains)
  for (ic in 1:sim$chains) { 
    if (sim$warmup2[ic] > 0) 
      lst[[ic]] <- sim$samples[[ic]][[n]][-(1:sim$warmup2[ic])][sim$permutation[[ic]]]
    else 
      lst[[ic]] <- sim$samples[[ic]][[n]][sim$permutation[[ic]]]
  } 
  do.call(c, lst)
}

get_samples <- function(n, sim, inc_warmup = TRUE) {
  # get chain samples
  # Args:
  #   n: parameter index (integer)
  #   sim: the sim list in stanfit 
  # 
  # Returns:
  #   a list of chains for the nth parameter; each chain is an 
  #   element of the list.  
  if (all(sim$warmup2 == 0)) inc_warmup <- TRUE # for the case warmup sample is discarded
  gcs <- function(s, inc_warmup, nw) {
    # Args:
    #   s: samples of all chains 
    #   nw: number of warmup 
    if (inc_warmup)  return(s[[n]])
    else return(s[[n]][-(1:nw)]) 
  } 
  ss <- mapply(gcs, sim$samples, inc_warmup, sim$warmup2, 
               SIMPLIFY = FALSE, USE.NAMES = FALSE) 
  ss 
} 

get_samples2 <- function(n, sim, inc_warmup = TRUE) {
  # serves the same purpose with get_samples, but with 
  # different implementation 
  # It seems that this one is fast. 
  if (all(sim$warmup2 == 0)) inc_warmup <- TRUE # for the case warmup sample is discarded
  npar <- length(n)
  lst <- vector("list", sim$chains)
  for (ic in 1:sim$chains) {
    lst[[ic]] <- 
      if (inc_warmup) sim$samples[[ic]][[n]] else sim$samples[[ic]][[n]][-(1:sim$warmup2[ic])]
  }
  lst
} 

par_traceplot <- function(sim, n, par_name, inc_warmup = TRUE, window = NULL, ...) {
  # same thin, n_save, warmup2 for all the chains
  thin <- sim$thin
  warmup2 <- sim$warmup2[1] 
  warmup <- sim$warmup 
  main <- paste("Trace of ", par_name) 
  chain_cols <- rstan_options("rstan_chain_cols")
  warmup_col <- rstan_options("rstan_warmup_bg_col") 

  start_i <- window[1] 
  window_size <- (window[2] - start_i) %/% thin
  id <- seq.int(start_i, by = thin, length.out = window_size)
  start_idx <- (if (warmup2 == 0) (start_i - warmup) else start_i) %/% thin
  if (start_idx < 1)  start_idx <- 1
  idx <- seq.int(start_idx, by = 1, length.out = window_size)

  yrange <- NULL 
  for (i in 1:sim$chains)
    yrange <- range(yrange, sim$samples[[i]][[n]][idx]) 

  if (inc_warmup) {
    plot(range(id), yrange, type = 'n', bty = 'l',
         xlab = 'Iterations', ylab = "", main = main, ...)
    rect(par("usr")[1], par("usr")[3], sim$warmup, par("usr")[4], 
         col = warmup_col, border = NA)
  } else {  
    plot(range(id), yrange, type = 'n', bty = 'l',
         xlab = 'Iterations (without warmup)', ylab = "", main = main, ...)
  } 

  for (i in 1:sim$chains)  
    lines(id, sim$samples[[i]][[n]][idx], 
          xlab = '', ylab = '', col = chain_cols[(i-1) %% 6 + 1], ...) 
} 

######

setGeneric(name = 'get_adaptation_info', 
           def = function(object, ...) { standardGeneric("get_adaptation_info")})

setMethod("get_adaptation_info", 
          definition = function(object) {
            if (object@mode == 1L) {
              cat("Stan model '", object@model_name, "' is of mode 'test_grad';\n",
                  "sampling is not conducted.\n", sep = '')
              return(invisible(NULL)) 
            } else if (object@mode == 2L) {
              cat("Stan model '", object@model_name, "' does not contain samples.\n", sep = '') 
              return(invisible(NULL)) 
            } 

            lai <- lapply(object@sim$samples, function(x) attr(x, "adaptation_info"))
            is_empty <- function(x) { 
              if (is.null(x)) return(TRUE) 
              if (is.character(x) && all(nchar(x) == 0)) return(TRUE)
              FALSE
            }
            if (all(sapply(lai, FUN = is_empty))) return(invisible(NULL))  
            return(lai) 
          }) 

setGeneric(name = "get_logposterior", 
           def = function(object, ...) { standardGeneric("get_logposterior")})


setMethod("get_logposterior", "stanfit",
          definition = function(object, inc_warmup = TRUE) {
            if (object@mode == 1L) {
              cat("Stan model '", object@model_name, "' is of mode 'test_grad';\n",
                  "sampling is not conducted.\n", sep = '')
              return(invisible(NULL)) 
            } else if (object@mode == 2L) {
              cat("Stan model '", object@model_name, "' does not contain samples.\n", sep = '') 
              return(invisible(NULL)) 
            } 

            llp <- lapply(object@sim$samples, function(x) x[['lp__']]) 
            if (inc_warmup) return(llp)
            if (object@sim$warmup2 == 0) return(llp)
            return(mapply(function(x, w) x[-(1:w)], 
                             llp, object@sim$warmup2,
                             SIMPLIFY = FALSE, USE.NAMES = FALSE)) 
          }) 

setGeneric(name = 'get_sampler_params', 
           def = function(object, ...) { standardGeneric("get_sampler_params")}) 

setMethod("get_sampler_params", 
          definition = function(object, inc_warmup = TRUE) {
            if (object@mode == 1L) {
              cat("Stan model '", object@model_name, "' is of mode 'test_grad';\n",
                  "sampling is not conducted.\n", sep = '')
              return(invisible(NULL)) 
            } else if (object@mode == 2L) {
              cat("Stan model '", object@model_name, "' does not contain samples.\n", sep = '') 
              return(invisible(NULL)) 
            } 

            ldf <- lapply(object@sim$samples, 
                          function(x) do.call(cbind, attr(x, "sampler_params")))   
            if (all(sapply(ldf, is.null))) return(invisible(NULL))  
            if (inc_warmup) return(invisible(ldf)) 
            return(mapply(function(x, w) x[-(1:w), , drop = FALSE], 
                             ldf, object@sim$warmup2, 
                             SIMPLIFY = FALSE, USE.NAMES = FALSE)) 
          }) 

setGeneric(name = 'get_elapsed_time',
           def = function(object, ...) { standardGeneric("get_elapsed_time")})

setMethod("get_elapsed_time",
          definition = function(object, inc_warmup = TRUE) {
            if (object@mode == 1L) {
              cat("Stan model '", object@model_name, "' is of mode 'test_grad';\n",
                  "sampling is not conducted.\n", sep = '')
              return(invisible(NULL))
            } else if (object@mode == 2L) {
              cat("Stan model '", object@model_name, "' does not contain samples.\n", sep = '')
              return(invisible(NULL))
            }

            ltime <- lapply(object@sim$samples,
                            function(x) attr(x, "elapsed_time"))
            t <- do.call(rbind, ltime)
            if (is.null(t)) return(t)
            cids <- sapply(object@stan_args, function(x) x$chain_id)
            rownames(t) <- paste0("chain:", cids)
            t
          })

setGeneric(name = 'get_posterior_mean', 
           def = function(object, ...) { standardGeneric("get_posterior_mean")}) 

setMethod("get_posterior_mean", signature = "stanfit",
          definition = function(object, pars) {
            if (object@mode == 1L) {
              cat("Stan model '", object@model_name, "' is of mode 'test_grad';\n",
                  "sampling is not conducted.\n", sep = '')
              return(invisible(NULL)) 
            } else if (object@mode == 2L) {
              cat("Stan model '", object@model_name, "' does not contain samples.\n", sep = '') 
              return(invisible(NULL)) 
            } 
            fnames <- flatnames(object@model_pars, object@par_dims, col_major = TRUE)
            if (!exists("posterior_mean_4all", envir = object@.MISC, inherits = FALSE)) {
              mean_pars <- lapply(object@sim$samples, function(x) attr(x, "mean_pars"))
              mean_lp__ <- lapply(object@sim$samples, function(x) attr(x, "mean_lp__"))
              m <- rbind(do.call(cbind, mean_pars), do.call(cbind, mean_lp__))
              name_allchains <- NULL
              if (ncol(m) > 1) {
                m <- cbind(m, apply(m, 1, mean))
                name_allchains <- "mean-all chains"
              } 
              cids <- sapply(object@stan_args, function(x) x$chain_id)
              colnames(m) <- c(paste0('mean-chain:', cids), name_allchains)
              rownames(m) <- fnames
              assign("posterior_mean_4all", m, envir = object@.MISC)
            }
            pars <- if (missing(pars)) object@model_pars else check_pars(c(object@model_pars, fnames), pars)
            pars <- remove_empty_pars(pars, object@sim$dims_oi)
            tidx <- pars_total_indexes(object@model_pars, object@par_dims, fnames, pars) 
            tidx <- lapply(tidx, function(x) attr(x, "row_major_idx"))
            return(object@.MISC$posterior_mean_4all[unlist(tidx), , drop = FALSE])
          })

setGeneric(name = "extract",
           def = function(object, ...) { standardGeneric("extract")}) 

setMethod("extract", signature = "stanfit",
          definition = function(object, pars, permuted = TRUE, 
                                inc_warmup = FALSE, include = TRUE) {
            # Extract the samples in different forms for different parameters. 
            #
            # Args:
            #   object: the object of "stanfit" class 
            #   pars: the names of parameters (including other quantiles) 
            #   permuted: if TRUE, the returned samples are permuted without
            #     warming up. And all the chains are merged. 
            #   inc_warmup: if TRUE, warmup samples are kept; otherwise, 
            #     discarded. If permuted is TRUE, inc_warmup is ignored. 
            #   include: if FALSE interpret pars as those to exclude
            #
            # Returns:
            #   If permuted is TRUE, return an array (matrix) of samples with each
            #   column being the samples for a parameter. 
            #   If permuted is FALSE, return array with dimensions
            #   (# of iter (with or w.o. warmup), # of chains, # of flat parameters). 

            if (object@mode == 1L) {
              cat("Stan model '", object@model_name, "' is of mode 'test_grad';\n",
                  "sampling is not conducted.\n", sep = '')
              return(invisible(NULL)) 
            } else if (object@mode == 2L) {
              cat("Stan model '", object@model_name, "' does not contain samples.\n", sep = '') 
              return(invisible(NULL)) 
            } 

            if(!include) pars <- setdiff(object@sim$pars_oi, pars)
            pars <- if (missing(pars)) object@sim$pars_oi else check_pars_second(object@sim, pars) 
            pars <- remove_empty_pars(pars, object@sim$dims_oi)
            tidx <- pars_total_indexes(object@sim$pars_oi, 
                                       object@sim$dims_oi, 
                                       object@sim$fnames_oi, 
                                       pars) 

            n_kept <- object@sim$n_save - object@sim$warmup2
            fun1 <- function(par_i) {
              # sss <- sapply(tidx[[par_i]], get_kept_samples2, object@sim)
              # if (is.list(sss))  sss <- do.call(c, sss)
              # the above two lines are slower than the following line of code
              sss <- do.call(cbind, lapply(tidx[[par_i]], get_kept_samples2, object@sim)) 
              dim(sss) <- c(sum(n_kept), object@sim$dims_oi[[par_i]]) 
              dimnames(sss) <- list(iterations = NULL)
              sss 
            } 
           
            if (permuted) {
              slist <- lapply(pars, fun1) 
              names(slist) <- pars 
              return(slist) 
            } 

            tidx <- unlist(tidx, use.names = FALSE) 
            tidxnames <- object@sim$fnames_oi[tidx] 
            sss <- lapply(tidx, get_samples2, object@sim, inc_warmup) 
            sss2 <- lapply(sss, function(x) do.call(c, x))  # concatenate samples from different chains
            sssf <- unlist(sss2, use.names = FALSE) 
  
            n2 <- object@sim$n_save[1]  ## assuming all the chains have equal iter 
            if (!inc_warmup) n2 <- n2 - object@sim$warmup2[1] 
            dim(sssf) <- c(n2, object@sim$chains, length(tidx)) 
            cids <- sapply(object@stan_args, function(x) x$chain_id)
            dimnames(sssf) <- list(iterations = NULL, chains = paste0("chain:", cids), parameters = tidxnames)
            sssf 
          })  

setMethod("summary", signature = "stanfit", 
          function(object, pars, 
                   probs = c(0.025, 0.25, 0.50, 0.75, 0.975), use_cache = TRUE, ...) { 
            # Summarize the samples (that is, compute the mean, SD, quantiles, for 
            # the samples in all chains and chains individually after removing
            # warmup samples, and n_eff and split R hat for all the kept samples.)
            # 
            # Returns: 
            #   A list with elements:
            #   summary: the summary for all the kept samples 
            #   c_summary: the summary for individual chains. 
            # 
            # Note: 
            #   This function is not straight in terms of implementation as it
            #   saves some standard summaries including n_eff and Rhat in the
            #   environment of the object. The summaries and created upon 
            #   the first time standard summary is called and resued later if possible. 
            #   In addition, the indexes complicate the implementation as internally
            #   we use column major indexes for vector/array parameters. But it might
            #   be better to use row major indexes for the output such as print.

            if (object@mode == 1L) {
              cat("Stan model '", object@model_name, "' is of mode 'test_grad';\n",
                  "sampling is not conducted.\n", sep = '')
              return(invisible(NULL)) 
            } else if (object@mode == 2L) {
              cat("Stan model '", object@model_name, "' does not contain samples.\n", sep = '') 
              return(invisible(NULL)) 
            } 

            if (isTRUE(all.equal(object@sim$n_save, object@sim$warmup2,
                check.attributes = FALSE, check.names = FALSE))) {
              cat("Stan model '", object@model_name, "' does not contain samples after warmup.\n", sep = '')
              return(invisible(NULL))
            }

            if (!exists("summary", envir = object@.MISC, inherits = FALSE) && use_cache) 
              assign("summary", summary_sim(object@sim), envir = object@.MISC)
           
            pars <- if (missing(pars)) object@sim$pars_oi else check_pars_second(object@sim, pars) 
            pars <- remove_empty_pars(pars, object@sim$dims_oi)
            if (missing(probs)) 
              probs <- c(0.025, 0.25, 0.50, 0.75, 0.975)  

            if (!use_cache) {
              # not using the cached (and not create cache, which takes time for too many pars)
              ss <-  summary_sim(object@sim, pars, probs) 
              s1 <- cbind(ss$msd[, 1, drop = FALSE], ss$sem, ss$msd[, 2, drop = FALSE], 
                          ss$quan, ss$ess, ss$rhat)
              colnames(s1) <- c("mean", "se_mean", "sd", colnames(ss$quan), 'n_eff', 'Rhat')
              s2 <- combine_msd_quan(ss$c_msd, ss$c_quan) 
              idx2 <- match(attr(ss, "row_major_idx"), attr(ss, "col_major_idx"))
              sf <- list(summary = s1[idx2, , drop = FALSE],
                         c_summary = s2[idx2, , , drop = FALSE])
              return(sf)
            } 
            m <- match(probs, default_summary_probs())
            if (any(is.na(m))) { # unordinary quantiles are requested 
              ss <-  summary_sim_quan(object@sim, pars, probs) 
              col_idx <- attr(ss, "col_major_idx") 
              ss$ess <- object@.MISC$summary$ess[col_idx, drop = FALSE] 
              ss$rhat <- object@.MISC$summary$rhat[col_idx, drop = FALSE] 
              ss$mean <- object@.MISC$summary$msd[col_idx, 1, drop = FALSE] 
              ss$sd <- object@.MISC$summary$msd[col_idx, 2, drop = FALSE] 
              ss$sem <- object@.MISC$summary$sem[col_idx]  
              s1 <- cbind(ss$mean, ss$sem, ss$sd, 
                          ss$quan, ss$ess, ss$rhat)
              colnames(s1) <- c("mean", "se_mean", "sd", colnames(ss$quan), 'n_eff', 'Rhat')
              s2 <- combine_msd_quan(object@.MISC$summary$c_msd[col_idx, , , drop = FALSE], ss$c_quan) 
              idx2 <- match(attr(ss, "row_major_idx"), col_idx)
              ss <- list(summary = s1[idx2, , drop = FALSE],
                         c_summary = s2[idx2, , , drop = FALSE])
              return(ss)
            }

            tidx <- pars_total_indexes(object@sim$pars_oi, 
                                       object@sim$dims_oi, 
                                       object@sim$fnames_oi, 
                                       pars) 
            tidx <- lapply(tidx, function(x) attr(x, "row_major_idx"))
            tidx <- unlist(tidx, use.names = FALSE)
            tidx_len <- length(tidx)

            ss <- object@.MISC$summary 
            qnames <- colnames(ss$quan)[m]
            
            if (!is.null(object@stan_args[[1]]$method) && 
                         object@stan_args[[1]]$method == "variational") {
              s1 <- cbind(ss$msd[tidx, 1, drop = FALSE],
                          ss$msd[tidx, 2, drop = FALSE], 
                          ss$quan[tidx, m, drop = FALSE])
              dim(s1) <- c(length(tidx), length(m) + 2L)
              colnames(s1) <- c("mean", "sd", qnames)
            }
            else {
              s1 <- cbind(ss$msd[tidx, 1, drop = FALSE], 
                          ss$sem[tidx, drop = FALSE], 
                          ss$msd[tidx, 2, drop = FALSE], 
                          ss$quan[tidx, m, drop = FALSE], 
                          ss$ess[tidx, drop = FALSE],
                          ss$rhat[tidx, drop = FALSE])
              dim(s1) <- c(length(tidx), length(m) + 5L) 
              colnames(s1) <- c("mean", "se_mean", "sd", qnames, 'n_eff', 'Rhat')
            }
            pars_names <- rownames(ss$msd)[tidx] 
            rownames(s1) <- pars_names 
            s2 <- combine_msd_quan(ss$c_msd[tidx, , , drop = FALSE], ss$c_quan[tidx, m, , drop = FALSE]) 
            # dim(s2) <- c(tidx_len, length(m) + 2, object@sim$chains)
            # dimnames(s2) <- list(parameter = pars_names, 
            #                      stats = c("mean", "sd", qnames), NULL) 
            ss <- list(summary = s1, c_summary = s2) 
            return(ss)
          })  

if (!isGeneric("traceplot")) {
  setGeneric(name = "traceplot",
             def = function(object, ...) { standardGeneric("traceplot") }) 
} 

if (!isGeneric("log_prob")) {
  setGeneric(name = "log_prob", 
             def = function(object, ...) { standardGeneric("log_prob") }) 
} 

if (!isGeneric("unconstrain_pars")) {
  setGeneric(name = "unconstrain_pars", 
             def = function(object, ...) { standardGeneric("unconstrain_pars") }) 
} 

setMethod("unconstrain_pars", signature = "stanfit", 
          function(object, pars) {
            # pars is a list as in specifying inits for a chain
            if (!is_sfinstance_valid(object)) 
              stop("the model object is not created or not valid")
            object@.MISC$stan_fit_instance$unconstrain_pars(pars)
          })

if (!isGeneric("constrain_pars")) {
  setGeneric(name = "constrain_pars", 
             def = function(object, ...) { standardGeneric("constrain_pars") }) 
} 

setMethod("constrain_pars", signature = "stanfit", 
          function(object, upars) {
            # upars is a vector on the unconstrained space (R^N*), 
            # where N* is the number of unconstrained parameters. 
            if (!is_sfinstance_valid(object)) 
              stop("the model object is not created or not valid")
            p <- object@.MISC$stan_fit_instance$constrain_pars(upars)
            idx_wo_lp <- which(object@model_pars != 'lp__')
            rstan_relist(p, create_skeleton(object@model_pars[idx_wo_lp], object@par_dims[idx_wo_lp]))
          })


setMethod("log_prob", signature = "stanfit", 
          function(object, upars, adjust_transform = TRUE, gradient = FALSE) {
            if (!is_sfinstance_valid(object)) 
              stop("the model object is not created or not valid")
            return(object@.MISC$stan_fit_instance$log_prob(upars, adjust_transform, gradient))
          }) 

if (!isGeneric("get_num_upars")) {
  setGeneric(name = "get_num_upars", 
             def = function(object, ...) { standardGeneric("get_num_upars") }) 
} 

setMethod("get_num_upars", signature = "stanfit", 
          function(object) {
            if (!is_sfinstance_valid(object)) 
              stop("the model object is not created or not valid")
            object@.MISC$stan_fit_instance$num_pars_unconstrained()
          })

if (!isGeneric("grad_log_prob")) {
  setGeneric(name = "grad_log_prob", 
             def = function(object, ...) { standardGeneric("grad_log_prob") }) 
} 

setMethod("grad_log_prob", signature = "stanfit", 
          function(object, upars, adjust_transform = TRUE) {
            if (!is_sfinstance_valid(object)) 
              stop("the model object is not created or not valid")
            object@.MISC$stan_fit_instance$grad_log_prob(upars, adjust_transform) 
          }) 

setMethod("traceplot", signature = "stanfit", 
          function(object, pars, inc_warmup = TRUE, ask = FALSE, 
                   nrow = 4, ncol = 2, window = NULL, include = TRUE, ...) { 
            # Args:
            #  nrow, defaults to 4
            #  ncol, defaults to 2 
            #  nrow and ncol are used to define mfrow for the whole plot area
            #  when there are many parameters. 
            #  window, for plotting only a window of the whole iterations
            #  default to NULL for all iterations
            #  include the elements of pars (FALSE -> exclude)

            if (object@mode == 1L) {
              cat("Stan model '", object@model_name, "' is of mode 'test_grad';\n",
                  "sampling is not conducted.\n", sep = '')
              return(invisible(NULL)) 
            } else if (object@mode == 2L) {
              cat("Stan model '", object@model_name, "' does not contain samples.\n", sep = '') 
              return(invisible(NULL)) 
            } 

            if(!include) pars <- setdiff(object@sim$pars_oi, pars)
            pars <- if (missing(pars)) object@sim$pars_oi else check_pars_second(object@sim, pars) 
            pars <- remove_empty_pars(pars, object@sim$dims_oi)
            tidx <- pars_total_indexes(object@sim$pars_oi, 
                                       object@sim$dims_oi, 
                                       object@sim$fnames_oi, 
                                       pars) 
            tidx <- lapply(tidx, function(x) attr(x, "row_major_idx"))
            tidx <- unlist(tidx, use.names = FALSE)
            mfrow_old <- par('mfrow')
            on.exit(par(mfrow = mfrow_old))
            par(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2, 2, 1) + 0.1, tck = -0.02)
            num_plots <- length(tidx) 
            if (num_plots %in% 2:nrow) par(mfrow = c(num_plots, 1)) 
            if (num_plots > nrow) par(mfrow = c(nrow, ncol)) 
    
            if (!is.null(window)) { 
              window <- sort(window)
              if (window[1] < 1) window[1] <- 1
              if (window[1] > object@sim$iter[1])
                stop("wrong specification of argument window", call. = FALSE)
              if (is.na(window[2]) || window[2] > object@sim$iter[1])
                window[2] <- object@sim$iter[1]
            } else { 
              window <- c(1, object@sim$iter[1]) 
            }
            if ((object@sim$warmup2 == 0 || !inc_warmup) && window[1] <= object@sim$warmup[1]) {
              window[1] <- object@sim$warmup[1] + 1
            }
            if (window[1] > window[2]) {
              stop("the given window does not include sample")
            } 
            if (window[1] > object@sim$warmup[1]) inc_warmup <- FALSE
            
            par_traceplot(object@sim, tidx[1], object@sim$fnames_oi[tidx[1]], 
                          inc_warmup = inc_warmup, window = window, ...)
            if (num_plots > nrow * ncol && ask) {
              ask_old <- devAskNewPage(ask = TRUE)
              on.exit(devAskNewPage(ask = ask_old), add = TRUE)
            }
            if (num_plots > 1) { 
              for (n in 2:num_plots)
                par_traceplot(object@sim, tidx[n], object@sim$fnames_oi[tidx[n]], 
                              inc_warmup = inc_warmup, window = window, ...)
            }
            invisible(NULL) 
          })  


is_sf_valid <- function(sf) {
  # Similar to is_sm_valid, this is only to test whether
  # the compiled DSO is loaded 
  return(is_sm_valid(sf@stanmodel)) 
} 

is_sfinstance_valid <- function(object) {
  # Args
  #  object: an instance of S4 class stanfit 
  exists("stan_fit_instance", envir = object@.MISC, inherits = FALSE) && 
  !(is_null_ptr(object@.MISC$stan_fit_instance@.xData$.pointer) ||
    is_null_ptr(object@.MISC$stan_fit_instance@.xData$.module) || 
    is_null_ptr(object@.MISC$stan_fit_instance@.xData$.cppclass)) 
} 


sflist2stanfit <- function(sflist) {
  # merge a list of stanfit objects into one
  # Args:
  #   sflist, a list of stanfit objects, each element of which 
  #   should have equal length of `iter`, `warmup`, and `thin`. 
  # Returns: 
  #   A new stanfit objects have all the chains in each element of sf_list. 
  #   The date would be where the new object is created. 
  # Note: 
  #   * method get_seed would not work well for this merged 
  #     stanfit object in that it only returns the seed used
  #     for the first object. But all the information is still there. 
  #   * When print function is called, the sampler info is obtained 
  #     only from the first chain. 
  #  
  sf_len <- length(sflist) 
  if (sf_len == 0) stop("'sflist' should have at least 1 element")
  if (!is.list(sflist) || 
      any(sapply(sflist, function(x) !is(x, "stanfit")))) {   
    stop("'sflist' must be a list of 'stanfit' objects")
  }

  non_zero_modes_idx <- which(sapply(sflist, function(x) x@mode) > 0)
  if (length(non_zero_modes_idx) > 0) { 
    stop("The following elements of 'sflist' do not contain samples: ",
         paste(non_zero_modes_idx, collapse = ', '), ".") 
  }   
  if (sf_len == 1) return(sflist[[1]])
  for (i in 2:sf_len) { 
    if (!identical(sflist[[i]]@sim$pars_oi, sflist[[1]]@sim$pars_oi) || 
        !identical(sflist[[i]]@sim$dims_oi, sflist[[1]]@sim$dims_oi))   
      stop("parameters in element ", i, " (stanfit object) are different from in element 1") 
    if (sflist[[i]]@sim$n_save[1] != sflist[[1]]@sim$n_save[1] ||  
        sflist[[i]]@sim$warmup2[1] != sflist[[1]]@sim$warmup2[1])
      stop("all 'stanfit' objects should have equal length of iterations and warmup") 
  } 
  n_chains = sum(sapply(sflist, function(x) x@sim$chains))  
  sim = list(samples = do.call(c, lapply(sflist, function(x) x@sim$samples)),  
             chains = n_chains, 
             iter = sflist[[1]]@sim$iter, 
             thin = sflist[[1]]@sim$thin, 
             warmup = sflist[[1]]@sim$warmup,
             n_save = rep(sflist[[1]]@sim$n_save[1], n_chains), 
             warmup2 = rep(sflist[[1]]@sim$warmup2[1], n_chains), 
             permutation = do.call(c, lapply(sflist, function(x) x@sim$permutation)), 
             pars_oi = sflist[[1]]@sim$pars_oi, 
             dims_oi = sflist[[1]]@sim$dims_oi, 
             fnames_oi = sflist[[1]]@sim$fnames_oi,
             n_flatnames = sflist[[1]]@sim$n_flatnames) 
  nfit <- new("stanfit", 
              model_name = sflist[[1]]@model_name, 
              model_pars = sflist[[1]]@model_pars,
              par_dims  = sflist[[1]]@par_dims, 
              mode = 0L, 
              sim = sim, 
              inits = do.call(c, lapply(sflist, function(x) x@inits)), 
              stan_args = do.call(c, lapply(sflist, function(x) x@stan_args)), 
              stanmodel = sflist[[1]]@stanmodel, 
              date = date(), 
              .MISC = new.env(parent = emptyenv())) 
  return(nfit)
} 
# sflist2stan(list(l1=ss1, l2=ss2))

names.stanfit <- function(x) dimnames(x)$parameters

`names<-.stanfit` <- function(x, value) {
  value <- as.character(value)
  len <- length(x@sim$fnames_oi)
  if(length(value) != len)
    stop(paste("parameter names must be of length", len))
  x@sim$fnames_oi <- value
  if(length(x@.MISC$summary)) {
    x@.MISC$summary <- rapply(x@.MISC$summary, f = function(y) {
      rownames(y) <- value
      return(y)
    }, how = "replace")
  }
  return(x)
}

as.array.stanfit <- function(x, ...) {
  if (x@mode != 0) return(numeric(0)) 
  out <- extract(x, permuted = FALSE, inc_warmup = FALSE, ...)
  # dimnames(out) <- dimnames(x)
  return(out)
} 
as.matrix.stanfit <- function(x, ...) {
  if (x@mode != 0) return(numeric(0)) 
  e <- extract(x, permuted = FALSE, inc_warmup = FALSE, ...) 
  out <- apply(e, 3, FUN = function(y) y)
  dimnames(out) <- dimnames(e)[-2]
  return(out)
}
 
as.data.frame.stanfit <- function(x, ...) {
  return( as.data.frame(as.matrix(x, ...)) )
}

dim.stanfit <- function(x) {
  if (x@mode != 0) return(numeric(0)) 
  c(x@sim$n_save[1] - x@sim$warmup2[1], x@sim$chains, x@sim$n_flatnames)
}

setGeneric("as.mcmc.list", function(object, ...) standardGeneric("as.mcmc.list"))

as.mcmc.list.stanfit <- function(object, pars, ...) {
  pars <- if (missing(pars)) object@sim$pars_oi else check_pars_second(object@sim, pars) 
  pars <- remove_empty_pars(pars, object@sim$dims_oi)
  tidx <- pars_total_indexes(object@sim$pars_oi, object@sim$dims_oi, object@sim$fnames_oi, pars)
  tidx <- lapply(tidx, function(x) attr(x, "row_major_idx"))
  tidx <- unlist(tidx, use.names = FALSE)

  lst <- vector("list", object@sim$chains)
  for (ic in 1:object@sim$chains) { 
    x <- do.call(cbind, object@sim$samples[[ic]])[,tidx,drop=FALSE]
    warmup2 <- object@sim$warmup2[ic] 
    if (warmup2 > 0) x <- x[-(1:warmup2), ]
    x <- as.matrix(x)
    end <- object@sim$iter
    thin <- object@sim$thin
    start <- end - (nrow(x) - 1) * thin
    class(x) <- 'mcmc'
    attr(x, "mcpar") <- c(start, end, thin)
    lst[[ic]] <- x 
  }
  class(lst) <- "mcmc.list"
  return(lst)
}

setMethod("as.mcmc.list", "stanfit", as.mcmc.list.stanfit)

As.mcmc.list <- function(object, pars, include = TRUE, ...) {
  if (missing(pars)) pars <- object@sim$pars_oi
  else if (!include) pars <- setdiff(object@sim$pars_oi, pars)
  pars <-  check_pars_second(object@sim, pars)
  return(as.mcmc.list.stanfit(object, pars = pars))
}

dimnames.stanfit <- function(x) {
  if (x@mode != 0) return(character(0)) 
  cids <- sapply(x@stan_args, function(x) x$chain_id)
  list(iterations = NULL, chains = paste0("chain:", cids), parameters = x@sim$fnames_oi) 
}
is.array.stanfit <- function(x)  return(x@mode == 0)

