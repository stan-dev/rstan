sbcFitFile <- function(save_progress, stanmodel, S) {
  file.path(save_progress, paste0(stanmodel@model_name, '-', S, '.rda'))
}

sbc <- function(stanmodel, data, M, ..., save_progress, load_incomplete=FALSE) {
  stopifnot(is(stanmodel, "stanmodel"))

  doSave <- !missing(save_progress)
  if (doSave && !dir.exists(save_progress)) {
    stop(paste0("Argument save_progress='", save_progress,
      "' provided but directory does not exist or is not a directory"))
  }

  # parameter names
  stan_code <- get_stancode(stanmodel)
  stan_code <- scan(what = character(), sep = "\n", quiet = TRUE, text = stan_code)
  pars_lines <- grep("[[:space:]]*(pars_)|(pars_\\[.*\\])[[:space:]]*=", stan_code, value = TRUE)
  if(length(pars_lines)==1){
    # if the parameter was defined and assigned in one line:
    pars_lines <- gsub(".*?=.*?(\\{|\\[)(.*?)(\\]|\\}).*","\\2",pars_lines)
    pars_names <- trimws(strsplit(pars_lines, split = ",")[[1]])
    pars_names <- unique(sub("^([a-z,A-Z,0-9,_]*)_.*", "\\1",
        pars_names))
     } else {
    pars_lines <- pars_lines[!grepl("^[[:space:]]*vector", pars_lines) &
        !grepl("^[[:space:]]*real", pars_lines)]
    pars_names <- trimws(sapply(strsplit(pars_lines, split = "=",
        fixed = TRUE), tail, n = 1))
    pars_names <- unique(sub("^([a-z,A-Z,0-9,_]*)_.*;", "\\1",
        pars_names))
     }
  noUnderscore <- grepl(";", pars_names, fixed=TRUE)
  if (any(noUnderscore)) {
    warning(paste("The following parameters were added to pars_ but did not",
      "have the expected underscore postfix:",
      paste(pars_names[noUnderscore], collapse = ' ')))
    return()
  }
  has_log_lik <- any(grepl("log_lik[[:space:]]*;[[:space:]]*", stan_code))

  if (!load_incomplete) {
    todo <- as.integer(seq(from = 0, to = .Machine$integer.max, length.out = M))
  } else {
    mn <- stanmodel@model_name
    runs <- dir(save_progress)
    runs <- runs[grepl(paste0("^", mn,'-(\\d+).rda$'), runs)]
    if (length(runs) == 0) {
      stop(paste("No completed runs found in", dir,
                 "matching regular expression", paste0("^", mn,'-(\\d+).rda$'),
                 "\nDid you use sbc(..., save_progress='/path/to/results')?"))
    }
    todo <- as.integer(sub(paste0(mn,'-(\\d+).rda'), "\\1", runs))
  }
  post <- parallel::mclapply(todo, FUN = function(S) {
    if (doSave) {
      file <- sbcFitFile(save_progress, stanmodel, S)
      if (file.exists(file)) return(TRUE)
    }
    out <- sampling(stanmodel, data, pars = c("ranks_", if (has_log_lik) "log_lik"), include = TRUE,
      chains = 1L, cores = 1L, seed = S, save_warmup = FALSE, thin = 1L, ...)
    out@stanmodel <- new("stanmodel")
    if (doSave) {
      save(out, file=file)
      return(TRUE)
    }
    out
  })
  if (doSave) {
    bad <- c()
    for (m in 1:length(todo)) {
      file <- sbcFitFile(save_progress, stanmodel, todo[m])
      got <- try(load(file), silent=TRUE)
      if (is(got, "try-error")) {
        bad <- c(bad, file)
        next
      }
      post[[m]] <- out
    }
    if (length(bad)) stop(paste("Remove corrupt files:",
      paste(bad, collapse=' '), "\nThen try again"))
  }
  bad <- sapply(post, FUN = function(x) x@mode != 0)
  if (any(bad)) {
    warning(sum(bad), " out of ", length(todo), " runs failed. Try decreasing 'init_r'")
    if(all(bad)) stop("cannot continue")
    post <- post[!bad]
  }

  noTwin <- is.na(match(pars_names, names(post[[1]]@par_dims)))
  if (any(noTwin)) {
    warning(paste("The following underscored priors did not have a matching",
      "parameter without the underscore:", paste(pars_names[noTwin], collapse = ' ')))
    pars_names <- NULL
  }
  pars_names <- try(flatnames(pars_names, post[[1]]@par_dims[pars_names]), silent = TRUE)
  if (!is.character(pars_names)) {
    warning("parameter names could not be calculated due to non-compliance with conventions; see help(sbc)")
    pars_names <- NULL
  }

  # divergences, etc.
  sampler_params <- simplify2array(sapply(post, FUN = get_sampler_params, inc_warmup = FALSE))

  # prior predictive distribution
  Y <- sapply(post, FUN = function(p) {
    means <- get_posterior_mean(p)[, 1]
    # will be just 'y_' if length 1, otherwise will have brackets
    means[grepl("y_$|y_\\[.*\\]", names(means))]
  })

  # realizations of true parameters from transformed data
  pars <- sapply(post, FUN = function(p) {
    means <- get_posterior_mean(p)[, 1]
    mark <- grepl("pars_\\[[[:digit:]]+\\]", names(means))
    return(means[mark])
  })
  if (length(pars) > 0L) {
    if (is.null(dim(pars))) pars <- matrix(pars, nrow=1)
    if (dim(pars)[1] != length(pars_names)) {
      warning("parameter names miscalculated due to non-compliance with conventions; see help(sbc)")
      pars_names <- NULL
    }
    if (is.null(dim(pars))) {
      pars <- t(as.matrix(pars))
    }
    rownames(pars) <- pars_names
  }

  # not actually ranks but rather unthinned binary values for draw > true
  ranks <- lapply(post, FUN = function(p) {
    r <- extract(p, pars = "ranks_", permuted = FALSE)[, 1, ]
    if (is.null(dim(r))) {
      r <- as.matrix(r)
    }
    colnames(r) <- pars_names
    r[] <- r > 0
    return(r)
  })

  if (has_log_lik) # high Pareto k values will be shown by the print method
    pareto_k <- sapply(post, FUN = function(x) suppressWarnings(loo(x))$diagnostics$pareto_k)

  out <- list(ranks = ranks, Y = Y, pars = pars, sampler_params = sampler_params,
              pareto_k = if (has_log_lik) pareto_k)
  class(out) <- "sbc"
  return(out)
}

  plot.sbc <- function(x, thin = 3, ...) {
    thinner <- seq(from = 1, to = nrow(x$ranks[[1]]), by = thin)
    u <- t(sapply(x$ranks, FUN = function(r) 1L + colSums(r[thinner, , drop = FALSE])))
    parameter <- as.factor(rep(colnames(u), each = nrow(u)))
    d <- data.frame(u = c(u), parameter)
    suppressWarnings(ggplot2::ggplot(d) +
      ggplot2::geom_freqpoly(ggplot2::aes(x = u), ...) +
      ggplot2::facet_wrap("parameter"))
  }

print.sbc <- function(x, ...) {
  divergences <- apply(x$sampler_params, MARGIN = 3, FUN = function(y) sum(y[,"divergent__"]))
  bad <- sum(divergences > 0L)
  cat(paste(bad, "chains had divergent transitions after warmup\n"))
  if (bad > 0L) cat(paste("there were a total of", sum(divergences),
                            "divergent transitions across all chains\n"))
  if (length(x$pareto_k)) {
    cut_pareto_k <- cut(c(x$pareto_k), breaks = c(-Inf, 0.5, 0.7, 1, Inf))
    cat("Aggregate Pareto k estimates:\n")
    print(prop.table(table(cut_pareto_k)))
  }
  return(invisible(NULL))
}
