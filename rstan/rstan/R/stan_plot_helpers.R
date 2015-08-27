.aesthetic <- function(which, ...) {
  dots <- list(...)
  which2 <- gsub("pt_", "", which)
  if (which2 %in% names(dots)) dots[[which2]]
  else .rstanvis_defaults[[which]]
}
.fill <- function(...) .aesthetic("fill", ...)
.color <- function(...) .aesthetic("color", ...)
.size <- function(...) .aesthetic("size", ...)
.alpha <- function(...) .aesthetic("alpha", ...)
.shape <- function(...) .aesthetic("shape", ...)
.pt_color <- function(...) .aesthetic("pt_color", ...)
.pt_size <- function(...) .aesthetic("pt_size", ...)

`%ifNULL%` <- function(x, replacement) {
  if (!is.null(x)) x
  else replacement
}

is.stanreg <- function(x) inherits(x, "stanreg")
is.stanfit <- function(x) inherits(x, "stanfit")

.get_stanreg_parnames <- function(x) {
  stopifnot(is.stanreg(x))
  pars <- rownames(x$stan_summary)
  sel <- grepl("mean_ppd|lp__|log-posterior", pars, ignore.case = TRUE)
  c(pars[!sel], pars[sel])
}

.check_object <- function(object, unconstrain = FALSE) {
  if (is.stanreg(object)) {
    if (object$algorithm != "sampling")
      stop("Plots not yet available for estimation methods other than MCMC (algorithm='sampling')", 
           call. = FALSE)
    if (unconstrain)
      stop("Option 'unconstrain' not yet available for stanreg objects.",
           call. = FALSE)
    .mode_check(object$stanfit)
  }
  else .mode_check(object)
}

.mode_check <- function(object) {
  stopifnot(is.stanfit(object))
  if (object@mode == 1L) {
    stop("Stan model '", object@model_name, "' is of mode 'test_grad';\n",
         "sampling is not conducted.\n",
         call. = FALSE)
  } else if (object@mode == 2L) {
    stop("Stan model '", object@model_name, "' does not contain samples.\n",
         call. = FALSE)
  }
  invisible(TRUE)
}


.reshape_sample <- function(x) {
  res <- lapply(seq_along(x), function(i) {
    data.frame(value = unlist(x[[i]], use.names = FALSE),
               parameter = rep(names(x[[i]]), each = length(x[[i]][[1L]])),
               chain = i)
  })
  res <- do.call(rbind, res)
  res$parameter <- as.character(res$parameter)
  res
}
.make_plot_data <- function(object, pars, include = TRUE,
                            inc_warmup = FALSE, unconstrain = FALSE) {
  
  window <- NULL
  if (is.stanreg(object)) {
    sim <- object$stanfit@sim
  }
  else sim <- object@sim
  if (!include) {
    if (missing(pars))
      stop("pars must be specified if include=FALSE.", call. = FALSE)
    pars <- setdiff(sim$pars_oi, pars)
  }
  pars <- if (missing(pars)) sim$pars_oi else check_pars_second(sim, pars)
  pars <- remove_empty_pars(pars, sim$dims_oi)
  if (unconstrain && "lp__" %in% pars) {
    if (length(pars) == 1L) stop("Can't unconstrain lp__", call. = FALSE)
    pars <- pars[-which(pars == "lp__")]
  }
  tidx <- pars_total_indexes(sim$pars_oi,
                             sim$dims_oi,
                             sim$fnames_oi,
                             pars)
  tidx <- lapply(tidx, function(x) attr(x, "row_major_idx"))
  tidx <- unlist(tidx, use.names = FALSE)
  num_plots <- length(tidx)
  
  if (!is.null(window)) {
    window <- sort(window)
    if (window[1] < 1) window[1] <- 1
    if (window[1] > sim$iter[1])
      stop("wrong specification of argument window", call. = FALSE)
    if (is.na(window[2]) || window[2] > sim$iter[1])
      window[2] <- sim$iter[1]
  } else {
    window <- c(1, sim$iter[1])
  }
  if ((sim$warmup2 == 0 || !inc_warmup) && window[1] <= sim$warmup[1]) {
    window[1] <- sim$warmup[1] + 1
  }
  if (window[1] > window[2]) {
    stop("the given window does not include sample")
  }
  if (window[1] > sim$warmup[1]) inc_warmup <- FALSE
  
  thin <- sim$thin
  warmup2 <- sim$warmup2[1]
  warmup <- sim$warmup
  start_i <- window[1]
  window_size <- (window[2] - start_i) %/% thin
  id <- seq.int(start_i, by = thin, length.out = window_size)
  start_idx <- (if (warmup2 == 0) (start_i - warmup) else start_i) %/% thin
  if (start_idx < 1)  start_idx <- 1
  idx <- seq.int(start_idx, by = 1, length.out = window_size)
  if (unconstrain) sim$samples <- .upars(object)
  samp_use <- lapply(sim$samples, function(chain) {
    out <- lapply(chain[tidx], function(x) x[idx])
    names(out) <- sim$fnames_oi[tidx]
    out
  })
  nchains <- length(samp_use)
  dat <- .reshape_sample(samp_use)
  dat$iteration <- idx
  dat$chain <- factor(dat$chain)
  lp <- which(dat$parameter == "lp__")
  if (!identical(lp, integer(0)))
    dat$parameter[lp] <- "log-posterior"
  if (is.stanreg(object))
    dat$parameter <- factor(dat$parameter, levels = sim$fnames_oi[tidx])
  list(samp = dat,
       nchains = nchains,
       nparams = length(tidx),
       warmup = warmup)
}

.upars <- function(object, permuted = FALSE, inc_warmup = TRUE) {
  if (is.stanreg(object)) object <- object$stanfit
  pars <- extract(object, permuted = permuted, inc_warmup = inc_warmup)
  nchains <- ncol(pars)
  pn <- dimnames(pars)$parameters
  param_names <- if ("lp__" %in% pn) pn[-which(pn == "lp__")] else pn
  skeleton <- get_inits(object)[[1L]]
  upars <- apply(pars, 1:2, FUN = function(theta) {
    unconstrain_pars(object, relist(theta, skeleton))
  })
  upars <- aperm(upars, c(3, 1, 2))
  out <- lapply(seq_len(nchains), function(chain) {
    x <- upars[chain,, ]
    plist <- lapply(seq_len(nrow(x)), function(param) x[param, ])
    names(plist) <- param_names
    plist
  })
  out
}

color_vector <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=50, c=50)[1:n]
}
color_vector_chain <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=80, c=50)[1:n]
}



# autocorrelation ---------------------------------------------------------
.ac_fun <- function(x, lag.max, partial = FALSE) {
  if (!partial)
    acf(x, lag.max = lag.max, plot = FALSE)$acf[,, 1L]
  else
    pacf(x, lag.max = lag.max, plot = FALSE)$acf[,, 1L]
}
.ac_plot_data <- function(dat, lags, partial = FALSE) {
  ch <- dat[, grep("chain", colnames(dat))]
  nc <- length(unique(ch))
  ac_list <- tapply(dat$value, INDEX = ch, FUN = .ac_fun, lag.max = lags,
                    partial = partial, simplify = FALSE)
  nl <- if (partial) lags else lags + 1
  ch <- factor(rep(1:nc, each = nl), labels = paste0("chain:", 1:nc))
  ll <- rep(seq(if (partial) 1 else 0, lags), nc)
  data.frame(chains = ch, ac = do.call(c, ac_list), lag = ll)
}
.ac_plot_data_multi <- function(dat, lags, partial = FALSE) {
  ch <- dat[, grep("chain", colnames(dat))]
  nc <- length(unique(ch))
  pa <- factor(dat[, grep("parameter", colnames(dat))])
  np <- length(unique(pa))
  ac_list <- tapply(dat$value, INDEX = list(ch, pa),
                    FUN = .ac_fun, lag.max = lags,
                    partial = partial, simplify = FALSE)
  nl <- if (partial) lags else lags + 1
  ch <- factor(rep(rep(1:nc, each = nl), np), labels = paste0("chain:", 1:nc))
  ll <- rep(seq(if (partial) 1 else 0, lags), nc * np)
  pp <- factor(rep(1:np, each = nc * nl), labels = levels(pa))
  data.frame(parameters = pp, chains = ch, ac = do.call(c, ac_list), lag = ll)
}



# rhat, ess, mcse ---------------------------------------------------------
.rhat_neff_mcse_hist <- function(which = c("rhat", "n_eff_ratio", "mcse_ratio"),
                                 object, pars=NULL, ...) {
  if (is.stanreg(object)) object <- object$stanfit
  thm <- .rstanvis_defaults$hist_theme

  if (which == "n_eff_ratio") {
    lp <- suppressWarnings(get_logposterior(object, inc_warmup = FALSE))
    SS <- prod(length(lp), length(lp[[1L]]))
  }
  
  if (!is.null(pars)) smry <- summary(object, pars = pars)$summary
  else smry <- summary(object)$summary
  
  xlab <- switch(which,
                 rhat = "Rhat statistic",
                 n_eff_ratio = "Effective sample size / Sample size",
                 mcse_ratio = "Monte Carlo SE / Posterior SD"
  )
  stat <- switch(which,
                 rhat = smry[, "Rhat"],
                 n_eff_ratio = smry[, "n_eff"] / SS,
                 mcse_ratio = smry[, "se_mean"] / smry[, "sd"]
  )
  df <- data.frame(stat)
  plot_labs <- labs(y = "", x = xlab)
  base <- ggplot(df, aes(x = stat))
  base +
    geom_histogram(fill = .fill(...), color = .color(...), ...) +
    plot_labs +
    thm
}




# nuts stuff --------------------------------------------------------------
.LP_NAME <- "log-posterior"
.LP_LAB <- "Log Posterior"
.METROP_LAB <- "Mean Metrop. Acceptance"
.STEPSIZE_LAB <- "Sampled Step Size"
.TREEDEPTH_LAB <- "Treedepth"
.NDIVERGENT_LAB <- "N Divergent"

.NUTS_VLINE_CLR <- "#222222"
.NUTS_FILL <- "#66a7e0"
.NUTS_CLR <- "#006dcc"
.NUTS_PT_CLR <- "#328ad6"
.NDIVERGENT_FILL <- "#ae0001"
.MAXTD_FILL <- "#eeba30"
.NDIVERGENT_CLR <-  "black"
.MAXTD_CLR <- "black"
.DIV_AND_MAXTD_SHAPE <- 21

.nuts_return <- function(graphs, ...) {
  gtable <- do.call(gridExtra::arrangeGrob, c(graphs, ...))
  plot(gtable)
  invisible(graphs)
}

.nuts_args_check <- function(...) {
  if ("par" %in% names(dots <- list(...)))
    stop("'par' argument should not be specified.", call. = FALSE)
}

.max_td <- function(x) {
  if (is.stanreg(x)) 
    x <- x$stanfit
  cntrl <- x@stan_args[[1L]]$control
  if (is.null(cntrl)) 11
  else {
    max_td <- cntrl$max_treedepth
    if (is.null(max_td)) 11
    else max_td  
  }
}

.sampler_params_post_warmup <- function(object, which = "stepsize__", as.df = FALSE) {
  if (is.stanreg(object))
    object <- object$stanfit
  sampler_params <- get_sampler_params(object)
  warmup_val <- floor(object@sim$warmup / object@sim$thin)
  sp <-lapply(1:length(sampler_params), function(i) {
    out <- sampler_params[[i]]
    out <- if (warmup_val == 0) out[, which] else out[-(1:warmup_val), which]
    names(out) <- (warmup_val + 1):(warmup_val + NROW(out))
    out
  })
  if (length(which) == 1L && as.df) .sampler_param_df(sp, warmup_val)
  else sp
}

.sampler_param_df <- function(sp, warmup_val) {
  sp_mat <- do.call("cbind", sp)
  colnames(sp_mat) <- paste0("chain:", 1:ncol(sp_mat))
  sp_mat <- cbind(iterations = (warmup_val+1):(warmup_val + nrow(sp_mat)), sp_mat)
  as.data.frame(sp_mat)
}

.reshape_df <- function(df) {
  iter_col <- grep("iteration", colnames(df))
  mdf <- reshape(df, direction = "long",
                 v.names = "value",
                 varying = list(colnames(df)[-iter_col]),
                 timevar = "variable",
                 times = colnames(df)[-iter_col],
                 idvar = colnames(df)[iter_col])
  rownames(mdf) <- 1:nrow(mdf)
  mdf$variable <- as.factor(mdf$variable)
  attr(mdf, "reshapeLong") <- NULL
  mdf
}

.sampler_param_vs_param <- function(p, sp, divergent = NULL, hit_max_td = NULL,
                                    p_lab, sp_lab, chain = 0, violin = FALSE,
                                    smoother = FALSE, ...) {
  xy_labs <- labs(y = if (missing(p_lab)) NULL else p_lab,
                  x = if (missing(sp_lab)) NULL else sp_lab)
  df <- data.frame(sp = do.call("c", sp), p = c(p))
  if (violin) df$sp <- as.factor(round(df$sp, 4))
  if (!is.null(divergent)) df$divergent <- do.call("c", divergent)
  if (!is.null(hit_max_td)) df$hit_max_td <- do.call("c", hit_max_td)
  
  base <- ggplot(df, aes_string(x = "sp", y = "p")) + xy_labs
  if (chain == 0) {
    if (violin)
      graph <- base + geom_violin(color = .NUTS_CLR, fill = .NUTS_FILL)
    else {
      graph <- base + geom_point(alpha = 0.5 * .alpha(...), color = .NUTS_PT_CLR,
                                 fill = .NUTS_FILL, shape = .shape(...), ...)
      if (smoother) graph <- graph + stat_smooth(se = FALSE)
      if (!is.null(divergent))
        graph <- graph + geom_point(data = subset(df, divergent == 1),
                                    aes_string(x = "sp", y = "p"),
                                    color = .NDIVERGENT_CLR, fill = .NDIVERGENT_FILL,
                                    alpha = 0.8, size = 3, shape = .DIV_AND_MAXTD_SHAPE)
      if (!is.null(hit_max_td))
        graph <- graph + geom_point(data = subset(df, hit_max_td == 1),
                                    aes_string(x = "sp", y = "p"),
                                    color = .MAXTD_CLR, fill = .MAXTD_FILL,
                                    alpha = 0.8, size = 3, shape = .DIV_AND_MAXTD_SHAPE)
    }
    return(graph)
  }
  chain_data <- data.frame(sp = sp[, chain], p = p[, chain])
  if (!is.null(divergent)) chain_data$div <- divergent[, chain]
  if (!is.null(hit_max_td)) chain_data$hit <- hit_max_td[, chain]
  chain_clr <- color_vector_chain(ncol(sp))[chain]
  chain_fill <- chain_clr
  if (violin) {
    chain_data$sp <- as.factor(round(chain_data$sp, 4))
    graph <- base +
      geom_violin(color = .NUTS_CLR, fill = .NUTS_FILL, ...) +
      geom_violin(data = chain_data, aes_string(x = "sp", y = "p"),
                  color = chain_clr, fill = chain_fill,
                  alpha = .alpha(...), ...)
    return(graph)
  }
  graph <- base + geom_point(alpha = 0.5 * .alpha(...), color = .NUTS_PT_CLR,
                             fill = .NUTS_FILL, shape = .shape(...), ...)
  if (smoother) graph <- graph + stat_smooth(se = FALSE)
  graph <- graph +
    geom_point(data = chain_data, aes_string(x = "sp", y = "p"),
               color = chain_fill, alpha = .alpha(...))
  if (smoother) graph <- graph +
    stat_smooth(data = chain_data, aes_string(x = "sp", y = "p"),
                color = chain_fill, se = FALSE)
  if (!is.null(divergent))
    graph <- graph + geom_point(data = chain_data[chain_data$div == 1,, drop=FALSE],
                                aes_string(x = "sp", y = "p"),
                                color = .NDIVERGENT_CLR, fill = .NDIVERGENT_FILL,
                                size = 3, shape = .DIV_AND_MAXTD_SHAPE)
  if (!is.null(hit_max_td))
    graph <- graph + geom_point(data = chain_data[chain_data$hit == 1,, drop=FALSE],
                                aes_string(x = "sp", y = "p"),
                                color = .MAXTD_CLR, fill = .MAXTD_FILL,
                                size = 3, shape = .DIV_AND_MAXTD_SHAPE)
  graph
}

.sampler_param_vs_sampler_param_violin <- function(df_x, df_y, lab_x, lab_y,
                                                   chain = 0) {
  
  xy_labs <- labs(y = lab_y, x = lab_x)
  df <- data.frame(x = do.call("c", df_x), y = do.call("c", df_y))
  df$x <- as.factor(df$x)
  
  base <- ggplot(df, aes_string("x","y")) + xy_labs
  graph <- base + geom_violin(color = .NUTS_CLR, fill = .NUTS_FILL)
  if (chain == 0) return(graph)
  chain_clr <- color_vector_chain(ncol(df_x))[chain]
  chain_fill <- chain_clr
  chain_data <- data.frame(x = as.factor(df_x[, chain]), y = df_y[, chain])
  graph + geom_violin(data = chain_data, aes_string("x","y"), color = chain_clr,
                      fill = chain_fill, alpha = 0.5)
}

.p_hist <- function(df, lab, chain = 0, ...) {
  mdf <- .reshape_df(df) # reshape2::melt(df, id.vars = grep("iteration", colnames(df), value = TRUE))
  base <- ggplot(mdf, aes_string(x = "value")) +
    geom_histogram(binwidth = diff(range(mdf$value))/30, fill = .NUTS_FILL,
                   color = .NUTS_CLR, size = .size(...), ...) +
    labs(x = if (missing(lab)) NULL else lab, y = "")
  if (chain == 0) {
    graph <- base +
      geom_vline(xintercept = mean(mdf$value), color = .NUTS_VLINE_CLR, size = .8) +
      geom_vline(xintercept = median(mdf$value), color = .NUTS_VLINE_CLR, lty = 2, size = 1)
    return(graph)
  }
  chain_data <- mdf[mdf$variable == paste0("chain:",chain), ]
  chain_clr <- color_vector_chain(ncol(df) - 1)[chain]
  chain_fill <- chain_clr
  base + geom_histogram(data = chain_data,
                        binwidth = diff(range(chain_data$value))/30,
                        fill = chain_fill, alpha = 0.5) +
    geom_vline(xintercept = mean(chain_data$value), color = chain_clr, size = .8) +
    geom_vline(xintercept = median(chain_data$value),
               color = chain_clr, lty = 2, size = 1)
}

.treedepth_ndivergent_hist <- function(df_td, df_nd, chain = 0,
                                       divergent = c("All", 0, 1), ...) {
  x_lab <- if (divergent == "All")
    "Treedepth" else paste0("Treedepth (N Divergent = ", divergent,")")
  plot_labs <- labs(x = x_lab, y = "")
  
  mdf_td <- .reshape_df(df_td) #reshape2::melt(df_td, id.vars = grep("iteration", colnames(df_td), value = TRUE))
  mdf_nd <- .reshape_df(df_nd) #reshape2::melt(df_nd, id.vars = grep("iteration", colnames(df_nd), value = TRUE))
  mdf <- cbind(mdf_td, div = mdf_nd$value)
  plot_data <- if (divergent == "All") mdf else mdf[mdf$div == divergent,,drop=FALSE]
  if (nrow(plot_data) == 0) return(NULL)
  
  graph <- ggplot(plot_data, aes_q(x = quote(factor(value))), na.rm = TRUE) +
    stat_bin(aes_q(y=quote(..count../sum(..count..))),
             width=1, fill = .NUTS_FILL, color = .NUTS_CLR, size = .size(...)) +
    plot_labs
  if (chain == 0) return(graph)
  chain_clr <- color_vector_chain(ncol(df_td) - 1)[chain]
  chain_fill <- chain_clr
  chain_data <- plot_data[plot_data$variable == paste0("chain:",chain),, drop=FALSE]
  graph + stat_bin(data = chain_data, aes_q(y=quote(..count../sum(..count..))),
                   fill = chain_fill, alpha = .alpha(...), width = 1)
}