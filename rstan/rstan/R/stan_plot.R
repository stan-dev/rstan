# suppress messages from ggplot2
quietgg <- function(gg) {
  suppressMessages(suppressWarnings(print(gg)))
  invisible(gg)
}

# traceplot ---------------------------------------------------------------
stan_trace <- function(object, pars, include = TRUE,
                       unconstrain = FALSE,
                       inc_warmup = FALSE,
                       nrow = NULL, ncol = NULL,
                       ..., window = NULL) {
  
  .check_object(object, unconstrain)
  plot_data <- .make_plot_data(object, pars, include, inc_warmup, unconstrain)
  
  thm <- .rstanvis_defaults$theme
  clrs <- rep_len(.rstanvis_defaults$chain_colors, plot_data$nchains)
  
  base <- ggplot(plot_data$samp,
                 aes_string(x = "iteration", y = "value", color = "chain"))
  if (inc_warmup) base <- base +
    annotate("rect", xmin = -Inf, xmax = plot_data$warmup,
             ymin = -Inf, ymax = Inf, fill = .rstanvis_defaults$grays[2L])
  
  graph <-
    base +
    geom_path(...) +
    scale_color_manual(values = clrs) +
    labs(x="",y="") +
    thm
  
  if (plot_data$nparams == 1)
    graph <- graph + ylab(unique(plot_data$samp$parameter))
  else
    graph <- graph + facet_wrap(~parameter, nrow = nrow, ncol = ncol, scales = "free")
  
  if (!is.null(window)) {
    if (!is.numeric(window) || length(window) != 2)
      stop("'window' should be a numeric vector of length 2.")
    graph <- graph + coord_cartesian(xlim = window)
  }
  graph
}


# scatterplot -------------------------------------------------------------
stan_scat <- function(object, pars, unconstrain = FALSE, inc_warmup = FALSE,
                      nrow = NULL, ncol = NULL, ...) {
  
  .check_object(object, unconstrain)
  thm <- .rstanvis_defaults$theme
  dots <- .add_aesthetics(list(...), c("fill", "pt_color", "pt_size", "alpha", "shape"))
  if (missing(pars) || length(pars) != 2L)
    stop("'pars' must contain exactly two parameter names", call. = FALSE)
#   ndivergent <- 
#     .sampler_params_post_warmup(object, "divergent__", as.df = TRUE)[, -1L]
#   treedepth <- 
#     .sampler_params_post_warmup(object, "treedepth__", as.df = TRUE)[, -1L]
#   max_td <- .max_td(object)
#   div <- unname(rowSums(ndivergent) == 1)
#   hit_max_td <- sapply(1:nrow(treedepth), function(i) any(treedepth[i,] == max_td))
  plot_data <- .make_plot_data(
    object, 
    pars = pars, 
    include = TRUE, 
    inc_warmup = inc_warmup, 
    unconstrain = unconstrain
  )
  p1 <- plot_data$samp$parameter == pars[1]
  val1 <- plot_data$samp[p1, "value"]
  val2 <- plot_data$samp[!p1, "value"]
  df <- data.frame(x = val1, y = val2)
#   nchains <- plot_data$nchains
#   sel <- seq_len(nrow(df) / nchains)
#   div <- df[div[sel], ]
#   td <- df[hit_max_td[sel], ]
  base <- ggplot(df, aes_string("x", "y"))
  graph <-
    base +
    do.call(geom_point, dots) +
#     geom_point(data = div, aes_string("x","y"), color = "red") +
#     geom_point(data = td, aes_string("x","y"), color = "yellow") +
    labs(x = pars[1], y = pars[2]) +
    thm
  
  graph
}


# histograms ---------------------------------------------------------------
stan_hist <- function(object, pars, include = TRUE,
                      unconstrain = FALSE,
                      inc_warmup = FALSE,
                      nrow = NULL, ncol = NULL,
                      ...) {

  .check_object(object, unconstrain)
  dots <- .add_aesthetics(list(...), c("fill", "color"))
  plot_data <- .make_plot_data(object, pars, include, inc_warmup, unconstrain)
  thm <- .rstanvis_defaults$hist_theme
  base <- ggplot(plot_data$samp, aes_string(x = "value", y = "..density.."))
  graph <-
    base +
    do.call(geom_histogram, dots) + 
    xlab("") +
    thm
  
  if (plot_data$nparams == 1)
    graph + xlab(unique(plot_data$samp$parameter))
  else
    graph + facet_wrap(~parameter, nrow = nrow, ncol = ncol, scales = "free")
}



# densities -----------------------------------------------------------
stan_dens <- function(object, pars, include = TRUE,
                      unconstrain = FALSE,
                      inc_warmup = FALSE,
                      nrow = NULL, ncol = NULL,
                      ...,
                      separate_chains = FALSE) {
  
  .check_object(object, unconstrain)
  plot_data <- .make_plot_data(object, pars, include, inc_warmup, unconstrain)
  clrs <- rep_len(.rstanvis_defaults$chain_colors, plot_data$nchains)
  thm <- .rstanvis_defaults$hist_theme
  base <- ggplot(plot_data$samp, aes_string(x = "value")) + xlab("")
  if (!separate_chains) {
    dots <- .add_aesthetics(list(...), c("fill", "color"))
    graph <-
      base +
      do.call(geom_density, dots) + 
      thm
  } else {
    dots <- .add_aesthetics(list(...), c("color", "alpha"))
    dots$mapping <- aes_string(fill = "chain")
    graph <-
      base +
      do.call(geom_density, dots) + 
      scale_fill_manual(values = clrs) +
      thm
  }
  if (plot_data$nparams == 1)
    graph + xlab(unique(plot_data$samp$parameter))
  else
    graph + facet_wrap(~parameter, nrow = nrow, ncol = ncol, scales = "free")
}



# autocorrelation ---------------------------------------------------------
stan_ac <- function(object, pars, include = TRUE,
                    unconstrain = FALSE,
                    inc_warmup = FALSE,
                    nrow = NULL, ncol = NULL,
                    ...,
                    separate_chains = FALSE,
                    lags = 25, partial = FALSE) {
  .check_object(object, unconstrain)
  plot_data <- .make_plot_data(object, pars, include, inc_warmup, unconstrain)
  clrs <- rep_len(.rstanvis_defaults$chain_colors, plot_data$nchains)
  thm <- .rstanvis_defaults$theme
  dots <- .add_aesthetics(list(...), c("size", "color", "fill"))

  dat_args <- list(dat = plot_data$samp, lags = lags, partial = partial)
  dat_fn <- ifelse(plot_data$nparams == 1, ".ac_plot_data",".ac_plot_data_multi")
  ac_dat <- do.call(dat_fn, dat_args)
  if (!separate_chains) {
    dots$position <- "dodge"
    dots$stat <- "summary"
    dots$fun.y <- "mean"
    y_lab <- paste("Avg.", if (partial) "partial", "autocorrelation")
    ac_labs <- labs(x = "Lag", y = y_lab)
    y_scale <- scale_y_continuous(breaks = seq(0, 1, 0.25))
    base <- ggplot(ac_dat, aes_string(x = "lag", y = "ac"))
    graph <- base +
      do.call(geom_bar, dots) + 
      y_scale +
      ac_labs +
      thm
    
    if (plot_data$nparams == 1) {
      y_lab <- ylab(paste0(y_lab, " (", pars,")"))
      return(graph + y_lab)
    }
    else return(graph + facet_wrap(~parameters, nrow = nrow, ncol = ncol,
                                   scales = "free_x"))
  }
  dots$position <- "identity"
  dots$stat <- "identity"
  ac_labs <- labs(x = "Lag", y = if (partial)
    "Partial autocorrelation" else "Autocorrelation")
  y_scale <- scale_y_continuous(breaks = seq(0, 1, 0.25),
                                labels = c("0","","0.5","",""))
  
  base <- ggplot(ac_dat, aes_string(x = "lag", y = "ac"))
  graph <- base +
    do.call(geom_bar, dots) +
    y_scale +
    ac_labs +
    thm
  
  if (plot_data$nparams == 1) {
    graph <- graph + facet_wrap(~chains, nrow = nrow, ncol = ncol)
    return(graph)
  } else { # nParams > 1
    graph <- graph + facet_grid(parameters ~ chains, scales = "free_x")
    return(graph)
  }
}



# parameter estimates -----------------------------------------------------
stan_plot <- function(object, pars, include = TRUE, unconstrain = FALSE,
                      ...) {
  
  inc_warmup <- FALSE
  .check_object(object, unconstrain)
  thm <- .rstanvis_defaults$multiparam_theme
  plot_data <- .make_plot_data(object, pars, include, inc_warmup, unconstrain)
  
  color_by_rhat <- FALSE # FIXME 
  dots <- list(...)
  defs <- list(point_est = "median", show_density = FALSE,
               show_outer_line = TRUE, ci_level = 0.8, outer_level = 0.95,
               fill_color = .rstanvis_defaults[["fill"]], 
               outline_color = .rstanvis_defaults[["color"]], 
               est_color = .rstanvis_defaults[["color"]])
  args <- names(defs)
  dotenv <- list()
  for (j in seq_along(args)) {
    if (args[j] %in% names(dots))
      dotenv[[args[j]]] <- dots[[args[j]]]
    else dotenv[[args[j]]] <- defs[[j]]
  }
  if (!(dotenv[["point_est"]] %in% c("mean", "median")))
    stop("Point estimate should be either 'mean' or 'median'", call. = FALSE)
  if (color_by_rhat)
    stop("'color_by_rhat' not yet available", call. = FALSE)
  if (dotenv[["ci_level"]] > dotenv[["outer_level"]])
    stop("'ci_level' should be less than 'outer_level'", call. = FALSE)
  ci_level <- dotenv[["ci_level"]]
  outer_level <- dotenv[["outer_level"]]
  message("ci_level: ", ci_level," (",100 * ci_level, "% intervals)")
  message("outer_level: ", outer_level," (",100 * outer_level, "% intervals)")
  outer_level <- dotenv[["outer_level"]]
  probs.use <- c(0.5 - outer_level / 2, 0.5 - ci_level / 2, 0.5,
                 0.5 + ci_level / 2, 0.5 + outer_level / 2)
  samp <- plot_data$samp
  nparams <- plot_data$nparams
  statmat <- as.matrix(aggregate(samp$value, by = list(parameter = samp$parameter),
                                 FUN = function(x,...) c(mean(x), quantile(x,...)),
                                 probs = probs.use))
  param_names <- rownames(statmat) <- statmat[, 1L]
  statmat <- apply(statmat[, -1L, drop=FALSE], 1:2, as.numeric)
  colnames(statmat) <- c("mean", "2.5%", "25%", "50%", "75%", "97.5%")
  y <- as.numeric(seq(plot_data$nparams, 1, by = -1))
  xlim.use <- c(min(statmat[,2L]), max(statmat[,6L]))
  xlim.use <- xlim.use + diff(xlim.use) * c(-.05, .05)
  xy.df <- data.frame(params = rownames(statmat), y, statmat)
  colnames(xy.df) <- c("params", "y", "mean", "ll", "l", "m", "h", "hh")
  if (dotenv[["point_est"]] == "mean") xy.df$m <- xy.df$mean
  
  p.base <- ggplot(xy.df)
  p.name <- scale_y_continuous(breaks = y, labels = param_names,
                               limits = c(0.5, nparams + 1))
  p.all <- p.base + xlim(xlim.use) + p.name + thm
  show_density <- dotenv[["show_density"]]
  outline_color <- dotenv[["outline_color"]] %ifNULL% .rstanvis_defaults[["color"]]
  fill_color <- dotenv[["fill_color"]]
  est_color <- dotenv[["est_color"]]
  if (dotenv[["show_outer_line"]] || show_density) {
    p.ci <- geom_segment(aes_string(x = "ll", xend = "hh", y = "y", yend = "y"),
                         color = outline_color)
    p.all <- p.all + p.ci
  }
  if (show_density) { # plot kernel density estimate
    npoint.den <- 512
    y.den <- x.den <- matrix(0, nrow = npoint.den, ncol = nparams)
    for(i in 1:nparams){
      d.temp <- density(samp[samp$parameter == param_names[i], "value"],
                        from = statmat[i,2L],
                        to = statmat[i,6L],
                        n = npoint.den)
      x.den[,i] <- d.temp$x
      y.max <- max(d.temp$y)
      y.den[,i] <- d.temp$y / y.max * 0.8 + y[i]
    }
    df.den <- data.frame(x = as.vector(x.den), y = as.vector(y.den),
                         name = rep(param_names, each = npoint.den))
    p.den <- geom_line(data = df.den, aes_string("x", "y", group = "name"),
                       color = outline_color)
    
    #shaded interval
    y.poly <- x.poly <- matrix(0, nrow = npoint.den + 2, ncol = nparams)
    for(i in 1:nparams){
      d.temp <- density(samp[samp$parameter == param_names[i], "value"],
                        from = statmat[i, 3L],
                        to = statmat[i, 5L],
                        n = npoint.den)
      x.poly[,i] <- c(d.temp$x[1L], as.vector(d.temp$x), d.temp$x[npoint.den])
      y.max <- max(d.temp$y)
      y.poly[,i] <- as.vector(c(0, as.vector(d.temp$y) / y.max * 0.8, 0) + y[i])
    }
    df.poly <- data.frame(x = as.vector(x.poly), y = as.vector(y.poly),
                          name = rep(param_names, each = npoint.den + 2))
    p.poly <- geom_polygon(data = df.poly, aes_string("x", "y", group = "name", fill = "y"))
    p.col <- scale_fill_gradient(low = fill_color, high = fill_color, guide = "none")
    
    #point estimate
    if (color_by_rhat) {
      rhat_colors <- dotenv[["rhat_colors"]]
      p.point <- geom_segment(aes_string(x = "m", xend = "m", y = "y", yend = "y + 0.25",
                                         color = "rhat_id"), size = 1.5)
      p.all + p.poly + p.den + p.col + p.point + rhat_colors #+ rhat_lgnd
    } else {
      p.point <- geom_segment(aes_string(x = "m", xend = "m", y = "y", yend = "y + 0.25"),
                              colour = est_color, size = 1.5)
      p.all + p.poly + p.den + p.col + p.point
    }
  } else {
    p.ci.2 <- geom_segment(aes_string(x = "l", xend = "h", y = "y", yend = "y"),
                           colour = fill_color, size = 2)
    if (color_by_rhat) {
      p.point <- geom_point(aes_string(x = "m", y = "y", fill = "rhat_id"),
                            color = "black", shape = 21, size = 4)
      p.all + p.ci.2 + p.point + rhat_colors # + rhat_lgnd
    } else {
      p.point <- geom_point(aes_string(x = "m", y = "y"), size = 4,
                            color = fill_color, fill = est_color, shape = 21)
      p.all + p.ci.2 + p.point
    }
  }
}



# rhat, ess, mcse ---------------------------------------------------------
stan_rhat <- function(object, pars, ...) {
  .check_object(object)
  .vb_check(object)
  if (missing(pars)) pars <- NULL
  .rhat_neff_mcse_hist(which = "rhat", object = object, pars = pars, ...)
}

stan_ess <- function(object, pars, ...) {
  .check_object(object)
  .vb_check(object)
  if (missing(pars)) pars <- NULL
  .rhat_neff_mcse_hist(which = "n_eff_ratio", object = object, pars = pars, ...)
}

stan_mcse <- function(object, pars, ...) {
  .check_object(object)
  .vb_check(object)
  if (missing(pars)) pars <- NULL
  .rhat_neff_mcse_hist(which = "mcse_ratio", object = object, pars = pars, ...)
}



# NUTS --------------------------------------------------------------------
stan_diag <- function(object, 
                      information = c("sample","stepsize","treedepth","divergence"),
                      chain = 0, ...) {
  .vb_check(object)
  nchains <- if (is.stanreg(object)) 
    ncol(object$stanfit) else ncol(object)
  if (!isTRUE(nchains > 1))
    stop("'stan_diag' requires more than one chain.", call. = FALSE)
  info <- match.arg(information)
  fn <- paste0("stan_", info)
  do.call(fn, list(object, chain, ...))
}

stan_stepsize <- function(object, chain = 0, ...) {
  .nuts_args_check(...)
  thm <- .rstanvis_defaults$theme
  stepsize <- .sampler_params_post_warmup(object, "stepsize__", as.df = TRUE)
  lp <- extract(if (is.stanreg(object)) object$stanfit else object,
                       pars = "lp__", permuted = FALSE)[,,1L]
  
  graphs <- list()
  graphs$stepsize_vs_lp <- .sampler_param_vs_param(p = lp, sp = stepsize[,-1L],
                                                   p_lab = .LP_LAB,
                                                   sp_lab = .STEPSIZE_LAB,
                                                   chain = chain, violin = TRUE)
  
  metrop <- .sampler_params_post_warmup(object, "accept_stat__", as.df = TRUE)
  graphs$stepsize_vs_metrop <-
    .sampler_param_vs_sampler_param_violin(round(stepsize[,-1L], 4), metrop[,-1L],
                                           lab_x = .STEPSIZE_LAB,
                                           lab_y = .METROP_LAB,
                                           chain = chain)
  graphs <- lapply(graphs, function(x) x + thm)
  .nuts_return(graphs, ...)
}

stan_sample <- function(object, chain = 0, ...) {
  .nuts_args_check(...)
  thm <- .rstanvis_defaults$theme
  hist_thm <- .rstanvis_defaults$hist_theme
  lp <- extract(if (is.stanreg(object)) object$stanfit else object,
                       pars = "lp__", permuted = FALSE)[,,1L]
  lp_df <- as.data.frame(cbind(iterations = 1:nrow(lp), lp))
  metrop <- .sampler_params_post_warmup(object, "accept_stat__", as.df = TRUE)
  graphs <- list()
  graphs$lp_hist <-
    .p_hist(lp_df, lab = .LP_LAB, chain = chain)
  graphs$metrop_hist <-
    .p_hist(metrop, lab = .METROP_LAB, chain = chain) + xlim(0,1)
  graphs <- lapply(graphs, function(x) x + thm)
  graphs$metrop_vs_lp <-
    .sampler_param_vs_param(p = lp, sp = metrop[,-1L], p_lab = .LP_LAB,
                            sp_lab = .METROP_LAB, chain = chain) +
    thm
  
  .nuts_return(graphs, ...)
}


stan_treedepth <- function(object, chain = 0, ...) {
  .nuts_args_check(...)
  thm <- .rstanvis_defaults$theme
  hist_thm <- .rstanvis_defaults$hist_theme
  lp <- extract(if (is.stanreg(object)) object$stanfit else object,
                       pars = "lp__", permuted = FALSE)[,,1L]
  treedepth <- .sampler_params_post_warmup(object, "treedepth__", as.df = TRUE)
  ndivergent <- .sampler_params_post_warmup(object, "divergent__", as.df = TRUE)
  metrop <- .sampler_params_post_warmup(object, "accept_stat__", as.df = TRUE)
  
  graphs <- graphs_nd <- list()
  
  graphs$treedepth_vs_lp <-
    .sampler_param_vs_param(p = lp, sp = treedepth[, -1L],
                            p_lab = .LP_LAB,
                            sp_lab = .TREEDEPTH_LAB,
                            chain = chain, violin = TRUE)
  graphs$treedepth_vs_metrop <-
    .sampler_param_vs_sampler_param_violin(treedepth[,-1L], metrop[,-1L],
                                           lab_x = .TREEDEPTH_LAB,
                                           lab_y = .METROP_LAB,
                                           chain = chain)
  
  graphs_nd$treedepth_ndivergent <-
    .treedepth_ndivergent_hist(treedepth, ndivergent, chain = chain,
                               divergent = "All")
  
  any_nd <- any(ndivergent[,-1L] != 0)
  if (any_nd) {
    graphs_nd$treedepth_ndivergent0 <-
      .treedepth_ndivergent_hist(treedepth, ndivergent, chain = chain,
                                 divergent = 0)
    graphs_nd$treedepth_ndivergent1 <-
      .treedepth_ndivergent_hist(treedepth, ndivergent, chain = chain,
                                 divergent = 1)
  }
  graphs <- lapply(graphs, function(x) x + thm)
  graphs_nd <- lapply(graphs_nd, function(x) x + hist_thm)
  graphs <- c(graphs, graphs_nd)
  .nuts_return(graphs, ...)
}


stan_divergence <- function(object, chain = 0, ...) {
  .nuts_args_check(...)
  thm <- .rstanvis_defaults$theme
  lp <- extract(if (is.stanreg(object)) object$stanfit else object,
                       pars = "lp__", permuted = FALSE)[,,1L]
  ndivergent <- .sampler_params_post_warmup(object, "divergent__", as.df = TRUE)
  metrop <- .sampler_params_post_warmup(object, "accept_stat__", as.df = TRUE)
  graphs <- list()
  graphs$ndivergent_vs_lp <-
    .sampler_param_vs_param(p = lp, sp = ndivergent[, -1L],
                            p_lab = .LP_LAB,
                            sp_lab = .NDIVERGENT_LAB,
                            chain = chain, violin = TRUE)
  graphs$ndivergent_vs_metrop <-
    .sampler_param_vs_sampler_param_violin(ndivergent[,-1L], metrop[,-1L],
                                           lab_x = .NDIVERGENT_LAB,
                                           lab_y = .METROP_LAB,
                                           chain = chain)
  graphs <- lapply(graphs, function(x) x + thm)
  .nuts_return(graphs, ...)
}


stan_par <- function(object, par, chain = 0, ...) {
  if (missing(par))
    stop("'par' must be specified", call. = FALSE)
  if (is.stanreg(object))
    object <- object$stanfit
  if (!isTRUE(ncol(object) > 1))
    stop("'stan_par' requires more than one chain.", call. = FALSE)
  thm <- .rstanvis_defaults$theme
  samp <- extract(object, pars = c("lp__", par), permuted = FALSE)
  par_sel <- which(dimnames(samp)$parameters == par)
  
  cntrl <- object@stan_args[[1L]]$control
  if (is.null(cntrl))
    max_td <- 11
  else {
    max_td <- cntrl$max_treedepth
    if (is.null(max_td))
      max_td <- 11
  }
  max_td <- .max_td(object)
  metrop <- .sampler_params_post_warmup(object, "accept_stat__", as.df = TRUE)[,-1L]
  stepsize <- .sampler_params_post_warmup(object, "stepsize__", as.df = TRUE)[,-1L]
  ndivergent <- .sampler_params_post_warmup(object, "divergent__", as.df = TRUE)[,-1L]
  treedepth <- .sampler_params_post_warmup(object, "treedepth__", as.df = TRUE)[,-1L]
  hit_max_td <- apply(treedepth, 2L, function(y) as.numeric(y == max_td))
  graphs <- list()
  par_samp <- samp[,, par_sel]
  lp <- samp[,, -par_sel]
  graphs[[paste0(par,"_vs_lp")]] <-
    .sampler_param_vs_param(sp = as.data.frame(lp),
                            p = par_samp,
                            divergent = ndivergent,
                            hit_max_td = as.data.frame(hit_max_td),
                            sp_lab = .LP_LAB, p_lab = par,
                            chain = chain)
  graphs[[paste0(par,"_vs_metrop")]] <-
    .sampler_param_vs_param(p = par_samp, sp = metrop,
                            divergent = ndivergent,
                            hit_max_td = as.data.frame(hit_max_td),
                            p_lab = par, sp_lab = .METROP_LAB,
                            chain = chain)
  graphs[[paste0(par,"_vs_stepsize")]] <-
    .sampler_param_vs_param(p = par_samp, sp = stepsize,
                            p_lab = par, sp_lab = .STEPSIZE_LAB,
                            chain = chain, violin = TRUE)
  
  graphs <- lapply(graphs, function(x) x + thm)
  .nuts_return(graphs, ...)
}
