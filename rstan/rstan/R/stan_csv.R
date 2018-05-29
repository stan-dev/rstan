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

paridx_fun <- function(names) {
  # Args:
  #   names: names (character vector) such as lp__, treedepth__, stepsize__,
  #          alpha, beta.1, 
  # Returns: 
  #   The indexes in the names that are parameters other than lp__,
  #   treedepth__, or stepsize__. The vector has attribute meta
  #   with the indexes of 'treedepth__', 'lp__', and 'stepsize__'
  #   if available. 
  
  sampler_param_names <- c('lp__', 'accept_stat__', 'treedepth__', 'stepsize__', 
                           'divergent__', 'n_leapfrog__', "energy__")
  metaidx <- match(sampler_param_names, names)
  names(metaidx) <- sampler_param_names
  paridx <- setdiff(seq_along(names), metaidx)
  attr(paridx, "meta") <- metaidx[!sapply(metaidx, is.na)]
  paridx
}

parse_stancsv_comments <- function(comments) {
  # Parse the comments in Stan CSV files to get information such as
  # iter, thin, seed, etc. This is specific to the CSV files
  # generated from Stan

  adapt_term_lineno <- which(grepl("Adaptation terminated", comments))[1]
  if (is.na(adapt_term_lineno)) adapt_term_lineno <- length(comments)
  time_lineno <- which(grepl("Elapsed Time", comments))
  has_time <- length(time_lineno) > 0
  len <- length(comments)
  if (length(adapt_term_lineno) < 1) 
    adapt_term_lineno <- len
  if (length(time_lineno) < 1)
    warning("line with \"Elapsed Time\" not found")

  if (adapt_term_lineno == len)
    adaptation_info <- ''
  else {
    if (has_time)
      adaptation_info <- paste(comments[(adapt_term_lineno+1):(time_lineno-1)], collapse = '\n')
    else
      adaptation_info <- paste(comments[(adapt_term_lineno+1):len], collapse = '\n')
  }
  if (has_time)
    time_info <- comments[time_lineno:len]
  else
    time_info <- ''
  if (adapt_term_lineno > 0) comments <- head(comments, adapt_term_lineno - 1L)

  has_eq <- sapply(comments, function(i) grepl('=', i))
  comments <- comments[has_eq] 
  comments <- gsub('^#+\\s*|\\s*|\\(Default\\)', '', comments)
  eq_pos <- regexpr("=", comments, fixed = TRUE)
  names0 <- substr(comments, 0, eq_pos - 1)
  values <- as.list(substring(comments, eq_pos + 1))
  
  id_idx <- which("id" == names0)
  if (length(id_idx) > 0) 
  names0[id_idx] <- "chain_id"
  
  compute_iter <- FALSE
  id_warmup <- which("num_warmup" == names0)
  if (length(id_warmup) > 0) {
    names0[id_warmup] <- "warmup"
    compute_iter <- TRUE
  }   
  
  id_numsamples <- which("num_samples" == names0)
  if (length(id_numsamples) > 0) {
    names0[id_numsamples] <- "iter"
  }
  names(values) <- names0;

  add_lst <- list(adaptation_info = adaptation_info,
                  has_time = has_time,
                  time_info = time_info)

  sampler_t <- NULL
  if (!is.null(values$algorithm) && is.null(values$sampler_t)) {
    if (values$algorithm == 'rwm' || values$algorithm == 'Metropolis')  
      sampler_t <- "Metropolis"
    else if (values$algorithm == 'hmc') {
       if (values$engine == 'static')  sampler_t <- "HMC"
       else {
         if (values$metric == 'unit_e') sampler_t <- "NUTS(unit_e)"
         else if (values$metric == 'diag_e') sampler_t <- "NUTS(diag_e)"
         else if (values$metric == 'dense_e') sampler_t <- "NUTS(dense_e)"
       } 
    } 
    add_lst <- c(add_lst, sampler_t = sampler_t)
  } 
  names1 <- intersect(c("thin", "iter", "warmup", "chain_id", "max_depth", 
                        "num_samples", "num_warmup", "id",
                        "max_treedepth", "save_warmup"), names0)
  names2 <- intersect(c("stepsize", "stepsize_jitter", "adapt_gamma", "adapt_kappa", 
                        "adapt_delta", "gamma", "kappa", "delta", "t0",
                        "adapt_t0"), names0) 
  for (z in names1) values[[z]] <- as.integer(values[[z]])
  for (z in names2) values[[z]] <- as.numeric(values[[z]])
  if (compute_iter) values[["iter"]] <- values[["iter"]] + values[["warmup"]]
  c(values, add_lst)  
}


read_stan_csv <- function(csvfiles, col_major = TRUE) {
  # Read the csv files saved from Stan (or RStan) to a stanfit object
  # Args:
  #   csvfiles: csv files fitted for the same model; each file contains 
  #     the sample of one chain 
  #   col_major: the order for array parameters. 
  # 

  if (length(csvfiles) < 1) 
    stop("csvfiles does not contain any CSV file name")

  g_skip <- 10
  
  ss_lst <- vector("list", length(csvfiles))
  cs_lst2 <- vector("list", length(csvfiles))

  for (i in seq_along(csvfiles)) {
    header <- read_csv_header(csvfiles[i])
    lineno <- attr(header, 'lineno')
    vnames <- strsplit(header, ",")[[1]]
    iter.count <- attr(header,"iter.count")
    variable.count <- length(vnames)
    df <- structure(replicate(variable.count,list(numeric(iter.count))),
                    names = vnames,
                    row.names = c(NA,-iter.count),
                    class = "data.frame")
    comments = character()
    con <- file(csvfiles[[i]],"rb")
    buffer.size <- min(ceiling(1000000/variable.count),iter.count)
    row.buffer <- matrix(ncol=variable.count,nrow=buffer.size)
    row <- 1
    buffer.pointer <- 1  
    while(length(char <- readBin(con,'int',size=1L)) > 0) {
      # back up 1 character, since we already looked at one to check for comment
      seek(con,origin="current",-1)
      if(char == 35){ #35 is '#'
        line <- readLines(con, n = 1)
        comments <- c(comments, line)
        next
      }
      if(char == 108){ #start of lp__ in header 108
        readLines(con, n = 1)
        next
      }
      row.buffer[buffer.pointer,] <- scan(con, nlines=1, sep="," ,quiet=TRUE)
      if(buffer.pointer == buffer.size){
        df[row:(row + buffer.size - 1), ] <- row.buffer
        row <- row + buffer.size
        buffer.pointer <- 0
      }
      buffer.pointer <- buffer.pointer + 1
      
    }
    if(buffer.pointer > 1){
      df[row:(row + buffer.pointer - 2), ] <- row.buffer[1:(buffer.pointer-1), ]
    }

    close(con)
    cs_lst2[[i]] <- parse_stancsv_comments(comments)
    ss_lst[[i]] <- df
  } 

  # use the first CSV file name as model name
  m_name <- sub("(_\\d+)*$", '', filename_rm_ext(basename(csvfiles[1])))

  sdate <- do.call(max, lapply(csvfiles, function(csv) file.info(csv)$mtime))
  sdate <- format(sdate, "%a %b %d %X %Y") # same format as date() 

  chains <- length(ss_lst)
  fnames <- names(ss_lst[[1]])
  n_save <- nrow(ss_lst[[1]])
  paridx <- paridx_fun(fnames)
  lp__idx <- attr(paridx, 'meta')["lp__"]
  par_fnames <- c(fnames[paridx], "lp__")
  pars_oi <- unique_par(par_fnames)
  dims_oi <- lapply(pars_oi, 
                    function(i) {
                      pat <- paste('^', i, '(\\.\\d+)*$', sep = '')
                      i_fnames <- par_fnames[grepl(pat, par_fnames)]
                      get_dims_from_fnames(i_fnames, i) 
                    })
  names(dims_oi) <- pars_oi
  midx <- if (!col_major) multi_idx_row2colm(dims_oi) else 1:length(par_fnames)
  if (chains > 1) {
    if (!all(sapply(ss_lst[-1], function(i) identical(names(i), fnames))))
      stop('the CSV files do not have same parameters')
    if (!all(sapply(ss_lst[-1], function(i) identical(length(i[[1]]), n_save)))) 
      stop('the number of iterations are not the same in all CSV files')
  } 
  mode <- 0L
  
  samples <- lapply(ss_lst, 
                    function(df) {
                      ss <- df[c(paridx, lp__idx)[midx]]
                      attr(ss, "sampler_params") <- df[setdiff(attr(paridx, 'meta'), lp__idx)] 
                      ss
                    })
  par_fnames <- par_fnames[midx]
  for (i in seq_along(samples)) {
    attr(samples[[i]], "adaptation_info") <- cs_lst2[[i]]$adaptation_info 
    attr(samples[[i]], "args") <- 
      list(sampler_t = cs_lst2[[i]]$sampler_t,
           chain_id = cs_lst2[[i]]$chain_id)
    if (cs_lst2[[i]]$has_time)
      attr(samples[[i]], "elapsed_time") <- get_time_from_csv(cs_lst2[[i]]$time_info)
  } 

  save_warmup <- sapply(cs_lst2, function(i) i$save_warmup)
  warmup <- sapply(cs_lst2, function(i) i$warmup)
  thin <- sapply(cs_lst2, function(i) i$thin)
  iter <- sapply(cs_lst2, function(i) i$iter)
  if (!all_int_eq(warmup) || !all_int_eq(thin) || !all_int_eq(iter)) 
    stop("not all iter/warmups/thin are the same in all CSV files")
  n_kept0 <- 1 + (iter - warmup - 1) %/% thin
  warmup2 <- 0
  if (max(save_warmup) == 0L) { # all equal to 0L
    n_kept <- n_save
  } else if (min(save_warmup) == 1L) { # all equals to 1L 
    warmup2 <- 1 + (warmup[1] - 1) %/% thin[1]
    n_kept <- n_save - warmup2 
  } 
  if (n_kept0[1] != n_kept) {
    warning("the number of iterations after warmup found (", n_kept, 
            ") does not match iter/warmup/thin from CSV comments (",
            paste(n_kept0, collapse = ','), ")")

    if (n_kept < 0) {
      warmup <- warmup + n_kept
      n_kept <- 0
      mode <- 2L
    }
    n_kept0 <- n_save
    iter <- n_save

    for (i in 1:length(cs_lst2)) {
      cs_lst2[[i]]$warmup <- warmup
      cs_lst2[[i]]$iter <- iter
    }
  }

  idx_kept <- if (warmup2 == 0) 1:n_kept else -(1:warmup2)
  for (i in seq_along(samples)) {
    m <- vapply(samples[[i]], function(x) mean(x[idx_kept]), numeric(1))
    attr(samples[[i]], "mean_pars") <- m[-length(m)]
    attr(samples[[i]], "mean_lp__") <- m["lp__"]
  }

  perm_lst <- lapply(1:chains, function(id) sample.int(n_kept))

  sim = list(samples = samples, 
             iter = iter[1], 
             thin = thin[1], 
             warmup = warmup[1], 
             chains = chains, 
             n_save = rep(n_save, chains),
             warmup2 = rep(warmup2, chains),
             permutation = perm_lst,
             pars_oi = pars_oi, 
             dims_oi = dims_oi,
             fnames_oi = dotfnames_to_sqrfnames(par_fnames), 
             n_flatnames = length(par_fnames))
  null_dso <- new("cxxdso", sig = list(character(0)), dso_saved = FALSE, dso_filename = character(0), 
                  modulename = character(0), system = R.version$system, cxxflags = character(0), 
                 .CXXDSOMISC = new.env(parent = emptyenv()))
  null_sm <- new("stanmodel", model_name = m_name, model_code = character(0), 
                 model_cpp = list(), dso = null_dso)

  nfit <- new("stanfit", 
              model_name = m_name,
              model_pars = pars_oi,
              par_dims = dims_oi, 
              mode = mode,
              sim = sim,
              inits = list(), 
              stan_args = cs_lst2,
              stanmodel = null_sm,
              date = sdate, # not the time of sampling
              .MISC = new.env(parent = emptyenv()))
  return(nfit)
}

read_one_stan_csv <- function(csvfile) {
  if (length(csvfile) != 1) 
    stop("'csvfile' must be of length 1")
  if (!file.exists(csvfile))
    stop("'csvfile' does not exist on the disk")
  
  mark <- 0L
  fields <- character()
  while(length(fields) == 0) {
    mark <- mark + 1L
    fields <- scan(csvfile, what = character(), sep = ",", 
                   comment.char = "#", nlines = mark, quiet = TRUE)
  }
  
  comments <- scan(csvfile, what = character(), sep = "\n", comment.char = "", 
                   nlines = mark + 2L, quiet = TRUE)
  comments <- gsub("#", "", comments, fixed = TRUE)
  comments <- gsub("(Default)", "", comments, fixed = TRUE)
  comments <- grep("=", comments, fixed = TRUE, value = TRUE)
  comments <- strsplit(comments, split = "=", fixed = TRUE)
  comments <- lapply(comments, FUN = trimws)
  comments <- sapply(comments, FUN = function(x) {
    y <- x[2]
    names(y) <- x[1]
    return(y)
  })
  
  method <- comments["algorithm"]
  if (method %in% c("meanfield", "fullrank")) {
    draws <- scan(csvfile, what = double(), sep = ",", comment.char = "", 
                  quiet = TRUE, skip = mark + 2L, 
                  nlines = mark + as.integer(comments["output_samples"]) + 3L)
    timings <- NULL
  }
  else { # sampling
    iter <- as.integer(comments["iter"])
    draws <- scan(csvfile, what = double(), sep = ",", comment.char = "", 
                  quiet = TRUE, skip = mark, nlines = mark + iter)
    timings <- scan(csvfile, what = character(), sep = "\n", comment.char = "", 
                    quiet = TRUE, skip = mark + iter)
  }
  draws <- matrix(draws, ncol = length(fields), byrow = TRUE)
  colnames(draws) <- fields
  draws <- as.data.frame(draws)
  
  attributes(draws)$comments <- comments
  attributes(draws)$timings <- timings
  return(draws)
}

if (!exists("trimws")) trimws <- function(x) gsub("^\\s+|\\s+$", "", x)
