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

###  deal with the optimization level for c++ compilation 

get_all_makefile_paths <- function() {
  out <- system2(file.path(Sys.getenv("R_HOME"), "bin", "R"),
                 args = "CMD SHLIB --dry-run", stdout = TRUE)
  out <- grep("SHLIB", out, value = TRUE)[1]
  makefiles <- strsplit(sub("SHLIB.*$", "", out), split = "-f ")[[1]][-1]
  makefiles <- gsub("'", "", makefiles)
  makefiles <- gsub("[[:space:]]*$", "", makefiles)
  makefiles <- gsub('\"', '', makefiles)
  makefiles <- makefiles[file.exists(makefiles)]
  return(makefiles)  
}

get_makefile_txt <- function() { 
  # get the all makefile content used for R CMD SHLIB 
  # The order or files to look for are following the code 
  # of function .SHLIB in packages tools/R/install.R 
  # 
  # Return: 
  #   a character vector in which each element is a line of 
  #   the makefile 
  makefiles <- get_all_makefile_paths()
  do.call(c, lapply(makefiles, function(f) readLines(f, warn = FALSE))) 
}

get_makefile_flags <- function(FLAGNAME, makefile_txt, headtotail = FALSE) { 
  # get current CXXFLAGS used for in R CMD SHLIB, which controls
  # how the generated C++ code for stan model are compiled 
  # at which optimization level. 
  # Args: 
  #   FLAGNAME: the name of the flag of interest such as CXXFLAGS 
  #   makefile_txt: a character vector each element of which is a line
  #     of the make file, typically obtained by get_makefile_txt() 
  # Return: 
  #   the whole line that defines CXXFLAGS 

  if (missing(makefile_txt)) 
    makefile_txt <- get_makefile_txt() 
  if (length(makefile_txt) == 0) return(as.character(NULL))

  lineno <- -1 
  iseq <- if (headtotail) 1:length(makefile_txt) else length(makefile_txt):1 
  for (i in iseq) { 
    pattern <- paste("^\\s*", FLAGNAME, "\\s*=", sep = '')
    if (!is.null(makefile_txt[i]) && grepl(pattern, makefile_txt[i])) {
      lineno <- i 
      break;
    }
  } 

  if (-1 != lineno) return(makefile_txt[lineno]) 
  paste(FLAGNAME, " = ", sep = '') 
} 

last_makefile <- function() {
  # Return the last makefile for 'R CMD SHLIB'. 
  # In essential, R CMD SHLIB uses GNU make.  
  # R CMD SHLIB uses a series of files as the whole makefile.  
  # This function return the last one, where we can set flags 
  # to overwrite what is set before. 
  # 
  makefiles <- get_all_makefile_paths()
  return(tail(makefiles, n = 1L))
} 

set_makefile_flags <- function(flags) { 
  # Set user-defined CXXFLAGS for R CMD SHLIB by 
  # set CXXFLAGS or others in the last makefile obtained by last_makefile 
  # Args: 
  #   flags: a named list with names CXXFLAGS and R_XTRA_CPPFLAGS 
  # 
  # Return: 
  #   TRUE if it is successful, otherwise, the function
  #   stops and reports an error. 
  # 
  lmf <- last_makefile() 
  homedotR <- dirname(lmf) 
  if (!file.exists(homedotR)) {
    if (!dir.create(homedotR, showWarnings = FALSE, recursive = TRUE))
      stop(paste("failed to create directory ", homedotR, sep = '')) 
  } 

  flagnames <- names(flags) 
  paste_bn <- function(...)  paste(..., sep = '\n') 
  flags <- lapply(flags, 
                  function(x) {
                    if (grepl('#set_by_rstan', x)) return(x)
                    paste(x, " #set_by_rstan")  
                  })

  # the Makevars file does not exist 
  if (!file.exists(lmf)) {
    if (file.access(homedotR, mode = 2) < 0 || 
        file.access(homedotR, mode = 1) < 0 ) 
      stop(paste("directory ", homedotR, " is not writable", sep = '')) 
    cat(paste("# created by rstan at ", date(), sep = ''), '\n', file = lmf) 
    cat(do.call(paste_bn, flags), file = lmf, append = TRUE) 
    return(invisible(NULL)) 
  } 

  if (file.access(lmf, mode = 2) < 0) 
    stop(paste(lmf, " is not writable", sep = '')) 

  lmf_txt <- readLines(lmf, warn = FALSE) 
  if (length(lmf_txt) < 1) {  # empty file 
    cat(do.call(paste_bn, flags), file = lmf, append = TRUE) 
    return(invisible(NULL))  
  } 

  # change the existing file 
  makefile_exta <- "" 
  for (fname in flagnames) { 
    found <- FALSE 
    for (i in length(lmf_txt):1) {
      pattern <- paste("^\\s*", fname, "\\s*=", sep = '')
      if (grepl(pattern, lmf_txt[i])) { 
        lmf_txt[i] <- flags[[fname]] 
        found <- TRUE 
        break
      }
    }
    if (!found)
      makefile_exta <- paste(makefile_exta, flags[[fname]], sep = '\n') 
  } 

  cat(paste(lmf_txt, collapse = '\n'), file = lmf) 
  cat(makefile_exta, file = lmf, append = TRUE) 
  invisible(NULL) 
} 

rm_last_makefile <- function(force = TRUE) {
  # remove the file given by last_makefile()
  # return: 
  #  0 for success, 1 for failure, invisibly. 
  #  -1 file does not exist
  lmf <- last_makefile() 
  msgs <- c(paste(lmf, " does not exist", sep = ''), 
            paste("removed file ", lmf, sep = ''), 
            paste("failed to remove file ", lmf, sep = '')) 
  ret <- -1 
  if (file.exists(lmf)) ret <- unlink(lmf, force = force) 
  cat(msgs[ret + 2], '\n') 
  invisible(ret) 
} 

get_cxxo_level <- function(str) {
  # obtain the optimization level from a string like -O3 
  # Return: 
  #   the optimization level after -O as a character 
  # 
  p <- regexpr('-O.', str, perl = TRUE) + 2 
  if (p == 1) return("")  # not found 
  substr(str, p, p)
} 

if_debug_defined <- function(str) {
  grepl("DDEBUG", str)
} 

if_ndebug_defined <- function(str) {
  grepl("DNDEBUG", str)
} 


rm_rstan_makefile_flags <- function() {
  # remove flags in $HOME/.R/Makevars with #rstan 
  lmf <- last_makefile() 
  if (!file.exists(lmf)) {
    message("no user Makevars file found; nothing need be done.")
    return(invisible(NULL))
  }
  if (file.access(lmf, mode = 2) < 0) 
    stop(paste(lmf, " is not writable", sep = '')) 
  lmf_txt <- readLines(lmf, warn = FALSE) 
  if (length(lmf_txt) < 1) return(invisible(NULL))
  for (i in length(lmf_txt):1) {
    if (grepl("#set_by_rstan", lmf_txt[i])) 
      lmf_txt[i] <- NA  
  } 
  lmf_txt <- lmf_txt[!is.na(lmf_txt)] 
  cat(paste(lmf_txt, collapse = '\n'), file = lmf) 
  message("compiler flags set by rstan are removed.") 
  invisible(NULL)
} 


#' Set environment variables
#' @param new a named vector of names and values to assign as 
#'  environment variables.
#' @param action should new values "replace", "prefix" or "suffix"
#'  existing variables with the same name.
.set_enviro_var <- function (envs, action = "replace") {
  if (length(envs) == 0) {
    return()
  }
  stopifnot(!is.null(names(envs)) && all(names(envs) != ""))
  stopifnot(is.character(action), length(action) == 1)
  action <- match.arg(action, c("replace", "prefix", "suffix"))
  envs <- envs[!duplicated(names(envs), fromLast = TRUE)]
  old <- Sys.getenv(names(envs), names = TRUE, unset = NA)
  set <- !is.na(envs)
  both_set <- set & !is.na(old)
  if (any(both_set)) {
    if (action == "prefix") {
      envs[both_set] <- paste(envs[both_set], old[both_set])
    }
    else if (action == "suffix") {
      envs[both_set] <- paste(old[both_set], envs[both_set])
    }
  }
  if (any(set)) 
    do.call("Sys.setenv", as.list(envs[set]))
  if (any(!set)) 
    Sys.unsetenv(names(envs)[!set])
  invisible(old)
}

#' Create a local version of makevars which is cleaned up after exiting
#'  the parent frame. 
#'  @note Has the same behavior as `withr::with_makevars` but closes the 
#'   temporary makevar file after the result object leaves scope.
#'  @param new a named vector of names and values to assign to the makevars values.
#'  @param path Path to Makevars. By default is `withr::makevars_user()`.
#'  @param assignment Assignment type to use.
#'  @param local_env The environment that when exited will cleanup the temp Makevars.
.local_makevars <- function(new, path = withr::makevars_user(),
                            assignment = c("=", ":=", "?=", "+="), local_env = parent.frame()) {
  RMU <- Sys.getenv("R_MAKEVARS_USER")
  assignment <- match.arg(assignment)
  makevars_file <- tempfile(pattern = "Makevars")
  setup_makevar <- function(state, RMU) {
    withr::set_makevars(new, path, state, assignment = assignment)
    old <- .set_enviro_var(envs = c(R_MAKEVARS_USER = state), action = "replace")
    withr::defer(.set_enviro_var(old), envir = local_env)
    invisible(old)
  }
  cleanup_makevar <- function(state, RMU) {
    unlink(state)
    if (RMU == "") {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKE_VARS_USER = RMU)
    }
  }
  return(withr::local_(setup_makevar, cleanup_makevar)(makevars_file, RMU))
}


#' Remove march flag from Makevars
.remove_march_makevars <- function() {
  makevar_files <- withr::makevars_user()
  cxx_flags <- grep("^CXX.*FLAGS", readLines(file.path(makevar_files)), value = TRUE)
  trimmed_flag <- trimws(gsub("-march[[:space:]]*=[[:space:]]*native", "",
    substr(cxx_flags, regexpr("-", cxx_flags), nchar(cxx_flags))))
  cxxflag_locations <- attr(regexpr("CXX.*FLAGS", cxx_flags), "match.length")
  cxx_flag_name <- substr(cxx_flags, 0, cxxflag_locations)
  names(trimmed_flag) <- cxx_flag_name
  no_march_makevar <- .local_makevars(trimmed_flag, local_env = parent.frame())
  return(no_march_makevar)
}
