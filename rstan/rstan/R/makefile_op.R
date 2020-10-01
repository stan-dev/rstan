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

#' Shim for tools::makevars_user()
#' @keywords internal
#' @export
makevars_user <- function() {
  if (getRversion() >= "3.3") {
    return(tools::makevars_user())
  }
  # Below is tools::makevars_user() from R 3.6.2
  m <- character()
  if (.Platform$OS.type == "windows") {
    if (!is.na(f <- Sys.getenv("R_MAKEVARS_USER", NA_character_))) {
      if (file.exists(f))
        m <- f
    }
    else if ((Sys.getenv("R_ARCH") == "/x64") && file.exists(f <- path.expand("~/.R/Makevars.win64")))
      m <- f
    else if (file.exists(f <- path.expand("~/.R/Makevars.win")))
      m <- f
    else if (file.exists(f <- path.expand("~/.R/Makevars")))
      m <- f
  }
  else {
    if (!is.na(f <- Sys.getenv("R_MAKEVARS_USER", NA_character_))) {
      if (file.exists(f))
        m <- f
    }
    else if (file.exists(f <- path.expand(paste0("~/.R/Makevars-",
            Sys.getenv("R_PLATFORM")))))
      m <- f
    else if (file.exists(f <- path.expand("~/.R/Makevars")))
      m <- f
  }
  m
}


#' Remove march flag from Makevars
.warn_march_makevars <- function() {
  makevar_files <- makevars_user()
  if (length(makevar_files)) {
    cxx_flags <- grep("^CXX.*FLAGS", readLines(file.path(makevar_files)), value = TRUE)
  } else {
    cxx_flags <- character(0)
  }
  if (length(cxx_flags) != 0) {
    has_march_native <- any(grepl("-march[[:space:]]*=[[:space:]]*native",
                                substr(cxx_flags, regexpr("-", cxx_flags), nchar(cxx_flags))))
    if (has_march_native && !rstan_options("disable_march_warning")) {
      warning(paste0("Detected -march=native in the Makevars file at '", makevar_files,
                     "'. Compiling with the -march=native flag on windows with Rtools",
                     " can cause crashes because of the compiler implementation. It's",
                     " recommended you remove -march=native from your Makevar file.",
                     " You can disable this warning by setting ",
                     'rstan_options(disable_march_warning = TRUE)'))
    }
  }
  return(NULL)
}
