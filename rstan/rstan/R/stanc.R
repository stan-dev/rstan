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

stanc <- function(file, model_code = '', model_name = "anon_model", verbose = FALSE, ...) {
  # Call stanc, written in C++ 
  # Args:
  # ..., to pass argument obfuscate_model_name 
  # 
  model_name2 <- deparse(substitute(model_code))  
  if (is.null(attr(model_code, "model_name2"))) 
    attr(model_code, "model_name2") <- model_name2 
  model_code <- get_model_strcode(file, model_code)  
  if (missing(model_name) || is.null(model_name)) 
    model_name <- attr(model_code, "model_name2") 
  if (verbose) 
    cat("\nTRANSLATING MODEL '", model_name, "' FROM Stan CODE TO C++ CODE NOW.\n", sep = '')
  SUCCESS_RC <- 0 
  EXCEPTION_RC <- -1
  PARSE_FAIL_RC <- -2 
  
  dotlst <- list(...) 
  omn_con <- "obfuscate_model_name" %in% names(dotlst) 

  obfuscate_name <- if (!omn_con) TRUE else as.logical(dotlst[["obfuscate_model_name"]]) 
  if (is.na(obfuscate_name))  obfuscate_name <- FALSE
  # model_name in C++, to avoid names that would be problematic in C++. 
  model_cppname <- legitimate_model_name(model_name, obfuscate_name = obfuscate_name) 
  r <- .Call("CPP_stanc261", model_code, model_cppname)
  # from the cpp code of stanc,
  # returned is a named list with element 'status', 'model_cppname', and 'cppcode' 
  r$model_name <- model_name  
  r$model_code <- model_code 
  if (is.null(r)) {
    stop(paste("failed to run stanc for model '", model_name,
               "' and no error message provided", sep = ''))
  } else if (r$status == PARSE_FAIL_RC) {
    stop(paste("failed to parse Stan model '", model_name,
               "' and no error message provided"), sep = '')
  } else if (r$status == EXCEPTION_RC) {
    error_msg <- paste("failed to parse Stan model '", model_name,
                       "' with error message:\n", r$msg, sep = '')
    error_msg_len <- min(8170, nchar(error_msg))
    warning.length <- getOption('warning.length')
    if (error_msg_len > warning.length) {
      options(warning.length = error_msg_len)
      on.exit(options(warning.length = warning.length), add = TRUE)
    }
    stop(error_msg)
  } 

  if (r$status == SUCCESS_RC && verbose)
    cat("successful in parsing the Stan model '", model_name, "'.\n", sep = '')

  r$status = !as.logical(r$status)
  return(r)
}


stan_version <- function() {
  .Call('CPP_stan_version')
}

rstudio_stanc <- function(filename) {
  output <- stanc(filename)
  message(filename, " is syntactically correct.")
  return(invisible(output))
}
