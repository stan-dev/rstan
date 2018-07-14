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

stanc <- function(file, model_code = '', model_name = "anon_model",
                  verbose = FALSE, obfuscate_model_name = TRUE,
                  allow_undefined = FALSE, 
                  isystem = c(if (!missing(file)) dirname(file), getwd())) {

  # Call stanc, written in C++ 
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
  
  # model_name in C++, to avoid names that would be problematic in C++. 
  model_cppname <- legitimate_model_name(model_name, obfuscate_name = obfuscate_model_name) 
  r <- .Call("CPP_stanc280", model_code, model_cppname, allow_undefined, isystem)
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
    lapply(r$msg, function(x) message(x))
    error_msg <- paste("failed to parse Stan model '", model_name,
                       "' due to the above error.", sep = '')
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
  output <- stanc_builder(filename)
  message(filename, " is syntactically correct.")
  return(invisible(output))
}

stanc_builder <- function(file, isystem = c(dirname(file), getwd()),
                          verbose = FALSE, obfuscate_model_name = FALSE, 
                          allow_undefined = FALSE) {
  stopifnot(is.character(file), length(file) == 1, file.exists(file))
  model_cppname <- sub("\\.stan$", "", basename(file))
  program <- readLines(file)
  includes <- grep("^[[:blank:]]*#include ", program)
  while(length(includes) > 0) {
    for (i in rev(includes)) {
      header <- sub("^[[:blank:]]*#include[[:blank:]]+", "", program[i])
      header <- gsub('\\"', '', header)
      header <- gsub("\\'", '', header)
      header <- sub("<", "", header, fixed = TRUE)
      header <- sub(">", "", header, fixed = TRUE)
      header <- sub("//.*$", "", header)
      header <- sub("/\\*.*$", "", header)
      header <- sub("#.*$", "", header)
      header <- sub("[[:blank:]]*$", "", header)
      files <- file.path(isystem, header)
      existent <- file.exists(files)
      if (any(existent))
        program <- append(program, values = readLines(files[which(existent)[1]]), 
                          after = i)
      else program <- append(program, values = readLines(header), after = 1)
      program[i] <- ""
    }
    includes <- grep("^[[:blank:]]*#include ", program)
  }    
  out <- stanc(model_code = paste(program, collapse = "\n"), 
               model_name = model_cppname, verbose = verbose, 
               obfuscate_model_name = obfuscate_model_name, 
               allow_undefined = allow_undefined)
  return(out)
}
