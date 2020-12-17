# This file is part of RStan
# Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020 Trustees of Columbia University
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

stan_version <- function() {
  .Call(CPP_stan_version)
}

rstudio_stanc <- function(filename) {
  output <- stanc(filename, allow_undefined = TRUE)
  if (!output$status) {
    msg <- output$errors[2]
    line <- as.integer(sub("^.*line ([[:digit:]]+),.*$", "\\1", msg))
    column <- as.integer(sub("^.*column ([[:digit:]]+),.*$", "\\1", msg))
    msg <- gsub("\n", "", msg, fixed = TRUE)
    msg <- sub("^.*:", "", msg)
    rstudioapi::sourceMarkers(name = "Stan",
                              markers = data.frame(type = "error",
                                                   file = filename,
                                                   line = line,
                                                   column = column, 
                                                   message = msg,
                                                   stringsAsFactors = FALSE),
                              basePath = dirname(filename),
                              autoSelect = "none")
  } else {
    message(filename, " is syntactically correct.")
    return(invisible(output))
  }
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

stanc <- function(file, model_code = '', model_name = "anon_model", 
                  verbose = FALSE, obfuscate_model_name = TRUE,
                  allow_undefined = FALSE,
                  isystem = c(if (!missing(file)) dirname(file), getwd())) {
  model_name2 <- deparse(substitute(model_code))
  if (is.null(attr(model_code, "model_name2")))
    attr(model_code, "model_name2") <- model_name2
  model_code <- get_model_strcode(file, model_code)
  if (missing(model_name) || is.null(model_name))
    model_name <- attr(model_code, "model_name2")
  
  model_code <- scan(text = model_code, what = character(), sep = "\n", quiet = TRUE)
  model_code <- gsub('#include /(.*$)', '#include "\\1"', model_code)
  has_pound <- any(grepl("#", model_code, fixed = TRUE))
  if (has_pound) {
    unprocessed <- tempfile(fileext = ".stan")
    processed <- tempfile(fileext = ".stan")
    on.exit(file.remove(unprocessed))
    writeLines(model_code, con = unprocessed)
    ARGS <- paste("-E -nostdinc -x c++ -P -C", 
                  paste("-I", isystem, " ", collapse = ""), 
                  "-o", processed, unprocessed, "-Wno-invalid-pp-token")
    CPP <- system2(file.path(R.home(component = "bin"), "R"), 
                   args = "CMD config CC", stdout = TRUE)
    pkgbuild::with_build_tools(system(paste(CPP, ARGS), 
                                      ignore.stdout = TRUE, ignore.stderr = TRUE),
                               required = rstan_options("required") && 
                                 identical(Sys.getenv("WINDOWS"), "TRUE") &&
                                !identical(Sys.getenv("R_PACKAGE_SOURCE"), "") )
    if (file.exists(processed)) {
      on.exit(file.remove(processed), add = TRUE)
      model_code <- paste(readLines(processed), collapse = "\n")
    }
  } else model_code <- paste(model_code, collapse = "\n")

  if (verbose)
    cat("\nTRANSLATING MODEL '", model_name, "' FROM Stan CODE TO C++ CODE NOW.\n", sep = '')
  model_cppname <- legitimate_model_name(model_name, obfuscate_name = obfuscate_model_name)
  
  stanc_ctx <- V8::v8()
  stanc_ctx$source("https://github.com/stan-dev/stanc3/releases/download/nightly/stanc.js")
  stopifnot(stanc_ctx$validate("stanc"))
  model_cppcode <- try(stanc_ctx$call("stanc", model_cppname, model_code, 
                                      as.array(ifelse(allow_undefined, "allow_undefined", ""))), 
                       silent = TRUE)
  if (inherits(model_cppcode, "try-error")) {
    stop("parser failed badly")
  } else if (length(model_cppcode$errors)) {
    model_cppcode$status <- FALSE
    return(model_cppcode)
  }
  return(list(status = TRUE, model_cppname = model_cppname, cppcode = model_cppcode$result,
              model_name = model_name, model_code = model_code))
}
