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

stanc_process <- function(file, model_code = '', model_name = "anon_model",
                          auto_format = FALSE,
                          isystem = c(dirname(file), getwd())) {
  model_name2 <- deparse(substitute(model_code))
  if (is.null(attr(model_code, "model_name2")))
    attr(model_code, "model_name2") <- model_name2

  model_code <- get_model_strcode(file, model_code)
  if (missing(model_name) || is.null(model_name))
    model_name <- attr(model_code, "model_name2")

  model_attr <- attributes(model_code)
  model_code <- scan(text = model_code, what = character(), sep = "\n", quiet = TRUE)

  # Remove trailing whitespaces
  model_code <- trimws(model_code, "r")

  includes <- grep("^[[:blank:]]*#include ", model_code)
  while(length(includes) > 0) {
    for (i in rev(includes)) {
      header <- sub("^[[:blank:]]*#include[[:blank:]]+", "", model_code[i])
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
        model_code <- append(model_code, values = readLines(files[which(existent)[1]]),
                          after = i)
      else model_code <- append(model_code, values = readLines(header), after = 1)
      model_code[i] <- ""
    }
    includes <- grep("^[[:blank:]]*#include ", model_code)
  }

  model_code <- gsub('#include /(.*$)', '#include "\\1"', model_code)
  has_pound <- any(grepl("#", model_code, fixed = TRUE))

  if (has_pound && isFALSE(auto_format)) {
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

  mostattributes(model_code) <- model_attr

  return(model_code)
}

stanc_builder <- function(file, isystem = c(dirname(file), getwd()),
                          verbose = FALSE, obfuscate_model_name = FALSE,
                          allow_undefined = isTRUE(getOption("stanc.allow_undefined", FALSE)),
                          allow_optimizations = isTRUE(getOption("stanc.allow_optimizations", FALSE)),
                          standalone_functions = isTRUE(getOption("stanc.standalone_functions", FALSE)),
                          use_opencl = isTRUE(getOption("stanc.use_opencl", FALSE)),
                          warn_pedantic = isTRUE(getOption("stanc.warn_pedantic", FALSE)),
                          warn_uninitialized = isTRUE(getOption("stanc.warn_uninitialized", FALSE))) {
  stopifnot(is.character(file), length(file) == 1, file.exists(file))
  model_name <- sub("\\.[^.]*$", "", filename_rm_ext(basename(file)))

  model_cppname <- legitimate_model_name(model_name, obfuscate_name = obfuscate_model_name)

  auto_format <- isTRUE(getOption("stanc.auto_format", FALSE))
  if (isTRUE(auto_format)) {
    model_code <- stanc_process(file = file,
                                model_name = model_name,
                                auto_format = TRUE,
                                isystem = isystem)

    stopifnot(stanc_ctx$validate("stanc"))
    formatted_code <- try(stanc_ctx$call("stanc", model_cppname,
                          model_code, as.array("auto-format")),
                          silent = TRUE)
    if (!inherits(formatted_code, "try-error") && !is.null(formatted_code$result)) {
      model_code <- formatted_code$result
    }
  }

  model_code <- stanc_process(file = file,
                              model_name = model_name,
                              auto_format = FALSE,
                              isystem = isystem)

  out <- stanc(model_code = model_code,
               model_name = model_name, verbose = verbose,
               obfuscate_model_name = obfuscate_model_name,
               allow_undefined = allow_undefined,
               allow_optimizations = allow_optimizations,
               standalone_functions = standalone_functions,
               use_opencl = use_opencl,
               warn_pedantic = warn_pedantic,
               warn_uninitialized = warn_uninitialized)
  return(out)
}

stanc <- function(file, model_code = '', model_name = "anon_model",
                  verbose = FALSE, obfuscate_model_name = TRUE,
                  allow_undefined = isTRUE(getOption("stanc.allow_undefined", FALSE)),
                  allow_optimizations = isTRUE(getOption("stanc.allow_optimizations", FALSE)),
                  standalone_functions = isTRUE(getOption("stanc.standalone_functions", FALSE)),
                  use_opencl = isTRUE(getOption("stanc.use_opencl", FALSE)),
                  warn_pedantic = isTRUE(getOption("stanc.warn_pedantic", FALSE)),
                  warn_uninitialized = isTRUE(getOption("stanc.warn_uninitialized", FALSE)),
                  isystem = c(if (!missing(file)) dirname(file), getwd())) {
  if (missing(file)) {
    file <- tempfile(fileext = ".stan")
    on.exit(file.remove(file))
    writeLines(model_code, con = file)
  } else if (isTRUE(!nzchar(file)) || is.null(file)) {
    stop("Empty or invalid filename!")
  }

  if (missing(model_name) && is.character(file) && length(file) == 1 && file.exists(file))
    model_name <- sub("\\.[^.]*$", "", filename_rm_ext(basename(file)))

  model_cppname <- legitimate_model_name(model_name, obfuscate_name = obfuscate_model_name)

  auto_format <- isTRUE(getOption("stanc.auto_format", FALSE))
  if (isTRUE(auto_format)) {
    model_code <- stanc_process(file = file,
                                model_code = model_code,
                                model_name = model_name,
                                auto_format = TRUE,
                                isystem = isystem)

    stopifnot(stanc_ctx$validate("stanc"))
    formatted_code <- try(stanc_ctx$call("stanc", model_cppname,
                          model_code, as.array("auto-format")),
                          silent = TRUE)
    if (!inherits(formatted_code, "try-error") && !is.null(formatted_code$result)) {
      model_code <- formatted_code$result
    }
  }

  model_code <- stanc_process(file = file,
                              model_code = model_code,
                              model_name = model_name,
                              auto_format = FALSE,
                              isystem = isystem)

  if (isTRUE(rstan_options("threads_per_chain") > 1L)) {
    Sys.setenv("STAN_NUM_THREADS" = rstan_options("threads_per_chain"))
  }

  if (verbose)
    cat("\nTRANSLATING MODEL '", model_name, "' FROM Stan CODE TO C++ CODE NOW.\n", sep = '')

  stopifnot(stanc_ctx$validate("stanc"))
  stanc_flags <- c("allow-undefined",
                   "O1",
                   "standalone-functions",
                   "use-opencl",
                   "warn-pedantic",
                   "warn-uninitialized")
  istanc_flags <- c(allow_undefined,
                    allow_optimizations,
                    standalone_functions,
                    use_opencl,
                    warn_pedantic,
                    warn_uninitialized)
  if (sum(istanc_flags) >= 1) {
    stanc_flags <- as.array(stanc_flags[istanc_flags])
  } else {
    stanc_flags <- as.array("")
  }
  model_cppcode <- try(stanc_ctx$call("stanc", model_cppname, model_code, stanc_flags),
                       silent = TRUE)
  if (inherits(model_cppcode, "try-error")) {
    stop("parser failed badly; maybe try installing the V8 package")
  } else if (length(model_cppcode$errors)) {
    model_cppcode$status <- FALSE
    stop(paste(model_cppcode$errors, collapse = "\n"))
  } else {
    model_cppcode$status <- TRUE
  }

  # Make sure that the model name is not NULL
  if (is.null(model_name)) model_name <- "anon_model"

  # Use model_name in locations_array__
  cppcode <- model_cppcode$result
  cppcode <- gsub(paste0(" (in ", shQuote('string'), ", line "),
                  paste0(" (in ", shQuote(model_name), ", line "),
                  cppcode, fixed = TRUE)

  if (isTRUE(rstan_options("threads_per_chain") > 1L)) {
    # Initialize Stan/math TBB arena and global control
    cppcode <- paste("#ifdef STAN_THREADS",
                     "#ifndef RSTAN_THREADING",
                     "#define RSTAN_THREADING",
                     "#include <stan/math/prim/core/init_threadpool_tbb.hpp>",
                     "auto tbb_init = stan::math::init_threadpool_tbb();",
                     "#endif",
                     "#endif",
                     cppcode,
                     sep = "\n")
  }

  # Define USE_STANC3 for StanHeaders 2.26
  cppcode <- paste("#ifndef USE_STANC3",
                   "#define USE_STANC3",
                   "#endif",
                   cppcode,
                   sep = "\n")

  # Remove leading space from keywords
  cppcode <- scan(text = cppcode, what = character(), sep = "\n", quiet = TRUE)
  cppcode <- gsub("^ class ", "class ", cppcode)
  cppcode <- gsub("^ namespace ", "namespace ", cppcode)
  cppcode <- gsub("^ private:", "private:", cppcode)
  cppcode <- gsub("^ public:", "public:", cppcode)
  cppcode <- paste(cppcode, collapse = "\n")

  return(list(status = model_cppcode$status,
              model_cppname = model_cppname, cppcode = cppcode,
              model_name = model_name, model_code = model_code))
}
