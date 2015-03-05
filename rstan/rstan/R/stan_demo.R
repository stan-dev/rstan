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

delete_stan_demo_folder <- function() {
  MODELS_HOME <- system.file('include', 'example-models-master', package = 'rstan')
  unlink(MODELS_HOME, recursive = TRUE) 
}

stan_demo <-
function(model = character(0), 
         method = c("sampling", "optimizing"), ...) {
  if(is.numeric(model)) {
    MODEL_NUM <- as.integer(model)
    model <- character(0)
  }
  else MODEL_NUM <- -1L
  MODELS_HOME <- system.file('include', 'example-models-master', package = 'rstan')
  if(MODELS_HOME == "") {
    MODELS_HOME <- file.path(tempdir(), "example-models-master")
    if(!file.exists(MODELS_HOME)) {
      writable <- file.access(system.file('include', package = 'rstan'), 
                              mode = 2) == 0
      choices <- c(if(writable) "Download to include directory of the rstan package",
                   "Download to the temporary directory", "Do not download")
      if(interactive()) choice <- menu(choices, title = "Do you want to download the example models?")
      else choice <- if(writable) 2L else 1L
      if(choice == 0 | choice == length(choices)) stop("cannot proceed without example models")
      if(choices[choice] == "Download to include directory of the rstan package") {
        FILE <- file.path(system.file('include', package = 'rstan'),
                          "example-models-master.zip")
      }
      else FILE <- file.path(tempdir(), "example-models-master.zip")
      if (!require(RCurl)) stop("cannot proceed without R package RCurl being installed")
      writeBin(getBinaryURL("https://github.com/stan-dev/example-models/archive/master.zip",
                            .opts = curlOptions(followlocation = TRUE)), 
               FILE)
      unzip(FILE, exdir = dirname(FILE))
      MODELS_HOME <- file.path(dirname(FILE), "example-models-master")
    }
  }
  WD <- getwd()
  on.exit(setwd(WD))
  setwd(MODELS_HOME)
  MODELS <- dir(MODELS_HOME, pattern = paste0(model, ".stan", "$"), 
                recursive = TRUE, full.names = FALSE)
  if(length(MODELS) == 0) {
    stop("'model' not found; leave 'model' unspecified to see all choices")
  }
  else if(length(MODELS) > 1) {
    if(MODEL_NUM %in%  1:length(MODELS)) {
      MODELS <- MODELS[MODEL_NUM]
    }
    else if(MODEL_NUM == 0) MODELS <- ""
    else MODELS <- select.list(MODELS)
    if(!nzchar(MODELS)) {
      return(dir(MODELS_HOME, pattern = paste0(model, ".stan", "$"), 
                 recursive = TRUE, full.names = FALSE))
    }
    model <- sub(".stan$", "", basename(MODELS))
  }
  MODEL_HOME <- dirname(MODELS)
  STAN_ENV <- new.env()
  if(file.exists(fp <- file.path(MODEL_HOME, paste0(model, ".data.R")))) {
    source(fp, local = STAN_ENV, verbose = FALSE, echo = TRUE)
  }
  method <- match.arg(method)
  dots <- list(...)
  if(is.null(dots$object)) dots$object <- stan_model(MODELS, model_name = model)
  dots$data <- STAN_ENV
  if(method == "sampling") return(do.call(sampling, args = dots))
  else return(do.call(optimizing, args = dots))
}
