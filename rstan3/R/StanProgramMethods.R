# This file is part of RStan
# Copyright (C) 2015 Jiqiang Guo and Benjamin Goodrich
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

StanProgram$methods(initialize = function(file, code, ...) {
  "Initializes an object of the StanProgram class"
  if (missing(file)) {
    .self$file <- NA_character_
    tf <- tempfile()
    writeLines(code, con = tf)
    file2 <- file.path(tempdir(), paste0(tools::md5sum(tf), ".stan"))
    if(!file.exists(file2)) file.rename(from = tf, to = file2)
    ret <- rstan::stanc(file2, ...)
  }
  else {
    ret <- rstan::stanc(file, ...)
    .self$file <<- file
    file2 <- file
  }
  .self$stan_code <<- scan(text = ret$model_code, what = character(), 
                           sep = "\n", quiet = TRUE)
  .self$cpp_code <<- scan(text = ret$cppcode, what = character(), 
                          sep = "\n", quiet = TRUE)
  .self$dso <<- rstan::stan_model(file2)@dso
  return(invisible(NULL))
})

StanProgram$methods(show = function() {
  "Shows the Stan code"
  cat(.self$stan_code, sep = "\n")
})

StanProgram$methods(expose = function() {
  "Brings user-defined Stan functions into the R environment
  Useful for unit-testing, posterior simulation, information criteria, etc."
  if(is.na(.self$file)) {
    file2 <- paste0(tempfile(), ".stan")
    writeLines(.self$stan_code, con = file2)
  }
  else file2 <- .self$file
  rstan::testify(file2)
})

StanProgram$methods(instantiate = function(data = list()) {
  "Specify the data to be conditioned on in the Stan program"
  stop("FIXME: Implement")
  # return an instance of StanProgramWithData
})

StanProgram$methods(save = function(filepath) {
  "Save this StanProgram object to the disk in serialized form"
  if(is.na(.self$file)) {
    tf <- tempfile()
    writeLines(.self$stan_code, con = tf)
    file2 <- file.path(tempdir(), paste0(tools::md5sum(tf), ".stan"))
  }
  else file2 <- .self$file
  if (missing(filepath)) filepath <- sub("stan$", "rda", file2)
  saveRDS(.self, file = filepath)
})

StanProgram$methods(help = help_from_instance)
