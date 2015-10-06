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

StanProgram$methods(initialize = function(file = file.choose(), code, auto_write = 
                                          rstan_options("auto_write"), ...) {
  "Initializes an object of the StanProgram class"
  if (!missing(code)) {
    tf <- tempfile("stanprogram")
    writeLines(code, con = tf)
    md5 <- tools::md5sum(tf)
    file <- file.path(tempdir(), paste0(md5, ".stan"))
    if(!file.exists(file)) file.rename(from = tf, to = file)
    ret <- rstan::stanc(file)
  }
  else {
    md5 <- tools::md5sum(file)
    ret <- rstan::stanc(file)
  }
  .self$stan_code <<- scan(text = ret$model_code, what = character(), 
                           sep = "\n", quiet = TRUE)
  .self$cpp_code <<- scan(text = ret$cppcode, what = character(), 
                          sep = "\n", quiet = TRUE)
  .self$dso <<- rstan::stan_model(file, auto_write = FALSE, ...)@dso
  if (auto_write) .self$save(file = sub("stan$", "rda", file))
  else .self$save(file = file.path(tempdir(), paste0(md5, ".rda")))
  return(invisible(NULL))
})

StanProgram$methods(show = function() {
  "Shows the Stan code"
  cat(.self$stan_code, sep = "\n")
})

StanProgram$methods(expose = function() {
  "Brings user-defined Stan functions into the R environment
  Useful for unit-testing, posterior simulation, information criteria, etc."
  rstan::testify(writeLines(.self$stan_code, con = tempfile("testify")))
})

StanProgram$methods(instantiate = function(data = list()) {
  "Specify the data to be conditioned on in the Stan program"
  stop("FIXME: Implement")
  # return an instance of StanProgramWithData
})

StanProgram$methods(save = function(file) {
  "Save this StanProgram object to the disk in serialized form"
  saveRDS(.self, file = file)
})

StanProgram$methods(identical = function(program) {
  "Test whether this StanProgram is identical to another"
  if (!is(program, "StanProgram")) return(FALSE)
  return(identical(.self$stan_code, program$stan_code))
})
  
StanProgram$methods(help = help_from_instance)
