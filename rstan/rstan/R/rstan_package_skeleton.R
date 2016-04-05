# This file is part of RStan
# Copyright (C) 2015, 2016 Jiqiang Guo and Benjamin Goodrich
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

rstan.package.skeleton <- function(name = "anRpackage", list = character(),
                                   environment = .GlobalEnv,
                                   path = ".", force = FALSE,
                                   code_files = character(),
                                   stan_files = character()) {
  
  if (length(stan_files) > 0 && !all(grepl("\\.stan$", stan_files))) {
    stop("all files named in 'stan_files' must end with a .stan extension")
  }
  
  mc <- match.call()
  mc$stan_files <- NULL
  mc[[1]] <- quote(utils::package.skeleton)
  eval(mc)
  
  if (R.version$major < 3 || 
      (R.version$major == 3 && R.version$minor < 2.2)) {
        
    warning("rstan.package.skeleton is only fully operational with R >= 3.2.2",
            "Follow the package skeleton of the rstanarm package on GitHub")
    return(invisible(NULL))
  }
  
  DIR <- file.path(path, name)
  download.file("https://raw.githubusercontent.com/stan-dev/rstanarm/master/cleanup",
                destfile = file.path(DIR, "cleanup"), quiet = TRUE)
  download.file("https://raw.githubusercontent.com/stan-dev/rstanarm/master/cleanup.win",
                destfile = file.path(DIR, "cleanup.win"), quiet = TRUE)
  cat("cleanup*", file = file.path(DIR, ".Rbuildignore"), sep = "\n")

  TOOLS <- file.path(DIR, "tools")
  dir.create(TOOLS)
  download.file("https://raw.githubusercontent.com/stan-dev/rstanarm/master/tools/make_cpp.R",
                destfile = file.path(TOOLS, "make_cpp.R"), quiet = TRUE)
  EXEC <- file.path(DIR, "exec")
  dir.create(EXEC)
  file.copy(stan_files, EXEC)
  
  INST <- file.path(DIR, "inst")
  dir.create(INST)
  CHUNKS <- file.path(DIR, "inst", "chunks")
  dir.create(CHUNKS)
  file.create(file.path(CHUNKS, "common_functions.stan"))
  file.create(file.path(CHUNKS, "license.stan"))
  
  SRC <- file.path(DIR, "src")
  dir.create(SRC, showWarnings = FALSE)
  download.file("https://raw.githubusercontent.com/stan-dev/rstanarm/master/src/Makevars",
                destfile = file.path(SRC, "Makevars"), quiet = TRUE)
  
  R <- file.path(DIR, "R")
  dir.create(R, showWarnings = FALSE)
  download.file("https://raw.githubusercontent.com/stan-dev/rstanarm/master/R/stanmodels.R",
                destfile = file.path(R, "stanmodels.R"), quiet = TRUE)
  cat('.onLoad <- function(libname, pkgname) { Rcpp::loadRcppModules() }', 
      file = file.path(R, "zzz.R"), sep = "\n", append = TRUE)

  if (length(stan_files) == 0) module_names <- "NAME"
  else module_names <- paste0("stan_fit4", sub("\\.stan$", "", basename(stan_files)), "_mod")
  cat("Depends: R (>= 3.0.2), Rcpp (>= 0.11.0)", 
      "Imports: rstan (>= 2.8.1)",
      "LinkingTo: StanHeaders (>= 2.9.0), rstan (>= 2.9.0), BH (>= 1.60.0), Rcpp (>= 0.11.0), RcppEigen",
      file = file.path(DIR, "DESCRIPTION"), sep = "\n", append = TRUE)
  cat("RcppModules: ", paste(module_names, collapse = ", "), "\n",
      file = file.path(DIR, "DESCRIPTION"), append = TRUE)
  cat("\n Stan specific notes:",
      "If you add any additional .stan files to the exec/ directory, ",
      "be sure to add an entry in the RcppModules: line of DESCRIPTION.",
      "Also, be sure to add useDynLib(mypackage, .registration = TRUE) to NAMESPACE.",
      "You can put into inst/chunks/common_functions.stan any function that is needed by any .stan file, ",
      "in which case any .stan file can have #include 'common_functions.stan' in its functions block.",
      "The precompiled stanmodel objects will appear in a named list called 'stanmodels'.",
      "The 'cleanup' and 'cleanup.win' scripts in the root of the directory must be made executable.",
      file = file.path(DIR, "Read-and-delete-me"), sep = "\n", append = TRUE)
  return(invisible(NULL))
}
