# It is easiest to install StanHeaders from CRAN via install.packages("StanHeaders")
# You only need this script if you want to install the develop (or some other) branch of StanHeaders,
# including its submodules. This requires the git2r and devtools packages

install_StanHeaders <- function(rstan_branch = "develop",
                                math_branch = "master",
                                library_branch = "master") {
  path_rstan <- tempfile(pattern = "git2r-")
  git2r::clone("http://github.com/stan-dev/rstan", path_rstan, branch = rstan_branch)
  git2r::clone("http://github.com/stan-dev/stan",
               file.path(path_rstan, "StanHeaders", "inst", "include", "upstream"),
               branch = library_branch)
  git2r::clone("http://github.com/stan-dev/math",
               file.path(path_rstan, "StanHeaders", "inst", "include", "mathlib"),
               branch = math_branch)
  # writeLines(c(".PHONY: static", readLines(file.path(path_rstan, "StanHeaders", "src", "Makevars.win")),
  #              "static: $(OBJECTS)", "\t@mkdir -p ../lib", "\t$(AR) -rs ../lib/libStanHeaders.a $(OBJECTS)"),
  #            con = file.path(path_rstan, "StanHeaders", "src", "Makevars"))
  utils::install.packages(file.path(path_rstan, "StanHeaders"), type = "source",
                          repos = NULL, INSTALL_opts = "--preclean")
}
