# It is easiest to install StanHeaders from CRAN via install.packages("StanHeaders")
# You only need this script if you want to install the develop (or some other) branch of StanHeaders,
# including its submodules. This requires the git2r and devtools packages

path_rstan <- tempfile(pattern = "git2r-")
git2r::clone("http://github.com/stan-dev/rstan", path_rstan, branch = "develop")
git2r::clone("http://github.com/stan-dev/stan", 
             file.path(path_rstan, "StanHeaders", "inst", "include", "upstream"), 
             branch = "develop") # may want to change this branch
git2r::clone("http://github.com/stan-dev/math", 
             file.path(path_rstan, "StanHeaders", "inst", "include", "mathlib"), 
             branch = "develop") # may want to change this branch
devtools::install(file.path(path_rstan, "StanHeaders"), args = "--preclean")
