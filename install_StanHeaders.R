# It is easiest to install StanHeaders from CRAN via install.packages("StanHeaders")
# You only need this script if you want to install the develop (or some other) branch of StanHeaders,
# including its submodules. This requires the git2r and devtools packages

path_rstan <- tempfile(pattern = "git2r-")
path_stan_dev <- "/home/rgiordan/Documents/git_repos/stan-dev/"

git2r::clone(file.path(path_stan_dev, "rstan"),
             path_rstan, branch = "add_hessian4")

git2r::clone(file.path(path_stan_dev, "stan"),
             file.path(path_rstan, "StanHeaders", "inst", "include", "upstream"),
             branch = "add_hessians2")

git2r::clone(file.path(path_stan_dev, "math"),
             file.path(path_rstan, "StanHeaders", "inst", "include", "mathlib"),
             branch = "develop")

devtools::install(file.path(path_rstan, "StanHeaders"), args = "--preclean")
