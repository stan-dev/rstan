# It is easiest to install StanHeaders from CRAN via install.packages("StanHeaders")
# You only need this script if you want to install the develop (or some other) branch of StanHeaders,
# including its submodules. This requires the git2r and devtools packages

install_StanHeaders <- function(branch = "develop") {
  path_rstan <- tempfile(pattern = "git2r-")

  git2r::clone("https://github.com/stan-dev/rstan", path_rstan, branch = branch)

  on.exit(setwd(getwd()))
  setwd(path_rstan)

  try(system("sh sh_b.sh --no-build-vignettes --no-manual"))
}
