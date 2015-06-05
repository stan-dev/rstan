stopifnot(require(rstan)) # temporarily rely on rstan 2.x to implement some things
rstan_options(auto_write = TRUE)
stopifnot(require(methods))
stopifnot(require(stats4))

source("../R/AllClass.R")
source("../R/StanProgramMethods.R")
