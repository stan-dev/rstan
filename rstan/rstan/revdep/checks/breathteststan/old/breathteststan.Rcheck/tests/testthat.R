library(testthat)
library(breathteststan)
options(warn = 2)
# Only one test per file to avoid hanging 32-bit compile
#test_check("breathteststan", filter = "stan_fit")
Sys.unsetenv("R_TESTS") # https://github.com/r-lib/testthat/issues/603
test_check("breathteststan")
