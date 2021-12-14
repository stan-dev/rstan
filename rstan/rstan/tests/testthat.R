library(testthat)
library(rstan)

# Have to test expose_stan_functions separately
test_check("rstan",filter = "^((?!expose_stan_functions).)*$",perl=T)
test_file("testthat/test-expose_stan_functions.R")
