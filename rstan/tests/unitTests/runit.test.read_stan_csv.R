test_read_stan_csv <- function() { 
  exfit <- read_stan_csv(dir(system.file('misc', package = 'rstan'), full.names = TRUE))
  checkEquals(exfit@model_pars, c("mu", "sigma", "z", "alpha", "lp__"), checkNames = FALSE)
}
