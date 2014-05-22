test_read_stan_csv <- function() { 
  exfit <- read_stan_csv(dir(system.file('misc', package = 'rstan'),
                             pattern='rstan_doc_ex_[[:digit:]].csv',
                             full.names = TRUE))
  checkEquals(exfit@model_pars, c("mu", "sigma", "z", "alpha", "lp__"), checkNames = FALSE)
}

test_read_stan_csv_incomplete_csv_0_samples <- function() {
  exfit <- read_stan_csv(dir(system.file('misc', package = 'rstan'),
                             pattern='rstan_doc_ex_incomplete_1.csv',
                             full.names = TRUE))
  checkEquals(exfit@model_pars, c("mu", "sigma", "z", "alpha", "lp__"), checkNames = FALSE)
}

test_read_stan_csv_incomplete_csv_few_samples <- function() {
  exfit <- read_stan_csv(dir(system.file('misc', package = 'rstan'),
                             pattern='rstan_doc_ex_incomplete_1.csv',
                             full.names = TRUE))
  checkEquals(exfit@model_pars, c("mu", "sigma", "z", "alpha", "lp__"), checkNames = FALSE)
}
