test_read_stan_csv <- function() { 
  exfit <- read_stan_csv(dir(system.file('misc', package = 'rstan'),
                             'rstan_doc_ex_?.csv',
                             full.names = TRUE))
  checkEquals(exfit@model_pars, c("mu", "sigma", "z", "alpha", "lp__"), checkNames = FALSE)
}

test_read_stan_csv_incomplete_csv() {
  exfit <- read_stan_csv(dir(system.file('misc', package = 'rstan'),
                             'rstan_doc_ex_incomplete.csv',
                             full.names = TRUE))
  checkEquals(exfit@model_pars, c("mu", "sigma", "z", "alpha", "lp__"), checkNames = FALSE)
}
