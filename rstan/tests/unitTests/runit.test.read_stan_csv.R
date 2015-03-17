test_read_stan_csv <- function() { 
  exfit <- read_stan_csv(dir(system.file('misc', package = 'rstan'),
                             pattern='rstan_doc_ex_[[:digit:]].csv',
                             full.names = TRUE))
  checkEquals(dim(as.array(exfit)), c(100, 4, 10))
  checkEquals(exfit@model_pars, c("mu", "sigma", "z", "alpha", "lp__"), checkNames = FALSE)
  t <- matrix(c(0.005308, 0.003964, 0.004776, 0.003744, 0.010337, 0.002867, 0.004711, 0.004117),
              byrow = TRUE, ncol = 2)
  colnames(t) <- c("warmup", "sample")
  rownames(t) <- paste0("chain:", 1:nrow(t))
  checkEquals(get_elapsed_time(exfit), t)
}

test_read_stan_csv_incomplete_csv_0_samples <- function() {
  exfit <- read_stan_csv(dir(system.file('misc', package = 'rstan'),
                             pattern='rstan_doc_ex_incomplete_1.csv',
                             full.names = TRUE))
  checkEquals(exfit@model_pars, c("mu", "sigma", "z", "alpha", "lp__"), checkNames = FALSE)
  checkTrue(is.null(get_elapsed_time(exfit)))
}

test_read_stan_csv_incomplete_csv_few_samples <- function() {
  exfit <- read_stan_csv(dir(system.file('misc', package = 'rstan'),
                             pattern='rstan_doc_ex_incomplete_1.csv',
                             full.names = TRUE))
  checkEquals(exfit@model_pars, c("mu", "sigma", "z", "alpha", "lp__"), checkNames = FALSE)
}
