test_get_elapsed_time <- function() {
  code <- 'parameters { real y; } model {y ~ normal(0,1); }'
  nchains = 5L
  fit <- stan(model_code = code, chains = nchains, iter = 30)
  checkEquals(dim(get_elapsed_time(fit)), c(nchains, 2L))
  nchains = 1L
  fit2 <- stan(fit = fit, chains = nchains, iter = 30)
  checkEquals(dim(get_elapsed_time(fit2)), c(nchains, 2L))
}
