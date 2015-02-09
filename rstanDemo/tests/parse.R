models <- rstanDemo::stan_demo(0)
stopifnot(all(sapply(models, FUN = function(f) rstan::stanc(f, verbose = TRUE)$status)))
