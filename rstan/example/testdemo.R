library(rstan)
mfile <- system.file("example-models/basic_distributions/binormal.stan", package = "rstanDemo")
cat(mfile, "\n")

# test directly
stan(file = mfile, test_grad = TRUE, chains = 1)

# test inside a function as rstanDemo would do
fun <- function(...) { 
  STAN_ENV <- new.env()
  sm <- stan_model(file = mfile)
  dots <- list(...)
  dots$object <- sm
  dots$data <- STAN_ENV
  a <- do.call(sampling, args = dots)
}

f <- fun(test_grad = TRUE, chains = 1)

