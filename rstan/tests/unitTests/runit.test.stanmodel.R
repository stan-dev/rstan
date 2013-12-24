
test_optimzing <- function() {
  m <- stan_model(model_code = 'parameters {real y;} model {y ~ normal(0,1);}')
  o <- optimizing(m, hessian = TRUE)
  checkEquals(o$par[1], 0, tolerance = 0.1, checkNames = FALSE)
  checkEquals(o$hessian[1,1], -1, tolerance = 0.1, checkNames = FALSE)
}
