
test_optimzing <- function() {
  m <- stan_model(model_code = 'parameters {real y;} model {y ~ normal(0,1);}')
  o <- optimizing(m, hessian = TRUE)
  checkEquals(o$par[1], 0, tolerance = 0.1, checkNames = FALSE)
  checkEquals(o$hessian[1,1], -1, tolerance = 0.1, checkNames = FALSE)

  mc <- '
    parameters {
      simplex[2] a; 
      real y;
    } 
    model {
      y ~ normal(0,1); 
      increment_log_prob(-square(a[1]) - square(a[2]));
    }
  ' 
  m2 <- stan_model(model_code = mc)
  o2 <- optimizing(m2, hessian = TRUE)
  checkEquals(o2$par[3], 0, tolerance = 0.1, checkNames = FALSE)
}
