
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
  set.seed(1287)
  o2 <- optimizing(m2, hessian = TRUE, seed = 4)
  set.seed(1287)
  o3 <- optimizing(m2, as_vector = FALSE, seed = 4)
  s <- list(a = 1:2, y = 3)
  checkEquals(o3$par, relist(o2$par, s))
  checkEquals(o2$par[3], 0, tolerance = 0.1, checkNames = FALSE)
}
