
test_optimizing <- function() {
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
      target += -square(a[1]) - square(a[2]);
    }
    generated quantities {
      matrix[4, 5] B;
      for (i in 1:4) for (j in 1:5) B[i,j] = 1;
    }
  ' 
  # TODO: add check the names of the returns
  m2 <- stan_model(model_code = mc)
  set.seed(1287)
  o2 <- optimizing(m2, hessian = TRUE, seed = 4)
  set.seed(1287)
  o3 <- optimizing(m2, as_vector = FALSE, seed = 4)
  s <- list(a = 1:2, y = 3)
  attr(s$a, "dim") <- 2
  checkEquals(o3$par, rstan:::rstan_relist(o2$par, s))
  checkEquals(o2$par[3], 0, tolerance = 0.1, checkNames = FALSE)
}
