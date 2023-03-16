test_that("expose_stan_functions works", {
  funmod <- "
    functions {
      real rtn_real(real x, vector y) {
        return x;
      }

      vector rtn_vector(real x, vector y) {
        return y;
      }
    }
  "
  expose_stan_functions(stanc(model_code = funmod))
  expect_equal(1.5, rtn_real(1.5, 1:4))
  expect_equal(1:4, rtn_vector(1.5, 1:4))
})
