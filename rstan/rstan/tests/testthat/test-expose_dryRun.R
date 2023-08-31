test_that("expose_functions works with dryRun", {
  funcode <- "
    functions {
      real rtn_real(real x) {
        return x;
      }
    }"

  expect_no_error(expose_stan_functions(stanc(model_code = funcode), dryRun = TRUE))
})
