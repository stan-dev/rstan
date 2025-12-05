test_that("long error messages work", {
  skip("Backwards compatibility")

  code <- "parameters { cov_matrix[3] y; }
           model { y ~ normal(0, 1); }"
  expect_error({r <- stanc(model_code = code)})
  warning_length <- getOption("warning.length")
  emsg <- geterrmessage()
  expect_true(nchar(emsg) < 1000)
  expect_equal(warning_length, getOption("warning.length"))
})
