test_that("Exposing standalone function works", {
  model_code <-
  '
  functions {
    real standard_normal_rng() {
      return normal_rng(0,1);
    }
  }
  '

  expose_stan_functions(stanc(model_code = model_code),
                        env = environment(), rebuild = T)
  PRNG <- get_rng(seed = 3)
  expect_equal(standard_normal_rng(PRNG), 0.8825979695)
})
