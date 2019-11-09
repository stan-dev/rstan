context("Bayesian fit with Stan")

test_that("A single record can be fitted", {
  skip_on_cran()
  library(breathtestcore)
  chains = 1
  student_t_df = 10
  dose = 100
  iter = 100
  sample_minutes = 15
  data = cleanup_data(simulate_breathtest_data(n_records = 1, seed = 100)$data)
  fit = stan_fit(data, dose = dose, student_t_df = student_t_df,
                 chains = chains, iter = iter  )
  expect_is(fit, "breathtestfit")
})

