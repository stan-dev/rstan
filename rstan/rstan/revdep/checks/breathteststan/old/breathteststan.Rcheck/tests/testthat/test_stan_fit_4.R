context("Bayesian fit with Stan")

test_that("Multiple chains return valid results similar to nlme", {
  skip_on_cran()
  skip_on_32bit()
  library(breathtestcore)
  chains = 2
  student_t_df = 10
  dose = 100
  iter = 500
  sample_minutes = 15
  data = cleanup_data(simulate_breathtest_data(n_records = 3, seed = 100)$data)
  fit = stan_fit(data, dose = dose, student_t_df = student_t_df,
                 chains = chains, iter = iter  )
  fit_nlme = nlme_fit(data, dose = dose)
  cf = coef(fit) %>%
    left_join(coef(fit_nlme), by = c("patient_id", "parameter", "method", "group")) %>%
    filter(method == "exp_beta") %>%
    mutate(rel_diff = 2*abs(value.x - value.y)/(value.x + value.y)) %>%
    select(parameter, rel_diff) %>%
    summarize(
      rel_diff = mean(rel_diff)
    )
  expect_lt(cf$rel_diff, 0.005)
})

