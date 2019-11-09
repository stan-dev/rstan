context("Bayesian fit of multiple groups")

test_that("Multiple records per patient return multiple groups (long version)", {
  skip_on_cran() # long
  skip_on_32bit() # nlme fit fails
  library(breathtestcore)
#  library(breathteststan)

  chains = 2
  student_t_df = 3  # Rough student distribution
  dose = 100
  iter = 500
  model = "breath_test_group_1"
  data("usz_13c", package = "breathtestcore")
  set.seed(4711)
  data = usz_13c %>%
    dplyr::filter( patient_id %in%
                     c("norm_001", "norm_002", "norm_003",
                       "pat_001", "pat_002", "pat_003")) %>%
    breathtestcore::cleanup_data()
  # fit nlme for comparison
  fit_nlme = breathtestcore::nlme_fit(data)
  expect_false(is.null(coef(fit_nlme)))
  ## The above fails on 32 bit
  # fit stan_group
  fit = stan_group_fit(data, dose = dose, student_t_df = student_t_df,
                       chains = chains, iter = iter, model = model  )

  cf = coef(fit)
  expect_equal(unique(cf$group), c("liquid_normal", "solid_normal", "solid_patient"))
  expect_gt(sigma(fit), 0.5)

  cf = coef(fit) %>%
    left_join(coef(fit_nlme), by = c("patient_id", "parameter", "method", "group"))  %>%
    filter(method == "exp_beta") %>%
    mutate(rel_diff = 2*abs(value.x - value.y)/(value.x + value.y))   %>%
    select(parameter, rel_diff) %>%
    summarize(
      rel_diff = mean(rel_diff)
    )
  expect_lt(cf$rel_diff, 0.04)
})


