context("Bayesian fit of multiple groups")

test_that("Multiple records per patient return multiple groups (CRAN version)", {
  skip_on_cran()
  skip_on_32bit()
  data("usz_13c", package = "breathtestcore")
  library(dplyr)

  set.seed(4711)
  data = usz_13c %>%
    dplyr::filter( patient_id %in%
                     c("norm_001", "norm_002", "norm_003")) %>%
    breathtestcore::cleanup_data()
  comment(data) = "comment"
  fit = stan_group_fit(data, iter = 300)
  expect_identical(names(fit), c("coef", "data", "stan_fit", "coef_chain"))
  expect_identical(comment(fit), "comment")
  expect_is(fit, "breathteststangroupfit")
})

