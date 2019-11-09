context("Coeffients by group")


test_that("Result with default parameters is tbl_df with required columns",{
  # This calls coef_by_group.breathtestfit, which include a post-hoc classic test
  # for contrasts.
  library(dplyr)
  library(breathtestcore)
  data("usz_13c", package = "breathtestcore")
  data = usz_13c %>%
    dplyr::filter( patient_id %in%  c("norm_001", "norm_002", "norm_004", "norm_007",
                                      "pat_001", "pat_002","pat_003")) %>%
    cleanup_data()
  comment(data) = "comment"

  fit = stan_fit(data, iter = 300, chains = 1)
  class(fit) = class(fit)[-1] # Remove class breathteststanfit
  cf = coef_by_group(fit) # S3 method
  expect_is(cf, "tbl_df")
  expect_is(cf, "coef_by_group")
  expect_identical(ncol(cf), 7L)
  expect_equal(names(cf), c("parameter", "method", "group", "estimate", "conf.low",
                 "conf.high", "diff_group"))
  expect_equal(comment(data), "comment")
  expect_identical(nrow(cf), 24L)
  expect_identical(sort(unique(cf$diff_group)), c("a", "b", "c"))
  expect_equal(unique(cf$group),
     c("liquid_normal", "solid_normal", "solid_patient"))
})

