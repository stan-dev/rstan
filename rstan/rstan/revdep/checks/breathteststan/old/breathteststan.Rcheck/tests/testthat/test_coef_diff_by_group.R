context("Coeffient differences by group")

test_that("Credible intervals are returned as coef_diff_by_group coefficients",{
  skip_on_cran()
  library(dplyr)
  library(breathtestcore)
  data("usz_13c", package = "breathtestcore")
  data = usz_13c %>%
    dplyr::filter( patient_id %in%
     c("norm_001", "norm_002", "norm_003", "norm_004", "pat_001", "pat_002","pat_003")) %>%
    cleanup_data()
  comment(data) = "comment"
  set.seed(4711)
  iter = 300
  chains = 1
  fit = stan_group_fit(data, iter = iter, chains = chains)
  expect_identical(names(fit), c("coef", "data", "stan_fit", "coef_chain"))
  cf_diff = coef_diff_by_group(fit)
  expect_is(cf_diff, "coef_diff_by_group_stan")
  expect_identical(nrow(cf_diff), 24L)
  expect_identical(names(cf_diff), c("parameter", "method", "groups",
                                     "estimate", "cred.low", "cred.high"))
  chain = attr(cf_diff, "chain")
  expect_is(chain, "data.frame")
#  expect_equal(nrow(chain), iter*chains*12)
  expect_identical(names(chain),
          c("key","value1","value2","group1","group2","diff"))
})

