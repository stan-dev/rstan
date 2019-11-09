context("S3 method to cast coef_diff by_group_stan class to array for mcmc plotting")


test_that("Result with default parameters is tbl_df with required columns",{
  skip_on_cran()
  library(breathtestcore)
  library(dplyr)
  data("usz_13c", package = "breathtestcore")

  data = usz_13c %>%
     dplyr::filter( patient_id %in%  c("norm_001", "norm_002", "norm_003",
                           "norm_004", "pat_001", "pat_002","pat_003")) %>%
     cleanup_data()
   fit = stan_group_fit(data, iter = 300, chains = 1)
   cf = coef_diff_by_group(fit)
   mt = as.matrix(cf)
   expect_equal(dim(mt), c(150, 3))
})
