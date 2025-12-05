test_that("getting diagnostics works", {
  skip("Backwards compatibility")

  ex_fit <- read_stan_csv("test_fit_diagnostics.csv")
  expect_equal(get_bfmi(ex_fit), 1.17558932146879)
  expect_equal(get_low_bfmi_chains(ex_fit), integer(0))

  expected_divergent <- logical(10)
  expected_divergent[c(2, 3)] <- TRUE
  expect_equal(get_divergent_iterations(ex_fit), expected_divergent)
  expect_equal(get_num_divergent(ex_fit), 2)

  expected_max_treedepth <- logical(10)
  expected_max_treedepth[c(1, 5)] <- TRUE
  expect_equal(get_max_treedepth_iterations(ex_fit), expected_max_treedepth)
  expect_equal(get_num_max_treedepth(ex_fit), 2)

  expect_equal(
    get_num_leapfrog_per_iteration(ex_fit),
    c(1, 3, 5, 5, 7, 5, 3, 1, 1, 7)
  )
})

