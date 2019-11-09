context("Step 5 - evaluate convergence")

load("testthat_objects/sample_stan_data.RData")
load("testthat_objects/sample_stan_fit.RData")

test_that("evaluate_convergence runs successfully", {
    expect_is(suppressWarnings(evaluate_convergence(Fit)), "logical")
})
