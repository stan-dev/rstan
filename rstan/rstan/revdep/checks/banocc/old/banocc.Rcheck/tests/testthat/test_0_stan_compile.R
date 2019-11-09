context("Step 0 - compiling bayesStanModel")

load("testthat_objects/banocc_model_test.RData")

test_that("banocc_model compiles to give an object of class 'stanmodel'", {
    expect_is(compiled_banocc_model, 'stanmodel')
})
