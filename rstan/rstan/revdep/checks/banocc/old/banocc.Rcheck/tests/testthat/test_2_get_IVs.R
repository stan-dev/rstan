context("Step 2 - obtain initial values")

load("testthat_objects/sample_stan_data.RData")

test_that("get_IVs returns a list", {
    expect_is(get_IVs(chains=3, data=Data), "list")
    expect_is(get_IVs(chains=1, data=Data), "list")
})

test_that("get_IVs returns a list of length chains", {
    expect_equal(length(get_IVs(chains=3, data=Data)), 3)
    expect_equal(length(get_IVs(chains=1, data=Data)), 1)
})

test_that("get_IVs elements have variable names", {
    test_names <- function(chains, data){
        sort(names(get_IVs(chains=chains, data=data)[[1]]))
    }
    expected_names <- sort(c("m", "O", "lambda"))
    expect_equal(test_names(3, Data), expected_names)
    expect_equal(test_names(1, Data), expected_names)
})

test_that("get_IVs gives error for 0 or negative chain numbers", {
    err_string <- "must be > 0"
    expect_error(get_IVs(chains=0, data=Data), err_string)
    expect_error(get_IVs(chains=-1, data=Data), err_string)
})

