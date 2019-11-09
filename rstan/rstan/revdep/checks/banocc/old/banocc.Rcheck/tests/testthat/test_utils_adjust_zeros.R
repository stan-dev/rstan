context("Utilities - adjust_zeros")

sample_data <- matrix(c(0, 0.8, 0.1, 0.1,
                        0.5, 0, 0, 0.5,
                        0.25, 0.25, 0.25, 0.25),
                      nrow=3, byrow=TRUE)

test_that("adjust_zeros returns a matrix or data frame", {
    test_return <- function(d, r){
        adjust_zeros(data=d, renorm=r)
    }
    expect_is(test_return(sample_data, TRUE),  "matrix")
    expect_is(test_return(sample_data, FALSE), "matrix")
    expect_is(test_return(as.data.frame(sample_data), TRUE),  "data.frame")
    expect_is(test_return(as.data.frame(sample_data), FALSE), "data.frame")
})

test_that("adjust_zeros replaces zero values when renorm=FALSE", {
    test_zeros <- function(z){
        adjust_zeros(data=sample_data, renorm=FALSE, zero_adj=z)
    }
    expect_zero_value <- function(z){
        d <- test_zeros(z)
        d_adj <- d[sample_data==0]
        eval(bquote(expect_equal(.(d_adj), rep(.(z), length(d_adj)))))
    }
    expect_zero_value(0)
    expect_zero_value(0.5)
    expect_zero_value(1e-4)
})

test_that("adjust_zeros gives row sums of 1 when renorm=TRUE", {
    test_rowsums <- function(z){
        apply(adjust_zeros(data=sample_data, renorm=TRUE, zero_adj=z), 1, sum)
    }
    expected_sums <- rep(1, nrow(sample_data))
    expect_equal(test_rowsums(0),    expected_sums)
    expect_equal(test_rowsums(0.5),  expected_sums)
    expect_equal(test_rowsums(1e-4), expected_sums)
})

test_that("adjust_zeros replaces zero values when renorm=TRUE", {
    test_zeros <- function(z){
        adjust_zeros(data=sample_data, renorm=TRUE, zero_adj=z)
    }
    expect_zero_value <- function(z){
        d <- test_zeros(z)
        expected_zero_vals <- unlist(apply(sample_data, 1, function(r){
            nz <- sum(r==0)
            rep(z/(1+nz*z), nz)
        }))
        eval(bquote(expect_equal(d[sample_data==0], .(expected_zero_vals))))
    }
    expect_zero_value(0)
    expect_zero_value(0.5)
    expect_zero_value(1e-4)
})
