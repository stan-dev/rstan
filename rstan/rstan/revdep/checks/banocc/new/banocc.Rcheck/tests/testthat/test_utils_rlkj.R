context("Utilities - rlkj")

d_vals <- c(2, 5)#seq(2, 200, 10)
eta_vals <- c(1, 50) #seq(1, 200, 10)

test_that("rlkj returns a matrix", {
    for (d in d_vals){
        for (eta in eta_vals){
            eval(bquote(expect_is(rlkj(d=.(d), eta=.(eta)), "matrix")))
            eval(bquote(
                expect_is(rlkj(d=.(d), eta=.(eta), cholesky=TRUE), "matrix")
                ))
        }
    }
})

test_that("rlkj returns square matrix", {
    test_square <- function(d, e, ch){
        r <- rlkj(d=d, eta=e, cholesky=ch)
        expect_equal(nrow(r), ncol(r))
    }
    for (d in d_vals){
        for (eta in eta_vals){
            test_square(d, eta, FALSE)
            test_square(d, eta, TRUE)
        }
    }
})

test_that("rlkj returns matrix with d rows", {
    test_nrow <- function(d, e, ch){
        eval(bquote(
            expect_equal(nrow(rlkj(d=.(d), eta=.(e), cholesky=.(ch))), .(d))))
    }
    for (d in d_vals){
        for (eta in eta_vals){
            test_nrow(d, eta, TRUE)
            test_nrow(d, eta, FALSE)
        }
    }
})

test_that("rlkj returns positive definite matrix if cholesky=FALSE", {
    test_pd <- function(d, e){
        eval(bquote(expect_true(all(eigen(rlkj(d=.(d), eta=.(e)))$values > 0))))
    }
    for (d in d_vals){
        for (eta in eta_vals){
            test_pd(d, eta)
        }
    }
})

test_that("rlkj returns symmetric matrix if cholesky=FALSE", {
    get_tri <- function(mat, type){
        if (type=="lower"){
            t(mat)[upper.tri(mat)]
        } else {
            mat[upper.tri(mat)]
        }
    }
    test_symmetric <- function(d, e){
        l <- rlkj(d=d, eta=e)
        expect_equal(get_tri(l, "lower"),
                     get_tri(l, "upper"))
    }
    for (d in d_vals){
        for (eta in eta_vals){
            test_symmetric(d, eta)
        }
    }
})

test_that("rlkj returns matrix with diagonal elts 1 if cholesky=FALSE", {
    test_diag <- function(d, e){
        eval(bquote(
            expect_equal(diag(rlkj(d=.(d), eta=.(e))), rep(1, .(d)))
            ))
    }
    for (d in d_vals){
        for (eta in eta_vals){
            test_diag(d, eta)
        }
    }
})

test_that("rlkj returns matrix with elts between -1 and 1 if cholesky=FALSE", {
    test_elts <- function(d, e){
        r <- rlkj(d=d, eta=e)
        r[upper.tri(r)]
    }
    for (d in d_vals){
        for (eta in eta_vals){
            expect_true(all(test_elts(d, eta) >= -1))
            expect_true(all(test_elts(d, eta) <= 1))
        }
    }
})

test_that("rlkj returns matrix with upper triangle of 0 if cholesky=TRUE", {
    get_uppertri <- function(mat){
        mat[upper.tri(mat)]
    }
    test_uppertri <- function(d, e){
        eval(bquote(
            expect_equal(get_uppertri(rlkj(d=.(d), eta=.(e), cholesky=TRUE)),
                         rep(0, choose(.(d), 2)))
            ))
    }
    for (d in d_vals){
        for (eta in eta_vals){
            test_uppertri(d, eta)
        }
    }
})

test_that("rlkj corr matrix is positive definite if cholesky=TRUE", {
    test_pd <- function(d, e){
        l <- rlkj(d=d, eta=e, cholesky=TRUE)
        eigen(l %*% t(l))$values
    }
    for (d in d_vals){
        for (eta in eta_vals){
            expect_true(all(test_pd(d, eta) > 0))
        }
    }
})

test_that("rlkj corr matrix is symmetric if cholesky=TRUE", {
    get_tri <- function(mat, type){
        if (type=="lower"){
            t(mat)[upper.tri(mat)]
        } else {
            mat[upper.tri(mat)]
        }
    }
    test_symmetric <- function(d, e){
        l <- rlkj(d=d, eta=e, cholesky=TRUE)
        r <- l %*% t(l)
        expect_equal(get_tri(r, "lower"), get_tri(r, "upper"))
    }
    for (d in d_vals){
        for (eta in eta_vals){
            test_symmetric(d, eta)
        }
    }
})

test_that("rlkj corr matrix has diagonal elements of 1 if cholesky=TRUE", {
    test_diag <- function(d, e){
        l <- rlkj(d=d, eta=e, cholesky=TRUE)
        diag(l %*% t(l))
    }
    for (d in d_vals){
        for (eta in eta_vals){
            expect_equal(test_diag(d, eta), rep(1, d))
        }
    }
})

test_that("rlkj corr matrix has elts between -1 and 1 if cholesky=TRUE", {
    test_elts <- function(d, e){
        l <- rlkj(d=d, eta=e, cholesky=TRUE)
        (l %*% t(l))[upper.tri(l)]
    }
    for (d in d_vals){
        for (eta in eta_vals){
            expect_true(all(test_elts(d, eta) >= -1))
            expect_true(all(test_elts(d, eta) <= 1))
        }
    }
})

test_that("rlkj gives error if eta < 1", {
    err_string <- "must be >= 1"
    for (d in d_vals){
        for (eta in c(-1, 0)){
            expect_error(rlkj(d=d, eta=eta, cholesky=TRUE),  err_string)
            expect_error(rlkj(d=d, eta=eta, cholesky=FALSE), err_string)
        }
    }
})

test_that("rlkj gives error if d < 2", {
    err_string <- "must be >= 2"
    for (d in c(-1, 0, 1)){
        for (eta in eta_vals){
            expect_error(rlkj(d=d, eta=eta, cholesky=TRUE),  err_string)
            expect_error(rlkj(d=d, eta=eta, cholesky=FALSE), err_string)
        }
    }
})
