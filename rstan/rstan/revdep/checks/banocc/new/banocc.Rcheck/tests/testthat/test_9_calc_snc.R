context("Step 7 - checking get_snc behavior")
load("testthat_objects/sample_stan_data.RData")
load("testthat_objects/sample_stan_fit.RData")
posterior_samples <- rstan::extract(Fit)

test_that("calc_snc only accepts numeric vector input", {
    test_vec <- function(s){
        calc_snc(sample_vec=s)
    }
    err_string <- "must be a numeric vector"
    expect_error(test_vec(list(1, 2, 3)),       err_string)
    expect_error(test_vec(matrix(1:4, ncol=2)), err_string)
    expect_error(test_vec(c("1", "2")),         err_string)
})

test_that("calc_snc returns NA if sd = 0", {
    expect_true(is.na(calc_snc(sample_vec=rep(0, 10))))
})

sample_vecs <- list(3:12, seq(-0.1, 0.1, 0.1), seq(-5, 4.9, 0.1))
snc_values  <- list(0,    1,                   0.59)

test_that("calc_snc returns a scalar", {
    for (i in seq_along(sample_vecs)){
        snc <- calc_snc(sample_vec=sample_vecs[[i]])
        eval(bquote(expect_is(.(snc), "numeric")))
        eval(bquote(expect_equal(length(.(snc)), 1)))
    }
})

test_that("calc_snc returns a value >= 0", {
    for (i in seq_along(sample_vecs)){
        eval(bquote(expect_true(calc_snc(sample_vec=.(sample_vecs[[i]])) >= 0)))
    }
})

test_that("calc_snc returns a value <= 1", {
    for (i in seq_along(sample_vecs)){
        eval(bquote(expect_true(calc_snc(sample_vec=.(sample_vecs[[i]])) <= 1)))
    }
})

test_that("get_snc returns a list", {
    test_is_list <- function(p){
        eval(bquote(expect_is(get_snc(posterior_samples=posterior_samples,
                                      parameter.names=.(p)), "list")))
    }
    test_is_list("m")
    test_is_list("S")
    test_is_list("lp__")
    test_is_list(c("m", "S", "lp__"))
})

test_that("get_snc returns list with parameter names", {
    get_names <- function(idx) sort(names(posterior_samples))[idx]
    test_names <- function(p){
        sort(names(get_snc(posterior_samples=posterior_samples,
                           parameter.names=p)))
    }
    expect_names <- function(idx){
        eval(bquote(expect_equal(test_names(p=.(get_names(idx))),
                                 .(get_names(idx)))))
    }
    idx <- seq_len(length(posterior_samples))
    expect_names(1)
    expect_names(idx)
})

test_that("get_snc elt dimensions are correct", {
    snc <- get_snc(posterior_samples=posterior_samples,
                   parameter.names=c("m", "S", "lp__"))
    expect_equal(dim(snc$lp__),     NULL)
    expect_equal(length(snc$lp__), 1)
    expect_equal(dim(snc$m),       NULL)
    expect_equal(length(snc$m),    Data$P)
    expect_equal(dim(snc$S),       c(Data$P, Data$P))
})

snc <- get_snc(posterior_samples=posterior_samples,
               parameter.names=c("m", "S", "lp__"))
test_that("get_snc matches eltwise calcn for vectors", {
    for (i in seq_len(Data$P)){
        eval(bquote(expect_equal(snc$m[.(i)],
                                 calc_snc(posterior_samples$m[, .(i)]))))
    }
})

test_that("get_snc matches eltwise calcn for matrices", {
    for (i in seq_len(Data$P)){
        for (k in seq_len(Data$P)){
            eval(bquote(
                expect_equal(snc$S[.(i), .(k)],
                             calc_snc(posterior_samples$S[, .(i), .(k)]))
                ))
        }
    }
})

test_that("get_snc matches eltwise calcn for scalars", {
    expect_equal(snc$lp__, calc_snc(as.vector(posterior_samples$lp__)))
})
