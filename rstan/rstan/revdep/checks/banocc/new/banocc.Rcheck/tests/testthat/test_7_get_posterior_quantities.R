context("Step 5 - testing get_posterior_quantiles")

load("testthat_objects/sample_stan_data.RData")
load("testthat_objects/sample_stan_fit.RData")
posterior_samples <- rstan::extract(Fit)

test_that("get_posterior_quantiles returns a list", {
    test_is_list <- function(probs, list, p.names){
        get_posterior_quantiles(posterior_samples=posterior_samples,
                                probs=probs, list=list,
                                parameter.names=p.names)
    }
    expect_is(test_is_list(0.5, FALSE, "m"),          "list")
    expect_is(test_is_list(0.5, TRUE, "m"),           "list")
    expect_is(test_is_list(0.5, FALSE, c("m", "S")),  "list")
    expect_is(test_is_list(0.5, TRUE, c("m", "S")),   "list")
    expect_is(test_is_list(0, FALSE, "m"),            "list")
    expect_is(test_is_list(0, TRUE, "m"),             "list")
    expect_is(test_is_list(1, FALSE, "m"),            "list")
    expect_is(test_is_list(1, TRUE, "m"),             "list")
    expect_is(test_is_list(0, FALSE, "lp__"),         "list")
    expect_is(test_is_list(0, TRUE,  "lp__"),         "list")
    expect_is(test_is_list(c(0, 0.5, 1), FALSE, "m"), "list")
    expect_is(test_is_list(c(0, 0.5, 1), TRUE, "m"),  "list")
})

test_that("get_posterior_quantiles returns with list parameter names", {
    test_names <- function(probs, list, idx){
        sort(names(get_posterior_quantiles(
            posterior_samples=posterior_samples,
            probs=probs, list=list,
            parameter.names=sort(names(posterior_samples))[idx])))
    }
    expected_names <- sort(names(posterior_samples))
    idx <- seq(1, length(posterior_samples))
    expect_equal(test_names(c(0, 0.5, 1), TRUE, 1),    expected_names[1])
    expect_equal(test_names(c(0, 0.5, 1), FALSE, 1),   expected_names[1])
    expect_equal(test_names(c(0, 0.5, 1), TRUE, idx),  expected_names)
    expect_equal(test_names(c(0, 0.5, 1), FALSE, idx), expected_names)
})

test_that("get_posterior_quantiles elt dimensions are correct when list=FALSE", {
    pq <- get_posterior_quantiles(posterior_samples=posterior_samples,
                                  probs=c(0, 0.5, 1), list=FALSE,
                                  parameter.names=c("m", "S", "lp__"))
    expect_equal(dim(pq$lp__),    NULL)
    expect_equal(length(pq$lp__), 3)
    expect_equal(dim(pq$m),       c(3, Data$P))
    expect_equal(dim(pq$S),       c(3, Data$P, Data$P))
})

test_that("get_posterior_quantiles elt lengths are correct when list=TRUE", {
    pq <- get_posterior_quantiles(posterior_samples=posterior_samples,
                                  probs=c(0, 0.5, 1), list=TRUE,
                                  parameter.names=c("m", "S", "lp__"))
    expect_equal(length(pq$lp__), 3)
    expect_equal(length(pq$m),    3)
    expect_equal(length(pq$S),    3)
})

test_that("get_posterior_quantiles quantile dimensions are correct when list=TRUE", {
    pq <- get_posterior_quantiles(posterior_samples=posterior_samples,
                                  probs=c(0, 0.5, 1), list=TRUE,
                                  parameter.names=c("m", "S", "lp__"))
    expect_equal(dim(pq$lp__[[1]]),    NULL)
    expect_equal(length(pq$lp__[[1]]), 1)
    expect_equal(dim(pq$m[[1]]),       NULL)
    expect_equal(length(pq$m[[1]]),    Data$P)
    expect_equal(dim(pq$S[[1]]),       c(Data$P, Data$P))
})

probs <- c(0, 0.5, 1)
pq <- get_posterior_quantiles(posterior_samples=posterior_samples,
                              probs=probs, list=FALSE,
                              parameter.names=c("lp__", "m", "S"))
test_that("get_posterior_quantiles matches eltwise calcn for vectors when list=FALSE",{
    for (i in seq_len(Data$P)){
        expect_equal(pq$m[, i],
                     quantile(posterior_samples$m[, i], probs=probs))
    }
})
test_that("get_posterior_quantiles matches eltwise calcn for matrices when list=FALSE", {
    for (i in seq_len(Data$P)){
        for (k in seq_len(Data$P)){
            q <- quantile(posterior_samples$S[, i, k], probs=probs)
            expect_equal(pq$S[, i, k], q)
        }
    }
})
test_that("get_posterior_quantiles matches eltwise calcn for scalars when list=FALSE", {
    q <- quantile(posterior_samples$lp__, probs=probs)
    expect_equal(pq$lp__, q)
})


pq <- get_posterior_quantiles(posterior_samples=posterior_samples,
                              probs=probs, list=TRUE,
                              parameter.names=c("lp__", "m", "S"))
test_that("get_posterior_quantiles matches eltwise calcn for vectors when list=TRUE", {
    for (i in seq_len(Data$P)){
        q <- unname(quantile(posterior_samples$m[, i], probs=probs))
        for (k in seq_along(probs)){
            expect_equal(pq$m[[k]][i], q[k])
        }
    }
})
test_that("get_posterior_quantiles matches eltwise calcn for matrices when list=TRUE", {
    for (i in seq_len(Data$P)){
        for (k in seq_len(Data$P)){
            q <- unname(quantile(posterior_samples$S[, i, k],
                                 probs=probs))
            for (l in seq_along(probs)){
                expect_equal(pq$S[[l]][i, k], q[l])
            }
        }
    }
})
test_that("get_posterior_quantiles matches eltwise calcn for scalars when list=TRUE", {
    q <- unname(quantile(posterior_samples$lp__, probs=probs))
    for (i in seq_along(probs)){
        expect_equal(pq$lp__[[i]][1], q[i])
    }
})

context("Step 5 - test get_posterior_estimates")

test_that("get_posterior_estimates returns a list", {
  test_is_list <- function(m, p.names){
      get_posterior_estimates(posterior_samples=posterior_samples,
                              estimate_method=m,
                              parameter.names=p.names)
  }
  expect_is(test_is_list("mean", "m"),                   "list")
  expect_is(test_is_list("mean", "lp__"),                "list")
  expect_is(test_is_list("mean", "S"),                   "list")
  expect_is(test_is_list("mean", c("m", "lp__", "S")),   "list")
  expect_is(test_is_list("median", "m"),                 "list")
  expect_is(test_is_list("median", "lp__"),              "list")
  expect_is(test_is_list("median", "S"),                 "list")
  expect_is(test_is_list("median", c("m", "lp__", "S")), "list")
})

test_that("get_posterior_estimates returns error if estimate_method not right", {
    test_error <- function(m){
        get_posterior_estimates(posterior_samples=posterior_samples,
                                estimate_method=m,
                                parameter.names=c("m", "S", "lp__"))
    }
    err_string <- "[Ii]nvalid estimate_method"
    expect_error(test_error("Mean"),   err_string)
    expect_error(test_error("Median"), err_string)
    expect_error(test_error("gobble"), err_string)
    expect_error(test_error(""),       err_string)
})

test_that("get_parameter_estimates names are parameter names", {
    test_names <- function(m, idx){
        sort(names(get_posterior_estimates(
            posterior_samples=posterior_samples,
            estimate_method=m,
            parameter.names=sort(names(posterior_samples))[idx])))
    }
    expected_names <- sort(names(posterior_samples))
    idx <- seq_along(posterior_samples)
    expect_equal(test_names("mean", 1),     expected_names[1])
    expect_equal(test_names("median", 1),   expected_names[1])
    expect_equal(test_names("mean", idx),   expected_names)
    expect_equal(test_names("median", idx), expected_names)
})

test_that("get_parameter_estimates median elts have dims matching parms", {
    pe <- get_posterior_estimates(posterior_samples=posterior_samples,
                                  estimate_method="median",
                                  parameter.names=c("lp__", "m", "S"))
    expect_equal(dim(pe$lp__),    NULL)
    expect_equal(length(pe$lp__), 1)
    expect_equal(dim(pe$m),       NULL)
    expect_equal(length(pe$m),    Data$P)
    expect_equal(dim(pe$S),       c(Data$P, Data$P))
})

test_that("get_parameter_estimates mean elts have dims matching parms", {
    pe <- get_posterior_estimates(posterior_samples=posterior_samples,
                                  estimate_method="median",
                                  parameter.names=c("lp__", "m", "S"))
    expect_equal(dim(pe$lp__),    NULL)
    expect_equal(length(pe$lp__), 1)
    expect_equal(dim(pe$m),       NULL)
    expect_equal(length(pe$m),    Data$P)
    expect_equal(dim(pe$S),       c(Data$P, Data$P))
})

pe <- get_posterior_estimates(posterior_samples=posterior_samples,
                              estimate_method="median",
                              parameter.names=c("lp__", "m", "S"))
test_that("get_parameter_estimates medians match eltwise for scalars", {
    expect_equal(pe$lp__, median(posterior_samples$lp__))
})
test_that("get_parameter_estimates medians match eltwise for vectors", {
    for (i in seq_len(Data$P)){
        expect_equal(pe$m[i], median(posterior_samples$m[, i]))
    }
})
test_that("get_parameter_estimates medians match eltwise for matrices", {
    for (i in seq_len(Data$P)){
        for (k in seq_len(Data$P)){
            expect_equal(pe$S[i, k], median(posterior_samples$S[, i, k]))
        }
    }
})

pe <- get_posterior_estimates(posterior_samples=posterior_samples,
                              estimate_method="mean",
                              parameter.names=c("lp__", "m", "S"))
test_that("get_parameter_estimates means match eltwise for scalars", {
    expect_equal(pe$lp__, mean(posterior_samples$lp__))
})
test_that("get_parameter_estimates means match eltwise for vectors", {

    for (i in seq_len(Data$P)){
        expect_equal(pe$m[i], mean(posterior_samples$m[, i]))
    }
})
test_that("get_parameter_estimates means match eltwise for matrices", {
    for (i in seq_len(Data$P)){
        for (k in seq_len(Data$P)){
            expect_equal(pe$S[i, k], mean(posterior_samples$S[, i, k]))
        }
    }
})
