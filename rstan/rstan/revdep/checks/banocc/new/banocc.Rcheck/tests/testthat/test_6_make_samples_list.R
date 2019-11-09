context("Step 4 - Making list of posterior samples")

load("testthat_objects/sample_stan_fit.RData")
load("testthat_objects/sample_stan_data.RData")

samples           <- rstan::extract(Fit, permuted=FALSE)
stan_args         <- attr(Fit, 'stan_args')
posterior_samples <- rstan::extract(Fit)
num_samples <- ceiling((stan_args[[1]]$iter - stan_args[[1]]$warmup) /
                       stan_args[[1]]$thin)


test_that("rstan::extract gives array with iter x chains x parms", {
    dims <- c("iterations", "chains", "parameters")
    expect_equal(names(dimnames(samples)), dims)
    expect_equal(dim(samples)[1], num_samples)
    expect_equal(dim(samples)[2], length(stan_args))
})

test_that("make_samples_list returns a list", {
    test_is_list <- function(t, cc){
        make_samples_list(samples=samples, thin=t,
                          concatenate.chains=cc)
    }
    expect_is(test_is_list(1, FALSE), "list")
    expect_is(test_is_list(1, TRUE),  "list")
    expect_is(test_is_list(3, FALSE), "list")
    expect_is(test_is_list(3, TRUE),  "list")
})

test_that("make_samples_list names are equal to those from rstan::extract with permuted=TRUE", {
    expected_names <- sort(names(posterior_samples))
    test_names <- function(t, cc){
        sort(names(make_samples_list(samples=samples,
                                     thin=t, concatenate.chains=cc)))
    }
    expect_equal(test_names(1, FALSE), expected_names)
    expect_equal(test_names(1, TRUE),  expected_names)
    expect_equal(test_names(3, FALSE), expected_names)
    expect_equal(test_names(3, TRUE),  expected_names)
})

test_that("make_samples_list gives warning if thin is greater than number of samples", {
    test_warn <- function(t, cc){
        make_samples_list(samples=samples, thin=t, concatenate.chains=cc)
    }
    warn_string <- "greater than number of samples"
    expect_warning(test_warn(stan_args[[1]]$iter, TRUE),  warn_string)
    expect_warning(test_warn(stan_args[[1]]$iter, FALSE), warn_string)
})

test_that("make_samples_list elements have correct dims when concatenate.chains=FALSE", {
    expected_parm_nums <- lapply(posterior_samples, function(ps){
        if (length(dim(ps)) > 1){
            prod(dim(ps)[-1])
        } else {
            1
        }
        })
    expect_dims_equal <- function(t){
        sl <- suppressWarnings(make_samples_list(samples=samples, thin=t,
                                                 concatenate.chains=FALSE))
        expected_dims <- c(ceiling(num_samples / t), length(stan_args))
        for (p in names(sl)){
            eval(bquote(expect_equal(dim(sl[[.(p)]]),
                         .(c(expected_dims, expected_parm_nums[[p]])))))
        }
    }
    expect_dims_equal(1)
    expect_dims_equal(num_samples + 1)
    expect_dims_equal(floor(num_samples/5))
})

test_that("make_samples_list elements have correct dims when concatenate.chains=TRUE", {
    expected_parm_dims <- lapply(posterior_samples, function(ps) dim(ps)[-1])
    expect_dims_equal <- function(t){
        sl <- suppressWarnings(make_samples_list(samples=samples, thin=t,
                                                 concatenate.chains=TRUE))
        expected_dims <- c(ceiling(num_samples / t) * length(stan_args))
        for (p in names(sl)){
            eval(bquote(expect_equal(dim(sl[[.(p)]]),
                         .(c(expected_dims, expected_parm_dims[[p]])))))
        }
    }
    expect_dims_equal(1)
    expect_dims_equal(num_samples + 1)
    expect_dims_equal(floor(num_samples/5))
})
