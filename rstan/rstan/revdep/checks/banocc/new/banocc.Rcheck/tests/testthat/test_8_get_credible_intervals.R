context("Step 6 - checking get_credible_intervals behavior")
load("testthat_objects/sample_stan_data.RData")
load("testthat_objects/sample_stan_fit.RData")
posterior_samples <- rstan::extract(Fit)

test_that("get_credible_intervals accepts only \"marginal.centered\" or \"marginal.hpd\" as type", {
    test_type <- function(t){
        get_credible_intervals(posterior_samples=posterior_samples, type=t)
    }
    err_string <- "'type' must be one of"
    expect_error(test_type("Marginal.centered"), err_string)
    expect_error(test_type("marginal.Centered"), err_string)
    expect_error(test_type("marginal.HPD"),      err_string)
    expect_error(test_type("Marginal.hpd"),      err_string)
    expect_error(test_type("marginal"),          err_string)
    expect_error(test_type("marginal.spin"),     err_string)
})

test_that("get_credible_intervals returns a list", {
    test_is_list <- function(list, t){
        get_credible_intervals(posterior_samples=posterior_samples,
                               list=list,
                               parameter.names=names(posterior_samples),
                               type=t)
    }
    expect_is(test_is_list(FALSE, "marginal.centered"), "list")
    expect_is(test_is_list(FALSE, "marginal.hpd"),      "list")
    expect_is(test_is_list(TRUE, "marginal.centered"),  "list")
    expect_is(test_is_list(TRUE, "marginal.hpd"),       "list")
})

test_that("get_credible_intervals names are parameter names", {
    test_names <- function(list, t, idx){
        sort(names(get_credible_intervals(
            posterior_samples=posterior_samples,
            list=list, type=t,
            parameter.names=sort(names(posterior_samples))[idx])))
    }
    expected_names <- sort(names(posterior_samples))
    idx <- seq_along(posterior_samples)
    expect_equal(test_names(FALSE, "marginal.centered", 1),   expected_names[1])
    expect_equal(test_names(FALSE, "marginal.centered", idx), expected_names)
    expect_equal(test_names(FALSE, "marginal.hpd", 1),        expected_names[1])
    expect_equal(test_names(FALSE, "marginal.hpd", idx),      expected_names)
    expect_equal(test_names(TRUE, "marginal.centered", 1),    expected_names[1])
    expect_equal(test_names(TRUE, "marginal.centered", idx),  expected_names)
    expect_equal(test_names(TRUE, "marginal.hpd", 1),         expected_names[1])
    expect_equal(test_names(TRUE, "marginal.hpd", idx),       expected_names)

})

test_that("get_credible_intervals elts have correct dims when list=FALSE", {
    expect_equal_dims <- function(t){
        ci <- get_credible_intervals(posterior_samples=posterior_samples,
                                     list=FALSE,
                                     parameter.names=names(posterior_samples),
                                     type=t)
        for (p in names(ci)){
            expected_dim <- c(2, dim(posterior_samples[[p]])[-1])
            if (length(expected_dim) > 1){
                eval(bquote(expect_equal(dim(ci[[.(p)]]), .(expected_dim))))
            } else {
                eval(bquote(expect_equal(length(ci[[.(p)]]), .(expected_dim))))
            }
        }
    }
    expect_equal_dims("marginal.centered")
    expect_equal_dims("marginal.hpd")
})

test_that("get_credible_intervals elts have correct dims when list=TRUE", {
    test_dims <- function(t, p.names){
        ci <- get_credible_intervals(posterior_samples=posterior_samples,
                                     list=TRUE, type=t,
                                     parameter.names=p.names)
        lapply(ci, function(k) lapply(k, dim))
    }
    test_lengths <- function(t, p.names){
        ci <- get_credible_intervals(posterior_samples=posterior_samples,
                                     list=TRUE, type=t,
                                     parameter.names=p.names)
        lapply(ci, function(k) lapply(k, length))
    }
    expect_equal_dims <- function(t){
        p.names <- c("lp__", "m", "S")
        expected_dims <- list(list(NULL, NULL),
                              list(NULL, NULL),
                              list(c(Data$P, Data$P), c(Data$P, Data$P)))
        expected_dims <- lapply(expected_dims, function(ed){
            names(ed) <- c("lower", "upper")
            return(ed)
            })
        expected_lens <- list(list(1, 1),
                              list(Data$P, Data$P),
                              list(Data$P^2, Data$P^2))
        expected_lens <- lapply(expected_lens, function(el){
            names(el) <- c("lower", "upper")
            return(el)
        })
        
        ci_dims <- test_dims(t, p.names)
        ci_lens <- test_lengths(t, p.names)

        for (i in seq_along(p.names)){
            eval(bquote(
                expect_equal(ci_dims[[.(p.names[i])]], .(expected_dims[[i]]))
                ))
            eval(bquote(
                expect_equal(ci_lens[[.(p.names[i])]], .(expected_lens[[i]]))
                ))
        }
        
    }
    expect_equal_dims("marginal.centered")
    expect_equal_dims("marginal.hpd")
})
