context("Step 3 - running BAnOCC")
load("testthat_objects/banocc_model_test.RData")
load("testthat_objects/sample_stan_data.RData")


test_that("get_gamma_param gives error if param <= 0", {
    err_string <- "must be > 0"
    expect_error(get_gamma_param(0, "a"),  err_string)
    expect_error(get_gamma_param(-1, "a"), err_string)
})

test_that("get_gamma_param returns param if param > 0", {
    expect_equal(get_gamma_param(pi, "a"), pi)
    expect_equal(get_gamma_param(2, "a"), 2)
})

n  <- rep(0, Data$P)
L  <- 10 * diag(Data$P)

init1    <- list(list(m = n, O=diag(Data$P), lamda=0.02))
init2    <- NULL
init_opt <- list(init1, init2)

use_matrix        <- list(TRUE, FALSE)

opt_idx <- expand.grid(init          = seq_along(init_opt),
                       use_matrix    = seq_along(use_matrix))

all_args <- lapply(seq_len(nrow(opt_idx) * 2), function(i){
    j <- ifelse(i %% 2 == 0, i / 2, (i + 1) / 2)
    idx <- c(opt_idx$init[j],
             opt_idx$use_matrix[j])
    args <- list(init=init_opt[[idx[1]]],
                 use_matrix=use_matrix[[idx[2]]])
    return(args)
})

sw_run_banocc <- function(use_matrix=FALSE, init=NULL){
    sink("banocc.out", type="output")
    if (use_matrix){
        data <- as.matrix(Data$C)
    } else {
        data <- as.data.frame(Data$C)
    }
    rb <- suppressWarnings(run_banocc(
        compiled_banocc_model=compiled_banocc_model, C=data, a=0.5,
        b=0.01, n=n, L=L,
        chains=1, iter=4, warmup=2, init=init,
        verbose=FALSE, num_level=0
        ))
    sink()
    return(rb)
}

on_cran <- !identical(Sys.getenv("NOT_CRAN"), "true")
on_bioc <- !identical(Sys.getenv("BBS_HOME"), "")
if (!on_cran && !on_bioc){
    rb <- lapply(seq_along(all_args), function(i){
        do.call(what=sw_run_banocc, args=all_args[[i]])
    })
}

test_that("run_banocc returns a list", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        expect_is(rb[[i]], "list")
    }
})

test_that("run_banocc takes a matrix", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        if (all_args[[i]]$use_matrix){
            expect_is(rb[[i]], "list")
        }
    }
})

test_that("run_banocc takes a data frame", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        if(!all_args[[i]]$use_matrix){
            expect_is(rb[[i]], "list")
        }
    }
})

rb_names <- c("Data", "Fit")

test_that("run_banocc returns list with correct names", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        expect_equal(sort(names(rb[[i]])), sort(rb_names))
    }
})

test_that("run_banocc Data elt is list", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        expect_is(rb[[i]]$Data, "list")
    }
})

data_names <- c("C", "N", "P", "n", "L", "a", "b")
test_that("run_banocc Data elt has correct names", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        expect_equal(sort(names(rb[[i]]$Data)), sort(data_names))
    }
})

test_that("run_banocc Fit elt is stanfit object", {
    skip_on_cran()
    skip_on_bioc()
    for (i in seq_along(rb)){
        expect_is(rb[[i]]$Fit, "stanfit")
    }
})
