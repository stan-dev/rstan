context("Step 4 - getting output from banocc")

load("testthat_objects/sample_stan_data.RData")
load("testthat_objects/sample_stan_fit.RData")

conf_alpha_opt    <- list(0.05, 0.5)
get_min_width_opt <- list(TRUE, FALSE)
calc_snc_opt      <- list(TRUE, FALSE)
eval_convergence  <- list(TRUE, FALSE)
use_stanfit       <- list(TRUE, FALSE)

opt_idx <- expand.grid(conf_alpha    = seq_along(conf_alpha_opt),
                       get_min_width = seq_along(get_min_width_opt),
                       calc_snc      = seq_along(calc_snc_opt),
                       eval_convergence = seq_along(eval_convergence),
                       use_stanfit = seq_along(use_stanfit))

all_args <- lapply(seq_len(nrow(opt_idx) * 2), function(i){
    j <- ifelse(i %% 2 == 0, i / 2, (i + 1) / 2)
    idx <- c(opt_idx$conf_alpha[j],
             opt_idx$get_min_width[j],
             opt_idx$calc_snc[j],
             opt_idx$eval_convergence[j],
             opt_idx$use_stanfit[j])
    args <- list(conf_alpha=conf_alpha_opt[[idx[1]]],
                 get_min_width=get_min_width_opt[[idx[2]]],
                 calc_snc=calc_snc_opt[[idx[3]]],
                 eval_convergence=eval_convergence[[idx[4]]],
                 use_stanfit=use_stanfit[[idx[5]]])
    return(args)
})

sw_get_b_output <- function(conf_alpha, get_min_width, calc_snc,
                            eval_convergence, use_stanfit=FALSE){
    if (use_stanfit){
        banoccfit <- Fit
    } else {
        banoccfit <- list(Fit=Fit, Data=Data)
    }
    gbo <- suppressWarnings(get_banocc_output(
        banoccfit=banoccfit, conf_alpha=conf_alpha,
        get_min_width=get_min_width, calc_snc=calc_snc,
        eval_convergence=eval_convergence, verbose=FALSE, num_level=0
        ))
    return(gbo)
}

gbo <- lapply(seq_along(all_args), function(i){
    do.call(what=sw_get_b_output, args=all_args[[i]])
})

test_that("get_banocc_output returns a list", {
    for (i in seq_along(gbo)){
        expect_is(gbo[[i]], "list")
    }
})

test_that("get_banocc_output takes a stanfit object", {
    for (i in seq_along(gbo)){
        if (all_args[[i]]$use_stanfit){
            expect_is(gbo[[i]], "list")
        }
    }
})

test_that("get_banocc_output takes a list", {
    for (i in seq_along(gbo)){
        if(!all_args[[i]]$use_stanfit){
            expect_is(gbo[[i]], "list")
        }
    }
})

gbo_names <- c("Fit", "CI.hpd", "Estimates.median")
extra_names <- c("Min.width", "SNC", "Data")

test_that("get_banocc_output returns list with correct names", {
    for (i in seq_along(gbo)){
        names_i <- gbo_names
        if (all_args[[i]]$get_min_width) names_i <- c(names_i, extra_names[1])
        if (all_args[[i]]$calc_snc) names_i <- c(names_i, extra_names[2])
        if (!all_args[[i]]$use_stanfit) names_i <- c(names_i, extra_names[3])
        expect_equal(sort(names(gbo[[i]])), sort(names_i))
    }
})

test_that("get_banocc_output Data elt is list or NULL", {
    for (i in seq_along(gbo)){
        if (!all_args[[i]]$use_stanfit){
            expect_is(gbo[[i]]$Data, "list")
        } else {
            expect_is(gbo[[i]]$Data, "NULL")
        }
    }
})

data_names <- c("C", "N", "P", "n", "L", "a", "b")
test_that("get_banocc_output Data elt has correct names", {
    for (i in seq_along(gbo)){
        if (!all_args[[i]]$use_stanfit){
            expect_equal(sort(names(gbo[[i]]$Data)), sort(data_names))
        }
    }
})

test_that("get_banocc_output Fit elt is stanfit object", {
    for (i in seq_along(gbo)){
        expect_is(gbo[[i]]$Fit, "stanfit")
    }
})

test_that("get_banocc_output CI.hpd elt is list with correct names", {
    ci_names <- c("lower", "upper")
    for (i in seq_along(gbo)){
        expect_is(gbo[[i]]$CI.hpd, "list")
        expect_equal(sort(names(gbo[[i]]$CI.hpd)), sort(ci_names))
    }
})

p <- ncol(Data$C)
test_that("get_banocc_output CI.hpd elt elements are pxp matrices", {
    for (i in seq_along(gbo)){
        for (k in seq_along(gbo[[i]]$CI.hpd)){
            expect_is(gbo[[i]]$CI.hpd[[k]], "matrix")
            expect_equal(dim(gbo[[i]]$CI.hpd[[k]]), c(p, p))
        }
    }
})

test_that("get_banocc_output CI.hpd elt elements have col and row names", {
    for (i in seq_along(gbo)){
        for (k in seq_along(gbo[[i]]$CI.hpd)){
            if (!all_args[[i]]$use_stanfit){
                ci <- gbo[[i]]$CI.hpd[[k]]
                expect_equal(colnames(ci), names(Data$C))
                expect_equal(rownames(ci), names(Data$C))
            }
        }
    }
})

test_that("get_banocc_output CI.hpd elt elements are between -1 and 1", {
    for (i in seq_along(gbo)){
        for (k in seq_along(gbo[[i]]$CI.hpd)){
            if (!all_args[[i]]$eval_convergence){
                expect_true(all(gbo[[i]]$CI.hpd[[k]] - 1 <= 1e-12))
                expect_true(all(gbo[[i]]$CI.hpd[[k]] + 1 >= 1e-12))
            } else {
                expect_true(all(is.na(gbo[[i]]$CI.hpd[[k]])))
            }
        }
    }
})

test_that("get_banocc_output Estimates.median elt is a pxp matrix", {
    for (i in seq_along(gbo)){
        expect_is(gbo[[i]]$Estimates.median, "matrix")
        expect_equal(dim(gbo[[i]]$Estimates.median), c(p, p))
    }
})

test_that("get_banocc_output Estimates.median has col and row names", {
    for (i in seq_along(gbo)){
        if (!all_args[[i]]$use_stanfit){
            est.med <- gbo[[i]]$Estimates.median
            expect_equal(colnames(est.med), names(Data$C))
            expect_equal(rownames(est.med), names(Data$C))
        }
    }
})

test_that("get_banocc_output Estimates.median elts are between -1 and 1", {
    for (i in seq_along(gbo)){
        est.med <- gbo[[i]]$Estimates.median
        if (!all_args[[i]]$eval_convergence){
            expect_true(all(est.med + 1 >= 1e-12))
            expect_true(all(est.med - 1 <= 1e-12))
        } else {
            expect_true(all(is.na(est.med)))
        }
    }
})

test_that("get_banocc_output Min.width elt is a pxp matrix", {
    for (i in seq_along(gbo)){
        if (all_args[[i]]$get_min_width){
            expect_is(gbo[[i]]$Min.width, "matrix")
            expect_equal(dim(gbo[[i]]$Min.width), c(p, p))
        }
    }
})

test_that("get_banocc_output Min.width has col and row names", {
    for (i in seq_along(gbo)){
        if (all_args[[i]]$get_min_width && !all_args[[i]]$use_stanfit){
            mw <- gbo[[i]]$Min.width
            expect_equal(colnames(mw), names(Data$C))
            expect_equal(rownames(mw), names(Data$C))
        }
    }
})

test_that("get_banocc_output Min.width elts are between 0 and 1", {
    for (i in seq_along(gbo)){
        if (all_args[[i]]$get_min_width){
            mw <- gbo[[i]]$Min.width
            if (!all_args[[i]]$eval_convergence){
                expect_true(all(mw - 1<= 1e-12))
                expect_true(all(mw    >= 1e-12))
            } else {
                expect_true(all(is.na(mw)))
            }
        }
    }
})

test_that("get_banocc_output SNC elt is a pxp matrix", {
    for (i in seq_along(gbo)){
        if (all_args[[i]]$calc_snc){
            expect_is(gbo[[i]]$SNC, "matrix")
            expect_equal(dim(gbo[[i]]$SNC), c(p, p))
        }
    }
})

test_that("get_banocc_output SNC has col and row names", {
    for (i in seq_along(gbo)){
        if (all_args[[i]]$calc_snc && !all_args[[i]]$use_stanfit){
            snc <- gbo[[i]]$SNC
            expect_equal(colnames(snc), names(Data$C))
            expect_equal(rownames(snc), names(Data$C))
        }
    }
})

test_that("get_banocc_output SNC elts are between 0 and 1", {
    for (i in seq_along(gbo)){
        if (all_args[[i]]$calc_snc){
            snc <- gbo[[i]]$SNC
            if (!all_args[[i]]$eval_convergence){
                if (any(!is.na(snc))){
                    snc <- na.omit(as.vector(snc))
                    expect_true(all(snc - 1 <= 1e-18))
                    expect_true(all(snc    >= -1e-18))
                }
            } else {
                expect_true(all(is.na(snc)))
            }
        }
    }
})

## Need to test that banocc takes: 1) A data frame or 2) a matrix with 3) compositional rows
## Need to test that C can have zeros
