
context("STAN");

test_that("call stan", {
    var.cov    <- c("sodium", "lvef", "any.vasodilator.use");
    var.resp   <- "y";
    var.trt    <- "trt";
    var.censor <- "censor";
    resptype   <- "survival";
    var.cov    <- c("sodium", "lvef", "any.vasodilator.use");
    var.estvar <- c("Estimate", "Variance");

    subgrp.effect <- bzGetSubgrpRaw(solvd.sub,
                                    var.resp   = var.resp,
                                    var.trt    = var.trt,
                                    var.cov    = var.cov,
                                    var.censor = var.censor,
                                    resptype   = resptype);

    rst.nse <- bzCallStan("nse", dat.sub=subgrp.effect,
                          var.estvar = var.estvar,
                          var.cov = var.cov,
                          par.pri = c(B=1000),
                          chains=4,
                          iter=40,
                          warmup=20, seed=1000);

    mus <- rst.nse$get.mus();

    expect_error(bzCallStan("nse", dat.sub=subgrp.effect,
                          var.estvar = "error",
                          var.cov = var.cov,
                          par.pri = c(B=1000),
                          chains=4,
                          iter=40,
                          warmup=20, seed=1000),
                 "Variables *");

    expect_error(bzCallStan("nse", dat.sub=subgrp.effect,
                            var.estvar = var.estvar,
                            var.cov = var.cov,
                            par.pri = c(ERR=1000),
                            chains=4,
                            iter=40,
                            warmup=20, seed=1000),
                 "Prior *");

    expect_error(bzCallStan("error", dat.sub=subgrp.effect,
                            var.estvar = var.estvar,
                            var.cov = var.cov,
                            par.pri = c(B=1000),
                            chains=4,
                            iter=40,
                            warmup=20, seed=1000),
                 "arg *");

    expect_is(rst.nse, "beanz.stan");
    expect_is(rst.nse$stan.rst, "stanfit");
    expect_equal(rst.nse$mdl, "No subgroup effect");
    expect_is(rst.nse$smps, "array");
    expect_equal(ncol(rst.nse$smps), 4);
    expect_equal(nrow(rst.nse$smps), 20);
    expect_is(rst.nse$dic, "numeric");
    expect_is(rst.nse$looic, "loo");
    expect_is(rst.nse$rhat, "numeric");
    expect_equal(nrow(mus), 80);
    expect_equal(ncol(mus), 8);
})
