
context("PRESENT");

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

rst.sr <- bzCallStan("sr",
                     dat.sub=subgrp.effect,
                     var.estvar = var.estvar,
                     var.cov = var.cov,
                     par.pri = c(B=1000, C=1000),
                     chains=4,
                     iter=40,
                     warmup=20, seed=1000);

rst.nse <- bzCallStan("nse",
                      dat.sub=subgrp.effect,
                      var.estvar = var.estvar,
                      var.cov = var.cov,
                      par.pri = c(B=1000),
                      chains=4,
                      iter=40,
                      warmup=20, seed=1000);

test_that("select subgroup", {
    mus <- array(1:100, dim=c(20,5));
    expect_identical(get.sel.subgrp(mus),       as.numeric(1:5));
    expect_identical(get.sel.subgrp(mus, "A"),  as.numeric(1:5));
    expect_identical(get.sel.subgrp(mus, 1:10), as.numeric(1:5));
    expect_identical(get.sel.subgrp(mus, 7:9),  as.numeric(1:5));
})

test_that("get all mus", {
    mus.1 <- get.all.mus(rst.sr);
    mus.2 <- get.all.mus(rst.sr, sel.grps = 1:2);
    mus.3 <- get.all.mus(rst.sr, sel.grps = 1:2, ref.stan.rst = rst.nse);
    mus.4 <- get.all.mus(rst.sr, sel.grps = 1:2,
                         ref.stan.rst = rst.nse, ref.sel.grps = 3:4);

    expect_equal(nrow(mus.1), 80);
    expect_equal(ncol(mus.1), 8);
    expect_equal(nrow(mus.2), 80);
    expect_equal(ncol(mus.2), 2);
    expect_identical(colnames(mus.2), c("Subgroup 1", "Subgroup 2"));
    expect_equal(nrow(mus.3), 80);
    expect_equal(ncol(mus.3), 3);
    expect_identical(colnames(mus.3), c("Subgroup 1", "Subgroup 2", "No subgroup effect(1)"));
    expect_equal(nrow(mus.4), 80);
    expect_equal(ncol(mus.4), 4);
    expect_identical(colnames(mus.4),
                     c("Subgroup 1", "Subgroup 2", "No subgroup effect(3)", "No subgroup effect(4)"));

})

test_that("beanz comparison", {

    bz.1 <- bzSummaryComp(rst.sr, sel.grps=NULL, cut=0, digits=3);
    bz.2 <- bzSummaryComp(rst.sr, sel.grps=c(1:3), cut=0, digits=3);
    bz.3 <- bzSummaryComp(rst.sr, sel.grps=c(1:3), cut=0.1, digits=3);

    expect_equal(nrow(bz.1), 28);
    expect_equal(nrow(bz.2), 3);
    expect_equal(colnames(bz.3)[9], "ProbLT0.1");
    expect_true(is.data.frame(bz.3));
})


test_that("beanz summary", {

    bz.1 <- bzSummary(rst.sr, sel.grps=NULL, cut=0, digits=3);
    bz.2 <- bzSummary(rst.sr, sel.grps=c(1:3), cut=0.1, digits=3);
    bz.3 <- bzSummary(rst.sr, sel.grps=c(1:3), ref.stan.rst = rst.nse);

    expect_equal(nrow(bz.1), 8);
    expect_equal(ncol(bz.1), 9);
    expect_equal(nrow(bz.2), 3);
    expect_equal(colnames(bz.2)[9], "ProbLT0.1");
    expect_true(is.data.frame(bz.2));
    expect_equal(nrow(bz.3), 4);
})



