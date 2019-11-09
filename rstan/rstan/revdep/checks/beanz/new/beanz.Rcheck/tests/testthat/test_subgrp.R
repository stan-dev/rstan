context("Subgroup")

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

test_that("get subgroups", {
    sg.1 <- r.get.subgroup(solvd.sub, var.cov);
    sg.2 <- r.get.subgroup(solvd.sub, var.cov[1:2]);

    expect_equal(nrow(sg.1$subgrp), 8);
    expect_equal(ncol(sg.1$subgrp), 3);
    expect_equal(nrow(sg.2$subgrp), 4);
    expect_equal(ncol(sg.2$subgrp), 2);
    expect_equal(sg.1$xgrp[15],     6);
    expect_equal(sg.2$xgrp[15],     3);
})

test_that("raw subgroup effect", {
    expect_error(bzGetSubgrpRaw(solvd.sub,
                                    var.resp   = "Error",
                                    var.trt    = var.trt,
                                    var.cov    = var.cov,
                                    var.censor = var.censor,
                                    resptype   = resptype),
                 "Variables *");

    expect_is(subgrp.effect, "data.frame");
    expect_equal(sum(subgrp.effect[,"N"]), nrow(solvd.sub));
    expect_equal(nrow(subgrp.effect), 8);
})

test_that("subgroup effect", {
    bz1 <- bzGetSubgrp(subgrp.effect,
                       var.ey = "Estimate",
                       var.variance = "Variance",
                       var.cov = var.cov);
    bz2 <- bzGetSubgrp(subgrp.effect,
                       var.ey = "Estimate",
                       var.variance = "Variance",
                       var.cov = var.cov[1:2]);

    expect_error(bzGetSubgrp(subgrp.effect,
                             var.ey = "Error",
                             var.variance = "Variance",
                             var.cov = var.cov),
                 "Variables *");

    expect_is(bz1, "data.frame");
    expect_equal(nrow(bz1), 8);
    expect_equal(nrow(bz2), 4);
})
