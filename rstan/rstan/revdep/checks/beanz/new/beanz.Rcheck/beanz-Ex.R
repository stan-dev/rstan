pkgname <- "beanz"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('beanz')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("bzCallStan")
### * bzCallStan

flush(stderr()); flush(stdout())

### Name: bzCallStan
### Title: Call STAN models
### Aliases: bzCallStan

### ** Examples

## Not run: 
##D var.cov    <- c("sodium", "lvef", "any.vasodilator.use");
##D var.resp   <- "y";
##D var.trt    <- "trt";
##D var.censor <- "censor";
##D resptype   <- "survival";
##D var.estvar <- c("Estimate", "Variance");
##D 
##D subgrp.effect <- bzGetSubgrpRaw(solvd.sub,
##D                                   var.resp   = var.resp,
##D                                   var.trt    = var.trt,
##D                                   var.cov    = var.cov,
##D                                   var.censor = var.censor,
##D                                   resptype   = resptype);
##D 
##D rst.nse    <- bzCallStan("nse", dat.sub=subgrp.effect,
##D                          var.estvar = var.estvar, var.cov = var.cov,
##D                          par.pri = c(B=1000),
##D                          chains=4, iter=600,
##D                          warmup=200, thin=2, seed=1000);
##D 
##D rst.sr     <- bzCallStan("sr", dat.sub=subgrp.effect,
##D                         var.estvar=var.estvar, var.cov = var.cov,
##D                         par.pri=c(B=1000, C=1000),
##D                         chains=4, iter=600,
##D                         warmup=200, thin=2, seed=1000);
## End(Not run)



cleanEx()
nameEx("bzComp")
### * bzComp

flush(stderr()); flush(stdout())

### Name: bzComp
### Title: Comparison of posterior treatment effects
### Aliases: bzComp bzSummaryComp bzPlotComp bzForestComp

### ** Examples

## Not run: 
##D var.cov    <- c("sodium", "lvef", "any.vasodilator.use");
##D var.resp   <- "y";
##D var.trt    <- "trt";
##D var.censor <- "censor";
##D resptype   <- "survival";
##D var.estvar <- c("Estimate", "Variance");
##D 
##D subgrp.effect <- bzGetSubgrpRaw(solvd.sub,
##D                              var.resp   = var.resp,
##D                              var.trt    = var.trt,
##D                              var.cov    = var.cov,
##D                              var.censor = var.censor,
##D                              resptype   = resptype);
##D 
##D rst.sr     <- bzCallStan("sr", dat.sub=subgrp.effect,
##D                          var.estvar=var.estvar, var.cov = var.cov,
##D                          par.pri=c(B=1000, C=1000),
##D                          chains=4, iter=500,
##D                          warmup=100, thin=2, seed=1000);
##D 
##D sel.grps <- c(1,4,5);
##D tbl.sub <- bzSummaryComp(rst.sr, sel.grps=sel.grps);
##D bzPlot(rst.sr, sel.grps = sel.grps);
##D bzForest(rst.sr, sel.grps = sel.grps);
## End(Not run)




cleanEx()
nameEx("bzGailSimon")
### * bzGailSimon

flush(stderr()); flush(stdout())

### Name: bzGailSimon
### Title: Gail-Simon Test
### Aliases: bzGailSimon

### ** Examples

## Not run: 
##D var.cov    <- c("sodium", "lvef", "any.vasodilator.use");
##D var.resp   <- "y";
##D var.trt    <- "trt";
##D var.censor <- "censor";
##D resptype   <- "survival";
##D subgrp.effect <- bzGetSubgrp(solvd.sub,
##D                                   var.resp   = var.resp,
##D                                   var.trt    = var.trt,
##D                                   var.cov    = var.cov,
##D                                   var.censor = var.censor,
##D                                   resptype   = resptype);
##D 
##D gs.pval <- bzGailSimon(subgrp.effect$Estimate,
##D                        subgrp.effect$Variance); 
## End(Not run)





cleanEx()
nameEx("bzGetSubgrpRaw")
### * bzGetSubgrpRaw

flush(stderr()); flush(stdout())

### Name: bzGetSubgrpRaw
### Title: Get subgroup treatment effect estimation and variance
### Aliases: bzGetSubgrpRaw

### ** Examples


## Not run: 
##D var.cov    <- c("sodium", "lvef", "any.vasodilator.use");
##D var.resp   <- "y";
##D var.trt    <- "trt";
##D var.censor <- "censor";
##D resptype   <- "survival";
##D subgrp.effect <- bzGetSubgrpRaw(solvd.sub,
##D                                   var.resp   = var.resp,
##D                                   var.trt    = var.trt,
##D                                   var.cov    = var.cov,
##D                                   var.censor = var.censor,
##D                                   resptype   = resptype);
## End(Not run)




cleanEx()
nameEx("bzPredSubgrp")
### * bzPredSubgrp

flush(stderr()); flush(stdout())

### Name: bzPredSubgrp
### Title: Predictive Distribution
### Aliases: bzPredSubgrp

### ** Examples

## Not run: 
##D var.cov    <- c("sodium", "lvef", "any.vasodilator.use");
##D var.resp   <- "y";
##D var.trt    <- "trt";
##D var.censor <- "censor";
##D resptype   <- "survival";
##D var.estvar <- c("Estimate", "Variance");
##D 
##D subgrp.effect <- bzGetSubgrp(solvd.sub,
##D                                   var.resp   = var.resp,
##D                                   var.trt    = var.trt,
##D                                   var.cov    = var.cov,
##D                                   var.censor = var.censor,
##D                                   resptype   = resptype);
##D 
##D rst.nse    <- bzCallStan("nse", dat.sub=subgrp.effect,
##D                          var.estvar = var.estvar, var.cov = var.cov,
##D                          par.pri = c(B=1000),
##D                          chains=4, iter=4000,
##D                          warmup=2000, thin=2, seed=1000);
##D 
##D pred.effect <- bzPredSubgrp(rst.nes,
##D                             dat.sub = solvd.sub,
##D                             var.estvar = var.estvar);
## End(Not run)



cleanEx()
nameEx("bzSummary")
### * bzSummary

flush(stderr()); flush(stdout())

### Name: bzSummary
### Title: Posterior subgroup treatment effects
### Aliases: bzSummary bzSummary bzPlot bzForest

### ** Examples

## Not run: 
##D sel.grps <- c(1,4,5);
##D tbl.sub <- bzSummary(rst.sr, ref.stan.rst=rst.nse, ref.sel.grps=1);
##D bzPlot(rst.sr, sel.grps = sel.grps, ref.stan.rst=rst.nse, ref.sel.grps=1);
##D bzForest(rst.sr, sel.grps = sel.grps, ref.stan.rst=rst.nse, ref.sel.grps=1);
## End(Not run)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
