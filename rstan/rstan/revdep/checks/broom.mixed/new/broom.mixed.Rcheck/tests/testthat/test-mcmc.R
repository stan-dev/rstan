## test tidy and glance methods from mcmc_tidiers
stopifnot(require("testthat"), require("broom.mixed"))

context("stan tidiers")

if (suppressPackageStartupMessages(require(rstan, quietly = TRUE))) {
  test_that("tidy returns indexes if requested on rstanarm fits", {

    # Make sure that (inst/)extdata/run_examples.R was run to generate rds
    rstan_example <- readRDS(system.file("extdata", "rstan_example.rds", package = "broom.mixed"))
    # check_tidy from helper-checkers
    td <- tidy(rstan_example)
    check_tidy(td, 18, 3, c("term", "estimate", "std.error"))

    td <- tidy(rstan_example, index = TRUE)
    check_tidy(td, 18, 4, c("term", "index", "estimate", "std.error"))

    td <- tidy(rstan_example, drop.pars = NULL)
    expect_equal(td[19, ][["term"]], "lp__")

    td <- tidy(rstan_example, conf.int = TRUE)
    check_tidy(td, 18, 5, c("term", "estimate", "std.error", "conf.low", "conf.high"))

    td <- tidy(rstan_example, rhat = TRUE)
    check_tidy(td, 18, 4, c("term", "estimate", "std.error", "rhat"))

    td <- tidy(rstan_example, ess = TRUE)
    check_tidy(td, 18, 4, c("term", "estimate", "std.error", "ess"))
  })
}

context("brms tidiers")

if (suppressPackageStartupMessages(require("brms", quietly = TRUE))) {
    load(system.file("extdata","brms_example.rda",
                     package="broom.mixed",
                     mustWork=TRUE))
    ## n.b. different S3 methods found depending on environment
    zz <-  tidy(brms_zip,effects="ran_vals")
    zz2 <- tidy(brms_zip)
    zz3 <- tidy(brms_multi)
    expect_warning(tidy(brms_multi_RE),"currently incorrect")
    suppressWarnings(zz4 <- tidy(brms_multi_RE))
    test_that("correct levels for models with zi/ranef",{
        expect_equal(zz[["level"]],
                     rep(c(paste("R",1:12,sep="-"),paste("VF",1:11,sep="-")),2))
    })

    test_that("component returned for brms models",
    {
        expect_equal(zz2[["component"]],
                     unlist(lapply(list(c(1,1),c(13,1),c(1,1)),
                                   rep,x=c("cond","zi"))))
    })

    test_that("component tags stripped from brms models",
    {
        expect_equal(c(table(zz2[["term"]])),
                     c(`(Intercept)` = 2L, minedno = 2L, `sd__(Intercept)` = 2L,
                       sppDESML = 1L, 
           `sppDESML:minedno` = 1L, sppDF = 1L, `sppDF:minedno` = 1L,
           sppDM = 1L, `sppDM:minedno` = 1L, sppECMA = 1L,
           `sppECMA:minedno` = 1L, sppECML = 1L, `sppECML:minedno` = 1L,
           sppPR = 1L, `sppPR:minedno` = 1L))
    })

    test_that("multi-component brms models",
    {
        check_tidy(zz3, 8, 9,
                   c("response", "effect", "component", "group",
                     "term", "estimate", "std.error",
                     "conf.low", "conf.high"))
    })
    
} ## if require("brms")

context("mcmc tidiers")

if (suppressPackageStartupMessages(require(coda, quietly = TRUE))) {
    data(line)
    x1 <- line[[1]]
    
    test_that("mcmc with ess",
    {

        expect_warning(td <- tidy(
                           x = x1,
                           conf.int = TRUE,
                           robust = TRUE,
                           rhat = TRUE,
                           index = TRUE,
                           ess = TRUE),
                       "ignoring 'rhat'")
        check_tidy(td, 3, 7,
                   c("term","index","estimate","std.error",
                     "conf.low","conf.high","ess"))
    })
}
