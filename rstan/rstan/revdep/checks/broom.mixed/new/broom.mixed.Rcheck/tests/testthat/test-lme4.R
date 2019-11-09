
stopifnot(require("testthat"), require("broom.mixed"))
## test tidy, augment, glance methods from lme4-tidiers.R

if (require(lme4, quietly = TRUE)) {
  load(system.file("extdata", "lme4_example.rda",
    package = "broom.mixed",
    mustWork = TRUE
    ))

context("lme4 models")

  d <- as.data.frame(ChickWeight)
  colnames(d) <- c("y", "x", "subj", "tx")
  fit <<- lmer(y ~ tx * x + (x | subj), data = d)

  test_that("tidy works on lme4 fits", {
    td <- tidy(fit)
    ## FIXME: fails if lmerTest has been loaded previously ...
    expect_equal(dim(td), c(12, 6))
    expect_equal(
      names(td),
      c(
        "effect", "group", "term", "estimate",
        "std.error", "statistic"
      )
    )
    expect_equal(
      td$term,
      c(
        "(Intercept)", "tx2", "tx3", "tx4", "x",
        "tx2:x", "tx3:x", "tx4:x",
        "sd__(Intercept)", "sd__x",
        "cor__(Intercept).x", "sd__Observation"
      )
    )
  })

  test_that("tidy/glance works on glmer fits", {
    gm <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
      cbpp, binomial,
      nAGQ = 0
    )
    ggm <- broom::glance(gm)
    expect_equal(names(ggm), c("sigma", "logLik", "AIC", "BIC", "deviance", "df.residual"))
    td <- tidy(gm)
    expect_equal(
      names(td),
      c(
        "effect", "group", "term", "estimate",
        "std.error", "statistic", "p.value"
      )
    )
    td_ran <- tidy(gm, "ran_pars")
    expect_equal(names(td_ran), c("effect", "group", "term", "estimate"))
  })

  test_that("glance includes deviance iff method='ML'", {
    expect(!("deviance" %in% names(glance(lmm0))))
    expect("REMLcrit" %in% names(glance(lmm0)))
    expect("deviance" %in% names(glance(lmm0ML)))
  })


  test_that("tidy works on non-linear fits", {
    startvec <- c(Asym = 200, xmid = 725, scal = 350)
    # use nAGQ = 0 to avoid warnings
    nm <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym | Tree,
      Orange,
      start = startvec, nAGQ = 0L
    )
    gnm <- broom::glance(nm)
    expect_equal(names(gnm), c("sigma", "logLik", "AIC", "BIC", "deviance", "df.residual"))
    td <- tidy(nm)
    expect_equal(
      names(td),
      c(
        "effect", "group", "term", "estimate",
        "std.error", "statistic"
      )
    )
    td_ran <- tidy(nm, "ran_pars")
    expect_equal(names(td_ran), c("effect", "group", "term", "estimate"))
  })

  test_that("scales works", {
    t1 <- tidy(fit, effects = "ran_pars")
    t2 <- tidy(fit, effects = "ran_pars", scales = "sdcor")
    expect_equal(t1$estimate, t2$estimate)
    expect_error(
      tidy(fit, effects = "ran_pars", scales = "varcov"),
      "unrecognized ran_pars scale"
    )
    t3 <- tidy(fit, effects = "ran_pars", scales = "vcov")
    expect_equal(
      t3$estimate[c(1, 2, 4)],
      t2$estimate[c(1, 2, 4)]^2
    )
    expect_error(
      tidy(fit, scales = "vcov"),
      "must be provided for each effect"
    )
  })

  test_that("tidy works with more than one RE grouping variable", {
    dd <- expand.grid(f = factor(1:10), g = factor(1:5), rep = 1:3)
    dd$y <- suppressMessages(simulate(~(1 | f) + (1 | g),
      newdata = dd,
      newparams = list(beta = 1, theta = c(1, 1)),
      family = poisson, seed = 101
    ))[[1]]
    gfit <- glmer(y ~ (1 | f) + (1 | g), data = dd, family = poisson)
    tnames <- as.character(tidy(gfit, effects = "ran_pars")$term)
    expect_equal(tnames, rep("sd__(Intercept)", 2))
  })

  test_that("augment works on lme4 fits with or without data", {
    au1 <- suppressWarnings(broom::augment(fit))
    au2 <- suppressWarnings(broom::augment(fit, d))
    ## FIXME: columns not ordered the same??
    expect_equal(au1, au2[names(au1)])
  })

  dNAs <<- d
  dNAs$y[c(1, 3, 5)] <- NA

  test_that("augment works on lme4 fits with NAs", {
    fitNAs <- lmer(y ~ tx * x + (x | subj), data = dNAs,
                     control=lmerControl(check.conv.grad=
                     .makeCC("warning", tol = 5e-2, relTol = NULL)))
    au <- suppressWarnings(broom::augment(fitNAs))
    expect_equal(nrow(au), sum(complete.cases(dNAs)))
  })

  test_that("augment works on lme4 fits with na.exclude", {
      fitNAs <- lmer(y ~ tx * x + (x | subj),
                     data = dNAs, na.action = "na.exclude",
                     control=lmerControl(check.conv.grad=
                     .makeCC("warning", tol = 5e-2, relTol = NULL)))

    # expect_error(suppressWarnings(augment(fitNAs)))
    au <- suppressWarnings(broom::augment(fitNAs, dNAs))

    # with na.exclude, should have NAs in the output where there were NAs in input
    expect_equal(nrow(au), nrow(dNAs))
    expect_equal(complete.cases(au), complete.cases(dNAs))
  })

  test_that("glance works on lme4 fits", {
    g <- broom::glance(fit)
    expect_equal(dim(g), c(1, 6))
  })

  test_that("ran_vals works", {
    td0 <- tidy(lmm0, "ran_vals")
    td1 <- tidy(lmm1, "ran_vals")
    expect_equal(dim(td0), c(18, 6))
    expect_equal(dim(td1), c(36, 6))
    if (packageVersion("lme4") >= "1.1.18") {
      td2 <- tidy(lmm2, "ran_vals")
      expect_equal(dim(td2), c(36, 6))
      expect_equal(names(td1), names(td2))
    }
  })
  test_that("confint preserves term names", {
    td3 <- tidy(lmm0, conf.int = TRUE, conf.method = "Wald", effects = "fixed")
    expect_equal(td3$term, c("(Intercept)", "Days"))
  })
}

test_that("tidy respects conf.level", {
     tmpf <- function(cl=0.95) {
         return(tidy(lmm0,conf.int=TRUE,conf.level=cl)[1,][["conf.low"]])
     }
     expect_equal(tmpf(),232.3019,tolerance=1e-4)
     expect_equal(tmpf(0.5),244.831,tolerance=1e-4)
})

test_that("effects='ran_pars' + conf.int works", {
    tt <- tidy(lmm0, effects="ran_pars", conf.int=TRUE, conf.method="profile",
               quiet=TRUE)[c("conf.low","conf.high")]
    tt0 <- structure(list(conf.low = c(26.007120448854, 27.8138472081303
), conf.high = c(52.9359835296834, 34.591049857869)), row.names = c(NA, 
-2L), class = c("tbl_df", "tbl", "data.frame"))
    tt0 <- structure(list(conf.low = c(26.00712, 27.81384),
                                conf.high = c(52.9359, 34.59104)),
                           row.names = c(NA, -2L),
                     class = c("tbl_df", "tbl", "data.frame"))
    ## ??? why do I need as.data.frame??
    ## otherwise [1] "Rows in x but not y: 2, 1. Rows in y but not x: 2, 1. "
    expect_equal(as.data.frame(tt0), as.data.frame(tt),
                 tolerance=1e-5)

})

test_that("augment returns a tibble", {
    ## GH 51
    expect_is(augment(fit), "tbl")
})

## KEEP THIS LAST to avoid screwing up S3 methods stack
if (require(lmerTest, quietly = TRUE)) {
  context("lmerTest")
  test_that("lmerTest results include p-values", {
    lmm1X <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
    expect("p.value" %in% names(tidy(lmm1X, effect = "fixed")),
           "no p value in lmerTest results")
  })
  detach("package:lmerTest")
}

