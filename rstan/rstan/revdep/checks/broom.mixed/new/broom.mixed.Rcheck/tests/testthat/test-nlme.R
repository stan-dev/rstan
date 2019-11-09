
## test tidy, augment, glance methods from nlme-tidiers.R
stopifnot(require("testthat"), require("broom.mixed"))

load(system.file("extdata", "nlme_example.rda", package = "broom.mixed",
                 mustWork=TRUE))

if (suppressPackageStartupMessages(require(nlme, quietly = TRUE))) {
  context("nlme models")

  d <- as.data.frame(ChickWeight)
  colnames(d) <- c("y", "x", "subj", "tx")
  fit <- lme(y ~ tx * x, random = ~x | subj, data = d)

  test_that("tidy works on nlme/lme fits", {
    td <- tidy(fit)
    expect_equal(
      names(td),
      c(
        "effect", "group", "term", "estimate",
        "std.error", "df", "statistic", "p.value"
      )
    )
  })

  test_that("tidy works on non-linear fits", {
    fm1 <- nlme(height ~ SSasymp(age, Asym, R0, lrc),
      data = Loblolly,
      fixed = Asym + R0 + lrc ~ 1,
      random = Asym ~ 1,
      start = c(Asym = 103, R0 = -8.5, lrc = -3.3)
    )
    td <- expect_warning(tidy(fm1), "ran_pars not yet implemented")
    expect_equal(
      names(td),
      c(
        "effect", "group", "term", "estimate",
        "std.error", "df", "statistic", "p.value"
      )
    )
    td_ran <- expect_warning(tidy(fm1, "ran_pars"))
    expect_is(td_ran, "tbl_df")
    expect_equal(ncol(td_ran), 0)
    expect_equal(nrow(td_ran), 0)
    td_fix <- tidy(fm1, "fixed")
    expect_equal(
      names(td_fix),
      c(
        "term", "estimate",
        "std.error", "df", "statistic", "p.value"
      )
    )
  })

  test_that("augment works on lme fits with or without data", {
    au1 <- augment(fit)
    au2 <- augment(fit, d)
    expect_equal(au1, au2)
    expect_equal(dim(au1), c(578, 7))
  })
  dNAs <- d
  dNAs$y[c(1, 3, 5)] <- NA

  test_that("augment works on lme fits with NAs and na.omit", {
    fitNAs <- lme(y ~ tx * x,
      random = ~x | subj, data = dNAs,
      na.action = "na.omit"
    )
    au <- augment(fitNAs)
    expect_equal(nrow(au), sum(complete.cases(dNAs)))
  })


  test_that("augment works on lme fits with na.omit", {
    fitNAs <- lme(y ~ tx * x,
      random = ~x | subj, data = dNAs,
      na.action = "na.exclude"
    )

    au <- augment(fitNAs, dNAs)

    # with na.exclude, should have NAs in the output where there were NAs in input
    expect_equal(nrow(au), nrow(dNAs))
    expect_equal(complete.cases(au), complete.cases(dNAs))
  })

  test_that("glance includes deviance iff method='ML'", {
    expect(!("deviance" %in% names(glance(lmm0))))
    expect("deviance" %in% names(glance(lmm0ML)))
  })

  ## FIXME: weak tests - only shows that no error occurs and
  ##  the right type is returned!
  test_that("glance works on nlme fits", {
    expect_is(glance(fit), "data.frame")
  })


  testFit <- function(fit, data = NULL) {
    test_that("Pinheiro/Bates fit works", {
      expect_is(tidy(fit, "fixed"), "data.frame")
      # TODO: Better idea than suppressWarnings to avoid "ran_pars" ??
      expect_is(suppressWarnings(tidy(fit)), "data.frame")
      expect_is(glance(fit), "data.frame")
      if (is.null(data)) {
        expect_is(augment(fit), "data.frame")
      } else {
        expect_is(augment(fit, data), "data.frame")
      }
    })
  }

  testFit(lme(score ~ Machine, data = Machines, random = ~1 | Worker))
  testFit(lme(score ~ Machine, data = Machines, random = ~1 | Worker, method = "ML"))
  testFit(lme(score ~ Machine, data = Machines, random = ~1 | Worker / Machine))
  testFit(lme(pixel ~ day + day^2, data = Pixel, random = list(Dog = ~day, Side = ~1)))
  testFit(lme(pixel ~ day + day^2 + Side,
    data = Pixel,
    random = list(Dog = ~day, Side = ~1)
  ))

  testFit(lme(yield ~ ordered(nitro) * Variety,
    data = Oats,
    random = ~1 / Block / Variety
  ))
  # There are cases where no data set is returned in the result
  # We can do nothing about this inconsistency but give a useful error message in augment
  fit <- nlme(conc ~ SSfol(Dose, Time, lKe, lKa, lCl),
    data = Theoph,
    random = pdDiag(lKe + lKa + lCl ~ 1)
  )
  test_that(
    "Fit without data in returned structure works when data are given", {
      expect_true(testFit(fit, Theoph))
    }
  )
  # When no data are passed, a meaningful message is issued
  expect_error(augment(fit), "explicit")

  context("gls models")
  
  test_that("basic gls tidying", {

      check_tidy(tidy(gls1), 3, 5,
                 c("term","estimate","std.error","statistic","p.value"))
      check_tidy(tidy(gls1, conf.int=TRUE), 3, 7,
                 c("term","estimate","std.error","statistic","p.value",
                   "conf.low","conf.high"))
    
  })      
}
