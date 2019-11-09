stopifnot(require("testthat"), require("broom.mixed"))
## test lmerTest

## HACK: need to make sure we find the right generic!
tidy <- broom.mixed:::tidy.merMod

if (require(lmerTest, quietly = TRUE)) {
  test_that("testing lmerTest p-values behind Douglas Bates' back", {
    lmm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
    td <- tidy(lmm1, "fixed")
    check_tidy(td, 2, 7, c(
      "effect", "term", "estimate",
      "std.error", "df", "statistic", "p.value"
    ))
    td_ran <- tidy(lmm1, "ran_pars")
    check_tidy(td_ran, 4, 4, c("effect", "group", "term", "estimate"))
    expect_false(all(is.na(td_ran$estimate)))
  })
}
