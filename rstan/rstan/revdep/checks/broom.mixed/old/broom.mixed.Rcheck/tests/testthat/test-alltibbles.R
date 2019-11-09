context("test all tidy methods return tibbles")
verbose <- FALSE
stopifnot(
  require("testthat"), require("broom.mixed"), require("stringr"),
  require("lme4"), require("glmmTMB")
)

## objects to SKIP (from RDA files)
skip_objects <- list(nlme_example.rda = c("lmm2"),
                  ## glance method not working for this case
                  brms_example.rda = c("brms_multi","brms_multi_RE")
                  )

skip_files <- c("efc.rds")

## need data available for augment() methods ...
data("sleepstudy", package = "lme4")
data("Salamanders", package = "glmmTMB")
efc <- readRDS(system.file("extdata","efc.rds",package="broom.mixed"))
ex <- list.files(system.file("extdata", package = "broom.mixed"),
  pattern = "\\.rd"
  )
ex <- setdiff(ex,skip_files)


## test tidy, augment, glance methods from lme4-tidiers.R

for (e in ex) {
    p <- stringr::str_extract(e, "^[^_]+")
    if (require(p,character.only=TRUE)) {
        f <- system.file("extdata", e, package = "broom.mixed")
        fn <- stringr::str_extract(f, "[^/]+$")
        if (verbose) cat(fn, "\n")
        if (grepl("\\.rds", e)) {
            x <- list()
            x[[1]] <- readRDS(f)
        } else {
            ## rda file
            L <- load(f)
            x <- mget(setdiff(L, skip_objects[[fn]]))
        }
        for (z in x) {
            testf <- function(fn_str, obj) {
                cc <- class(obj)[1]
                if (sum(grepl(cc, methods(fn_str))) > 0) {
                    if (verbose) cat(sprintf("found method %s for %s\n", fn_str, cc))
                    return(expect_is(get(fn_str)(obj), "tbl_df"))
                } else {
                    return(TRUE)
                }
            }
            testf("glance", z)
            testf("augment", z)
            testf("tidy", z)
        } ## loop over objects
    } ## if package available
}
