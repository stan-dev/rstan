library(bayesplot)
context("PPC: misc. functions")

source(test_path("data-for-ppc-tests.R"))

# melt_yrep ---------------------------------------------------------------
expect_molten_yrep <- function(yrep) {
  y <- rnorm(ncol(yrep))
  yrep <- validate_yrep(yrep, y)

  x <- melt_yrep(yrep)
  expect_equal(ncol(x), 4)
  expect_equal(nrow(x), prod(dim(yrep)))

  rep_nums <- rep(seq_len(nrow(yrep)), length(y))
  obs_nums <- sort(rep(seq_len(length(y)), nrow(yrep)))

  expect_identical(colnames(x), c("y_id", "rep_id", "rep_label", "value"))
  expect_equal(x$y_id, obs_nums)
  expect_equal(x$rep_id, rep_nums)

  expect_s3_class(x, "data.frame")
  expect_s3_class(x$rep_label, "factor")
  expect_type(x$rep_id, "integer")
  expect_type(x$y_id, "integer")
  expect_type(x$value, "double")
}

test_that("melt_yrep returns correct structure", {
  expect_molten_yrep(yrep)
  expect_molten_yrep(yrep2)

  load(test_path("data-for-binomial.rda"))
  expect_molten_yrep(Ey)
  expect_molten_yrep(validate_yrep(yrep, y))
})


# melt_and_stack ----------------------------------------------------------
test_that("melt_and_stack returns correct structure", {
  molten_yrep <- melt_yrep(yrep)
  d <- melt_and_stack(y, yrep)
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), nrow(molten_yrep) + length(y))

  sorted_names <- sort(c(colnames(molten_yrep), c("is_y", "is_y_label")))
  expect_equal(sort(colnames(d)), sorted_names)
})


# ppc_group_data ----------------------------------------------------------

d <- ppc_group_data(y, yrep, group)
d_stat <- ppc_group_data(y, yrep, group, stat = "mean")

test_that("ppc_group_data returns correct structure", {
  expect_identical(colnames(d), c("group", "variable", "value"))
  expect_s3_class(d, c("grouped_df", "tbl_df", "tbl", "data.frame"))

  expect_identical(colnames(d_stat), colnames(d))
  expect_s3_class(d, c("grouped_df", "tbl_df", "tbl", "data.frame"))

  nr <- length(unique(d$variable)) * length(unique(group))
  expect_equal(nrow(d_stat), nr)
})

test_that("ppc_group_data with stat returns correct values for y", {
  for (lev in levels(group)) {
    mean_y_group <- with(d_stat, value[group == lev & variable == "y"])
    expect_equal(mean_y_group, mean(y[group == lev]),
                 info = paste("group =", lev))
  }
})

test_that("ppc_group_data with stat returns correct values for yrep", {
  for (lev in levels(group)) {
    for (j in 1:nrow(yrep)) {
      var <- paste0("yrep_", j)
      mean_yrep_group <- with(d_stat, value[group == lev & variable == var])
      expect_equal(mean_yrep_group, mean(yrep[j, group == lev]),
                   info = paste("group =", lev, "|", "rep =", j))
    }
  }
})


# is_whole_number, all_counts --------------------------------------------
test_that("is_whole_number works correctly", {
  expect_equal(is_whole_number(c(1L, 2, 3/3, 4/5)),
               c(rep(TRUE, 3), FALSE))
  expect_true(!is_whole_number("1"))
})
test_that("all_counts works correctly", {
  expect_true(all_counts(1))
  expect_true(all_counts(0:5))
  expect_true(all_counts(matrix(rpois(10, 1), 2, 5)))
  expect_false(all_counts(rnorm(5)))
  expect_false(all_counts(c("1", "2")))
  expect_false(all_counts(c(1, 1.5)))
  expect_false(all_counts(c(-1, 2)))
})

