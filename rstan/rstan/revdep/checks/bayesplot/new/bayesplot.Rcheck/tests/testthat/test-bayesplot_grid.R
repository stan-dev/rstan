library(bayesplot)
context("bayesplot_grid")

y <- example_y_data()
yrep <- example_yrep_draws()[1:25, ]
gr <- gridExtra::arrangeGrob(ppc_stat(y, yrep, binwidth = 1))
p1 <- ppc_scatter_avg(y, yrep)
p2 <- ppc_stat(y, yrep, binwidth = 1)

test_that("as_bayesplot_grid works", {
  expect_s3_class(as_bayesplot_grid(gr), "bayesplot_grid")
  expect_s3_class(as_bayesplot_grid(gr), "gtable")
})

test_that("bayesplot_grid throws correct errors", {
  expect_error(bayesplot_grid(xlim = 2),
               "No plots specified")
  expect_error(bayesplot_grid(gr, plots = list(p1, p2)),
               "'...' and 'plots' can't both be specified")
  expect_error(bayesplot_grid(plots = gr),
               "'plots' must be a list of ggplot objects")
  expect_error(bayesplot_grid(gr),
               "objects in '...' must be ggplot objects.")
  expect_error(bayesplot_grid(p1, p2, titles = c("plot1")),
               "length(titles) == length(plots) is not TRUE", fixed = TRUE)
  expect_error(bayesplot_grid(p1, p2, subtitles = c("plot1")),
               "length(subtitles) == length(plots) is not TRUE", fixed = TRUE)
})

test_that("bayesplot_grid works", {
  expect_message(
    a <- bayesplot_grid(p1, p2, xlim = c(-200, 200), ylim = c(0, 200)),
    "Adding another scale for 'y'"
  )
  expect_silent(
    b <- bayesplot_grid(plots = list(p1, p2),
                        titles = c("plot1", "plot2"),
                        subtitles = c("plot1_sub", "plot2_sub"),
                        legends = FALSE)
  )

  expect_s3_class(a, "bayesplot_grid")
  expect_s3_class(b, "bayesplot_grid")
  expect_equal(length(a$grobs), 2)
  expect_equal(length(b$grobs), 2)
})
