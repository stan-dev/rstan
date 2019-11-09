library(bayesplot)
context("PPC: predictive errors")

source(test_path("data-for-ppc-tests.R"))

test_that("ppc_error_hist and ppc_error_scatter return ggplot object", {
  expect_gg(ppc_error_hist(y, yrep[1:5, ]))
  expect_gg(ppc_error_scatter(y, yrep[1:5, ]))

  expect_gg(ppc_error_hist(y, yrep[1,, drop = FALSE]))
  expect_gg(ppc_error_scatter(y, yrep[1,, drop = FALSE]))

  expect_gg(ppc_error_hist(y2, yrep2))
  expect_gg(ppc_error_scatter(y2, yrep2))
})

test_that("ppc_error_hist_grouped returns ggplot object", {
  expect_gg(ppc_error_hist_grouped(y, yrep[1:5, ], group))
  expect_gg(ppc_error_hist_grouped(y, yrep[1,, drop = FALSE], group,
                                   freq = FALSE, binwidth = 1))
  expect_error(ppc_error_hist_grouped(y2, yrep2, group2),
               "'group' must have more than one unique value")
})

test_that("ppc_error_scatter_avg returns ggplot2 object", {
  expect_gg(ppc_error_scatter_avg(y, yrep))
  expect_gg(ppc_error_scatter_avg(y, yrep[1:5, ]))
})

test_that("ppc_error_scatter_avg same as ppc_error_scatter if nrow(yrep) = 1", {
  expect_equal(ppc_error_scatter_avg(y2, yrep2),
               ppc_error_scatter(y2, yrep2))
  expect_equal(ppc_error_scatter_avg(y, yrep[1,, drop=FALSE]),
               ppc_error_scatter(y, yrep[1,, drop = FALSE]))
})

test_that("ppc_error_scatter_avg_vs_x returns ggplot2 object", {
  expect_gg(ppc_error_scatter_avg_vs_x(y, yrep, x = rnorm(length(y))))
  expect_gg(ppc_error_scatter_avg_vs_x(y, yrep[1:5, ], x = rnorm(length(y))))
})

test_that("ppc_error_binned returns ggplot object", {
  load(test_path("data-for-binomial.rda"))
  expect_gg(ppc_error_binned(y, Ey))
  expect_gg(ppc_error_binned(y[1:5], Ey[, 1:5]))
  expect_gg(ppc_error_binned(rep(y, 2), cbind(Ey, Ey)))
})
