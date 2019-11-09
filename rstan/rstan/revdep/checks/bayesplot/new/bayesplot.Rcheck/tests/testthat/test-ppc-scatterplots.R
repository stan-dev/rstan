library(bayesplot)
context("PPC: scatterplots")

source(test_path("data-for-ppc-tests.R"))

test_that("ppc_scatter returns ggplot object", {
  expect_gg(ppc_scatter(y, yrep[1,, drop = FALSE]))
  expect_gg(ppc_scatter(y, yrep[1:3, ]))
  expect_gg(ppc_scatter(y2, yrep2))
})

test_that("ppc_scatter_avg returns ggplot object", {
  expect_gg(ppc_scatter_avg(y, yrep))
  expect_gg(ppc_scatter_avg(y2, yrep2))
})

test_that("ppc_scatter_avg same as ppc_scatter if nrow(yrep) = 1", {
  expect_equal(ppc_scatter_avg(y2, yrep2),
               ppc_scatter(y2, yrep2))
  expect_equal(ppc_scatter_avg(y, yrep[1,, drop=FALSE]),
               ppc_scatter(y, yrep[1,, drop = FALSE]))
})

test_that("ppc_scatter_avg_grouped returns a ggplot object", {
  expect_gg(ppc_scatter_avg_grouped(y, yrep, group))
  expect_gg(ppc_scatter_avg_grouped(y, yrep, as.numeric(group)))
  expect_gg(ppc_scatter_avg_grouped(y, yrep, as.integer(group)))

  expect_error(ppc_scatter_avg_grouped(y2, yrep2, group2),
               "'group' must have more than one unique value")
})
