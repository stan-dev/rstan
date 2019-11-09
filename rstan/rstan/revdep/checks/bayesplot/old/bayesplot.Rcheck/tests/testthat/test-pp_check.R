library(bayesplot)
context("PPC: pp_check generic and default method")

test_that("default pp_check method works", {
  y <- example_y_data()
  yrep <- example_yrep_draws()
  g <- example_group_data()

  expect_equal(
    pp_check(y, yrep[1:50, ], ppc_dens_overlay),
    ppc_dens_overlay(y, yrep[1:50, ])
  )
  expect_equal(
    pp_check(y, yrep, fun = "stat_grouped", group = g, stat = "median"),
    ppc_stat_grouped(y, yrep, group = g, stat = "median")
  )
})

test_that("pp_check method can be defined", {
  pp_check.foo <- function(object, ..., type = c("multiple", "overlaid")) {
    y <- object[["y"]]
    yrep <- object[["yrep"]]
    switch(match.arg(type),
           multiple = ppc_hist(y, yrep[1:min(8, nrow(yrep)),, drop = FALSE]),
           overlaid = ppc_dens_overlay(y, yrep)
    )
  }

  x <- structure(
    list(y = rnorm(50), yrep = matrix(rnorm(500), ncol = 50)),
    class = "foo"
  )
  expect_gg(pp_check(x))
  expect_gg(pp_check(x, type = "overlaid"))
})
