library(bayesplot)
library(ggplot2)
context("Convenience functions (for ggplot objects)")


# abline_01, vline_ and hline_ ------------------------------------------
test_that("abline_01 returns the correct object", {
  expect_equal(
    abline_01(color = "green", linetype = 2),
    geom_abline(intercept = 0, slope = 1, color = "green", linetype = 2, na.rm = TRUE)
  )
})
test_that("vline_* and hline_* return correct objects", {
  expect_equal(
    vline_0(color = "red"),
    geom_vline(xintercept = 0, color = "red", na.rm = TRUE)
  )
  expect_equal(
    hline_0(size = 2, linetype = 3),
    geom_hline(yintercept = 0, size = 2, linetype = 3, na.rm = TRUE)
  )
  expect_equal(
    vline_at(c(3,4), na.rm = FALSE),
    geom_vline(xintercept = c(3,4))
  )
  expect_equal(
    hline_at(c(3,4), na.rm = FALSE),
    geom_hline(yintercept = c(3,4))
  )
})
test_that("vline_at with 'fun' works", {
  x <- example_mcmc_draws(chains = 1)
  expect_equal(vline_at(x, colMeans),
               geom_vline(xintercept = colMeans(x), na.rm = TRUE))
})
test_that("calc_v (internal function) works", {
  a <- 1:4
  expect_identical(calc_v(a, "mean"), 2.5)
  expect_identical(calc_v(a, median), 2.5)
  expect_equal(calc_v(c(a, NA), mean), NA_real_)
  expect_identical(calc_v(c(a, NA), min, list(na.rm = TRUE)), 1L)
})

# lbub --------------------------------------------------------------------
test_that("lbub works", {
  f1 <- lbub(p = 0.5)
  f2 <- lbub(p = 0.5, med = FALSE)

  expect_type(f1, "closure")
  expect_type(f2, "closure")
  expect_identical(
    f1(1:50),
    setNames(c(13.25, 25.5, 37.75), c("25%", "50%", "75%"))
  )
  expect_identical(
    f2(1:50),
    setNames(c(13.25, 37.75), c("25%", "75%"))
  )
})

# plot and facet backgrounds ----------------------------------------------
test_that("grid_lines returns correct theme object", {
  thm <- theme_default() + grid_lines(size = 1.5, color = "purple")
  expect_equal(thm$panel.grid.major, element_line(size = 1.5, color = "purple"))
  expect_equal(thm$panel.grid.minor, element_line(size = 0.75, color = "purple"))
})
test_that("panel_bg returns correct theme object", {
  bg1 <- panel_bg()
  bg2 <- panel_bg(fill = "blue", linetype = 2)

  expect_identical(bg1, theme(panel.background = element_rect()))
  expect_identical(bg2, theme(panel.background = element_rect(fill = "blue", linetype = 2)))
  expect_identical(panel_bg(on = FALSE), theme(panel.background = element_blank()))
})
test_that("plot_bg returns correct theme object", {
  bg1 <- plot_bg()
  bg2 <- plot_bg(fill = "blue", linetype = 2)

  expect_identical(bg1, theme(plot.background = element_rect()))
  expect_identical(bg2, theme(plot.background = element_rect(fill = "blue", linetype = 2)))
  expect_identical(plot_bg(on = FALSE), theme(plot.background = element_blank()))
})
test_that("facet_bg returns correct theme object", {
  bg1 <- facet_bg()
  bg2 <- facet_bg(fill = "blue", linetype = 2)

  expect_identical(bg1, theme(strip.background = element_rect()))
  expect_identical(bg2, theme(strip.background = element_rect(fill = "blue", linetype = 2)))
  expect_identical(facet_bg(on = FALSE), theme(strip.background = element_blank()))
})

# legend position and text ------------------------------------------------
test_that("legend_none returns correct theme object", {
  none <- legend_none()
  expect_s3_class(none, "theme")
  expect_equivalent(none, list(legend.position = "none"))
  expect_false(attr(none, "complete"))
})
test_that("legend_move returns correct theme object", {
  left <- legend_move("left")
  expect_s3_class(left, "theme")
  expect_equivalent(left, list(legend.position = "left"))
  expect_false(attr(left, "complete"))

  pos <- legend_move(c(0.25, 0.5))
  expect_s3_class(pos, "theme")
  expect_equivalent(pos, list(legend.position = c(0.25, 0.5)))
  expect_false(attr(pos, "complete"))
})
test_that("legend_text returns correct theme object", {
  expect_equal(
    legend_text(size = 16, color = "purple"),
    theme(legend.text = element_text(color = "purple", size = 16))
  )
})

# axis and facet text --------------------------------------------------
test_that("xaxis_text returns correct theme object", {
  expect_identical(xaxis_text(FALSE), theme(axis.text.x = element_blank()))
  expect_equal(
    xaxis_text(face = "bold", angle = 30),
    theme(axis.text.x = element_text(face = "bold", angle = 30))
  )
})
test_that("yaxis_text returns correct theme object", {
  expect_identical(yaxis_text(FALSE), theme(axis.text.y = element_blank()))
  expect_equivalent(
    yaxis_text(face = "bold", angle = 30),
    theme(axis.text.y = element_text(face = "bold", angle = 30))
  )
})
test_that("facet_text returns correct theme object", {
  expect_identical(facet_text(FALSE), theme(strip.text = element_blank()))
  expect_equal(
    facet_text(size = 12, color = "blue"),
    theme(strip.text = element_text(color = "blue", size = 12))
  )
})

# axis titles -------------------------------------------------------------
test_that("xaxis_title returns correct theme object", {
  expect_identical(xaxis_title(FALSE), xlab(NULL))
  expect_equal(
    xaxis_title(face = "bold", angle = 30),
    theme(axis.title.x = element_text(face = "bold", angle = 30))
  )
})
test_that("yaxis_title returns correct theme object", {
  expect_identical(yaxis_title(FALSE), ylab(NULL))
  expect_equal(
    yaxis_title(face = "bold", angle = 30),
    theme(axis.title.y = element_text(face = "bold", angle = 30))
  )
})

# tick marks --------------------------------------------------
test_that("xaxis_ticks returns correct theme object", {
  expect_identical(xaxis_ticks(FALSE), theme(axis.ticks.x = element_blank()))
  expect_equal(
    xaxis_ticks(size = 0.5, color = "red"),
    theme(axis.ticks.x = element_line(size = 0.5, color = "red"))
  )
})
test_that("yaxis_ticks returns correct theme object", {
  expect_identical(yaxis_ticks(FALSE), theme(axis.ticks.y = element_blank()))
  expect_equal(
    yaxis_ticks(size = 0.5, color = "red"),
    theme(axis.ticks.y = element_line(size = 0.5, color = "red"))
  )
})


# overlay functions -------------------------------------------------------
test_that("overlay_function returns the correct object", {
  expect_error(overlay_function(), 'argument "fun" is missing')
  expect_equal(
    overlay_function(fun = "dnorm"),
    stat_function(fun = "dnorm", inherit.aes = FALSE)
  )
})


