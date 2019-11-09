context("Code Style")

test_that("Code style is in line with INWT style conventions", {
  lintr::expect_lint_free(linters = INWTUtils::selectLinters())
})
