context("Handling missing values")

# Create example data
suppressWarnings(RNGversion("3.5.0"))
set.seed(10143)
n <- 75
x1 <- rnorm(n, sd = 2)
x2 <- rnorm(n, sd = 1)
x3 <- rnorm(n, sd = 3)
x4 <- rnorm(n, sd = 0.5)

y <- 0.4 + 0.3 * x1 + 0.3 * x1 * x3 + 0.4 * x1 ^ 2 * x2 ^ 3 + rnorm(n, sd = 0.3)

yUncertainty <- rexp(n, 10) * 0.01
data <- data.frame(x1, x2, x3, x4, y)
formula <- y ~ x1 + x2 + x3


test_that("Relevant variables in complete data are not changed by missings handling", {
  expect_equal(BMSC:::handleMissingData(data, formula, yUncertainty)$data,
               data[, c("x1", "x2", "x3", "y")])
  expect_equal(BMSC:::handleMissingData(data,
                                                           formula,
                                                           yUncertainty)$yUncertainty,
               yUncertainty)
  expect_message(BMSC:::handleMissingData(data, formula, yUncertainty)$data, NA)
  expect_warning(BMSC:::handleMissingData(data, formula, yUncertainty)$data, NA)
  
  data$x4[1:5] <- NA
  expect_equal(BMSC:::handleMissingData(data, formula, yUncertainty)$data,
               data[, c("x1", "x2", "x3", "y")])
})


test_that("A message is returned if some rows are missing", {
  posMisX1 <- sample(1:nrow(data), 5)
  posMisX3 <- sample(1:nrow(data), 3)
  
  data$x1[posMisX1] <- NA
  data$x3[posMisX3] <- NA
  
  expect_equal(BMSC:::handleMissingData(data, formula, yUncertainty)$data,
               (data[-c(posMisX1, posMisX3), c("x1", "x2", "x3", "y")]))
  expect_equal(BMSC:::handleMissingData(data,
                                                           formula,
                                                           yUncertainty)$yUncertainty,
               yUncertainty[-c(posMisX1, posMisX3)])
  expect_message(BMSC:::handleMissingData(data, formula, yUncertainty)$data,
                 "7 rows \\(9.33% of the data\\) were excluded due to missing values.")
  expect_error(BMSC:::handleMissingData(data, formula, yUncertainty)$data, NA)
})


test_that("An error occures if there is no complete row in the data", {
  data$x2 <- NA
  expect_error(BMSC:::handleMissingData(data, formula, yUncertainty)$data,
               "There are no rows without missing values in the data.")
  
})
