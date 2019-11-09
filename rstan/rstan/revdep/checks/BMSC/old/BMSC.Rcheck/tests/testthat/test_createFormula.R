context("Formula creation")

test_that("createFormula returns correct input - usual cases", {
  input1 <- as.formula("dep ~ x1 + x2 + x3")
  input2 <- "dep ~ 1 + x1 + x2"
  input3 <- "dep ~ 1 + x1^2 + x2*x3 + I(x2^4)"
  input4 <- "dep ~ 1 + log(x1) + sqrt(x2)"
  expect_equal(createFormula(input1), input1)
  expect_equal(createFormula(input2), as.formula("dep ~ x1 + x2"))
  expect_equal(createFormula(input3, 2, 2), createFormula(input1, 2, 2))
  expect_equal(createFormula(input4), createFormula(input2))
})

test_that("createFormula returns correct input - no or only intercept", {
  input1 <- as.formula("dep ~ 1") # Intercept only
  input2 <- "dep ~ -1 + x1" # No intercept
  input3 <- "dep ~ 0 + x1" # No intercept
  input4 <- "dep ~ x1 + x2 - 1" # No intercept
  expect_equal(createFormula(input1), input1)
  expect_equal(createFormula(input2, 2, 1), as.formula("dep ~ x1 + I(x1^2) - 1"))
  expect_equal(createFormula(input3, 1, 2), as.formula("dep ~ x1 - 1"))
  expect_equal(createFormula(input4, 2, 3),
               as.formula(paste0("dep ~ x1 + x2 + I(x1^2) + I(x2^2) + x1:x2 + ",
                                 "I(x2^2):x1 + I(x1^2):x2 + I(x1^2):I(x2^2) ",
                                 "- 1")))
})

test_that("createFormula/tryAsFormula returns correct errors", {
  expect_error(createFormula("someNonsense"),
               '"someNonsense" cannot be turned into a formula.')
  expect_error(tryAsFormula("some Stuff"),
               '"some Stuff" cannot be turned into a formula.')
  expect_error(createFormula("y ~ x1 + x2", maxExponent = 0),
               "maxExponent and interactionDepth must be positive integers.")
  expect_error(createFormula("y ~ x1 + x2", interactionDepth = -1),
               "maxExponent and interactionDepth must be positive integers.")
  expect_error(createFormula("~ x1 + x2"), "Formula is not complete.")
  expect_error(createFormula("~ something"), "Formula is not complete.")
  expect_error(createFormula("Y~x"), NA)
})

test_that("computing models with createFormula", {
  mod1 <- lm(createFormula("mpg ~ cyl + disp",
                           maxExponent = 2,
                           interactionDepth = 2), mtcars)
  expect_equal(names(mod1$coefficients),
               c("(Intercept)", "cyl", "disp", "I(cyl^2)", "I(disp^2)",
                 "cyl:disp", "cyl:I(disp^2)", "disp:I(cyl^2)",
                 "I(cyl^2):I(disp^2)"))
  
  mod2 <- lm(createFormula("mpg ~ cyl + disp - 1",
                           maxExponent = 2,
                           interactionDepth = 2), mtcars)
  expect_true(all(!grepl("intercept", names(mod2$coefficients)))) # No intercept
})
