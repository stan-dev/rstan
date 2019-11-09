context("Polynomials")

test_that("addPowToVars", {
  vars <- c("x1", "someName")
  expect_equal(addPowToVars(vars = vars, power = 1), vars)
  expect_equal(addPowToVars(vars = vars, power = 2),
               c("I(x1^2)", "I(someName^2)"))
  expect_equal(addPowToVars(vars = vars[1], power = 2), "I(x1^2)")
})

test_that("makePoly", {
  vars <- c("just", "some_", "names")
  res1 <- makePoly(vars, maxExponent = 1)
  res2 <- makePoly(vars, maxExponent = 2)
  expect_equal(res1, vars)
  expect_equal(makePoly(vars, maxExponent = 3),
               c(makePoly(vars, maxExponent = 2),
                 "I(just^3)", "I(some_^3)", "I(names^3)"))
})
