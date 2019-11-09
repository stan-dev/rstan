context("Adding interactions")

test_that("extractVarname", {
  x <- c("x1", "I(x2^3)", "I(x_5^2)", "var_2[3]", "I(var_2[3]^10)")
  expect_equal(extractVarname(x), c("x1", "x2", "x_5", "var_2[3]", "var_2[3]"))
})

test_that("sortAndPaste", {
  expect_equal(sortAndPaste(c("b", "a", "a")), "a:a:b")
  expect_equal(sortAndPaste(c("onylOne")), "onylOne")
})

test_that("addInteractionToVars", {
  vars1 <- c("var1", "var2")
  vars2 <- c("x_1", "I(x_1^2)", "x_2", "I(x_2^2)")
  expect_equal(addInteractionToVars(1, vars1), vars1)
  expect_equal(addInteractionToVars(2, vars1), "var1:var2")
  expect_equal(addInteractionToVars(3, vars1), character(0))
  expect_equal(addInteractionToVars(1, vars2), vars2)
  expect_equal(addInteractionToVars(2, vars2),
               c("x_1:x_2", "I(x_2^2):x_1", "I(x_1^2):x_2", "I(x_1^2):I(x_2^2)"))
  expect_equal(addInteractionToVars(3, vars2), character(0))
  expect_equal(addInteractionToVars(3, c(vars2, "oneMore")),
               c("oneMore:x_1:x_2", "I(x_2^2):oneMore:x_1",
                 "I(x_1^2):oneMore:x_2", "I(x_1^2):I(x_2^2):oneMore"))
})

test_that("makeInteractions", {
  
  vars1 <- c("var1", "var2")
  vars2 <- c("x_1", "I(x_1^2)", "x_2", "I(x_2^2)")
  vars3 <- c(vars2, "one_more")
  
  expect_equal(sort(makeInteractions(vars3, 3)),
               sort(c(addInteractionToVars(1, vars3),
                      addInteractionToVars(2, vars3),
                      addInteractionToVars(3, vars3))))
  
  expect_equal(makeInteractions(vars1, 1), vars1)
  expect_equal(sort(makeInteractions(vars1, 2)), sort(c(vars1, ("var1:var2"))))
  expect_equal(makeInteractions(vars1, 3), makeInteractions(vars1, 2))
  
  expect_equal(makeInteractions(vars2, 1), vars2)
  expect_equal(sort(makeInteractions(vars2, 2)),
               c("I(x_1^2)", "I(x_1^2):I(x_2^2)", "I(x_1^2):x_2", "I(x_2^2)",
                 "I(x_2^2):x_1", "x_1", "x_1:x_2", "x_2"))
  expect_equal(makeInteractions(vars2, 2), makeInteractions(vars2, 3))
  
  expect_equal(makeInteractions(vars3, 1), vars3)
  expect_equal(sort(makeInteractions(vars3, 3)),
               sort(c(
                 # 1-way
                 "I(x_1^2)", "I(x_2^2)", "one_more", "x_1", "x_2",
                 # 2-way
                 "I(x_1^2):I(x_2^2)", "I(x_1^2):one_more", "I(x_1^2):x_2",
                 "I(x_2^2):one_more", "I(x_2^2):x_1",
                 "one_more:x_1", "one_more:x_2",
                 "x_1:x_2",
                 # 3-way
                 "I(x_1^2):I(x_2^2):one_more", "I(x_1^2):one_more:x_2",
                 "I(x_2^2):one_more:x_1",
                 "one_more:x_1:x_2")))
  expect_equal(makeInteractions(vars3, 3), makeInteractions(vars3, 4))
})
