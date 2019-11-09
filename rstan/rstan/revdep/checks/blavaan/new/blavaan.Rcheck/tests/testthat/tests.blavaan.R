test_that("blavaan arguments", {
  x1 <- rnorm(100)
  x2 <- rnorm(100)
  y1 <- 0.5 + 2*x1 + rnorm(100)
  Data <- data.frame(y1 = y1, x1 = x1, x2 = x2)

  model <- ' y1 ~ x1 '

  ## auto convergence in stan
  expect_error(bsem(model, data=Data, fixed.x=TRUE, target="stan", convergence="auto"))

  ## seed length != # chains
  expect_error(bsem(model, data=Data, fixed.x=TRUE, seed=1))

  ## supply ordinals
  expect_error(bsem(model, data=Data, fixed.x=TRUE, ordered=c("y1", "x1", adapt=2, burnin=2, sample=2)))

  ## unknown cp
  expect_error(bsem(model, data=Data, ov.cp="blah", fixed.x=TRUE))

  ## cp/std.lv clash
  expect_error(bsem(model, data=Data, fixed.x=TRUE, std.lv=TRUE, cp="fa"))
  
  model2 <- ' y1 ~ b1*x1 + b2*x2
              b1 + b2 == 0 '

  ## equality constraint with multiple variables on lhs
  expect_error(bsem(model2, data=Data, fixed.x=TRUE))

  ## do.fit=FALSE
  fit <- bsem(model, data=Data, fixed.x=TRUE, adapt=2,
              burnin=2, sample=2, do.fit=FALSE)
  expect_s4_class(fit, "blavaan")

  ## named variable that clashes
  names(Data)[1] <- "lambda"
  model2 <- ' lambda ~ b1*x1 + b2*x2 '
  expect_error(bsem(model2, data=Data))
              
})
