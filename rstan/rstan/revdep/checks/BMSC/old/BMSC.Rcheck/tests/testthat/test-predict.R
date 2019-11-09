context("Predict method")

# Simulate data
suppressWarnings(RNGversion("3.5.0"))
set.seed(6245)
n <- 50
x1 <- rnorm(n)
x2 <- 2 * rnorm(n)
y <- 0.4 + 0.3 * x1 + 3 * x2 + rnorm(n, sd = 0.3)
dat <- data.frame(x1, x2, y)

datNew <- data.frame(x1 = rnorm(30),
                     x2 = rnorm(30))

# Insert NAs to into new data
datNewMis <- datNew
datNewMis$x1[sample(1:nrow(datNewMis), 5)] <- NA
datNewMis$x2[sample(1:nrow(datNewMis), 2)] <- NA

# Compute models (with intercept)
models <- constrSelEst(formula = y ~ x1 + x2,
                       data = dat, chains = 1, iterations = 500)

# Without intercept, but interactions and exponents > 1
modelsNoInt <- constrSelEst(formula = y ~ x1 + x2,
                            maxExponent = 2,
                            interactionDepth = 2,
                            intercept = FALSE,
                            data = dat,
                            scale = TRUE, chains = 1,
                            iterations = 500)


test_that("predict just works", {
  
  pred <- predict(object = models[[2]], newdata = datNew)
  expect_is(pred, "numeric")
  expect_length(pred, nrow(datNew))
  
  pred <- predict(object = modelsNoInt[[4]], newdata = dat)
  expect_is(pred, "numeric")
  expect_length(pred, nrow(dat))
  
})


test_that("Dependent variable in newdata is ignored", {
  datNewWithDep <- datNew
  datNewWithDep$y <- rnorm(30)
  expect_equal(predict(object = models[[2]], newdata = datNew),
               predict(object = models[[2]], newdata = datNewWithDep))
})


test_that("predict inserts NA values at correct positions", {
  
  pred <- predict(object = models[[length(models)]], newdata = datNewMis)
  expect_true(all_equal(which(apply(datNewMis, 1, function(x) any(is.na(x)))),
                        which(is.na(pred))))
  
  pred <- predict(object = modelsNoInt[[length(models)]], newdata = datNewMis)
  expect_true(all_equal(which(apply(datNewMis, 1, function(x) any(is.na(x)))),
                        which(is.na(pred))))
  
  pred <- predict(object = modelsNoInt[[2]], newdata = datNewMis)
  expect_true(all_equal(which(apply(datNewMis, 1, function(x) any(is.na(x)))),
                        which(is.na(pred))))
})


test_that("predict inserts NA values at correct positions (not all variables in model)", {
  
  pred <- predict(object = models[[1]], newdata = datNewMis)
  expect_true(all_equal(which(apply(datNewMis[, "x2", drop = FALSE], 1, function(x) any(is.na(x)))),
                        which(is.na(pred))))
  
  pred <- predict(object = modelsNoInt[[1]], newdata = datNewMis)
  expect_true(all_equal(which(apply(datNewMis[, "x2", drop = FALSE], 1, function(x) any(is.na(x)))),
                        which(is.na(pred))))
})
