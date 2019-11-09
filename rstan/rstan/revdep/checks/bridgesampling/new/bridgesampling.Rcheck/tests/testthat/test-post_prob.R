
context('post_prob with lists')


test_that("post_prob works with lists and with NAs.", {
  bridge_o <- structure(list(logml = c(4291.14352476047, 4293.29076119542,
4291.96372581169, 4293.02187182362, NA, NA, 4290.9761730488,
4293.32075269401, 4293.5762219227, 4294.02761288449), niter = c(104,
16, 52, 8, 1000, 1000, 167, 16, 21, 44), method = "normal", repetitions = 10), .Names = c("logml",
"niter", "method", "repetitions"), class = "bridge_list")

  H0L <- structure(list(logml = c(-20.8088381186739, -20.8072772698116,
-20.808454454621, -20.8083419072281, -20.8087870541247, -20.8084887398113,
-20.8086023582344, -20.8079083169745, -20.8083048489095, -20.8090050811436
), niter = c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4), method = "normal",
    repetitions = 10), .Names = c("logml", "niter", "method",
"repetitions"), class = "bridge_list")

  H1L <- structure(list(logml = c(-17.961665507006, -17.9611290723151,
-17.9607509604499, -17.9608629535992, -17.9602093576442, -17.9600223300432,
-17.9610157118017, -17.9615557696561, -17.9608437034849, -17.9606743200309
), niter = c(4, 4, 4, 4, 4, 4, 4, 4, 3, 4), method = "normal",
    repetitions = 10), .Names = c("logml", "niter", "method",
"repetitions"), class = "bridge_list")

  H0 <- structure(list(logml = -20.8084543022433, niter = 4, method = "normal"),
                .Names = c("logml", "niter", "method"), class = "bridge")

  expect_is(post_prob(H1L, H0L), "matrix")
  expect_warning(post_prob(H1L, H0L, H0), "recycled")
  expect_warning(post_prob(H1L, H0L, 4), "ignored")
  expect_warning(post_prob(H0, H0L, 4), "ignored")
  expect_warning(post_prob(H1L, H0L, bridge_o), "NA")
  expect_error(post_prob(H1L, 4, 5, 6), "one object")
  expect_error(post_prob(H0, 4, 5, 6), "one object")
})
