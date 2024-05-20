test_that("ess works", {
  x11 <- scan(text =
  "1.35753917 -0.85304564  0.63602474  1.73601154  1.50308977 -1.24181974
   1.15386190  0.32766624 -1.29115560  1.92305046  0.64939076 -0.91140710
   1.38623164 -0.96530474 -0.31625826 -0.68192232 -0.86569785  0.96665068
   0.96252216 -0.05462611  0.34276742  1.47861666 -1.00667409 -0.59733478
   0.37617140  0.35661588 -0.26311752 -0.49101728  0.16869817 -1.29999245
   0.73691138 -0.31461659 -0.93353768  0.16547296"
  )
  x12 <- rep(1, length(x11))
  x1 <- list(x = x11, y = x12)

  x21 <- scan(text =
    "1.36448595 -0.82685756  0.63388523  1.73215725  1.50707870 -1.24387243
     1.16280682  0.34093785 -1.29084492  1.92193932  0.65552925 -0.90864304
     1.39626293 -0.95800944 -0.33344170 -0.69712612 -0.84638531  0.97138811
     0.96125515 -0.05611147  0.34584569  1.48678021 -1.00881337 -0.59505457
     0.39640764  0.38449461 -0.27761440 -0.49926682  0.16758883 -1.31532165
     0.74995372 -0.31061169 -0.93410850  0.16593008"
  )
  x2 <- list(x = x21, y = rep(1, length(x21)))

  lst <- list(samples = list(c1 = x1, c2 = x2),
              n_save = rep(34, 2),
              permutation = NULL,
              warmup2 = rep(17, 2),
              chains = 2,
              n_flatnames = 2)

  # Just make sure this can run without segfault first.
  ess <- rstan:::rstan_ess(lst, 1)
  ess <- rstan:::rstan_ess(lst, 2)

  # This is NaN because variable rho_hat_odd in function effective_sample_size
  # in ../../src/chains.cpp is NaN and then max_t=1.
  expect_true(is.nan(ess))
})

test_that("ess and rhat work", {
  c1 <- read.csv("blocker1.csv", comment.char = "#", header = TRUE)[, -(1:4)]
  # c1 <- do.call(cbind, c1)
  c2 <- read.csv("blocker2.csv", comment.char = "#", header = TRUE)[, -(1:4)]
  # c2 <- do.call(cbind, c2)
  lst <- list(samples = list(c1 = c1, c2 = c2),
              n_save = c(nrow(c1), nrow(c2)),
              permutation = NULL,
              warmup2 = rep(0, 2), chains = 2, n_flatnames = ncol(c1))
  ess <- rstan:::rstan_ess(lst, 1)
  # cat("ess=", ess, "\n")
  expect_equal(ess, 466.0988, tolerance = 0.001)
  rhat <- rstan:::rstan_splitrhat(lst, 1)
  # cat("rhat=", rhat, "\n")
  expect_equal(rhat, 1.007182, tolerance = 0.001)
  ess2 <- rstan:::rstan_ess(lst, 5)
  # cat("ess=", ess2, "\n")
  expect_equal(ess2, 518.0513, tolerance = 0.001)
  rhat2 <- rstan:::rstan_splitrhat(lst, 5)
  # cat("rhat=", rhat2, "\n")
  expect_equal(rhat2, 1.003782, tolerance = 0.001);
})

