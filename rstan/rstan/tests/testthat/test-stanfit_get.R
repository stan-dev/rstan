test_that("get_elapsed_time works", {
  skip("Backwards compatibility")

  code <- "
    parameters {
      real y;
    }
    model {
      y ~ normal(0,1);
    }
  "
  nchains <- 5L
  sm <- stan_model(model_code = code)
  fit <- sampling(sm, chains = nchains, iter = 30)
  expect_equal(dim(get_elapsed_time(fit)), c(nchains, 2L))
  nchains <- 1L
  fit2 <- sampling(sm, chains = nchains, iter = 30)
  expect_equal(dim(get_elapsed_time(fit2)), c(nchains, 2L))
})

test_that("get_adaptation_info works", {
  skip("Backwards compatibility")

  code <- "
    parameters {
      real y;
    }
    model {
      y ~ normal(0,1);
    }
  "
  fit <- stan(model_code = code, chains = 5L, iter = 30L)
  info <- get_adaptation_info(fit)
  expect_true(all(!is.na(info)))
  expect_equal(length(info), 5L)
})
