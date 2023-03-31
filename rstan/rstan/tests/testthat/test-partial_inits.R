test_that("partial_inits works", {
  skip("Backwards compatibility")

  stan_model_code <- "
    data {
      int<lower=0> N;
      array[N] real;
    }

    parameters {
      real mu;
      real alpha;
      simplex[3] delta;
      real eta;
      real<lower=0> sigma;
    }

    transformed parameters {
      real beta;
      beta = -alpha;
    }

    model {
      mu ~ normal(0, 10);
      y ~ normal(mu, 1);
      alpha ~ normal(0, 1);
      eta ~ normal(0, 1);
      sigma ~ exponential(2);
    }

    generated quantities {
      real gamma;
      gamma = mu + alpha;
    }
  "

  rr <- stan_model(
    model_code = stan_model_code,
    model_name = "m1",
    verbose = FALSE
  )

  y <- rnorm(20)
  dat <- list(N = length(y), y = y)
  f <- sampling(rr, data = dat, init = 0, iter = 10)
  i1 <- get_inits(f)
  # Check 1-20.
  for (i in seq_along(i1)) {
    expect_equal(i1[[i]]$mu, 0)
    expect_equal(i1[[i]]$alpha, 0)
    expect_equal(i1[[i]]$eta, 0)
    expect_equal(i1[[i]]$sigma, 1)
    expect_equal(i1[[i]]$delta, as.array(rep(1, 3) / 3))
  }

  f3 <- sampling(
    rr,
    data = dat,
    iter = 10,
    chains = 1,
    init = list(list(mu = 2)), seed = 3, thin = 1,
    enable_random_init = TRUE
  )
  i3 <- get_inits(f3)
  # Check 21-29.
  for (i in seq_along(i3)) {
    expect_equal(i3[[i]]$mu, 2)
    # This assumes that it is not possible that the following is initialized to
    # the expect randomly.
    expect_false(identical(i3[[i]]$alpha, 0))
    expect_false(identical(i3[[i]]$alpha, 1))
    expect_false(identical(i3[[i]]$eta, 0))
    expect_false(identical(i3[[i]]$eta, 1))
    expect_false(identical(i3[[i]]$sigma, 1))
    expect_false(identical(i3[[i]]$sigma, 2))
    expect_false(identical(i3[[i]]$delta, as.array(rep(1, 3) / 3)))
    expect_false(identical(i3[[i]]$delta, as.array(c(.1, .2, .7))))
  }

  # Allow partial inits, but mu is not specified correctly.
  f4 <- sampling(rr,
    data = dat, iter = 10, chains = 1,
    init = list(list(mu = rep(2, 2))), seed = 3, thin = 1,
    enable_random_init = TRUE
  )
  emsg <- geterrmessage()

  # Check 30-31.
  # checkTrue(grepl('.*mismatch.*', emsg))
  # checkTrue(grepl('.*().*(2).*', emsg))
  expect_true(grepl("Error", emsg))

  f5 <- sampling(rr,
    data = dat, iter = 10, chains = 1,
    seed = 3, thin = 1, init_r = 2,
    enable_random_init = TRUE
  )

  f6 <- sampling(rr,
    data = dat, iter = 10, chains = 1,
    seed = 3, thin = 1, init_r = 4,
    enable_random_init = TRUE
  )

  f7 <- sampling(rr,
    data = dat, iter = 10, chains = 1,
    seed = 4, thin = 1, init_r = 4,
    enable_random_init = TRUE
  )

  i5 <- get_inits(f5)[[1]]
  i6 <- get_inits(f6)[[1]]
  i7 <- get_inits(f7)[[1]]

  # Check 32-34
  expect_equal(i5$alpha * 2, i6$alpha)
  expect_equal(log(i5$sigma) * 2, log(i6$sigma))
  expect_false(identical(i7$alpha, i6$alpha))

  f8 <- sampling(rr,
    data = dat, iter = 10, chains = 1,
    init = list(list(mu = 2, sigma = 2)),
    enable_random_init = TRUE
  )
  i8 <- get_inits(f8)[[1]]
  # Check 35.
  expect_equal(i8$sigma, 2)

  f9 <- sampling(rr,
    data = dat, iter = 10, chains = 2, seed = 3,
    init = list(list(mu = 2), list(mu = 1, sigma = 2)),
    enable_random_init = TRUE
  )
  f10 <- sampling(rr,
    data = dat, iter = 10, chains = 2, seed = 3,
    init = list(list(mu = 2, eta = 1), list(mu = 1, sigma = 2, eta = 1)),
    enable_random_init = TRUE
  )
  i9 <- get_inits(f9)
  i10 <- get_inits(f10)

  # Check 36, this check is based on the implementation in which
  # though say mu is specified by the user, a random init is still
  # generated but not used.
  expect_equal(i9[[1]]$delta, i10[[1]]$delta)
})
