
test_that("paridx_fun works", {
  names <- c("alpha", "beta.1", "beta.2", "lp__", "treedepth__", "stepsize__")
  paridx <- rstan:::paridx_fun(names)
  meta <- 4:6
  names(meta) <- c("lp__", "treedepth__", "stepsize__")
  exp <- c(1, 2, 3); attr(exp, 'meta') <- meta
  expect_equal(paridx, exp)
  paridx2 <- rstan:::paridx_fun(names[-5])
  meta <- 4:5; names(meta) <- c("lp__", "stepsize__")
  attr(exp, 'meta') <- meta
  expect_equal(paridx2, exp)
  paridx3 <- rstan:::paridx_fun(names[-(4:6)])
  meta <- integer(0); names(meta) <- character(0)
  attr(exp, 'meta') <- meta
  expect_equal(paridx3, exp)
})

test_that("parse_stancsv works", {
  comments <- c("# stan_version_major = 1",
                "# stan_version_minor = 3",
                "# stan_version_patch = 0",
                "# model = dogs",
                "# method = sample (Default)",
                "#   sample",
                "#     num_samples = 1000 (Default)",
                "#     num_warmup = 1000 (Default)",
                "#     save_warmup = 0 (Default)",
                "#     thin = 1 (Default)",
                "#     adapt",
                "#       engaged = 1 (Default)",
                "#       gamma = 0.050000000000000003 (Default)",
                "#       delta = 0.65000000000000002 (Default)",
                "#       kappa = 0.75 (Default)",
                "#       t0 = 10 (Default)",
                "#     algorithm = hmc (Default)",
                "#       hmc",
                "#         engine = nuts (Default)",
                "#           nuts",
                "#             max_depth = 10 (Default)",
                "#         metric = diag_e (Default)",
                "#         stepsize = 1 (Default)",
                "#         stepsize_jitter = 0 (Default)",
                "# id = 0 (Default)",
                "# data = dogs.data.R",
                "# init = 2 (Default)",
                "# random",
                "#   seed = 3086139456",
                "# output",
                "#   file = samples.csv (Default)",
                "#   append_sample = 0 (Default)",
                "#   diagnostic_file =  (Default)",
                "#   append_diagnostic = 0 (Default)",
                "#   refresh = 100 (Default)",
                "# Adaptation terminated",
                "# Step size = 0.876835",
                "# Diagonal elements of inverse mass matrix:",
                "# 0.0110854, 0.0230723",
                "#  Elapsed Time: 5.92159 seconds (Warm-up)",
                "#                5.68831 seconds (Sampling)",
                "#                11.6099 seconds (Total)")
  lst <- rstan:::parse_stancsv_comments(comments)
  expect_equal(lst$chain_id, 0)
  existence <- c("seed", "chain_id", "iter", "warmup", "thin",
                 "save_warmup", "stepsize", "time_info", "has_time",
                 "adaptation_info") %in% names(lst)
  expect_true(all(existence))
  expect_equal(lst$seed, "3086139456")
  expect_equal(lst$stepsize, 1)
  expect_equal(lst$sampler_t, "NUTS(diag_e)")
  expect_equal(lst$has_time, TRUE)
})
