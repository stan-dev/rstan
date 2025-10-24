test_that("rstan monitor snapshot is stable", {
  excode <- '
  transformed data {
    array[20] real y;
    y[1] = 0.5796;  y[2] = 0.2276;   y[3]  = -0.2959;
    y[4] = -0.3742; y[5] = 0.3885;   y[6]  = -2.1585;
    y[7] = 0.7111;  y[8] = 1.4424;   y[9]  = 2.5430;
    y[10] = 0.3746; y[11] = 0.4773;  y[12] = 0.1803;
    y[13] = 0.5215; y[14] = -1.6044; y[15] = -0.6703;
    y[16] = 0.9459; y[17] = -0.382;  y[18] = 0.7619;
    y[19] = 0.1006; y[20] = -1.7461;
  }
  parameters {
    real mu;
    real<lower=0, upper=10> sigma;
    matrix[3,2] z;
    real<lower=0> alpha;
  }
  model {
    y ~ normal(mu, sigma);
    for (i in 1:3)
      z[i] ~ normal(0, 1);
    alpha ~ exponential(2);
  }
'

  exfit <- suppressWarnings(stan(
    model_code = excode,
    iter = 500,
    seed = 0,
    cores = 1
  ))

  expect_snapshot_value(
    monitor(exfit) |> as.data.frame() |> round(6),
    style = "json2"
  )
})
