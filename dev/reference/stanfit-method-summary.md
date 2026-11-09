# Summary method for stanfit objects

Summarize the distributions of estimated parameters and derived
quantities using the posterior draws.

## Usage

``` r
# S4 method for class 'stanfit'
summary(object, pars, probs = c(0.025, 0.25, 0.50, 0.75, 0.975),
  use_cache = TRUE, ...)
```

## Arguments

- object:

  An instance of class
  [`stanfit`](https://mc-stan.org/rstan/dev/reference/stanfit-class.md).

- pars:

  A character vector of parameter names. Defaults to all parameters as
  well as the log-posterior (`lp__`).

- probs:

  A numeric vector of
  [`quantile`](https://rdrr.io/r/stats/quantile.html)s of interest. The
  default is `c(0.025,0.25,0.5,0.75,0.975)`.

- use_cache:

  Logical, defaulting to `TRUE`. When `use_cache=TRUE` the summary
  quantities for all parameters are computed and cached for future use.
  Setting `use_cache=FALSE` can be used to avoid performing the summary
  computations for all parameters if `pars` is given as some specific
  parameters.

- ...:

  Currently unused.

## Value

The `summary` method returns a named list with elements `summary` and
`c_summary`, which contain summaries for for all chains merged and
individual chains, respectively. Included in the summaries are
quantiles, means, standard deviations (sd), effective sample sizes
(n_eff), and split Rhats (the potential scale reduction derived from all
chains after splitting each chain in half and treating the halves as
chains). For the summary of all chains merged, Monte Carlo standard
errors (se_mean) are also reported.

## See also

- [`monitor`](https://mc-stan.org/rstan/dev/reference/monitor.md), which
  computes similar summaries but accepts an array of MCMC draws as its
  input rather than a `stanfit` object.

- The RStan vignettes for more example usage.

## Examples

``` r
# \dontrun{
ecode <- '
  parameters {
    array[2] real<lower=0> y;
  } 
  model {
    y ~ exponential(1);
  }
'
fit <- stan(model_code = ecode)
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 5e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.05 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.005 seconds (Warm-up)
#> Chain 1:                0.005 seconds (Sampling)
#> Chain 1:                0.01 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 1e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.01 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.005 seconds (Warm-up)
#> Chain 2:                0.005 seconds (Sampling)
#> Chain 2:                0.01 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 1e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.01 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.005 seconds (Warm-up)
#> Chain 3:                0.006 seconds (Sampling)
#> Chain 3:                0.011 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 1e-06 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.01 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 0.006 seconds (Warm-up)
#> Chain 4:                0.005 seconds (Sampling)
#> Chain 4:                0.011 seconds (Total)
#> Chain 4: 
s <- summary(fit, probs = c(0.1, 0.9))
s$summary  # all chaines merged
#>            mean    se_mean        sd        10%       90%     n_eff      Rhat
#> y[1]  1.0147886 0.02021218 1.0144742  0.1057559  2.355153 2519.1611 1.0001895
#> y[2]  0.9780169 0.01959189 0.9753713  0.1125876  2.239543 2478.4896 0.9994361
#> lp__ -3.1545421 0.03678723 1.1550466 -4.6083260 -2.117495  985.8365 1.0005713
s$c_summary  # individual chains
#> , , chains = chain:1
#> 
#>          stats
#> parameter      mean        sd        10%       90%
#>      y[1]  1.028724 0.9628692  0.1191143  2.380378
#>      y[2]  0.947648 0.9555451  0.1080397  2.127725
#>      lp__ -3.086648 0.9806253 -4.4762992 -2.109054
#> 
#> , , chains = chain:2
#> 
#>          stats
#> parameter       mean       sd        10%       90%
#>      y[1]  0.9984413 1.003530  0.1066804  2.278915
#>      y[2]  0.9747653 0.971302  0.1209582  2.263833
#>      lp__ -3.1312997 1.207422 -4.6207371 -2.109312
#> 
#> , , chains = chain:3
#> 
#>          stats
#> parameter      mean       sd         10%       90%
#>      y[1]  0.986253 1.054671  0.09536197  2.304043
#>      y[2]  0.996872 0.938427  0.11428575  2.266265
#>      lp__ -3.179568 1.166695 -4.68105504 -2.142655
#> 
#> , , chains = chain:4
#> 
#>          stats
#> parameter       mean       sd         10%       90%
#>      y[1]  1.0457358 1.034874  0.09856248  2.424349
#>      y[2]  0.9927822 1.034218  0.10953289  2.265077
#>      lp__ -3.2206526 1.244873 -4.70781301 -2.128995
#> 
# }
```
