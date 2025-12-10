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
#> Chain 1: Gradient evaluation took 1.8e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.18 seconds.
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
#> Chain 1:  Elapsed Time: 0.006 seconds (Warm-up)
#> Chain 1:                0.005 seconds (Sampling)
#> Chain 1:                0.011 seconds (Total)
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
#> Chain 2:  Elapsed Time: 0.006 seconds (Warm-up)
#> Chain 2:                0.006 seconds (Sampling)
#> Chain 2:                0.012 seconds (Total)
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
#> Chain 3:  Elapsed Time: 0.006 seconds (Warm-up)
#> Chain 3:                0.005 seconds (Sampling)
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
#>            mean    se_mean        sd        10%       90%     n_eff     Rhat
#> y[1]  0.9595214 0.01819913 0.9581889  0.1054167  2.159105 2772.0503 1.000208
#> y[2]  0.9816504 0.01882405 0.9677815  0.1013469  2.232863 2643.1905 1.000009
#> lp__ -3.1484041 0.03789293 1.1233914 -4.6473975 -2.117413  878.9128 1.001904
s$c_summary  # individual chains
#> , , chains = chain:1
#> 
#>          stats
#> parameter       mean        sd        10%       90%
#>      y[1]  0.9986302 0.9678030  0.1078815  2.232374
#>      y[2]  0.9767451 0.9441381  0.1024321  2.237703
#>      lp__ -3.1136490 1.0916640 -4.6080260 -2.112247
#> 
#> , , chains = chain:2
#> 
#>          stats
#> parameter       mean        sd        10%       90%
#>      y[1]  0.9039287 0.8824154  0.1171884  2.018244
#>      y[2]  0.9613771 0.9842030  0.1082879  2.286797
#>      lp__ -3.1312855 1.1349815 -4.4844423 -2.117125
#> 
#> , , chains = chain:3
#> 
#>          stats
#> parameter       mean        sd         10%       90%
#>      y[1]  0.9667037 0.9991372  0.09625611  2.198628
#>      y[2]  1.0474490 0.9950988  0.11646897  2.343141
#>      lp__ -3.1101334 1.0805761 -4.52121870 -2.104319
#> 
#> , , chains = chain:4
#> 
#>          stats
#> parameter       mean        sd         10%       90%
#>      y[1]  0.9688230 0.9782246  0.10038614  2.197503
#>      y[2]  0.9410304 0.9447339  0.06700093  2.115310
#>      lp__ -3.2385485 1.1803231 -4.83945038 -2.149419
#> 
# }
```
