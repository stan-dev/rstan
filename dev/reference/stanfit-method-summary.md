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
#> Chain 1: Gradient evaluation took 3e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.03 seconds.
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
#> Chain 2:                0.006 seconds (Sampling)
#> Chain 2:                0.011 seconds (Total)
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
#> Chain 4:  Elapsed Time: 0.005 seconds (Warm-up)
#> Chain 4:                0.005 seconds (Sampling)
#> Chain 4:                0.01 seconds (Total)
#> Chain 4: 
s <- summary(fit, probs = c(0.1, 0.9))
s$summary  # all chaines merged
#>            mean    se_mean        sd         10%       90%     n_eff     Rhat
#> y[1]  0.9993487 0.02079046 1.0159236  0.10042471  2.330060 2387.7791 1.000669
#> y[2]  0.9685756 0.02036467 0.9653798  0.09135474  2.225154 2247.2001 1.000072
#> lp__ -3.2001021 0.03787490 1.1664838 -4.74101666 -2.131617  948.5376 1.002617
s$c_summary  # individual chains
#> , , chains = chain:1
#> 
#>          stats
#> parameter       mean        sd        10%       90%
#>      y[1]  1.0694339 1.0612259  0.0953882  2.499088
#>      y[2]  0.9829245 0.9651658  0.0706803  2.308414
#>      lp__ -3.2950690 1.2429727 -4.9694643 -2.150289
#> 
#> , , chains = chain:2
#> 
#>          stats
#> parameter       mean        sd         10%       90%
#>      y[1]  0.9308465 0.9684359  0.07974216  2.204396
#>      y[2]  0.9072754 0.9450972  0.08752071  2.160447
#>      lp__ -3.2605769 1.2392431 -4.83487860 -2.127506
#> 
#> , , chains = chain:3
#> 
#>          stats
#> parameter       mean        sd        10%       90%
#>      y[1]  0.9809423 1.0068234  0.1058002  2.310432
#>      y[2]  0.9677877 0.9153224  0.1031050  2.126079
#>      lp__ -3.1030783 1.0605475 -4.4806192 -2.114384
#> 
#> , , chains = chain:4
#> 
#>          stats
#> parameter      mean       sd        10%       90%
#>      y[1]  1.016172 1.021521  0.1174006  2.347030
#>      y[2]  1.016315 1.030443  0.1001645  2.259636
#>      lp__ -3.141684 1.102647 -4.6362905 -2.128575
#> 
# }
```
