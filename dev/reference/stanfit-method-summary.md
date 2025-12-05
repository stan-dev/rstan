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
#> Chain 1:  Elapsed Time: 0.006 seconds (Warm-up)
#> Chain 1:                0.006 seconds (Sampling)
#> Chain 1:                0.012 seconds (Total)
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
#> Chain 2:                0.005 seconds (Sampling)
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
#>           mean    se_mean        sd        10%       90%    n_eff      Rhat
#> y[1]  1.002613 0.01835231 0.9807767  0.1136825  2.364517 2856.005 0.9997273
#> y[2]  1.021541 0.02185917 1.0419919  0.1211635  2.324770 2272.278 1.0012315
#> lp__ -3.122314 0.03506663 1.1292094 -4.5617387 -2.116804 1036.957 1.0009056
s$c_summary  # individual chains
#> , , chains = chain:1
#> 
#>          stats
#> parameter       mean        sd         10%       90%
#>      y[1]  0.9846447 0.8879088  0.14156529  2.286021
#>      y[2]  0.9794331 0.9806301  0.09451098  2.250496
#>      lp__ -3.0964163 1.1156421 -4.56092285 -2.108688
#> 
#> , , chains = chain:2
#> 
#>          stats
#> parameter      mean        sd        10%       90%
#>      y[1]  1.044494 1.0251767  0.1287719  2.503673
#>      y[2]  1.019454 0.9513040  0.1459870  2.237624
#>      lp__ -3.047046 0.9963802 -4.2953768 -2.121551
#> 
#> , , chains = chain:3
#> 
#>          stats
#> parameter       mean        sd         10%       90%
#>      y[1]  1.0008474 0.9777855  0.09455245  2.391080
#>      y[2]  0.9965311 1.0605187  0.12263733  2.252994
#>      lp__ -3.1305049 1.1071382 -4.57698837 -2.100179
#> 
#> , , chains = chain:4
#> 
#>          stats
#> parameter       mean       sd        10%       90%
#>      y[1]  0.9804651 1.025950  0.1037339  2.260858
#>      y[2]  1.0907465 1.160920  0.1076011  2.603155
#>      lp__ -3.2152872 1.275162 -4.8367761 -2.129478
#> 
# }
```
