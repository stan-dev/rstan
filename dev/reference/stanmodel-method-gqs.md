# Draw samples of generated quantities from a Stan model

Draw samples from the generated quantities block of a
[`stanmodel`](https://mc-stan.org/rstan/dev/reference/stanmodel-class.md).

## Usage

``` r
# S4 method for class 'stanmodel'
gqs(object, data = list(), draws, 
    seed = sample.int(.Machine$integer.max, size = 1L))
```

## Methods

- `object`:

  `signature(object = "stanmodel")` Evaluate the generated quantities
  block of a Stan program by supplying `data` and the `draws` output
  from a previous Stan program.

## Arguments

- object:

  An object of class
  [`stanmodel`](https://mc-stan.org/rstan/dev/reference/stanmodel-class.md).

- data:

  A named `list` or `environment` providing the data for the model or a
  character vector for all the names of objects used as data. See the
  **Passing data to Stan** section in
  [`stan`](https://mc-stan.org/rstan/dev/reference/stan.md).

- draws:

  A matrix of posterior draws, typically created by calling `as.matrix`
  on a
  [`stanfit`](https://mc-stan.org/rstan/dev/reference/stanfit-class.md).

- seed:

  The seed for random number generation. The default is generated from 1
  to the maximum integer supported by R on the machine. When a seed is
  specified by a number, `as.integer` will be applied to it. If
  `as.integer` produces `NA`, the seed is generated randomly. The seed
  can also be specified as a character string of digits, such as
  `"12345"`, which is converted to integer.

## Value

An object of S4 class
[`stanmodel`](https://mc-stan.org/rstan/dev/reference/stanmodel-class.md)
representing the fitted results.

## See also

[`stanmodel`](https://mc-stan.org/rstan/dev/reference/stanmodel-class.md),
[`stanfit`](https://mc-stan.org/rstan/dev/reference/stanfit-class.md),
[`stan`](https://mc-stan.org/rstan/dev/reference/stan.md)

## Examples

``` r
# \dontrun{
m <- stan_model(model_code = 'parameters {real y;} model {y ~ normal(0,1);}')
f <- sampling(m, iter = 300)
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.02 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 300 [  0%]  (Warmup)
#> Chain 1: Iteration:  30 / 300 [ 10%]  (Warmup)
#> Chain 1: Iteration:  60 / 300 [ 20%]  (Warmup)
#> Chain 1: Iteration:  90 / 300 [ 30%]  (Warmup)
#> Chain 1: Iteration: 120 / 300 [ 40%]  (Warmup)
#> Chain 1: Iteration: 150 / 300 [ 50%]  (Warmup)
#> Chain 1: Iteration: 151 / 300 [ 50%]  (Sampling)
#> Chain 1: Iteration: 180 / 300 [ 60%]  (Sampling)
#> Chain 1: Iteration: 210 / 300 [ 70%]  (Sampling)
#> Chain 1: Iteration: 240 / 300 [ 80%]  (Sampling)
#> Chain 1: Iteration: 270 / 300 [ 90%]  (Sampling)
#> Chain 1: Iteration: 300 / 300 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0 seconds (Warm-up)
#> Chain 1:                0 seconds (Sampling)
#> Chain 1:                0 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 1e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.01 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:   1 / 300 [  0%]  (Warmup)
#> Chain 2: Iteration:  30 / 300 [ 10%]  (Warmup)
#> Chain 2: Iteration:  60 / 300 [ 20%]  (Warmup)
#> Chain 2: Iteration:  90 / 300 [ 30%]  (Warmup)
#> Chain 2: Iteration: 120 / 300 [ 40%]  (Warmup)
#> Chain 2: Iteration: 150 / 300 [ 50%]  (Warmup)
#> Chain 2: Iteration: 151 / 300 [ 50%]  (Sampling)
#> Chain 2: Iteration: 180 / 300 [ 60%]  (Sampling)
#> Chain 2: Iteration: 210 / 300 [ 70%]  (Sampling)
#> Chain 2: Iteration: 240 / 300 [ 80%]  (Sampling)
#> Chain 2: Iteration: 270 / 300 [ 90%]  (Sampling)
#> Chain 2: Iteration: 300 / 300 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0 seconds (Warm-up)
#> Chain 2:                0 seconds (Sampling)
#> Chain 2:                0 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 0 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:   1 / 300 [  0%]  (Warmup)
#> Chain 3: Iteration:  30 / 300 [ 10%]  (Warmup)
#> Chain 3: Iteration:  60 / 300 [ 20%]  (Warmup)
#> Chain 3: Iteration:  90 / 300 [ 30%]  (Warmup)
#> Chain 3: Iteration: 120 / 300 [ 40%]  (Warmup)
#> Chain 3: Iteration: 150 / 300 [ 50%]  (Warmup)
#> Chain 3: Iteration: 151 / 300 [ 50%]  (Sampling)
#> Chain 3: Iteration: 180 / 300 [ 60%]  (Sampling)
#> Chain 3: Iteration: 210 / 300 [ 70%]  (Sampling)
#> Chain 3: Iteration: 240 / 300 [ 80%]  (Sampling)
#> Chain 3: Iteration: 270 / 300 [ 90%]  (Sampling)
#> Chain 3: Iteration: 300 / 300 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0 seconds (Warm-up)
#> Chain 3:                0 seconds (Sampling)
#> Chain 3:                0 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 0 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:   1 / 300 [  0%]  (Warmup)
#> Chain 4: Iteration:  30 / 300 [ 10%]  (Warmup)
#> Chain 4: Iteration:  60 / 300 [ 20%]  (Warmup)
#> Chain 4: Iteration:  90 / 300 [ 30%]  (Warmup)
#> Chain 4: Iteration: 120 / 300 [ 40%]  (Warmup)
#> Chain 4: Iteration: 150 / 300 [ 50%]  (Warmup)
#> Chain 4: Iteration: 151 / 300 [ 50%]  (Sampling)
#> Chain 4: Iteration: 180 / 300 [ 60%]  (Sampling)
#> Chain 4: Iteration: 210 / 300 [ 70%]  (Sampling)
#> Chain 4: Iteration: 240 / 300 [ 80%]  (Sampling)
#> Chain 4: Iteration: 270 / 300 [ 90%]  (Sampling)
#> Chain 4: Iteration: 300 / 300 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 0 seconds (Warm-up)
#> Chain 4:                0 seconds (Sampling)
#> Chain 4:                0 seconds (Total)
#> Chain 4: 
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
mc <-
'
parameters {real y;}
generated quantities {real y_rep = normal_rng(y, 1);}
'
m2 <- stan_model(model_code = mc)
f2 <- gqs(m2, draws = as.matrix(f))
f2
#> Inference for Stan model: anon_model.
#> 1 chains, each with iter=600; warmup=0; thin=1; 
#> post-warmup draws per chain=600, total post-warmup draws=600.
#> 
#>       mean se_mean   sd  2.5%   25%  50%  75% 97.5% n_eff Rhat
#> y_rep 0.15    0.08 1.43 -2.57 -0.79 0.08 1.07  2.83   352    1
#> 
#> Samples were drawn using  at Sun Feb 15 14:42:09 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
# }
```
