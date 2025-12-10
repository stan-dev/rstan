# Compute summaries of MCMC draws and monitor convergence

Similar to the `print` method for `stanfit` objects, but `monitor` takes
an array of simulations as its argument rather than a `stanfit` object.
For a 3-D array (iterations \* chains \* parameters) of MCMC draws,
`monitor` computes means, standard deviations, quantiles, Monte Carlo
standard errors, split Rhats, and effective sample sizes. By default,
half of the iterations are considered warmup and are excluded.

## Usage

``` r
monitor(sims, warmup = floor(dim(sims)[1]/2), 
        probs = c(0.025, 0.25, 0.5, 0.75, 0.975), 
        digits_summary = 1, print = TRUE, ...)
# S3 method for class 'simsummary'
print(x, digits = 3, se = FALSE, ...)
# S3 method for class 'simsummary'
x[i, j, drop = if (missing(i)) TRUE else length(j) == 1]
```

## Arguments

- sims:

  A 3-D array (iterations \* chains \* parameters) of MCMC simulations
  from any MCMC algorithm.

- warmup:

  The number of warmup iterations to be excluded when computing the
  summaries. The default is half of the total number of iterations. If
  `sims` doesn't include the warmup iterations then `warmup` should be
  set to zero.

- probs:

  A numeric vector specifying quantiles of interest. The defaults is
  `c(0.025,0.25,0.5,0.75,0.975)`.

- digits_summary:

  The number of significant digits to use when printing the summary,
  defaulting to 1. Applies to the quantities other than the effective
  sample size, which is always rounded to the nearest integer.

- print:

  Logical, indicating whether to print the summary after the
  computations are performed.

- ...:

  Additional arguments passed to the underlying `print` method.

- x:

  An object of class `simsummary` created by `monitor`

- digits:

  An integer scalar defaulting to 3 for the number of decimal places to
  print

- se:

  A logical scalar defaulting to `FALSE` indicating whether to print the
  estimated standard errors of the estimates

- i:

  A vector indicating which rows of the object created by `monitor` to
  select

- j:

  A vector indicating which columns of the object crated by `monitor` to
  select

- drop:

  A logical scalar indicating whether the resulting object should return
  a vector where possible

## Value

A 2-D array with rows corresponding to parameters and columns to the
summary statistics that can be printed and subset.

## References

The Stan Development Team *Stan Modeling Language User's Guide and
Reference Manual*. <https://mc-stan.org>.

## See also

S4 class
[`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.md) and
particularly its
[`print`](https://mc-stan.org/rstan/reference/print.stanfit.md) method.

## Examples

``` r
csvfiles <- dir(system.file('misc', package = 'rstan'),
                pattern = 'rstan_doc_ex_[0-9].csv', full.names = TRUE)
fit <- read_stan_csv(csvfiles)
# The following is just for the purpose of giving an example
# since print can be used for a stanfit object.
monitor(extract(fit, permuted = FALSE, inc_warmup = TRUE))
#> Inference for the input samples (4 chains: each with iter = 200; warmup = 100):
#> 
#>           Q5   Q50   Q95  Mean  SD  Rhat Bulk_ESS Tail_ESS
#> mu      -0.3   0.1   0.4   0.1 0.2  1.01      348      229
#> sigma    0.9   1.1   1.5   1.2 0.2  1.01      227      178
#> z[1,1]  -1.6   0.0   1.4   0.0 0.9  1.01      308      301
#> z[2,1]  -1.4   0.1   1.7   0.1 1.0  1.01      329      246
#> z[3,1]  -1.9  -0.1   1.7  -0.1 1.1  1.00      455      294
#> z[1,2]  -1.7   0.0   1.5   0.0 1.0  1.01      278      138
#> z[2,2]  -1.5   0.1   1.5   0.0 0.9  1.01      438      251
#> z[3,2]  -1.6   0.1   1.9   0.1 1.0  1.02      344      124
#> alpha    0.0   0.4   1.6   0.5 0.5  1.01      354      116
#> lp__   -21.5 -17.2 -14.6 -17.5 2.3  1.03      117      208
#> 
#> For each parameter, Bulk_ESS and Tail_ESS are crude measures of 
#> effective sample size for bulk and tail quantities respectively (an ESS > 100 
#> per chain is considered good), and Rhat is the potential scale reduction 
#> factor on rank normalized split chains (at convergence, Rhat <= 1.05).
```
