# Print a summary for a fitted model represented by a `stanfit` object

Print basic information regarding the fitted model and a summary for the
parameters of interest estimated by the samples included in a `stanfit`
object.

## Usage

``` r
# S3 method for class 'stanfit'
print(x, pars = x@sim$pars_oi, 
      probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
      digits_summary = 2, include = TRUE, ...)
```

## Arguments

- x:

  An object of S4 class `stanfit`.

- pars:

  A character vector of parameter names. The default is all parameters
  for which samples are saved. If `include = FALSE`, then the specified
  parameters are excluded from the printed summary.

- probs:

  A numeric vector of quantiles of interest. The default is
  `c(0.025,0.25,0.5,0.75,0.975)`.

- digits_summary:

  The number of significant digits to use when printing the summary,
  defaulting to 2. Applies to the quantities other than the effective
  sample size, which is always rounded to the nearest integer.

- include:

  Logical scalar (defaulting to `TRUE`) indicating whether to include or
  exclude the parameters named by the `pars` argument.

- ...:

  Additional arguments passed to the `summary` method for `stanfit`
  objects.

## Details

The information regarding the fitted model includes the number of
iterations, the number of chains, the total number of saved iterations,
the estimation algorithm used, and the timestamp indicating when
sampling finished.

The parameter summaries computed include means, standard deviations
(sd), quantiles, Monte Carlo standard errors (se_mean), split Rhats, and
effective sample sizes (n_eff). The summaries are computed after
dropping the warmup iterations and merging together the draws from all
chains.

In addition to the model parameters, summaries for the log-posterior
(`lp__`) are also reported.

## See also

S4 class
[`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.md) and
particularly its method
[`summary`](https://mc-stan.org/rstan/reference/stanfit-method-summary.md),
which is used to obtain the values that are printed.
