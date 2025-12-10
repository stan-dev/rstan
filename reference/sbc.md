# Simulation Based Calibration (sbc)

Check whether a model is well-calibrated with respect to the prior
distribution and hence possibly amenable to obtaining a posterior
distribution conditional on observed data.

## Usage

``` r
sbc(stanmodel, data, M, ..., save_progress, load_incomplete=FALSE)
  # S3 method for class 'sbc'
plot(x, thin = 3, ...)
  # S3 method for class 'sbc'
print(x, ...)
```

## Arguments

- stanmodel:

  An object of
  [`stanmodel-class`](https://mc-stan.org/rstan/reference/stanmodel-class.md)
  that is first created by calling the
  [`stan_model`](https://mc-stan.org/rstan/reference/stan_model.md)
  function

- data:

  A named `list` or `environment` providing the data for the model, or a
  character vector for all the names of objects to use as data. This is
  the same format as in
  [`stan`](https://mc-stan.org/rstan/reference/stan.md) or
  [`sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.md).

- M:

  The number of times to condition on draws from the prior predictive
  distribution

- ...:

  Additional arguments that are passed to
  [`sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.md),
  such as `refresh = 0` when calling `sbc`. For the `plot` and `print`
  methods, the additional arguments are not used.

- x:

  An object produced by `sbc`

- thin:

  An integer vector of length one indicating the thinning interval when
  plotting, which defaults to 3

- save_progress:

  If a directory is provided, stanfit objects are saved to disk making
  it easy to resume a partial `sbc` run after interruption.

- load_incomplete:

  When `save_progress` is used, load whatever runs have been saved to
  disk and ignore argument `M`.

## Details

This function assumes adherence to the following conventions in the
underlying Stan program:

1.  Realizations of the unknown parameters are drawn in the
    `transformed data` block of the Stan program and are postfixed with
    an underscore, such as `theta_`. These are considered the “true”
    parameters being estimated by the corresponding symbol declared in
    the `parameters` block, which should have the same name except for
    the trailing underscore, such as `theta`.

2.  The realizations of the unknown parameters are then conditioned on
    when drawing from the prior predictive distribution, also in the
    `transformed data` block. There is no restriction on the symbol name
    that holds the realizations from the prior predictive distribution
    but for clarity, it should not end with a trailing underscore.

3.  The realizations of the unknown parameters should be copied into a
    `vector` in the `generated quantities` block named `pars_`.

4.  The realizations from the prior predictive distribution should be
    copied into an object (of the same type) in the
    `generated quantities` block named `y_`. Technically, this step is
    optional and could be omitted to conserve RAM, but inspecting the
    realizations from the prior predictive distribution is a good way to
    judge whether the priors are reasonable.

5.  The `generated quantities` block must contain an integer array named
    `ranks_` whose only values are zero or one, depending on whether the
    realization of a parameter from the posterior distribution exceeds
    the corresponding “true” realization, such as `theta > theta_;`.
    These are not actually "ranks" but can be used afterwards to
    reconstruct (thinned) ranks.

6.  The `generated quantities` block may contain a vector named
    `log_lik` whose values are the contribution to the log-likelihood by
    each observation. This is optional but facilitates calculating
    Pareto k shape parameters to judge whether the posterior
    distribution is sensitive to particular observations.

Although the user can pass additional arguments to
[`sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.md)
through the ..., the following arguments are hard-coded and should not
be passed through the ...:

1.  `pars = "ranks_"` because nothing else needs to be stored for each
    posterior draw

2.  `include = TRUE` to ensure that `"ranks_"` is included rather than
    excluded

3.  `chains = 1` because only one chain is run for each integer less
    than `M`

4.  `seed` because a sequence of seeds is used across the `M` runs to
    preserve independence across runs

5.  `save_warmup = FALSE` because the warmup realizations are not
    relevant

6.  `thin = 1` because thinning can and should be done after the Markov
    Chain is finished, as is done by the `thin` argument to the `plot`
    method in order to make the histograms consist of approximately
    independent realizations

Other arguments will take the default values used by
[`sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.md)
unless passed through the .... Specifying `refresh = 0` is recommended
to avoid printing a lot of intermediate progress reports to the screen.
It may be necessary to pass a list to the `control` argument of
[`sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.md)
with elements `adapt_delta` and / or `max_treedepth` in order to obtain
adequate results.

Ideally, users would want to see the absence of divergent transitions
(which is shown by the `print` method) and other warnings, plus an
approximately uniform histogram of the ranks for each parameter (which
are shown by the `plot` method). See the vignette for more details.

## Value

The `sbc` function outputs a list of S3 class `"sbc"`, which contains
the following elements:

1.  `ranks` A list of `M` matrices, each with number of rows equal to
    the number of saved iterations and number of columns equal to the
    number of unknown parameters. These matrices contain the
    realizations of the `ranks_` object from the `generated quantities`
    block of the Stan program.

2.  `Y` If present, a matrix of realizations from the prior predictive
    distribution whose rows are equal to the number of observations and
    whose columns are equal to `M`, which are taken from the `y_` object
    in the `generated quantities` block of the Stan program.

3.  `pars` A matrix of realizations from the prior distribution whose
    rows are equal to the number of parameters and whose columns are
    equal to `M`, which are taken from the `pars_` object in the
    `generated quantities` block of the Stan program.

4.  `pareto_k` A matrix of Pareto k shape parameter estimates or `NULL`
    if there is no `log_lik` symbol in the `generated quantities` block
    of the Stan program

5.  `sampler_params` A three-dimensional array that results from
    combining calls to
    [`get_sampler_params`](https://mc-stan.org/rstan/reference/stanfit-class.md)
    for each of the `M` runs. The resulting matrix has rows equal to the
    number of post-warmup iterations, columns equal to six, and `M`
    floors. The columns are named `"accept_stat__"`, `"stepsize__"`,
    `"treedepth__"`, `"n_leapfrog__"`, `"divergent__"`, and
    `"energy__"`. The most important of which is `"divergent__"`, which
    should be all zeros and perhaps `"treedepth__"`, which should only
    rarely get up to the value of `max_treedepth` passed as an element
    of the `control` list to
    [`sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.md)
    or otherwise defaults to \\10\\.

The `print` method outputs the number of divergent transitions and
returns `NULL` invisibly. The `plot` method returns a
[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html) object
with histograms whose appearance can be further customized.

## References

The Stan Development Team *Stan Modeling Language User's Guide and
Reference Manual*. <https://mc-stan.org>.

Talts, S., Betancourt, M., Simpson, D., Vehtari, A., and Gelman, A.
(2018). Validating Bayesian Inference Algorithms with Simulation-Based
Calibration. arXiv preprint arXiv:1804.06788.
<https://arxiv.org/abs/1804.06788>

## See also

[`stan_model`](https://mc-stan.org/rstan/reference/stan_model.md) and
[`sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.md)

## Examples

``` r
scode <- "
data {
  int<lower = 1> N;
  real<lower = 0> a;
  real<lower = 0> b;
}
transformed data { // these adhere to the conventions above
  real pi_ = beta_rng(a, b);
  int y = binomial_rng(N, pi_);
}
parameters {
  real<lower = 0, upper = 1> pi;
}
model {
  target += beta_lpdf(pi | a, b);
  target += binomial_lpmf(y | N, pi);
}
generated quantities { // these adhere to the conventions above
  int y_ = y;
  vector[1] pars_;
  int ranks_[1] = {pi > pi_};
  vector[N] log_lik;
  pars_[1] = pi_;
  for (n in 1:y) log_lik[n] = bernoulli_lpmf(1 | pi);
  for (n in (y + 1):N) log_lik[n] = bernoulli_lpmf(0 | pi);
}
"
```
