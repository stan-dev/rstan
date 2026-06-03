# Obtain a point estimate by maximizing the joint posterior

Obtain a point estimate by maximizing the joint posterior from the model
defined by class `stanmodel`.

## Usage

    <!-- %% optimizing(object, \dots)   -->
      # S4 method for class 'stanmodel'
    optimizing(object, data = list(),
        seed = sample.int(.Machine$integer.max, 1), init = 'random',
        check_data = TRUE, sample_file = NULL,
        algorithm = c("LBFGS", "BFGS", "Newton"),
        verbose = FALSE, hessian = FALSE, as_vector = TRUE,
        draws = 0, constrained = TRUE, importance_resampling = FALSE, ...)

## Methods

- optimizing:

  `signature(object = "stanmodel")` Call Stan's optimization methods to
  obtain a point estimate for the model defined by S4 class `stanmodel`
  given the data, initial values, etc.

## Arguments

- object:

  An object of class
  [`stanmodel`](https://mc-stan.org/rstan/dev/reference/stanmodel-class.md).

- data:

  A named `list` or `environment` providing the data for the model or a
  character vector for all the names of objects used as data. See the
  **Passing data to Stan** section in
  [`stan`](https://mc-stan.org/rstan/dev/reference/stan.md).

- seed:

  The seed for random number generation. The default is generated from 1
  to the maximum integer supported by R on the machine. Even if multiple
  chains are used, only one seed is needed, with other chains having
  seeds derived from that of the first chain to avoid dependent samples.
  When a seed is specified by a number, `as.integer` will be applied to
  it. If `as.integer` produces `NA`, the seed is generated randomly. The
  seed can also be specified as a character string of digits, such as
  `"12345"`, which is converted to integer.

- init:

  Initial values specification. See the detailed documentation for the
  `init` argument in
  [`stan`](https://mc-stan.org/rstan/dev/reference/stan.md) with one
  exception. If specifying inits using a list then only a single named
  list of values should be provided. For example, to initialize a
  parameter `alpha` to `value1` and `beta` to `value2` you can specify
  `list(alpha = value1, beta = value2)`.

- check_data:

  Logical, defaulting to `TRUE`. If `TRUE` the data will be
  preprocessed; otherwise not. See the **Passing data to Stan** section
  in [`stan`](https://mc-stan.org/rstan/dev/reference/stan.md).

- sample_file:

  A character string of file name for specifying where to write samples
  for *all* parameters and other saved quantities. If not provided,
  files are not created. When the folder specified is not writable,
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) is used.

- algorithm:

  One of `"Newton"`, `"BFGS"`, and `"LBFGS"` (the default) indicating
  which optimization algorithm to use.

- verbose:

  `TRUE` or `FALSE` (the default): flag indicating whether to print
  intermediate output from Stan on the console, which might be helpful
  for model debugging.

- hessian:

  `TRUE` or `FALSE` (the default): flag indicating whether to calculate
  the Hessian (via numeric differentiation of the gradient function in
  the unconstrained parameter space).

- as_vector:

  `TRUE` (the default) or `FALSE`: flag indicating whether a vector is
  used to store the point estimate found. A list can be used instead by
  specifying it to be `FALSE`.

- draws:

  A non-negative integer (that defaults to zero) indicating how many
  times to draw from a multivariate normal distribution whose parameters
  are the mean vector and the inverse negative Hessian in the
  unconstrained space. If `draws > 0` and `importance_resampling=TRUE`
  then `log_p` and `log_g` will be computed and returned (see
  description in the **Value** section).

- constrained:

  A logical scalar indicating, if `draws > 0`, whether the draws should
  be transformed to the constrained space defined in the parameters
  block of the Stan program. Defaults to `TRUE`.

- importance_resampling:

  A logical scalar (defaulting to `FALSE`) indicating whether to do
  importance resampling to compute diagnostics on the draws from the
  normal approximation to the posterior distribution. If `TRUE` and
  `draws > 0` then `log_p` and `log_g` will be computed and returned
  (see description in the **Value** section).

- ...:

  Other optional parameters:

  - `iter` (`integer`), the maximum number of iterations, defaulting to
    2000.

  - `save_iterations` (logical), a flag indicating whether to save the
    iterations, defaulting to `FALSE`.

  - `refresh` (`integer`), the number of interations between screen
    updates, defaulting to 100.

  - `init_alpha` (`double`), for BFGS and LBFGS, the line search step
    size for first iteration, defaulting to 0.001.

  - `tol_obj` (`double`), for BFGS and LBFGS, the convergence tolerance
    on changes in objective function value, defaulting to 1e-12.

  - `tol_rel_obj` (`double`), for BFGS and LBFGS, the convergence
    tolerance on relative changes in objective function value,
    defaulting to 1e4.

  - `tol_grad` (`double`), for BFGS and LBFGS, the convergence tolerance
    on the norm of the gradient, defaulting to 1e-8.

  - `tol_rel_grad` (`double`), for BFGS and LBFGS, the convergence
    tolerance on the relative norm of the gradient, defaulting to 1e7.

  - `tol_param` (`double`), for BFGS and LBFGS, the convergence
    tolerance on changes in parameter value, defaulting to 1e-8.

  - `history_size` (`integer`), for LBFGS, the number of update vectors
    to use in Hessian approximations, defaulting to 5.

  Refer to the manuals for both CmdStan and Stan for more details.

## Value

A list with components:

- par:

  The point estimate found. Its form (vector or list) is determined by
  the `as_vector` argument.

- value:

  The value of the log-posterior (up to an additive constant, the
  `"lp__"` in Stan) corresponding to `par`.

- return_code:

  The value of the return code from the optimizer; anything that is not
  zero is problematic.

- hessian:

  The Hessian matrix if `hessian` is `TRUE`

- theta_tilde:

  If `draws > 0`, the matrix of parameter draws in the constrained or
  unconstrained space, depending on the value of the `constrained`
  argument.

- log_p:

  If `draws > 0` and `importance_resampling=TRUE`, a vector of length
  `draws` that contains the value of the log-posterior evaluated at each
  row of `theta_tilde`.

- log_g:

  If `draws > 0`, a vector of length `draws` that contains the value of
  the logarithm of the multivariate normal density evaluated at each row
  of `theta_tilde`.

If the optimization is not completed for reasons such as feeding wrong
data, it returns `NULL`.

## See also

[`stanmodel`](https://mc-stan.org/rstan/dev/reference/stanmodel-class.md)

## Examples

``` r
# \dontrun{
m <- stan_model(model_code = 'parameters {real y;} model {y ~ normal(0,1);}')
#> recompiling to avoid crashing R session
f <- optimizing(m, hessian = TRUE)
#> Error in dimnames(x) <- dn: length of 'dimnames' [2] not equal to array extent
# }
```
