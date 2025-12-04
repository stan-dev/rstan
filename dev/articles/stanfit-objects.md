# Accessing the contents of a stanfit object

This vignette demonstrates how to access most of data stored in a
stanfit object. A stanfit object (an object of class `"stanfit"`)
contains the output derived from fitting a Stan model using Markov chain
Monte Carlo or one of Stan’s variational approximations (meanfield or
full-rank). Throughout the document we’ll use the stanfit object
obtained from fitting the Eight Schools example model:

``` r
library(rstan)
fit <- stan_demo("eight_schools", refresh = 0)
```

    Error in stanc(file = file, model_code = model_code, model_name = model_name, : 0
    Syntax error in 'string', line 3, column 9 to column 10, parsing error:
       -------------------------------------------------
         1:  data {
         2:    int<lower=0> J;          // number of schools
         3:    real y[J];               // estimated treatment effects
                      ^
         4:    real<lower=0> sigma[J];  // s.e. of effect estimates
         5:  }
       -------------------------------------------------

    Ill-formed declaration. ";" expected after variable declaration.
    It looks like you are trying to use the old array syntax.
    Please use the new syntax:
    array[J] real y;

``` r
class(fit)
```

    Error: object 'fit' not found

## Posterior draws

There are several functions that can be used to access the draws from
the posterior distribution stored in a stanfit object. These are
`extract`, `as.matrix`, `as.data.frame`, and `as.array`, each of which
returns the draws in a different format.

  

#### extract()

The `extract` function (with its default arguments) returns a list with
named components corresponding to the model parameters.

``` r
list_of_draws <- extract(fit)
```

    Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'extract': object 'fit' not found

``` r
print(names(list_of_draws))
```

    Error: object 'list_of_draws' not found

In this model the parameters `mu` and `tau` are scalars and `theta` is a
vector with eight elements. This means that the draws for `mu` and `tau`
will be vectors (with length equal to the number of post-warmup
iterations times the number of chains) and the draws for `theta` will be
a matrix, with each column corresponding to one of the eight components:

``` r
head(list_of_draws$mu)
```

    Error: object 'list_of_draws' not found

``` r
head(list_of_draws$tau)
```

    Error: object 'list_of_draws' not found

``` r
head(list_of_draws$theta)
```

    Error: object 'list_of_draws' not found

  

#### as.matrix(), as.data.frame(), as.array()

The `as.matrix`, `as.data.frame`, and `as.array` functions can also be
used to retrieve the posterior draws from a stanfit object:

``` r
matrix_of_draws <- as.matrix(fit)
```

    Error: object 'fit' not found

``` r
print(colnames(matrix_of_draws))
```

    Error: object 'matrix_of_draws' not found

``` r
df_of_draws <- as.data.frame(fit)
```

    Error: object 'fit' not found

``` r
print(colnames(df_of_draws))
```

    Error: object 'df_of_draws' not found

``` r
array_of_draws <- as.array(fit)
```

    Error: object 'fit' not found

``` r
print(dimnames(array_of_draws))
```

    Error: object 'array_of_draws' not found

The `as.matrix` and `as.data.frame` methods essentially return the same
thing except in matrix and data frame form, respectively. The `as.array`
method returns the draws from each chain separately and so has an
additional dimension:

``` r
print(dim(matrix_of_draws))
```

    Error: object 'matrix_of_draws' not found

``` r
print(dim(df_of_draws))
```

    Error: object 'df_of_draws' not found

``` r
print(dim(array_of_draws))
```

    Error: object 'array_of_draws' not found

By default all of the functions for retrieving the posterior draws
return the draws for *all* parameters (and generated quantities). The
optional argument `pars` (a character vector) can be used if only a
subset of the parameters is desired, for example:

``` r
mu_and_theta1 <- as.matrix(fit, pars = c("mu", "theta[1]"))
```

    Error: object 'fit' not found

``` r
head(mu_and_theta1)
```

    Error: object 'mu_and_theta1' not found

  

## Posterior summary statistics and convergence diagnostics

Summary statistics are obtained using the `summary` function. The object
returned is a list with two components:

``` r
fit_summary <- summary(fit)
```

    Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'summary': object 'fit' not found

``` r
print(names(fit_summary))
```

    Error: object 'fit_summary' not found

In `fit_summary$summary` all chains are merged whereas
`fit_summary$c_summary` contains summaries for each chain individually.
Typically we want the summary for all chains merged, which is what we’ll
focus on here.

The summary is a matrix with rows corresponding to parameters and
columns to the various summary quantities. These include the posterior
mean, the posterior standard deviation, and various quantiles computed
from the draws. The `probs` argument can be used to specify which
quantiles to compute and `pars` can be used to specify a subset of
parameters to include in the summary.

For models fit using MCMC, also included in the summary are the Monte
Carlo standard error (`se_mean`), the effective sample size (`n_eff`),
and the R-hat statistic (`Rhat`).

``` r
print(fit_summary$summary)
```

    Error: object 'fit_summary' not found

If, for example, we wanted the only quantiles included to be 10% and
90%, and for only the parameters included to be `mu` and `tau`, we would
specify that like this:

``` r
mu_tau_summary <- summary(fit, pars = c("mu", "tau"), probs = c(0.1, 0.9))$summary
```

    Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'summary': object 'fit' not found

``` r
print(mu_tau_summary)
```

    Error: object 'mu_tau_summary' not found

Since `mu_tau_summary` is a matrix we can pull out columns using their
names:

``` r
mu_tau_80pct <- mu_tau_summary[, c("10%", "90%")]
```

    Error: object 'mu_tau_summary' not found

``` r
print(mu_tau_80pct)
```

    Error: object 'mu_tau_80pct' not found

  

## Sampler diagnostics

For models fit using MCMC the stanfit object will also contain the
values of parameters used for the sampler. The `get_sampler_params`
function can be used to access this information.

The object returned by `get_sampler_params` is a list with one component
(a matrix) per chain. Each of the matrices has number of columns
corresponding to the number of sampler parameters and the column names
provide the parameter names. The optional argument inc_warmup
(defaulting to `TRUE`) indicates whether to include the warmup period.

``` r
sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
```

    Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'get_sampler_params': object 'fit' not found

``` r
sampler_params_chain1 <- sampler_params[[1]]
```

    Error: object 'sampler_params' not found

``` r
colnames(sampler_params_chain1)
```

    Error: object 'sampler_params_chain1' not found

To do things like calculate the average value of `accept_stat__` for
each chain (or the maximum value of `treedepth__` for each chain if
using the NUTS algorithm, etc.) the `sapply` function is useful as it
will apply the same function to each component of `sampler_params`:

``` r
mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
```

    Error: object 'sampler_params' not found

``` r
print(mean_accept_stat_by_chain)
```

    Error: object 'mean_accept_stat_by_chain' not found

``` r
max_treedepth_by_chain <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
```

    Error: object 'sampler_params' not found

``` r
print(max_treedepth_by_chain)
```

    Error: object 'max_treedepth_by_chain' not found

  

## Model code

The Stan program itself is also stored in the stanfit object and can be
accessed using `get_stancode`:

``` r
code <- get_stancode(fit)
```

    Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'get_stancode': object 'fit' not found

The object `code` is a single string and is not very intelligible when
printed:

``` r
print(code)
```

    Error: object 'code' not found

A readable version can be printed using `cat`:

``` r
cat(code)
```

    Error: object 'code' not found

  

## Initial values

The `get_inits` function returns initial values as a list with one
component per chain. Each component is itself a (named) list containing
the initial values for each parameter for the corresponding chain:

``` r
inits <- get_inits(fit)
```

    Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'get_inits': object 'fit' not found

``` r
inits_chain1 <- inits[[1]]
```

    Error: object 'inits' not found

``` r
print(inits_chain1)
```

    Error: object 'inits_chain1' not found

  

## (P)RNG seed

The `get_seed` function returns the (P)RNG seed as an integer:

``` r
print(get_seed(fit))
```

    Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'get_seed': object 'fit' not found

  

## Warmup and sampling times

The `get_elapsed_time` function returns a matrix with the warmup and
sampling times for each chain:

``` r
print(get_elapsed_time(fit))
```

    Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'get_elapsed_time': object 'fit' not found
