# Draw samples from a Stan model

Draw samples from the model defined by class `stanmodel`.

## Usage

    <!-- %% sampling(object, \dots)   -->
      # S4 method for class 'stanmodel'
    sampling(object, data = list(), pars = NA,
        chains = 4, iter = 2000, warmup = floor(iter/2), thin = 1,
        seed = sample.int(.Machine$integer.max, 1),
        init = 'random', check_data = TRUE,
        sample_file = NULL, diagnostic_file = NULL, verbose = FALSE,
        algorithm = c("NUTS", "HMC", "Fixed_param"),
        control = NULL, include = TRUE,
        cores = getOption("mc.cores", 1L),
        open_progress = interactive() && !isatty(stdout()) &&
                        !identical(Sys.getenv("RSTUDIO"), "1"),
        show_messages = TRUE, ...)

## Methods

- `sampling`:

  `signature(object = "stanmodel")` Call a sampler (NUTS, HMC, or
  Fixed_param depending on parameters) to draw samples from the model
  defined by S4 class `stanmodel` given the data, initial values, etc.

## Arguments

- object:

  An object of class
  [`stanmodel`](https://mc-stan.org/rstan/dev/reference/stanmodel-class.md).

- data:

  A named `list` or `environment` providing the data for the model or a
  character vector for all the names of objects used as data. See the
  **Passing data to Stan** section in
  [`stan`](https://mc-stan.org/rstan/dev/reference/stan.md).

- pars:

  A vector of character strings specifying parameters of interest. The
  default is `NA` indicating all parameters in the model. If
  `include = TRUE`, only samples for parameters named in `pars` are
  stored in the fitted results. Conversely, if `include = FALSE`,
  samples for all parameters *except* those named in `pars` are stored
  in the fitted results.

- chains:

  A positive integer specifying the number of Markov chains. The default
  is 4.

- iter:

  A positive integer specifying the number of iterations for each chain
  (including warmup). The default is 2000.

- warmup:

  A positive integer specifying the number of warmup (aka burnin)
  iterations per chain. If step-size adaptation is on (which it is by
  default), this also controls the number of iterations for which
  adaptation is run (and hence these warmup samples should not be used
  for inference). The number of warmup iterations should be smaller than
  `iter` and the default is `iter/2`.

- thin:

  A positive integer specifying the period for saving samples. The
  default is 1, which is usually the recommended value.

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
  init argument in
  [`stan`](https://mc-stan.org/rstan/dev/reference/stan.md).

- check_data:

  Logical, defaulting to `TRUE`. If `TRUE` the data will be
  preprocessed; otherwise not. See the **Passing data to Stan** section
  in [`stan`](https://mc-stan.org/rstan/dev/reference/stan.md).

- sample_file:

  An optional character string providing the name of a file. If
  specified the draws for *all* parameters and other saved quantities
  will be written to the file. If not provided, files are not created.
  When the folder specified is not writable,
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) is used. When
  there are multiple chains, an underscore and chain number are appended
  to the file name prior to the `.csv` extension.

- diagnostic_file:

  An optional character string providing the name of a file. If
  specified the diagnostics data for *all* parameters will be written to
  the file. If not provided, files are not created. When the folder
  specified is not writable,
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) is used. When
  there are multiple chains, an underscore and chain number are appended
  to the file name prior to the `.csv` extension.

- verbose:

  `TRUE` or `FALSE`: flag indicating whether to print intermediate
  output from Stan on the console, which might be helpful for model
  debugging.

- algorithm:

  One of sampling algorithms that are implemented in Stan. Current
  options are `"NUTS"` (No-U-Turn sampler, Hoffman and Gelman 2011,
  Betancourt 2017), `"HMC"` (static HMC), or `"Fixed_param"`. The
  default and preferred algorithm is `"NUTS"`.

- control:

  A named `list` of parameters to control the sampler's behavior. See
  the details in the documentation for the `control` argument in
  [`stan`](https://mc-stan.org/rstan/dev/reference/stan.md).

- include:

  Logical scalar defaulting to `TRUE` indicating whether to include or
  exclude the parameters given by the `pars` argument. If `FALSE`, only
  entire multidimensional parameters can be excluded, rather than
  particular elements of them.

- cores:

  Number of cores to use when executing the chains in parallel, which
  defaults to 1 but we recommend setting the `mc.cores` option to be as
  many processors as the hardware and RAM allow (up to the number of
  chains).

- open_progress:

  Logical scalar that only takes effect if `cores > 1` but is
  recommended to be `TRUE` in interactive use so that the progress of
  the chains will be redirected to a file that is automatically opened
  for inspection. For very short runs, the user might prefer `FALSE`.

- show_messages:

  Either a logical scalar (defaulting to `TRUE`) indicating whether to
  print the summary of Informational Messages to the screen after a
  chain is finished or a character string naming a path where the
  summary is stored. Setting to `FALSE` is not recommended unless you
  are very sure that the model is correct up to numerical error.

- ...:

  Additional arguments can be `chain_id`, `init_r`, `test_grad`,
  `append_samples`, `refresh`, `enable_random_init`. See the
  documentation in
  [`stan`](https://mc-stan.org/rstan/dev/reference/stan.md).

## Value

An object of S4 class `stanfit` representing the fitted results. Slot
`mode` for this object indicates if the sampling is done or not.

## See also

[`stanmodel`](https://mc-stan.org/rstan/dev/reference/stanmodel-class.md),
[`stanfit`](https://mc-stan.org/rstan/dev/reference/stanfit-class.md),
[`stan`](https://mc-stan.org/rstan/dev/reference/stan.md)

## Examples

``` r
# \dontrun{
m <- stan_model(model_code = 'parameters {real y;} model {y ~ normal(0,1);}')
#> recompiling to avoid crashing R session
f <- sampling(m, iter = 100)
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.02 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: There aren't enough warmup iterations to fit the
#> Chain 1:          three stages of adaptation as currently configured.
#> Chain 1:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 1:          the given number of warmup iterations:
#> Chain 1:            init_buffer = 7
#> Chain 1:            adapt_window = 38
#> Chain 1:            term_buffer = 5
#> Chain 1: 
#> Chain 1: Iteration:  1 / 100 [  1%]  (Warmup)
#> Chain 1: Iteration: 10 / 100 [ 10%]  (Warmup)
#> Chain 1: Iteration: 20 / 100 [ 20%]  (Warmup)
#> Chain 1: Iteration: 30 / 100 [ 30%]  (Warmup)
#> Chain 1: Iteration: 40 / 100 [ 40%]  (Warmup)
#> Chain 1: Iteration: 50 / 100 [ 50%]  (Warmup)
#> Chain 1: Iteration: 51 / 100 [ 51%]  (Sampling)
#> Chain 1: Iteration: 60 / 100 [ 60%]  (Sampling)
#> Chain 1: Iteration: 70 / 100 [ 70%]  (Sampling)
#> Chain 1: Iteration: 80 / 100 [ 80%]  (Sampling)
#> Chain 1: Iteration: 90 / 100 [ 90%]  (Sampling)
#> Chain 1: Iteration: 100 / 100 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0 seconds (Warm-up)
#> Chain 1:                0 seconds (Sampling)
#> Chain 1:                0 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: WARNING: There aren't enough warmup iterations to fit the
#> Chain 2:          three stages of adaptation as currently configured.
#> Chain 2:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 2:          the given number of warmup iterations:
#> Chain 2:            init_buffer = 7
#> Chain 2:            adapt_window = 38
#> Chain 2:            term_buffer = 5
#> Chain 2: 
#> Chain 2: Iteration:  1 / 100 [  1%]  (Warmup)
#> Chain 2: Iteration: 10 / 100 [ 10%]  (Warmup)
#> Chain 2: Iteration: 20 / 100 [ 20%]  (Warmup)
#> Chain 2: Iteration: 30 / 100 [ 30%]  (Warmup)
#> Chain 2: Iteration: 40 / 100 [ 40%]  (Warmup)
#> Chain 2: Iteration: 50 / 100 [ 50%]  (Warmup)
#> Chain 2: Iteration: 51 / 100 [ 51%]  (Sampling)
#> Chain 2: Iteration: 60 / 100 [ 60%]  (Sampling)
#> Chain 2: Iteration: 70 / 100 [ 70%]  (Sampling)
#> Chain 2: Iteration: 80 / 100 [ 80%]  (Sampling)
#> Chain 2: Iteration: 90 / 100 [ 90%]  (Sampling)
#> Chain 2: Iteration: 100 / 100 [100%]  (Sampling)
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
#> Chain 3: WARNING: There aren't enough warmup iterations to fit the
#> Chain 3:          three stages of adaptation as currently configured.
#> Chain 3:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 3:          the given number of warmup iterations:
#> Chain 3:            init_buffer = 7
#> Chain 3:            adapt_window = 38
#> Chain 3:            term_buffer = 5
#> Chain 3: 
#> Chain 3: Iteration:  1 / 100 [  1%]  (Warmup)
#> Chain 3: Iteration: 10 / 100 [ 10%]  (Warmup)
#> Chain 3: Iteration: 20 / 100 [ 20%]  (Warmup)
#> Chain 3: Iteration: 30 / 100 [ 30%]  (Warmup)
#> Chain 3: Iteration: 40 / 100 [ 40%]  (Warmup)
#> Chain 3: Iteration: 50 / 100 [ 50%]  (Warmup)
#> Chain 3: Iteration: 51 / 100 [ 51%]  (Sampling)
#> Chain 3: Iteration: 60 / 100 [ 60%]  (Sampling)
#> Chain 3: Iteration: 70 / 100 [ 70%]  (Sampling)
#> Chain 3: Iteration: 80 / 100 [ 80%]  (Sampling)
#> Chain 3: Iteration: 90 / 100 [ 90%]  (Sampling)
#> Chain 3: Iteration: 100 / 100 [100%]  (Sampling)
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
#> Chain 4: WARNING: There aren't enough warmup iterations to fit the
#> Chain 4:          three stages of adaptation as currently configured.
#> Chain 4:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 4:          the given number of warmup iterations:
#> Chain 4:            init_buffer = 7
#> Chain 4:            adapt_window = 38
#> Chain 4:            term_buffer = 5
#> Chain 4: 
#> Chain 4: Iteration:  1 / 100 [  1%]  (Warmup)
#> Chain 4: Iteration: 10 / 100 [ 10%]  (Warmup)
#> Chain 4: Iteration: 20 / 100 [ 20%]  (Warmup)
#> Chain 4: Iteration: 30 / 100 [ 30%]  (Warmup)
#> Chain 4: Iteration: 40 / 100 [ 40%]  (Warmup)
#> Chain 4: Iteration: 50 / 100 [ 50%]  (Warmup)
#> Chain 4: Iteration: 51 / 100 [ 51%]  (Sampling)
#> Chain 4: Iteration: 60 / 100 [ 60%]  (Sampling)
#> Chain 4: Iteration: 70 / 100 [ 70%]  (Sampling)
#> Chain 4: Iteration: 80 / 100 [ 80%]  (Sampling)
#> Chain 4: Iteration: 90 / 100 [ 90%]  (Sampling)
#> Chain 4: Iteration: 100 / 100 [100%]  (Sampling)
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
# }
```
