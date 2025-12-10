# Set and read options used in RStan

Set and read options used in RStan. Some settings as options can be
controlled by the user.

## Usage

``` r
rstan_options(...)
```

## Arguments

- ...:

  Arguments of the form `opt = val` set option `opt` to value `val`.
  Arguments of the form `opt` set the function to return option `opt`'s
  value. Each argument must be a character string.

## Details

The available options are:

1.  `plot_rhat_breaks`: The cut off points for Rhat for which we would
    indicate using a different color. This is a numeric vector,
    defaulting to `c(1.1, 1.2, 1.5, 2)`. The value for this option will
    be sorted in ascending order, so for example
    `plot_rhat_breaks = c(1.2, 1.5)` is equivalent to
    `plot_rhat_breaks = c(1.5, 1.2)`.

2.  `plot_rhat_cols`: A vector of the same length as `plot_rhat_breaks`
    that indicates the colors for the breaks.

3.  `plot_rhat_nan_col`: The color for Rhat when it is `Inf` or `NaN`.

4.  `plot_rhat_large_col`: The color for Rhat when it is larger than the
    largest value of `plot_rhat_breaks`.

5.  `rstan_alert_col`: The color used in method `plot` of S4 class
    [`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.md) to
    show that the vector/array parameters are truncated.

6.  `rstan_chain_cols`: The colors used in methods `plot` and
    `traceplot` of S4 class
    [`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.md)
    for coloring different chains.

7.  `rstan_warmup_bg_col`: The background color for the warmup area in
    the traceplots.

8.  `boost_lib`: The path for the Boost C++ library used to compile Stan
    models. This option is valid for the whole R session if not changed
    again.

9.  `eigen_lib`: The path for the Eigen C++ library used to compile Stan
    models. This option is valid for the whole R session if not changed
    again.

10. `auto_write`: A logical scalar (defaulting to `FALSE`) that controls
    whether a compiled instance of a
    [`stanmodel-class`](https://mc-stan.org/rstan/reference/stanmodel-class.md)
    is written to the hard disk in the same directory as the `.stan`
    program.

11. `threads_per_chain`: A positive integer (defaulting to `1`). If the
    model was compiled with threading support, the number of threads to
    use in parallelized sections \_within\_ an MCMC chain (e.g., when
    using the Stan functions \`reduce_sum()\` or \`map_rect()\`). The
    actual number of CPU cores used is \`chains \* threads_per_chain\`
    where \`chains\` is the number of parallel chains. For an example of
    using threading, see the Stan case study \[Reduce Sum: A Minimal
    Example\](https://mc-stan.org/users/documentation/case-studies/reduce_sum_tutorial.html).

## Value

The values as a `list` for existing options and `NA` for non-existent
options. When only one option is specified, its old value is returned.
