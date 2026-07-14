# Check HMC diagnostics after sampling

These functions print summaries of important HMC diagnostics or extract
those diagnostics from a `stanfit` object. See the **Details** section,
below.

## Usage

``` r
check_hmc_diagnostics(object)
check_divergences(object)
check_treedepth(object)
check_energy(object)

get_divergent_iterations(object)
get_max_treedepth_iterations(object)
get_num_leapfrog_per_iteration(object)

get_num_divergent(object)
get_num_max_treedepth(object)

get_bfmi(object)
get_low_bfmi_chains(object)
```

## Arguments

- object:

  A stanfit object.

## Details

The `check_hmc_diagnostics` function calls the other `check_*` functions
internally and prints an overall summary, but the other functions can
also be called directly:

- `check_divergences` prints the number (and percentage) of iterations
  that ended with a divergence,

- `check_treedepth` prints the number (and percentage) of iterations
  that saturated the max treedepth,

- `check_energy` prints E-BFMI values for each chain for which E-BFMI is
  less than 0.2.

The `get_*` functions are for programmatic access to the diagnostics.

- `get_divergent_iterations` and `get_max_treedepth_iterations` return a
  logical vector indicating problems for individual iterations,

- `get_num_divergent` and `get_num_max_treedepth` return the number of
  offending interations,

- `get_num_leapfrog_per_iteration` returns an integer vector with the
  number of leapfrog evalutions for each iteration,

- `get_bfmi` returns per-chain E-BFMI values and `get_low_bfmi_chains`
  returns the indices of chains with low E-BFMI.

The following are several of many resources that provide more
information on these diagnostics:

- Brief explanations of some of the problems detected by these
  diagnostics can be found in the [*Brief Guide to Stan's
  Warnings*](https://mc-stan.org/misc/warnings.html).

- Betancourt (2017) provides much more depth on these diagnostics as
  well as a conceptual introduction to Hamiltonian Monte Carlo in
  general.

- Gabry et al. (2018) and the bayesplot package
  [vignettes](https://mc-stan.org/bayesplot/articles/) demonstrate
  various visualizations of these diagnostics that can be made in R.

## References

The Stan Development Team *Stan Modeling Language User's Guide and
Reference Manual*. <https://mc-stan.org/>.

Betancourt, M. (2017). A conceptual introduction to Hamiltonian Monte
Carlo. <https://arxiv.org/abs/1701.02434>.

Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and Gelman, A.
(2018). Visualization in Bayesian workflow. *Journal of the Royal
Statistical Society Series A*, accepted for publication. arXiv preprint:
https://arxiv.org/abs/1709.01449.

## Examples

``` r
# \dontrun{
schools <- stan_demo("eight_schools")
#> 
#> > J <- 8
#> 
#> > y <- c(28, 8, -3, 7, -1, 1, 18, 12)
#> 
#> > sigma <- c(15, 10, 16, 11, 9, 11, 10, 18)
#> 
#> > tau <- 25
#> 
#> SAMPLING FOR MODEL 'eight_schools' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 9e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.09 seconds.
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
#> Chain 1:  Elapsed Time: 0.036 seconds (Warm-up)
#> Chain 1:                0.035 seconds (Sampling)
#> Chain 1:                0.071 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'eight_schools' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.02 seconds.
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
#> Chain 2:  Elapsed Time: 0.037 seconds (Warm-up)
#> Chain 2:                0.022 seconds (Sampling)
#> Chain 2:                0.059 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'eight_schools' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 2e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.02 seconds.
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
#> Chain 3:  Elapsed Time: 0.028 seconds (Warm-up)
#> Chain 3:                0.02 seconds (Sampling)
#> Chain 3:                0.048 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'eight_schools' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 3e-06 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.03 seconds.
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
#> Chain 4:  Elapsed Time: 0.034 seconds (Warm-up)
#> Chain 4:                0.022 seconds (Sampling)
#> Chain 4:                0.056 seconds (Total)
#> Chain 4: 
#> Warning: There were 55 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
check_hmc_diagnostics(schools)
#> 
#> Divergences:
#> 55 of 4000 iterations ended with a divergence (1.375%).
#> Try increasing 'adapt_delta' to remove the divergences.
#> 
#> Tree depth:
#> 0 of 4000 iterations saturated the maximum tree depth of 10.
#> 
#> Energy:
#> E-BFMI indicated no pathological behavior.
check_divergences(schools)
#> 55 of 4000 iterations ended with a divergence (1.375%).
#> Try increasing 'adapt_delta' to remove the divergences.
check_treedepth(schools)
#> 0 of 4000 iterations saturated the maximum tree depth of 10.
check_energy(schools)
#> E-BFMI indicated no pathological behavior.
# }
```
