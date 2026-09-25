# Merge a list of stanfit objects into one

This function takes a list of `stanfit` objects and returns a
consolidated `stanfit` object. The `stanfit` objects to be merged need
to have the same configuration of iteration, warmup, and thin, besides
being from the same model. This could facilitate some parallel usage of
RStan. For example, if we call
[`stan`](https://mc-stan.org/rstan/dev/reference/stan.md) by parallel
and it returns a list of `stanfit` objects, this function can be used to
create one `stanfit` object from the list.

## Usage

``` r
sflist2stanfit(sflist)
```

## Arguments

- sflist:

  A list of `stanfit` objects.

## Value

An S4 object of `stanfit` consolidated from all the input `stanfit`
objects.

## Note

This function should be called in rare circumstances because
[`sampling`](https://mc-stan.org/rstan/dev/reference/stanmodel-method-sampling.md)
has a `cores` argument that allows multiple chains to be executed in
parallel. However, if you need to depart from that, the best practice is
to use `sflist2stanfit` on a list of `stanfit` objects created with the
same `seed` but different `chain_id` (see example below). Using the same
seed but different chain_id can make sure the random number generations
for all chains are not correlated.

This function would do some check to see if the `stanfit` objects in the
input list can be merged. But the check is not sufficient. So generally,
it is the user's responsibility to make sure the input is correct so
that the merging makes sense.

The date in the new `stanfit` object is when it is merged.

`get_seed` function for the new consolidated `stanfit` object only
returns the seed used in the first chain of the new object.

The sampler such as NUTS2 that is displayed in the printout by `print`
is the sampler used for the first chain. The `print` method assumes the
samplers are the same for all chains.

The included `stanmodel` object, which includes the compiled model, in
the new `stanfit` object is from the first element of the input list.

## References

The Stan Development Team *Stan Modeling Language User's Guide and
Reference Manual*. <https://mc-stan.org/>.

## See also

[`stan`](https://mc-stan.org/rstan/dev/reference/stan.md)

## Examples

``` r
# \dontrun{
library(rstan)
scode <- "
data {
  int<lower=1> N;
} 
parameters {
  array[N] real y1;
  array[N] real y2; 
} 
model {
  y1 ~ normal(0, 1);
  y2 ~ double_exponential(0, 2);
} 
"
seed <- 123 # or any other integer 
foo_data <- list(N = 2)
foo <- stan(model_code = scode, data = foo_data, chains = 1, iter = 1)
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 3e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.03 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: No variance estimation is
#> Chain 1:          performed for num_warmup < 20
#> Chain 1: 
#> Chain 1: Iteration: 1 / 1 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0 seconds (Warm-up)
#> Chain 1:                0 seconds (Sampling)
#> Chain 1:                0 seconds (Total)
#> Chain 1: 
f1 <- stan(fit = foo, data = foo_data, chains = 1, seed = seed, chain_id = 1) 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.02 seconds.
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
#> Chain 1:  Elapsed Time: 0.012 seconds (Warm-up)
#> Chain 1:                0.009 seconds (Sampling)
#> Chain 1:                0.021 seconds (Total)
#> Chain 1: 
f2 <- stan(fit = foo, data = foo_data, chains = 2, seed = seed, chain_id = 2:3) 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
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
#> Chain 2:  Elapsed Time: 0.012 seconds (Warm-up)
#> Chain 2:                0.009 seconds (Sampling)
#> Chain 2:                0.021 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
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
#> Chain 3:  Elapsed Time: 0.011 seconds (Warm-up)
#> Chain 3:                0.013 seconds (Sampling)
#> Chain 3:                0.024 seconds (Total)
#> Chain 3: 
f12 <- sflist2stanfit(list(f1, f2)) 

## parallel stan call for unix-like OS
library(parallel)

if (.Platform$OS.type == "unix") {
sflist1 <- 
  mclapply(1:4, mc.cores = 2, 
           function(i) stan(fit = foo, data = foo_data, seed = seed, 
                      chains = 1, chain_id = i, refresh = -1))
f3 <- sflist2stanfit(sflist1)
}
if (.Platform$OS.type == "windows") { # also works on non-Windows
CL <- makeCluster(2)
clusterExport(cl = CL, c("foo_data", "foo", "seed")) 
sflist1 <- parLapply(CL, 1:4, fun = function(cid) {  
  require(rstan)
  stan(fit = foo, data = foo_data, chains = 1, 
       iter = 2000, seed = seed, chain_id = cid)
})

fit <- sflist2stanfit(sflist1)
print(fit)
stopCluster(CL)
} # end example for Windows 
# }
```
