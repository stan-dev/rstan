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

    Trying to compile a simple C file

    Warning: There were 2 divergent transitions after warmup. See
    https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
    to find out why this is a problem and how to eliminate them.

    Warning: Examine the pairs() plot to diagnose sampling problems

``` r

class(fit)
```

    [1] "stanfit"
    attr(,"package")
    [1] "rstan"

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
print(names(list_of_draws))
```

    [1] "mu"    "tau"   "eta"   "theta" "lp__" 

In this model the parameters `mu` and `tau` are scalars and `theta` is a
vector with eight elements. This means that the draws for `mu` and `tau`
will be vectors (with length equal to the number of post-warmup
iterations times the number of chains) and the draws for `theta` will be
a matrix, with each column corresponding to one of the eight components:

``` r

head(list_of_draws$mu)
```

    [1]  6.803688 18.797813 14.921268  8.115805  8.610047  8.347953

``` r

head(list_of_draws$tau)
```

    [1] 0.2784982 8.8754590 0.2824127 4.2320069 7.1922227 4.5074509

``` r

head(list_of_draws$theta)
```

              
    iterations      [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
          [1,]  6.952096  6.943472  6.528042  6.934846  6.844753  6.871275
          [2,] 14.014660 15.231572 10.702913 20.790946  4.815897 19.034372
          [3,] 15.154398 15.365242 15.144703 15.364579 14.613567 15.207572
          [4,] 10.722967  7.793154  7.128528  3.511368  1.982854 -1.011026
          [5,]  3.298002  8.039675 14.572189 12.970073  9.836153  4.629030
          [6,] 10.882316  6.583171 12.250103 10.470163  6.353657  8.142503
              
    iterations      [,7]      [,8]
          [1,]  6.849246  6.953424
          [2,] 16.024938 17.861339
          [3,] 15.083187 15.102879
          [4,]  4.710894 10.397039
          [5,] 12.937909  5.541763
          [6,] 15.591744  6.124221

  

#### as.matrix(), as.data.frame(), as.array()

The `as.matrix`, `as.data.frame`, and `as.array` functions can also be
used to retrieve the posterior draws from a stanfit object:

``` r

matrix_of_draws <- as.matrix(fit)
print(colnames(matrix_of_draws))
```

     [1] "mu"       "tau"      "eta[1]"   "eta[2]"   "eta[3]"   "eta[4]"  
     [7] "eta[5]"   "eta[6]"   "eta[7]"   "eta[8]"   "theta[1]" "theta[2]"
    [13] "theta[3]" "theta[4]" "theta[5]" "theta[6]" "theta[7]" "theta[8]"
    [19] "lp__"    

``` r

df_of_draws <- as.data.frame(fit)
print(colnames(df_of_draws))
```

     [1] "mu"       "tau"      "eta[1]"   "eta[2]"   "eta[3]"   "eta[4]"  
     [7] "eta[5]"   "eta[6]"   "eta[7]"   "eta[8]"   "theta[1]" "theta[2]"
    [13] "theta[3]" "theta[4]" "theta[5]" "theta[6]" "theta[7]" "theta[8]"
    [19] "lp__"    

``` r

array_of_draws <- as.array(fit)
print(dimnames(array_of_draws))
```

    $iterations
    NULL

    $chains
    [1] "chain:1" "chain:2" "chain:3" "chain:4"

    $parameters
     [1] "mu"       "tau"      "eta[1]"   "eta[2]"   "eta[3]"   "eta[4]"  
     [7] "eta[5]"   "eta[6]"   "eta[7]"   "eta[8]"   "theta[1]" "theta[2]"
    [13] "theta[3]" "theta[4]" "theta[5]" "theta[6]" "theta[7]" "theta[8]"
    [19] "lp__"    

The `as.matrix` and `as.data.frame` methods essentially return the same
thing except in matrix and data frame form, respectively. The `as.array`
method returns the draws from each chain separately and so has an
additional dimension:

``` r

print(dim(matrix_of_draws))
print(dim(df_of_draws))
print(dim(array_of_draws))
```

    [1] 4000   19
    [1] 4000   19
    [1] 1000    4   19

By default all of the functions for retrieving the posterior draws
return the draws for *all* parameters (and generated quantities). The
optional argument `pars` (a character vector) can be used if only a
subset of the parameters is desired, for example:

``` r

mu_and_theta1 <- as.matrix(fit, pars = c("mu", "theta[1]"))
head(mu_and_theta1)
```

              parameters
    iterations        mu  theta[1]
          [1,]  8.131062 12.390647
          [2,]  6.737244  6.665255
          [3,]  3.571911  5.889686
          [4,] 12.730690 14.184127
          [5,]  5.641158  9.128763
          [6,]  9.728532  8.057433

  

## Posterior summary statistics and convergence diagnostics

Summary statistics are obtained using the `summary` function. The object
returned is a list with two components:

``` r

fit_summary <- summary(fit)
print(names(fit_summary))
```

    [1] "summary"   "c_summary"

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

                      mean    se_mean        sd        2.5%         25%
    mu         8.065056952 0.11178768 5.0959280  -1.7659388   4.7136998
    tau        6.661004311 0.17565043 5.6786017   0.2508909   2.5358807
    eta[1]     0.398848839 0.01568783 0.9299691  -1.4683919  -0.2000774
    eta[2]    -0.009003272 0.01518145 0.8689278  -1.7311259  -0.5764054
    eta[3]    -0.204610425 0.01377734 0.9303807  -1.9952128  -0.8313132
    eta[4]    -0.030810503 0.01351348 0.8980300  -1.7898506  -0.6482907
    eta[5]    -0.372222600 0.01624670 0.8735441  -2.0663086  -0.9292707
    eta[6]    -0.242308519 0.01401813 0.8903689  -1.9415121  -0.8324735
    eta[7]     0.338707139 0.01483374 0.8900252  -1.4786836  -0.2394864
    eta[8]     0.046451882 0.01432622 0.9647513  -1.8591361  -0.5882879
    theta[1]  11.618253016 0.18128230 8.5610044  -2.2303289   6.1177681
    theta[2]   7.848480954 0.09770967 6.3250393  -4.8414097   3.9999836
    theta[3]   6.097747710 0.13467292 7.8609953 -11.7070513   1.9275837
    theta[4]   7.519868683 0.17504005 6.8598769  -5.9574799   3.3818802
    theta[5]   5.045244407 0.13733539 6.3362539  -8.5661291   1.2262338
    theta[6]   6.010320596 0.11813458 6.8401198  -9.4506485   2.1056825
    theta[7]  10.768271055 0.11334758 6.7777191  -0.8663478   6.1667709
    theta[8]   8.390164917 0.13708187 8.0275767  -7.8679573   3.6793160
    lp__     -39.580160912 0.07599636 2.6621323 -45.6441363 -41.1754482
                      50%         75%      97.5%    n_eff      Rhat
    mu         7.85789418  11.3323203  19.068724 2078.063 1.0011083
    tau        5.33314162   9.2189657  21.062401 1045.164 1.0049045
    eta[1]     0.41473553   1.0147012   2.200372 3514.075 0.9997701
    eta[2]    -0.01751924   0.5601695   1.747928 3275.975 1.0001652
    eta[3]    -0.21159681   0.4242515   1.597250 4560.270 0.9999808
    eta[4]    -0.03594472   0.5828865   1.744649 4416.184 1.0003755
    eta[5]    -0.39806140   0.1805950   1.369440 2890.943 1.0014219
    eta[6]    -0.27255740   0.3337835   1.512549 4034.221 0.9999748
    eta[7]     0.35192654   0.9237267   2.089386 3600.005 0.9994109
    eta[8]     0.05999758   0.6941355   1.924542 4534.900 0.9994169
    theta[1]  10.41177233  15.8675550  32.546850 2230.173 1.0022281
    theta[2]   7.83252313  11.7355360  20.268807 4190.360 0.9998955
    theta[3]   6.63534430  10.9284052  20.295178 3407.176 1.0022850
    theta[4]   7.67270397  11.7352723  20.483821 1535.882 1.0026176
    theta[5]   5.47982886   9.2374763  16.379214 2128.630 1.0009484
    theta[6]   6.44389725  10.6123435  18.271347 3352.535 1.0026729
    theta[7]  10.30617582  14.6007699  26.608410 3575.548 1.0000497
    theta[8]   8.17248840  13.0347634  25.349553 3429.327 1.0002239
    lp__     -39.29970561 -37.7013595 -35.188197 1227.082 1.0018561

If, for example, we wanted the only quantiles included to be 10% and
90%, and for only the parameters included to be `mu` and `tau`, we would
specify that like this:

``` r

mu_tau_summary <- summary(fit, pars = c("mu", "tau"), probs = c(0.1, 0.9))$summary
print(mu_tau_summary)
```

            mean   se_mean       sd       10%      90%    n_eff     Rhat
    mu  8.065057 0.1117877 5.095928 1.9157182 14.44746 2078.063 1.001108
    tau 6.661004 0.1756504 5.678602 0.9572312 13.75330 1045.164 1.004905

Since `mu_tau_summary` is a matrix we can pull out columns using their
names:

``` r

mu_tau_80pct <- mu_tau_summary[, c("10%", "90%")]
print(mu_tau_80pct)
```

              10%      90%
    mu  1.9157182 14.44746
    tau 0.9572312 13.75330

  

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
sampler_params_chain1 <- sampler_params[[1]]
colnames(sampler_params_chain1)
```

    [1] "accept_stat__" "stepsize__"    "treedepth__"   "n_leapfrog__" 
    [5] "divergent__"   "energy__"     

To do things like calculate the average value of `accept_stat__` for
each chain (or the maximum value of `treedepth__` for each chain if
using the NUTS algorithm, etc.) the `sapply` function is useful as it
will apply the same function to each component of `sampler_params`:

``` r

mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
print(mean_accept_stat_by_chain)
```

    [1] 0.8420072 0.9405478 0.7536999 0.9326877

``` r

max_treedepth_by_chain <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
print(max_treedepth_by_chain)
```

    [1] 4 5 4 4

  

## Model code

The Stan program itself is also stored in the stanfit object and can be
accessed using `get_stancode`:

``` r

code <- get_stancode(fit)
```

The object `code` is a single string and is not very intelligible when
printed:

``` r

print(code)
```

    [1] "data {\n  int<lower=0> J;          // number of schools\n  array[J] real y;               // estimated treatment effects\n  array[J] real<lower=0> sigma;  // s.e. of effect estimates\n}\nparameters {\n  real mu;\n  real<lower=0> tau;\n  vector[J] eta;\n}\ntransformed parameters {\n  vector[J] theta;\n  theta = mu + tau * eta;\n}\nmodel {\n  target += normal_lpdf(eta | 0, 1);\n  target += normal_lpdf(y | theta, sigma);\n}"
    attr(,"model_name2")
    [1] "schools"

A readable version can be printed using `cat`:

``` r

cat(code)
```

    data {
      int<lower=0> J;          // number of schools
      array[J] real y;               // estimated treatment effects
      array[J] real<lower=0> sigma;  // s.e. of effect estimates
    }
    parameters {
      real mu;
      real<lower=0> tau;
      vector[J] eta;
    }
    transformed parameters {
      vector[J] theta;
      theta = mu + tau * eta;
    }
    model {
      target += normal_lpdf(eta | 0, 1);
      target += normal_lpdf(y | theta, sigma);
    }

  

## Initial values

The `get_inits` function returns initial values as a list with one
component per chain. Each component is itself a (named) list containing
the initial values for each parameter for the corresponding chain:

``` r

inits <- get_inits(fit)
inits_chain1 <- inits[[1]]
print(inits_chain1)
```

    $mu
    [1] -1.805984

    $tau
    [1] 3.709496

    $eta
    [1] -0.8315136  0.7793653  0.3112285 -1.2241347  1.9311699 -1.7715798 -1.0577400
    [8]  1.3319387

    $theta
    [1] -4.8904808  1.0850681 -0.6514836 -6.3469073  5.3576824 -8.3776523 -5.7296666
    [8]  3.1348370

  

## (P)RNG seed

The `get_seed` function returns the (P)RNG seed as an integer:

``` r

print(get_seed(fit))
```

    [1] 346871191

  

## Warmup and sampling times

The `get_elapsed_time` function returns a matrix with the warmup and
sampling times for each chain:

``` r

print(get_elapsed_time(fit))
```

            warmup sample
    chain:1  0.013  0.010
    chain:2  0.013  0.016
    chain:3  0.013  0.008
    chain:4  0.014  0.016
