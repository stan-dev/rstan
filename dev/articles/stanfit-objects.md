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

    [1]  5.22198096 10.25431910 11.48222576  0.07918646 -3.30819970 13.58606549

``` r

head(list_of_draws$tau)
```

    [1]  5.996095  5.909190  2.802256 11.933985  5.368431  6.944701

``` r

head(list_of_draws$theta)
```

              
    iterations      [,1]      [,2]      [,3]       [,4]       [,5]      [,6]
          [1,]  2.283421  8.578929  3.168522  6.7439015  4.6031011  6.907852
          [2,] 11.359170  6.478249 12.446666  9.9424217  6.7237383  2.968103
          [3,] 12.249217 10.539249  8.504100 12.0364981  8.5697172  7.527797
          [4,] -8.860330 10.688263  6.029505  3.8369223 -0.8734694 10.821918
          [5,] -4.906261 -3.894452 -4.229853  0.4010875 -4.9593446 -7.039061
          [6,] 23.468212  5.095051  1.616788 13.2977748  5.4846604  9.664369
              
    iterations      [,7]       [,8]
          [1,]  3.594139 -0.4996312
          [2,] 12.851494 18.5077162
          [3,] 12.471018 13.9139741
          [4,]  2.494871  8.9744022
          [5,] -5.959557  1.5743520
          [6,] 15.388620 13.4485342

  

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
          [1,] 16.808565 19.377383
          [2,]  9.085071  5.254745
          [3,]  9.490775 12.740426
          [4,]  7.024909  2.712334
          [5,] 10.325071 -1.006234
          [6,] 12.216922 13.034301

  

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
    mu         8.095883193 0.10639687 5.0640173  -1.4324607   4.9183331
    tau        6.293738721 0.14506881 5.3735697   0.2015424   2.3575233
    eta[1]     0.381015597 0.01459744 0.9507598  -1.5796632  -0.2234058
    eta[2]     0.009007562 0.01417983 0.9034728  -1.7880656  -0.5770494
    eta[3]    -0.187210058 0.01452286 0.9459167  -2.0682485  -0.8217128
    eta[4]    -0.035375600 0.01382681 0.8657225  -1.7561734  -0.5939173
    eta[5]    -0.361676495 0.01413001 0.8717491  -2.0227518  -0.9591457
    eta[6]    -0.214510464 0.01409324 0.8956887  -1.9366989  -0.8332017
    eta[7]     0.350821227 0.01497281 0.9004367  -1.4300924  -0.2491383
    eta[8]     0.052958067 0.01345953 0.9322137  -1.7588616  -0.5763021
    theta[1]  11.322773164 0.15505640 8.2776572  -1.9203819   6.0332449
    theta[2]   8.004101819 0.09630841 6.2422112  -4.2607888   4.0382535
    theta[3]   6.299784074 0.11559691 7.5421778 -10.6624969   2.3119588
    theta[4]   7.767177199 0.08539246 6.3793292  -5.0272230   3.8436346
    theta[5]   5.247303714 0.09140613 6.3314446  -8.9518058   1.5579495
    theta[6]   6.336858566 0.09616920 6.8222195  -8.2020483   2.3757087
    theta[7]  10.608966215 0.11328368 6.8471616  -1.9536568   6.0474640
    theta[8]   8.422207428 0.12279133 7.7450219  -6.6975366   3.8857270
    lp__     -39.674936458 0.07277398 2.6759517 -45.6028614 -41.2705096
                      50%         75%      97.5%    n_eff      Rhat
    mu         7.90163063  11.3031398  18.378286 2265.336 1.0007683
    tau        5.08593552   8.8503482  19.143525 1372.075 1.0005621
    eta[1]     0.39903034   1.0292366   2.190002 4242.172 0.9997237
    eta[2]     0.01004062   0.5826053   1.824950 4059.646 0.9994544
    eta[3]    -0.19934059   0.4211110   1.755174 4242.298 1.0009679
    eta[4]    -0.04264697   0.5060588   1.721352 3920.250 1.0017374
    eta[5]    -0.39280179   0.2076239   1.452278 3806.259 0.9996797
    eta[6]    -0.21183422   0.3730205   1.590079 4039.175 1.0004399
    eta[7]     0.37472434   0.9473720   2.146286 3616.594 0.9999667
    eta[8]     0.04892119   0.6660938   1.915073 4797.014 1.0002494
    theta[1]  10.27992160  15.3176698  30.916918 2849.938 0.9996924
    theta[2]   7.89781756  11.7876773  21.066178 4200.960 0.9998463
    theta[3]   6.76232287  11.1205656  19.795872 4256.975 1.0004553
    theta[4]   7.67299298  11.8187479  20.456410 5580.986 0.9994345
    theta[5]   5.75700790   9.5487355  16.429356 4797.942 0.9996267
    theta[6]   6.80506899  10.7812451  18.395870 5032.449 0.9999338
    theta[7]  10.09710317  14.5408069  26.214924 3653.310 0.9994559
    theta[8]   8.28648790  12.7220356  25.181417 3978.413 0.9997589
    lp__     -39.42810189 -37.7446463 -35.271979 1352.086 1.0019371

If, for example, we wanted the only quantiles included to be 10% and
90%, and for only the parameters included to be `mu` and `tau`, we would
specify that like this:

``` r

mu_tau_summary <- summary(fit, pars = c("mu", "tau"), probs = c(0.1, 0.9))$summary
print(mu_tau_summary)
```

            mean   se_mean       sd       10%      90%    n_eff     Rhat
    mu  8.095883 0.1063969 5.064017 1.8547907 14.36239 2265.336 1.000768
    tau 6.293739 0.1450688 5.373570 0.9052414 13.12377 1372.075 1.000562

Since `mu_tau_summary` is a matrix we can pull out columns using their
names:

``` r

mu_tau_80pct <- mu_tau_summary[, c("10%", "90%")]
print(mu_tau_80pct)
```

              10%      90%
    mu  1.8547907 14.36239
    tau 0.9052414 13.12377

  

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

    [1] 0.8228008 0.9521825 0.9379850 0.8681285

``` r

max_treedepth_by_chain <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
print(max_treedepth_by_chain)
```

    [1] 4 4 4 4

  

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
    chain:1  0.015  0.012
    chain:2  0.015  0.018
    chain:3  0.018  0.018
    chain:4  0.015  0.014
