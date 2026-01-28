# Extract samples from a fitted Stan model

Extract samples from a fitted model represented by an instance of class
[`stanfit`](https://mc-stan.org/rstan/dev/reference/stanfit-class.md).

## Usage

``` r
# S4 method for class 'stanfit'
extract(object, pars, permuted = TRUE, inc_warmup = FALSE, 
  include = TRUE)
```

## Methods

- extract:

  `signature(object = "stanfit")` Extract samples from a fitted model
  represented by an instance of class
  [`stanfit`](https://mc-stan.org/rstan/dev/reference/stanfit-class.md).

## Arguments

- object:

  An object of class
  [`stanfit`](https://mc-stan.org/rstan/dev/reference/stanfit-class.md).

- pars:

  An optional character vector providing the parameter names (or other
  quantity names) of interest. If not specified, all parameters and
  other quantities are used. The log-posterior with name `lp__` is also
  included by default.

- permuted:

  A logical scalar indicating whether the draws after the *warmup*
  period in each chain should be *permuted* and *merged*. If `FALSE`,
  the original order is kept. For each `stanfit` object, the permutation
  is fixed (i.e., extracting samples a second time will give the same
  sequence of iterations).

- inc_warmup:

  A logical scalar indicating whether to include the warmup draws. This
  argument is only relevant if `permuted` is `FALSE`.

- include:

  A logical scalar indicating whether the parameters named in `pars`
  should be included (`TRUE`) or excluded (`FALSE`).

## Value

When `permuted = TRUE`, this function returns a named list, every
element of which is an array representing samples for a parameter with
all chains merged together.

When `permuted = FALSE`, an array is returned; the first dimension is
for the iterations, the second for the number of chains, the third for
the parameters. Vectors and arrays are expanded to one parameter (a
scalar) per cell, with names indicating the third dimension. See the
examples (with comments) below. The
[`monitor`](https://mc-stan.org/rstan/dev/reference/monitor.md) function
can be applied to the returned array to obtain a summary (similar to the
`print` method for
[`stanfit`](https://mc-stan.org/rstan/dev/reference/stanfit-class.md)
objects).

## See also

S4 class
[`stanfit`](https://mc-stan.org/rstan/dev/reference/stanfit-class.md),
[`as.array.stanfit`](https://mc-stan.org/rstan/dev/reference/stanfit2array-method.md),
and [`monitor`](https://mc-stan.org/rstan/dev/reference/monitor.md)

## Examples

``` r
# \dontrun{
ex_model_code <- '
  parameters {
    array[2, 3] real alpha;
    array[2] real beta; 
  } 
  model {
    for (i in 1:2) for (j in 1:3) 
      alpha[i, j] ~ normal(0, 1); 
    for (i in 1:2) 
      beta ~ normal(0, 2); 
  } 
'

## fit the model 
fit <- stan(model_code = ex_model_code, chains = 4) 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 4e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.04 seconds.
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
#> Chain 1:  Elapsed Time: 0.013 seconds (Warm-up)
#> Chain 1:                0.013 seconds (Sampling)
#> Chain 1:                0.026 seconds (Total)
#> Chain 1: 
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
#> Chain 2:  Elapsed Time: 0.013 seconds (Warm-up)
#> Chain 2:                0.013 seconds (Sampling)
#> Chain 2:                0.026 seconds (Total)
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
#> Chain 3:  Elapsed Time: 0.013 seconds (Warm-up)
#> Chain 3:                0.013 seconds (Sampling)
#> Chain 3:                0.026 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 2e-06 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.02 seconds.
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
#> Chain 4:  Elapsed Time: 0.013 seconds (Warm-up)
#> Chain 4:                0.014 seconds (Sampling)
#> Chain 4:                0.027 seconds (Total)
#> Chain 4: 

## extract alpha and beta with 'permuted = TRUE' 
fit_ss <- extract(fit, permuted = TRUE) # fit_ss is a list 
## list fit_ss should have elements with name 'alpha', 'beta', 'lp__'
alpha <- fit_ss$alpha  
beta <- fit_ss$beta 
## or extract alpha by just specifying pars = 'alpha' 
alpha2 <- extract(fit, pars = 'alpha', permuted = TRUE)$alpha 
print(identical(alpha, alpha2)) 
#> [1] TRUE

## or extract alpha by excluding beta and lp__
alpha3 <- extract(fit, pars = c('beta', 'lp__'), 
                  permuted = TRUE, include = FALSE)$alpha
print(identical(alpha, alpha3))
#> [1] TRUE

## get the samples for alpha[1,1] and beta[2] 
alpha_11 <- alpha[, 1, 1] 
beta_2 <- beta[, 2] 

## extract samples with 'permuted = FALSE' 
fit_ss2 <- extract(fit, permuted = FALSE) # fit_ss2 is an array  

## the dimensions of fit_ss2 should be  
## "# of iterations * # of chains * # of parameters"
dim(fit_ss2) 
#> [1] 1000    4    9

## since the third dimension of `fit_ss2` indicates 
## parameters, the names should be 
##  alpha[1,1], alpha[2,1], alpha[1,2], alpha[2,2], 
##  alpha[1,3], alpha[2,3], beta[1], beta[2], and lp__ 
## `lp__` (the log-posterior) is always included 
## in the samples.  
dimnames(fit_ss2) 
#> $iterations
#> NULL
#> 
#> $chains
#> [1] "chain:1" "chain:2" "chain:3" "chain:4"
#> 
#> $parameters
#> [1] "alpha[1,1]" "alpha[2,1]" "alpha[1,2]" "alpha[2,2]" "alpha[1,3]"
#> [6] "alpha[2,3]" "beta[1]"    "beta[2]"    "lp__"      
#> 
# }

# Create a stanfit object from reading CSV files of samples (saved in rstan
# package) generated by funtion stan for demonstration purpose from model as follows. 
# 
excode <- '
  transformed data {
    array[20] real y;
    y[1] <- 0.5796;  y[2]  <- 0.2276;   y[3] <- -0.2959; 
    y[4] <- -0.3742; y[5]  <- 0.3885;   y[6] <- -2.1585;
    y[7] <- 0.7111;  y[8]  <- 1.4424;   y[9] <- 2.5430; 
    y[10] <- 0.3746; y[11] <- 0.4773;   y[12] <- 0.1803; 
    y[13] <- 0.5215; y[14] <- -1.6044;  y[15] <- -0.6703; 
    y[16] <- 0.9459; y[17] <- -0.382;   y[18] <- 0.7619;
    y[19] <- 0.1006; y[20] <- -1.7461;
  }
  parameters {
    real mu;
    real<lower=0, upper=10> sigma;
    vector[2] z[3];
    real<lower=0> alpha;
  } 
  model {
    y ~ normal(mu, sigma);
    for (i in 1:3) 
      z[i] ~ normal(0, 1);
    alpha ~ exponential(2);
  } 
'
# exfit <- stan(model_code = excode, save_dso = FALSE, iter = 200, 
#               sample_file = "rstan_doc_ex.csv")
# 
exfit <- read_stan_csv(dir(system.file('misc', package = 'rstan'),
                       pattern='rstan_doc_ex_[[:digit:]].csv',
                       full.names = TRUE))

ee1 <- extract(exfit, permuted = TRUE)
print(names(ee1))
#> [1] "mu"    "sigma" "z"     "alpha" "lp__" 

for (name in names(ee1)) {
  cat(name, "\n")
  print(dim(ee1[[name]]))
}
#> mu 
#> [1] 400
#> sigma 
#> [1] 400
#> z 
#> [1] 400   3   2
#> alpha 
#> [1] 400
#> lp__ 
#> [1] 400

ee2 <- extract(exfit, permuted = FALSE)
print(dim(ee2))
#> [1] 100   4  10
print(dimnames(ee2))
#> $iterations
#> NULL
#> 
#> $chains
#> [1] "chain:1" "chain:2" "chain:3" "chain:4"
#> 
#> $parameters
#>  [1] "mu"     "sigma"  "z[1,1]" "z[2,1]" "z[3,1]" "z[1,2]" "z[2,2]" "z[3,2]"
#>  [9] "alpha"  "lp__"  
#> 
```
