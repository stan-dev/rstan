# Demonstrate examples included in Stan

Stan includes a variety of examples and most of the BUGS example models
that are translated into Stan modeling language. One example is chosen
from a list created from matching user input and gets fitted in the
demonstration.

## Usage

``` r
stan_demo(model = character(0), 
            method = c("sampling", "optimizing", "meanfield", "fullrank"), ...)
```

## Arguments

- model:

  A character string for model name to specify which model will be used
  for demonstration. The default is an empty string, which prompts the
  user to select one the available models. If `model = 0`, a character
  vector of all models is returned without any user intervention. If
  `model = i` where `i > 0`, then the ith available model is chosen
  without user intervention, which is useful for testing.

- method:

  Whether to call
  [`sampling`](https://mc-stan.org/rstan/dev/reference/stanmodel-method-sampling.md)
  (the default),
  [`optimizing`](https://mc-stan.org/rstan/dev/reference/stanmodel-method-optimizing.md),
  or one of the variants of
  [`vb`](https://mc-stan.org/rstan/dev/reference/stanmodel-method-vb.md)
  for the demonstration

- ...:

  Further arguments passed to `method`.

## Value

An S4 object of `stanfit`, unless `model = 0`, in which case a character
vector of paths to available models is returned.

## References

The Stan Development Team *Stan Modeling Language User's Guide and
Reference Manual*. <https://mc-stan.org/>.

## See also

[`sampling`](https://mc-stan.org/rstan/dev/reference/stanmodel-method-sampling.md),
[`optimizing`](https://mc-stan.org/rstan/dev/reference/stanmodel-method-optimizing.md)

## Examples

``` r
  # \dontrun{
     dogsfit <- stan_demo("dogs") # run the dogs model
#> 
#> > n_dogs <- 31
#> 
#> > n_trials <- 25
#> 
#> > y <- structure(c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
#> +     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
#> +     1, 0, 1, 1, 1, 0, 1, 1, 1,  .... [TRUNCATED] 
#> 
#> SAMPLING FOR MODEL 'dogs' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000203 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.03 seconds.
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
#> Chain 1:  Elapsed Time: 2.025 seconds (Warm-up)
#> Chain 1:                2.032 seconds (Sampling)
#> Chain 1:                4.057 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'dogs' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.00019 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 1.9 seconds.
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
#> Chain 2:  Elapsed Time: 2.201 seconds (Warm-up)
#> Chain 2:                2.045 seconds (Sampling)
#> Chain 2:                4.246 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'dogs' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 0.000192 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 1.92 seconds.
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
#> Chain 3:  Elapsed Time: 2.156 seconds (Warm-up)
#> Chain 3:                1.951 seconds (Sampling)
#> Chain 3:                4.107 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'dogs' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 0.000203 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 2.03 seconds.
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
#> Chain 4:  Elapsed Time: 2.102 seconds (Warm-up)
#> Chain 4:                2.002 seconds (Sampling)
#> Chain 4:                4.104 seconds (Total)
#> Chain 4: 
     fit1 <- stan_demo(1) # run model_names[1]
#> 
#> > "score1" <- c(0.446, 0.334, -0.233, -0.304, -0.256, 
#> +     0.515, -0.652, 0.271, -0.556, 0.955, 0.307, -0.434, -0.272, 
#> +     -0.404, 0.415, 0.39, - .... [TRUNCATED] 
#> 
#> > "party" <- c(1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 
#> +     0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 
#> +     0, 0, 0, 1, 1, 0, 1, 0, 1 .... [TRUNCATED] 
#> 
#> > "x" <- c(0.62073, 0.50803, 0.38243, 0.29769, 0.32751, 
#> +     0.53766, 0.20003, 0.52211, 0.31267, 0.64143, 0.69113, 0.4386, 
#> +     0.30166, 0.25836,  .... [TRUNCATED] 
#> 
#> > "N" <- 357
#> 
#> SAMPLING FOR MODEL 'ideo_interactions' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 6e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.06 seconds.
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
#> Chain 1:  Elapsed Time: 0.11 seconds (Warm-up)
#> Chain 1:                0.126 seconds (Sampling)
#> Chain 1:                0.236 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'ideo_interactions' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 4e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.04 seconds.
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
#> Chain 2:  Elapsed Time: 0.109 seconds (Warm-up)
#> Chain 2:                0.126 seconds (Sampling)
#> Chain 2:                0.235 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'ideo_interactions' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 4e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.04 seconds.
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
#> Chain 3:  Elapsed Time: 0.11 seconds (Warm-up)
#> Chain 3:                0.11 seconds (Sampling)
#> Chain 3:                0.22 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'ideo_interactions' NOW (CHAIN 4).
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
#> Chain 4:  Elapsed Time: 0.11 seconds (Warm-up)
#> Chain 4:                0.121 seconds (Sampling)
#> Chain 4:                0.231 seconds (Total)
#> Chain 4: 
  # }
```
