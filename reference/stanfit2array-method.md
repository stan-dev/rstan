# Create array, matrix, or data.frame objects from samples in a `stanfit` object

The samples (without warmup) included in a
[`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.md) object
can be coerced to an `array`, `matrix`, or `data.frame`. Methods are
also provided for checking and setting names and dimnames.

## Usage

``` r
# S3 method for class 'stanfit'
as.array(x, ...) 
  # S3 method for class 'stanfit'
as.matrix(x, ...)
  # S3 method for class 'stanfit'
as.data.frame(x, ...)
  # S3 method for class 'stanfit'
is.array(x) 
  # S3 method for class 'stanfit'
dim(x)
  # S3 method for class 'stanfit'
dimnames(x)
  # S3 method for class 'stanfit'
names(x)
  # S3 method for class 'stanfit'
names(x) <- value
```

## Arguments

- x:

  An object of S4 class
  [`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.md).

- ...:

  Additional parameters that can be passed to
  [`extract`](https://mc-stan.org/rstan/reference/stanfit-method-extract.md)
  for extracting samples from `x`. For now, `pars` is the only
  additional parameter supported.

- value:

  For the `names` replacement method, a character vector to set/replace
  the parameter names in `x`.

## Details

`as.array` and `as.matrix` can be applied to a `stanfit` object to
coerce the samples without warmup to `array` or `matrix`. The
`as.data.frame` method first calls `as.matrix` and then coerces this
matrix to a `data.frame`.

The array has three named dimensions: iterations, chains, parameters.
For `as.matrix`, all chains are combined, leaving a matrix of iterations
by parameters.

## Value

`as.array`, `as.matrix`, and `as.data.frame` return an array, matrix,
and data.frame, respectively.

`dim` and `dimnames` return the dim and dimnames of the array object
that could be created, while `names` returns the third element of the
`dimnames`, which are the names of the margins of the posterior
distribution. The `names` assignment method allows for assigning more
interpretable names to them.

`is.array` returns `TRUE` for `stanfit` objects that include samples;
otherwise `FALSE`.

When the `stanfit` object does not contain samples, empty objects are
returned from `as.array`, `as.matrix`, `as.data.frame`, `dim`,
`dimnames`, and `names`.

## See also

S4 class
[`stanfit`](https://mc-stan.org/rstan/reference/stanfit-class.md) and
its method
[`extract`](https://mc-stan.org/rstan/reference/stanfit-method-extract.md)

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
      beta[i] ~ normal(0, 2); 
    # beta ~ normal(0, 2) // vectorized version
  } 
'

## fit the model 
fit <- stan(model_code = ex_model_code, chains = 4) 
#> Error in stanc(file = file, model_code = model_code, model_name = model_name,     verbose = verbose, obfuscate_model_name = obfuscate_model_name,     allow_undefined = allow_undefined, allow_optimizations = allow_optimizations,     standalone_functions = standalone_functions, use_opencl = use_opencl,     warn_pedantic = warn_pedantic, warn_uninitialized = warn_uninitialized,     isystem = isystem): 0
#> Syntax error in 'string', line 10, column 4 to column 5, lexing error:
#>    -------------------------------------------------
#>      8:      for (i in 1:2)
#>      9:        beta[i] ~ normal(0, 2);
#>     10:      # beta ~ normal(0, 2) // vectorized version
#>              ^
#>     11:    }
#>    -------------------------------------------------
#> 
#> Invalid character found.

dim(fit)
#> [1] 100   4  10
dimnames(fit)
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
is.array(fit) 
#> [1] TRUE
a <- as.array(fit)
m <- as.matrix(fit)
d <- as.data.frame(fit)
# }
```
