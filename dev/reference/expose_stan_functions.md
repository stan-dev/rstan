# Expose user-defined Stan functions to R for testing and simulation

The Stan modeling language allows users to define their own functions in
a `functions` block at the top of a Stan program. The
`expose_stan_functions` utility function uses
[`sourceCpp`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html) to export
those user-defined functions to the specified environment for testing
inside R or for doing posterior predictive simulations in R rather than
in the `generated quantities` block of a Stan program.

## Usage

``` r
expose_stan_functions(stanmodel, includes = NULL, 
                        show_compiler_warnings = FALSE, ...)
  get_rng(seed = 0L)
  get_stream()
```

## Arguments

- stanmodel:

  A
  [`stanmodel`](https://mc-stan.org/rstan/dev/reference/stanmodel-class.md)
  object, a
  [`stanfit`](https://mc-stan.org/rstan/dev/reference/stanfit-class.md)
  object, a list produced by
  [`stanc`](https://mc-stan.org/rstan/dev/reference/stanc.md) or the
  path to a Stan program (`.stan` file). In any of these cases, the
  underlying Stan program should contain a non-empty `functions` block.

- includes:

  If not `NULL` (the default), then a character vector of length one
  (possibly containing one or more `"\n"`) of the form
  `'#include "/full/path/to/my_header.hpp"'`, which will be inserted
  into the C++ code in the model's namespace and can be used to provide
  definitions of functions that are declared but not defined in
  `stanmodel`

- show_compiler_warnings:

  Logical scalar defaulting to `FALSE` that controls whether compiler
  warnings, which can be numerous and have never been relevant, are
  shown

- seed:

  An integer vector of length one indicating the state of Stan's
  pseudo-random number generator

- ...:

  Further arguments passed to
  [`sourceCpp`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html).

## Details

The `expose_stan_functions` function requires as much compliance with
the C++14 standard as is implemented in the RTools toolchain for
Windows. On Windows, you will likely need to specify
`CXX14 = g++ -std=c++1y` in the file whose path is
[`normalizePath`](https://rdrr.io/r/base/normalizePath.html)`("~/.R/Makevars")`
in order for `expose_stan_functions` to work. Outside of Windows, the
necessary compiler flags are set programatically, which is likely to
suffice.

There are a few special types of user-defined Stan functions for which
some additional details are relevant:

### (P)RNG functions

If a user-defined Stan function ends in `_rng`, then it can use the
Boost pseudo-random number generator used by Stan. When exposing such
functions to R, `base_rng__` and `pstream__` arguments will be added to
the [`formals`](https://rdrr.io/r/base/formals.html). The `base_rng__`
argument should be passed the result of a call to `get_rng` (perhaps
specifying its `seed` argument for reproducibility) and the `pstream__`
should be passed the result of a call to `get_stream`, which can be used
to see the result of `print` and `reject` calls in the user-defined Stan
functions. These arguments default to `get_stream()` and `get_rng()`
respectively.

### LP functions

If a user-defined Stan function ends in `_lp`, then it can modify the
log-probability used by Stan to evaluate Metropolis proposals or as an
objective function for optimization. When exposing such functions to R,
a `lp__` argument will be added to the
[`formals`](https://rdrr.io/r/base/formals.html). This `lp__` argument
defaults to zero, but a [`double`](https://rdrr.io/r/base/double.html)
precision scalar may be passed to this argument when the function is
called from R. Such a user-defined Stan function can terminate with
`return target();` or can execute `print(target());` to verify that the
calculation is correct.

## Value

The names of the new functions in `env` are returned invisibly.

## See also

[`sourceCpp`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html) and the
section in the Stan User Manual on user-defined functions

## Examples

``` r
# \dontrun{
model_code <-
  '
  functions {
    real standard_normal_rng() {
      return normal_rng(0,1);
   }
  }
'
expose_stan_functions(stanc(model_code = model_code))
standard_normal_rng()
#> [1] 0
PRNG <- get_rng(seed = 3)
o <- get_stream()
standard_normal_rng(PRNG, o)
#> [1] 1.165035
# }
```
