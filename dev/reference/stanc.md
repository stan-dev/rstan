# Translate Stan model specification to C++ code

Translate a model specification in Stan code to C++ code, which can then
be compiled and loaded for sampling.

## Usage

``` r
stanc(
    file, model_code = '', model_name = "anon_model", verbose = FALSE,
    obfuscate_model_name = TRUE,
    allow_undefined = isTRUE(getOption("stanc.allow_undefined", FALSE)),
    allow_optimizations = isTRUE(getOption("stanc.allow_optimizations", FALSE)),
    standalone_functions = isTRUE(getOption("stanc.standalone_functions", FALSE)),
    use_opencl = isTRUE(getOption("stanc.use_opencl", FALSE)),
    warn_pedantic = isTRUE(getOption("stanc.warn_pedantic", FALSE)),
    warn_uninitialized = isTRUE(getOption("stanc.warn_uninitialized", FALSE)),
    isystem = c(if (!missing(file)) dirname(file), getwd()))
  stanc_builder(
    file, isystem = c(dirname(file), getwd()),
    verbose = FALSE, obfuscate_model_name = FALSE,
    allow_undefined = isTRUE(getOption("stanc.allow_undefined", FALSE)),
    allow_optimizations = isTRUE(getOption("stanc.allow_optimizations", FALSE)),
    standalone_functions = isTRUE(getOption("stanc.standalone_functions", FALSE)),
    use_opencl = isTRUE(getOption("stanc.use_opencl", FALSE)),
    warn_pedantic = isTRUE(getOption("stanc.warn_pedantic", FALSE)),
    warn_uninitialized = isTRUE(getOption("stanc.warn_uninitialized", FALSE)))
```

## Arguments

- file:

  A character string or a connection that R supports specifying the Stan
  model specification in Stan's modeling language.

- model_code:

  Either a character string containing a Stan model specification or the
  name of a character string object in the workspace. This parameter is
  used only if parameter `file` is not specified, so it defaults to the
  empty string.

- model_name:

  A character string naming the model. The default is `"anon_model"`.
  However, the model name will be derived from `file` or `model_code`
  (if `model_code` is the name of a character string object) if
  `model_name` is not specified.

- verbose:

  Logical, defaulting to `FALSE`. If `TRUE` more intermediate
  information is printed during the translation procedure.

- obfuscate_model_name:

  Logical, defaulting to `TRUE`, indicating whether to use a
  randomly-generated character string for the name of the C++ class.
  This prevents name clashes when compiling multiple models in the same
  R session.

- isystem:

  A character vector naming a path to look for file paths in `file` that
  are to be included within the Stan program named by `file`. See the
  Details section below.

- allow_undefined:

  A logical scalar defaulting to `FALSE` indicating whether to allow
  Stan functions to be declared but not defined in `file` or
  `model_code`. If `TRUE`, then it is the caller's responsibility to
  provide a function definition in another header file or linked shared
  object.

- allow_optimizations:

  A logical scalar defaulting to `FALSE` indicating whether to allow
  level-1 compiler optimization to the Stan code.

- standalone_functions:

  A logical scalar defaulting to `FALSE` indicating whether to generate
  the standalone functions C++ code.

- use_opencl:

  A logical scalar defaulting to `FALSE` indicating whether to try to
  use `matrix_cl` signatures.

- warn_pedantic:

  A logical scalar defaulting to `FALSE` indicating whether to emit
  warnings about common mistakes in Stan programs.

- warn_uninitialized:

  A logical scalar defaulting to `FALSE` indicating whether emit
  warnings about uninitialized variables.

## Details

The `stanc_builder` function supports the standard C++ convention of
specifying something like `#include "my_includes.txt"` on an entire line
within the file named by the `file` argument. In other words,
`stanc_builder` would look for `"my_includes.txt"` in (or under) the
directories named by the `isystem` argument and — if found — insert its
contents verbatim at that position before calling `stanc` on the
resulting `model_code`. This mechanism reduces the need to copy common
chunks of code across Stan programs. It is possible to include such
files recursively.

Note that line numbers referred to in parser warnings or errors refer to
the postprocessed Stan program rather than `file`. In the case of a
parser error, the postprocessed Stan program will be printed after the
error message. Line numbers referred to in messages while Stan is
executing also refer to the postprocessed Stan program which can be
obtained by calling
[`get_stancode`](https://mc-stan.org/rstan/dev/reference/stanfit-class.md).

## Value

A list with named entries:

1.  `model_name` Character string for the model name.

2.  `model_code` Character string for the model's Stan specification.

3.  `cppcode` Character string for the model's C++ code.

4.  `status` Logical indicating success/failure (always `TRUE`) of
    translating the Stan code.

## Note

Unlike R, in which variable identifiers may contain dots (e.g. `a.1`),
Stan prohibits dots from occurring in variable identifiers. Furthermore,
C++ reserved words and Stan reserved words may not be used for variable
names; see the Stan User's Guide for a complete list.

## References

The Stan Development Team *Stan Modeling Language User's Guide and
Reference Manual*. <https://mc-stan.org/>.

The Stan Development Team *CmdStan Interface User's Guide*.
<https://mc-stan.org>.

## See also

[`stan_model`](https://mc-stan.org/rstan/dev/reference/stan_model.md)
and [`stan`](https://mc-stan.org/rstan/dev/reference/stan.md)

## Examples

``` r
stanmodelcode <- "
data {
  int<lower=0> N;
  array[N] real y;
}

parameters {
  real mu;
}

model {
  mu ~ normal(0, 10);
  y ~ normal(mu, 1);
}
"

r <- stanc(model_code = stanmodelcode, model_name = "normal1")
str(r)
#> List of 5
#>  $ status       : logi TRUE
#>  $ model_cppname: chr "model1cc36f0f1c33_normal1"
#>  $ cppcode      : chr "#ifndef USE_STANC3\n#define USE_STANC3\n#endif\n// Code generated by stanc v2.38.0-103-g543515f\n#include <stan"| __truncated__
#>  $ model_name   : chr "normal1"
#>  $ model_code   : chr "data {\n  int<lower=0> N;\n  array[N] real y;\n}\nparameters {\n  real mu;\n}\nmodel {\n  mu ~ normal(0, 10);\n"| __truncated__
#>   ..- attr(*, "model_name2")= chr "file1cc31b2535b6"
```
