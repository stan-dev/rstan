# Class representing model compiled from C++

A `stanmodel` object represents the model compiled from C++ code. The
`sampling` method defined in this class may be used to draw samples from
the model and `optimizing` method is for obtaining a point estimate by
maximizing the log-posterior.

## Objects from the Class

Instances of `stanmodel` are usually created by calling function
`stan_model` or function `stan`.

## Slots

- `model_name`::

  The model name, an object of type `character`.

- `model_code`::

  The Stan model specification, an object of type `character`.

- `model_cpp`::

  Object of type `list` that includes the C++ code for the model.

- `mk_cppmodule`::

  A function to return a RCpp module. This function will be called in
  function `sampling` and `optimzing` with one argument (the instance of
  `stanmodel` itself).

- `dso`::

  Object of S4 class `cxxdso`. The container for the dynamic shared
  objects compiled from the C++ code of the model, returned from
  function `cxxfunction` in package inline.

## Methods

- `show`:

  `signature(object = "stanmodel")`: print the Stan model specification.

- `vb`:

  `signature(object = "stanmodel")`: use the variational Bayes
  algorithms.

- `sampling`:

  `signature(object = "stanmodel")`: draw samples for the model (see
  [`sampling`](https://mc-stan.org/rstan/dev/reference/stanmodel-method-sampling.md)).

- `optimizing`:

  `signature(object = "stanmodel")`: obtain a point estimate by
  maximizing the posterior (see
  [`optimizing`](https://mc-stan.org/rstan/dev/reference/stanmodel-method-optimizing.md)).

- `get_cppcode`:

  `signature(object = "stanmodel")`: returns the C++ code for the model
  as a character string. This is part of the C++ code that is compiled
  to the dynamic shared object for the model.

- `get_stancode`:

  `signature(object = "stanmodel")`: returns the Stan code for the model
  as a character string

- `get_cxxflags`:

  `signature(object = "stanmodel")`: return the `CXXFLAGS` used for
  compiling the model. The returned string is like `CXXFLAGS = -O3`.

## Note

Objects of class `stanmodel` can be saved for use across R sessions only
if `save_dso = TRUE` is set during calling functions that create
`stanmodel` objects (e.g., `stan` and `stan_model`).

Even if `save_dso = TRUE`, the model cannot be loaded on a platform
(operating system, 32 bits or 64 bits, etc.) that differs from the one
on which it was compiled.

## See also

[`stan_model`](https://mc-stan.org/rstan/dev/reference/stan_model.md),
[`stanc`](https://mc-stan.org/rstan/dev/reference/stanc.md)
[`sampling`](https://mc-stan.org/rstan/dev/reference/stanmodel-method-sampling.md),
[`optimizing`](https://mc-stan.org/rstan/dev/reference/stanmodel-method-optimizing.md),
[`vb`](https://mc-stan.org/rstan/dev/reference/stanmodel-method-vb.md)

## Examples

``` r
showClass("stanmodel")
#> Class "stanmodel" [package "rstan"]
#> 
#> Slots:
#>                                                                        
#> Name:    model_name   model_code    model_cpp mk_cppmodule          dso
#> Class:    character    character         list     function       cxxdso
```
