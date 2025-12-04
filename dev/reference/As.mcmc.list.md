# Create an mcmc.list from a stanfit object

Create an [`mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html)
(coda) from a `stanfit` object.

## Usage

``` r
As.mcmc.list(object, pars, include = TRUE, ...)
```

## Arguments

- object:

  object of class `"stanfit"`

- pars:

  optional character vector of parameters to include

- include:

  logical scalar indicating whether to include (the default) or exclude
  the parameters named in `pars`

- ...:

  unused

## Value

An object of class
[`mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html).
