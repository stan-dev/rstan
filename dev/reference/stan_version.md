# Obtain the version of Stan

The stan version is in form of `major.minor.patch`; the first version is
1.0.0, indicating major version 1, minor version 0 and patch level 0.
Functionality only changes with minor versions and backward
compatibility will only be affected by major versions.

## Usage

``` r
stan_version()
```

## Value

A character string giving the version of Stan used in this version of
RStan.

## References

The Stan Development Team *Stan Modeling Language User's Guide and
Reference Manual*. <https://mc-stan.org/>.

## See also

[`stan`](https://mc-stan.org/rstan/dev/reference/stan.md) and
[`stan_model`](https://mc-stan.org/rstan/dev/reference/stan_model.md)

## Examples

``` r
  stan_version() 
#> [1] "2.37.0"
```
