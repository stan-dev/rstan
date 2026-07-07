# Obtain the full path of file `Makeconf`

Obtain the full path of file `Makeconf`, in which, for example the flags
for compiling C/C++ code are configured.

## Usage

``` r
makeconf_path()
```

## Value

An character string for the full path of file `Makeconf`.

## Details

The configuration for compiling shared objects using `R CMD SHLIB` are
set in file `Makeconf`. To change how the C++ code is compiled, modify
this file. For RStan, package inline compiles the C++ code using
`R CMD SHLIB`. To speed up compiled Stan models, increase the
optimization level to `-O3` defined in property `CXXFLAGS` in the file
`Makeconf`. This file may also be modified to specify alternative C++
compilers, such as clang++ or later versions of g++.

## See also

[`stan`](https://mc-stan.org/rstan/dev/reference/stan.md)

## Examples

``` r
makeconf_path() 
#> [1] "/opt/R/4.6.1/lib/R/etc/Makeconf"
```
