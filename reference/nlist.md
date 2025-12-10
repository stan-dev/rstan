# Created named lists

Create a named list using specified names or, if names are omitted,
using the names of the objects in the list. The code
\]`list(a = a, b = b)` becomes `nlist(a,b)` and `list(a = a, b = 2)`
becomes `nlist(a, b = 2)`, etc. This is convenient when creating the
list of data to pass to Stan.

## Usage

``` r
nlist(...)
```

## Arguments

- ...:

  The objects to include in the list.

## Value

A named list.

## Examples

``` r
# All variables already defined
x <- 1
y <- 2
nlist(x, y)
#> $x
#> [1] 1
#> 
#> $y
#> [1] 2
#> 

# Define some variables in the call and take the rest from the environment
nlist(x, y, z = 3)
#> $x
#> [1] 1
#> 
#> $y
#> [1] 2
#> 
#> $z
#> [1] 3
#> 
```
