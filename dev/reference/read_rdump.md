# Read data in an R dump file to a list

Create an R list from an R dump file

## Usage

``` r
read_rdump(f, keep.source = FALSE, ...)
```

## Arguments

- f:

  A character string providing the dump file name.

- keep.source:

  logical: should the source formatting be retained when echoing
  expressions, if possible?

- ...:

  passed to `source`

## Details

The R dump file can be read directly by R function `source`, which by
default would read the data into the user's workspace (the global
environment). This function instead read the data to a list, making it
convenient to prepare data for the `stan` model-fitting function.

## Value

A list containing all the data defined in the dump file with keys
corresponding to variable names.

## See also

[`stan_rdump`](https://mc-stan.org/rstan/dev/reference/stan_rdump.md);
[`dump`](https://rdrr.io/r/base/dump.html)

## Examples

``` r
x <- 1; y <- 1:10; z <- array(1:10, dim = c(2,5)) 
stan_rdump(ls(pattern = '^[xyz]'), file.path(tempdir(), "xyz.Rdump"))
l <- read_rdump(file.path(tempdir(), 'xyz.Rdump'))
print(l)
#> $x
#> [1] 1
#> 
#> $y
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
#> $z
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    3    5    7    9
#> [2,]    2    4    6    8   10
#> 
unlink(file.path(tempdir(), "xyz.Rdump"))
```
