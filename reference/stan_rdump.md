# Dump the data for a Stan model to R dump file in the limited format that Stan can read.

This function takes a vector of names of R objects and outputs text
representations of the objects to a file or connection. The file created
by `stan_rdump` is typically used as data input of the Stan package
(<https://mc-stan.org/>) or
[`source`](https://rdrr.io/r/base/source.html)d into another R session.
The usage of this function is very similar to `dump` in R.

## Usage

``` r
stan_rdump(list, file = "", append = FALSE, 
          envir = parent.frame(),
          width = options("width")$width, 
          quiet = FALSE)
```

## Arguments

- list:

  A vector of character string: the names of one or more R objects to be
  dumped. See the note below.

- file:

  Either a character string naming a file or a
  [connection](https://rdrr.io/r/base/connections.html). `""` indicates
  output to the console.

- append:

  Logical: if `TRUE` and `file` is a character string, output will be
  appended to `file`; otherwise, it will overwrite the contents of
  `file`.

- envir:

  The environment to search for objects.

- width:

  The width for maximum characters on a line. The output is broken into
  lines with `width.`

- quiet:

  Whether to suppress warning messages that would appear when a variable
  is not found or not supported for dumping (not being numeric or it
  would not be converted to numeric) or a variable name is not allowed
  in Stan.

## Value

An invisible character vector containing the names of the objects that
were dumped.

## Note

`stan_rdump` only dumps numeric data, which first can be a scalar,
vector, matrix, or (multidimensional) array. Additional types supported
are `logical` (`TRUE` and `FALSE`), `factor`, `data.frame` and a
specially structured `list`.

The conversion for logical variables is to map `TRUE` to 1 and `FALSE`
to 0. For `factor` variable, function `as.integer` is used to do the
conversion (If we want to transform a factor `f` to approximately its
original numeric values, see the help of function `factor` and do the
transformation before calling `stan_rdump`). In the case of
`data.frame`, function `data.matrix` is applied to the data frame before
dumping. See the notes in
[stan](https://mc-stan.org/rstan/reference/stan.md) for the specially
structured `list`, which will be converted to `array` before dumping.

`stan_rdump` will check whether the names of objects are legal variable
names in Stan. If an illegal name is found, data will be dumped with a
warning. However, passing the name checking does not necessarily mean
that the name is legal. More details regarding rules of variable names
in Stan can be found in Stan's manual.

If objects with specified names are not found, a warning will be issued.

## References

The Stan Development Team *Stan Modeling Language User's Guide and
Reference Manual*. <https://mc-stan.org>.

## See also

[`dump`](https://rdrr.io/r/base/dump.html)

## Examples

``` r
# set variables in global environment
a <- 17.5
b <- c(1,2,3)
# write variables a and b to file ab.data.R in temporary directory
stan_rdump(c('a','b'), file.path(tempdir(), "ab.data.R"))
unlink(file.path(tempdir(), "ab.data.R"))

x <- 1; y <- 1:10; z <- array(1:10, dim = c(2,5)) 
stan_rdump(ls(pattern = '^[xyz]'), file.path(tempdir(), "xyz.Rdump"))
cat(paste(readLines(file.path(tempdir(), "xyz.Rdump")), collapse = '\n'), '\n')
#> x <- 1
#> y <- 
#> c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#> z <- 
#> structure(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
#> .Dim = c(2, 5)) 
unlink(file.path(tempdir(), "xyz.Rdump"))
```
