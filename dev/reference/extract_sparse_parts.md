# Extract the compressed representation of a sparse matrix

Create a list of vectors that represents a sparse matrix.

## Usage

``` r
extract_sparse_parts(A)
```

## Arguments

- A:

  A [`matrix`](https://rdrr.io/r/base/matrix.html) or
  [`Matrix`](https://rdrr.io/pkg/Matrix/man/Matrix.html).

## Details

The Stan Math Library has a function called `csr_matrix_times_vector`,
which inputs a matrix in compressed row storage form and a dense vector
and returns their product without fillin. To use the
`csr_matrix_times_vector` function with a large sparse matrix, it is
optimal in terms of memory to simply pass the three vectors that
characterize the compressed row storage form of the matrix to the `data`
block of the Stan program. The `extract_sparse_parts` function provides
a convenient means of obtaining these vectors.

## Value

A named list with components

1.  `w` A numeric vector containing the non-zero elements of `A`.

2.  `v` An integer vector containing the column indices of the non-zero
    elements of `A`.

3.  `u` An integer vector indicating where in `w` a given row's non-zero
    values start.

## Examples

``` r
  A <- rbind(
    c(19L, 27L,  0L,  0L),
    c( 0L,  0L,  0L,  0L),
    c( 0L,  0L,  0L, 52L),
    c(81L,  0L, 95L, 33L)
  )
  str(extract_sparse_parts(A))
#> Loading required namespace: Matrix
#> List of 3
#>  $ w: num [1:6] 19 27 52 81 95 33
#>  $ v: int [1:6] 1 2 4 1 3 4
#>  $ u: int [1:5] 1 3 3 4 7
```
