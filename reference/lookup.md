# Look up the Stan function that corresponds to a R function or name.

This function helps to map between R functions and Stan functions.

## Usage

``` r
lookup(FUN, ReturnType = character())
```

## Arguments

- FUN:

  A character string naming a R function or a R function for which the
  (near) equivalent Stan function is sought. If no matching R function
  is found, `FUN` is reinterpreted as a
  [`regexp`](https://rdrr.io/r/base/regex.html) and matches are sought.

- ReturnType:

  A character string of positive length naming a valid return type for a
  Stan function: `int`, `int[]`, `matrix`, `real`, `real[,]`, `real[]`,
  `row_vector`, `T[]`, `vector`, or `void`. If `"ANY"` is passed, then
  the entire [`data.frame`](https://rdrr.io/r/base/data.frame.html) is
  returned and can be inspected with the
  [`View`](https://rdrr.io/r/utils/View.html) function, for example.

## Value

Ordinarily, a data.frame with rows equal to the number of partial
matches and four columns:

1.  `StanFunction` Character string for the Stan function's name.

2.  `Arguments` Character string indicating the arguments to that Stan
    function.

3.  `ReturnType` Character string indicating the return type of that
    Stan function.

4.  `Page` Integer indicating the page of the Stan reference manual
    where that Stan function is defined.

If there are no matching Stan functions, a character string indicating
so is returned.

## References

The Stan Development Team *Stan Modeling Language User's Guide and
Reference Manual*. <https://mc-stan.org/>.

The Stan Development Team *CmdStan Interface User's Guide*.
<https://mc-stan.org>.

## Examples

``` r
lookup(dnorm)      # Stan equivalents for the normal PDF (in log form)
#>           StanFunction
#> 415 normal_id_glm_lpdf
#> 418         normal_log
#> 419        normal_lpdf
#> 553    std_normal_lpdf
#>                                                                        Arguments
#> 415 (real, matrix, real, vector, T);(vector, row_vector, vector, vector, vector)
#> 418                                     (real, real, T);(vector, vector, vector)
#> 419                                     (real, real, T);(vector, vector, vector)
#> 553                                                                 (T);(vector)
#>     ReturnType
#> 415     T;real
#> 418     T;real
#> 419     T;real
#> 553     T;real
lookup("foo")      # fails
#> [1] "no matching Stan functions"
lookup("Student")  # succeeds even though there is no such R function
#>                      StanFunction
#> 374 multi_student_t_cholesky_lpdf
#> 375  multi_student_t_cholesky_rng
#> 376           multi_student_t_log
#> 377          multi_student_t_lpdf
#> 378           multi_student_t_rng
#> 557            student_t_ccdf_log
#> 558                 student_t_cdf
#> 559             student_t_cdf_log
#> 560               student_t_lccdf
#> 561                student_t_lcdf
#> 562                 student_t_log
#> 563                student_t_lpdf
#> 564                 student_t_rng
#>                                                     Arguments
#> 374                            (vector, real, vector, matrix)
#> 375 (real, vector, matrix);(real, array[] row_vector, matrix)
#> 376                            (vector, real, vector, matrix)
#> 377                            (vector, real, vector, matrix)
#> 378 (real, vector, matrix);(real, array[] row_vector, matrix)
#> 557    (real, real, real, T);(vector, vector, vector, vector)
#> 558    (real, real, real, T);(vector, vector, vector, vector)
#> 559    (real, real, real, T);(vector, vector, vector, vector)
#> 560    (real, real, real, T);(vector, vector, vector, vector)
#> 561    (real, real, real, T);(vector, vector, vector, vector)
#> 562    (real, real, real, T);(vector, vector, vector, vector)
#> 563    (real, real, real, T);(vector, vector, vector, vector)
#> 564                             (int, int, T);(int, int, int)
#>                ReturnType
#> 374                  real
#> 375 vector;array[] vector
#> 376                  real
#> 377                  real
#> 378 vector;array[] vector
#> 557                T;real
#> 558                T;real
#> 559                T;real
#> 560                T;real
#> 561                T;real
#> 562                T;real
#> 563                T;real
#> 564                T;real
lookup("^poisson") # every Stan function that starts with poisson
#>             StanFunction
#> 459     poisson_ccdf_log
#> 460          poisson_cdf
#> 461      poisson_cdf_log
#> 462        poisson_lccdf
#> 463         poisson_lcdf
#> 464          poisson_log
#> 465 poisson_log_glm_lpmf
#> 466      poisson_log_log
#> 467     poisson_log_lpmf
#> 468      poisson_log_rng
#> 469         poisson_lpmf
#> 470          poisson_rng
#>                                                     Arguments ReturnType
#> 459                                    (int, T);(int, vector)     T;real
#> 460                                    (int, T);(int, vector)     T;real
#> 461                                    (int, T);(int, vector)     T;real
#> 462                                    (int, T);(int, vector)     T;real
#> 463                                    (int, T);(int, vector)     T;real
#> 464                                    (int, T);(int, vector)     T;real
#> 465 (int, matrix, real, vector);(int, matrix, vector, vector)  real;real
#> 466                                    (int, T);(int, vector)     T;real
#> 467                                    (int, T);(int, vector)     T;real
#> 468                                                (T);(real)      T;int
#> 469                                    (int, T);(int, vector)     T;real
#> 470                                                (T);(real)      T;int
```
