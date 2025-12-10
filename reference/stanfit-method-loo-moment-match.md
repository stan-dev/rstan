# Moment matching for efficient approximate leave-one-out cross-validation (LOO)

A
[`loo_moment_match`](https://mc-stan.org/loo/reference/loo_moment_match.html)
method that is customized for stanfit objects. The `loo_moment_match`
method for stanfit objects —a wrapper around the
[`loo_moment_match`](https://mc-stan.org/loo/reference/loo_moment_match.html)
(loo package)— updates a loo object using moment matching (Paananen et
al., 2020).

## Usage

``` r
# S3 method for class 'stanfit'
loo_moment_match(x,
    loo = loo,
    ...)
```

## Arguments

- x:

  An object of S4 class `stanfit`.

- loo:

  A loo object that is modified.

- ...:

  Further arguments.

## Value

The `loo_moment_match()` methods return an updated `loo` object.

## References

Paananen, T., Piironen, J., Buerkner, P.-C., Vehtari, A. (2020).
Implicitly Adaptive Importance Sampling. [preprint
arXiv:1906.08850](https://arxiv.org/abs/1906.08850)

## See also

[`loo()`](https://mc-stan.org/rstan/reference/stanfit-method-loo.md),
[`loo_moment_match()`](https://mc-stan.org/rstan/reference/stanfit-method-loo.md)
